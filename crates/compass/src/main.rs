mod cli;
mod model_debug;
mod penalties;

use std::{
    collections::BTreeSet,
    path::{Path, PathBuf},
};

use clap::Parser;
use cli::{Commands, CompassCli, ModelDebugMode};
use gsmm::{
    mat::parse_mat_model,
    model::{Model, ModelConfig},
};
use polars::prelude::*;
use tracing::info;
use tracing_subscriber::{fmt, prelude::*, EnvFilter};

fn main() {
    let args = CompassCli::parse();

    tracing_subscriber::registry()
        .with(fmt::layer().with_file(true).with_line_number(true))
        .with(EnvFilter::from_default_env())
        .init();

    match args.command {
        Commands::Demo {
            input_data,
            gene_column,
            model,
            output,
        } => {
            info!("Input data: {:?}", input_data);
            let data = read_expression_data(&input_data, gene_column);
            let model = load_model(model);
            al_gore_rhythm(
                data,
                model,
                PaddingMode::Zero,
                GeneResolutionMode::PreferPrimary,
                output,
            );
        }
        Commands::ModelDebug { model, mode } => {
            let model = load_model(model);
            match mode {
                ModelDebugMode::ListReactions => model_debug::list_reactions(model),
                ModelDebugMode::ReactionInfo { reaction } => {
                    model_debug::reaction_info(model, reaction)
                }
            }
        }
    }
    // Read the expression file.

    // Preprocess the data
}

pub struct ExpressionData {
    data: LazyFrame,
    gene_column: PlSmallStr,
    data_columns: Vec<PlSmallStr>,
}

/// Read expression data from `path` into a polars DataFrame.
/// Aggregates columns by gene_column
/// Assumes the file is a variety of CSV:
///  - If the file extension is ".tsv" then assume the separator is tabs instead
/// # Panics
/// If anything goes wrong lol.
pub fn read_expression_data(path: &Path, gene_column: Option<String>) -> ExpressionData {
    let file_extension = path.extension().unwrap().to_str().unwrap();
    let sep = match file_extension {
        "tsv" => b'\t',
        _ => b',',
    };
    let data = LazyCsvReader::new(path)
        .with_separator(sep)
        .finish()
        .unwrap();
    let schema = data.clone().first().collect_schema().unwrap();
    let gene_col = match gene_column {
        Some(col) => col.into(),
        None => {
            let (col, dtype) = schema.get_at_index(0).expect("No columns in the schema");
            if let DataType::String = dtype {
                col.clone()
            } else {
                panic!("The first column {col:?} is dtype {dtype:?} rather than String. Expected a list of gene symbols.");
            }
        }
    };
    info!("Column with gene symbols detected as: {}", gene_col);
    // Treating all columns of type f64 as data columns
    // TODO: Could use a more robust/configurable method to detect data columns
    let data_columns = schema
        .iter()
        .filter_map(|(s, dtype)| {
            if matches!(dtype, DataType::Float64) {
                Some(s.clone())
            } else {
                None
            }
        })
        .collect::<Vec<_>>();
    info!("Detected {} data columns", data_columns.len());
    // For all the data columns, sum the values by gene
    let data = data
        .group_by([col(gene_col.clone())])
        .agg([cols(data_columns.clone()).sum()]);
    //warn!("Truncating data columns to 2 for debugging");
    //let data_columns: Vec<_> = data_columns.into_iter().take(2).collect();
    //let data = data.select([col(gene_col.clone()), cols(data_columns.clone())]);
    ExpressionData {
        data,
        gene_column: gene_col,
        data_columns,
    }
}

pub fn load_model(model: cli::MetabolicModelConfig) -> Model {
    // TODO: remove hardcoded paths
    let path = match model.model {
        cli::MetabolicModel::Recon1Mat => "../compass/Resources/Metabolic Models/RECON1_mat",
        cli::MetabolicModel::Recon2Mat => "../compass/Resources/Metabolic Models/RECON2_mat",
    };
    let path = PathBuf::from(path);
    info!(
        "Loading model from: {}",
        path.file_name().unwrap().to_str().unwrap()
    );
    // TODO: stop hardcoding mouse species
    parse_mat_model(
        ModelConfig {
            model_name: model.model.name().to_owned(),
            species: model.species,
            isoform_summing: model.isoform_summing,
        },
        &path,
    )
}

#[derive(Debug, Clone, Copy)]
pub enum PaddingMode {
    Zero,
    NaN,
}

pub enum GeneResolutionMode {
    /// Sum gene expression for all symbols of a gene
    SumAll,
    /// Look up gene expression by:
    /// 1. The primary symbol of the gene
    /// If 1 is not found, then:
    /// 2. The mean of all alternative symbols
    /// Returns NaN if no symbols are found
    /// Note this is what the python code does.
    PreferPrimary,
}

pub fn al_gore_rhythm(
    data: ExpressionData,
    model: Model,
    padding_mode: PaddingMode,
    gene_mode: GeneResolutionMode,
    output: PathBuf,
) {
    let df = gene_expr(&data, &model, padding_mode, gene_mode).cache();

    // TODO: Other metabolic and_functions?
    // For now, default to AND=mean and OR=sum. Which should be doable as a matrix multiplication, no?
    let mut rxn_expr_cols = Vec::new();
    let rxn_names = model
        .reactions()
        .iter()
        .map(|rxn| rxn.name().to_string())
        .collect::<Series>()
        .with_name("_reaction_name".into())
        .into();
    rxn_expr_cols.push(rxn_names);
    // Loop columns first, as the data columns are the memory-expensive parts.
    // TODO: Rayon here?
    for sample in data.data_columns.clone() {
        // TODO: With polars, should I be cloning the lazy frame like this?
        // Then I can do this with only one column in memory at a time
        let col_df = df.clone().select([col(sample.clone())]).collect().unwrap();
        assert!(col_df.width() == 1);
        // The column height is limited by the number of genes in the model
        // So should be fine to just collect it all at once
        let col = col_df.column(&sample).unwrap().f64().unwrap();
        // Handling nulls here again technically.
        let data = col.iter().map(|f| f.unwrap_or(0.0)).collect::<Vec<_>>();
        info!(
            "Sample {} has {} or {} genes",
            sample,
            data.len(),
            col.len()
        );
        let eval = penalties::GeneRuleEval {
            gene_expr: &data,
            or_op: penalties::GeneOrSum,
            and_op: penalties::GeneAndMean,
        };
        let rxn_expr: Series = model
            .reactions()
            .iter()
            .map(|rxn| {
                let res = if let Some(rule) = rxn.rule() {
                    eval.evaluate(rule)
                } else {
                    // Hmm, is this the correct way to fill things without a rule?
                    0.0
                };
                if rxn.name() == "r0739" && sample == "Ob-DHA-e_S154_L007_R1_001" {
                    eval.evaluate_debug(rxn.rule().unwrap(), true);
                }
                res
            })
            .map(|f| (f + 1.0).log2())
            .map(|f| 1.0 / (1.0 + f))
            .collect();
        rxn_expr_cols.push(rxn_expr.with_name(sample).into());
    }
    let mut rxn_expr_df = DataFrame::new(rxn_expr_cols).unwrap();
    println!("{:#?}", rxn_expr_df);
    let penalty_file = std::fs::File::create(output).unwrap();
    CsvWriter::new(penalty_file)
        .finish(&mut rxn_expr_df)
        .unwrap();

    /*let tmp_df = data
        .data
        .clone()
        .join(
            gene_df,
            [col(data.gene_column).str().to_lowercase()],
            [col(data.gene_column)],
            JoinArgs::new(JoinType::Inner),
        )
        .group_by([col(GENE_INDEX_KEY)])
        .agg([col(GENE_NAME).sum(), cols()])
        .sort([GENE_INDEX_KEY], Default::default());
    println!("{:#?}", tmp_df.collect().unwrap());*/
    // Correction:
    // Then join again on the gene_df, padding all columns with zeros for missing genes.
    // Replace gene symbol with canonical one from the model?
    // Yeah for debugging keep around list of them.
}

fn gene_expr(
    data: &ExpressionData,
    model: &Model,
    padding_mode: PaddingMode,
    gene_mode: GeneResolutionMode,
) -> LazyFrame {
    println!("{:#?}", data.data.clone().first().collect().unwrap());
    // This is a somewhat ugly way to construct things, but whatever
    // The goal here is to transform the data such that only genes from the metabolic model are present
    // And they any missing genes are padded with zeros (or NaNs, if configured that way)
    // Then the gene evaluation can be done by a simple indexing operation.

    // Create a full list of the symbols the model can process along with the associated gene index
    // Join the gene index frame with the expression data
    // Then sort and aggregate the data by gene index so they are in the same order as the model.
    let mut symbol_vec = Vec::new();
    let mut is_primary_symbol = Vec::new();
    let mut index_vec = Vec::new();
    println!("Number of genes {}", model.genes().len());
    for (i, gene) in model.genes().iter().enumerate() {
        if let Some(name) = gene.name() {
            for symbol in gene.get_symbols() {
                if name == symbol.as_str() {
                    is_primary_symbol.push(Some(true));
                } else {
                    is_primary_symbol.push(Some(false));
                }
                symbol_vec.push(Some(symbol.to_lowercase()));
                index_vec.push(i as u32);
            }
        } else {
            symbol_vec.push(None);
            is_primary_symbol.push(None);
            index_vec.push(i as u32);
        }
    }

    /*let mut seen_indices =
        BTreeSet::from_iter((0..model.genes().len()).into_iter().map(|x| x as u32));
    for index in index_vec.iter() {
        seen_indices.remove(index);
    }
    for index in seen_indices {
        // We need the index to have the correct number of genes around
        symbol_vec.push(None);
        is_primary_symbol.push(None);
        index_vec.push(index);
    }*/

    // TODO: Alternative naming scheme for reserved keys?
    const GENE_INDEX_KEY: &str = "_gene_index";
    const GENE_PRIMARY_KEY: &str = "_gene_primary";
    let gene_df = DataFrame::new(vec![
        Column::new(data.gene_column.clone().into(), symbol_vec),
        Column::new(GENE_INDEX_KEY.into(), index_vec.clone()),
        Column::new(GENE_PRIMARY_KEY.into(), is_primary_symbol),
    ])
    .unwrap();
    println!("Gene df: {gene_df:#?}");
    let gene_df = gene_df.lazy();
    /*let debug_gene_df = gene_df.clone().collect().unwrap();
    println!("{debug_gene_df:#?}");
    println!(
        "{:#?}",
        debug_gene_df
            .lazy()
            .filter(col(GENE_PRIMARY_KEY))
            .collect()
            .unwrap()
    );*/

    // Turn the gene column into a list of strings
    // Keep around the regular one for the join
    const GENE_LIST_KEY: &str = "_gene_list";
    let df = data.data.clone().select([
        col(data.gene_column.clone())
            .cast(DataType::List(Box::new(DataType::String)))
            .alias(GENE_LIST_KEY),
        col("*"),
    ]);

    /*let debug_stuff = df
        .clone()
        .filter(col(data.gene_column.clone()).eq(lit("cyp2c29")))
        //.filter(col(GENE_INDEX_KEY).eq(lit("CYP2C29")))
        .select([col("Ob-DHA-e_S154_L007_R1_001")])
        .collect()
        .unwrap();
    info!("Debug {:#?}", debug_stuff);*/
    let df = df.join(
        gene_df,
        [col(data.gene_column.clone()).str().to_lowercase()],
        [col(data.gene_column.clone()).str().to_lowercase()],
        JoinArgs::new(JoinType::Right),
    );
    // Here we
    // 1. Select the gene symbols that appear in the model
    // 2. Sort by the model's gene index
    // 3. Sum the data columns by the gene index
    let null_filler = match padding_mode {
        PaddingMode::Zero => 0.0,
        PaddingMode::NaN => f64::NAN,
    };
    let expr = match gene_mode {
        GeneResolutionMode::SumAll => df
            .sort([GENE_INDEX_KEY], Default::default())
            .group_by([col(GENE_INDEX_KEY)])
            .agg([
                cols(data.data_columns.clone()).fill_null(null_filler).sum(),
                col(GENE_LIST_KEY).flatten(),
            ]),
        GeneResolutionMode::PreferPrimary => {
            // First, add the secondary expression to the data frame
            // This is the mean of all alternative symbols
            let primary_expr = df
                .clone()
                .filter(col(GENE_PRIMARY_KEY))
                .sort([GENE_INDEX_KEY], Default::default())
                .group_by([col(GENE_INDEX_KEY)])
                .agg([
                    when(cols(data.data_columns.clone()).count().gt(0))
                        .then(cols(data.data_columns.clone()).sum())
                        .otherwise(lit(NULL)),
                    col(GENE_LIST_KEY).flatten(),
                ]);

            let secondary_expr = df
                .clone()
                .filter(not(col(GENE_PRIMARY_KEY)))
                .sort([GENE_INDEX_KEY], Default::default())
                .group_by([col(GENE_INDEX_KEY)])
                .agg([
                    // Drop NaNs to ensure they aren't included in the mean
                    cols(data.data_columns.clone()).drop_nans().mean(),
                    col(GENE_LIST_KEY).flatten(),
                ]);

            println!(
                "Primary expr: {:#?}",
                primary_expr.clone().collect().unwrap().shape()
            );

            println!(
                "Secondary expr: {:#?}",
                secondary_expr.clone().collect().unwrap().shape()
            );

            let index_df = DataFrame::new(vec![Column::new(GENE_INDEX_KEY.into(), index_vec)])
                .unwrap()
                .lazy();

            let combined_expr = primary_expr
                .join(
                    secondary_expr,
                    [col(GENE_INDEX_KEY)],
                    [col(GENE_INDEX_KEY)],
                    JoinArgs::new(JoinType::Full),
                )
                // Need to padd gene index s.t. all genes are present
                .join(
                    index_df,
                    [col(GENE_INDEX_KEY)],
                    [col(GENE_INDEX_KEY)],
                    JoinArgs::new(JoinType::Right),
                )
                /*.select([
                    col("*").exclude_dtype([DataType::Float64]),
                    dtype_cols([DataType::Float64]).fill_null(null_filler),
                ])*/
                .sort([GENE_INDEX_KEY], Default::default());

            let combined_expr = combined_expr.select(
                data.data_columns
                    .iter()
                    .map(|col_name| {
                        col(col_name.clone())
                            .fill_null(col(format!("{}_right", col_name)))
                            .alias(col_name.clone())
                    })
                    .collect::<Vec<_>>(),
            );

            combined_expr

            /*let secondary_expr = df
                .clone()
                .filter(not(col(GENE_PRIMARY_KEY)))
                .sort([GENE_INDEX_KEY], Default::default())
                .group_by([col(GENE_INDEX_KEY)])
                .agg([
                    cols(data.data_columns.clone())
                        .fill_null(null_filler)
                        .mean(),
                    col(GENE_LIST_KEY).flatten(),
                ]);
            // Then we need to fill nulls in the primary expression with the secondary expression

            let primary_expr = df
            .clone()
            .select(when(col(GENE_PRIMARY_KEY))
                .then(col("*"))
                .otherwise()
                )
            .filter(col(GENE_PRIMARY_KEY))
            .sort([GENE_INDEX_KEY], Default::default()).fill_null(secondary_expr.);
            todo!();*/
        }
    };
    println!("Expr: {:#?}", expr.clone().collect().unwrap());
    expr
}
