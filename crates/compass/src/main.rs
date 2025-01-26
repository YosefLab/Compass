mod cli;
mod penalties;

use std::path::{Path, PathBuf};

use clap::Parser;
use cli::{Commands, CompassCli};
use gsmm::{
    mat::parse_mat_model,
    model::{Model, Species},
};
use polars::prelude::*;
use tracing::info;
use tracing_subscriber::{fmt, prelude::*, EnvFilter};

fn main() {
    let args = CompassCli::parse();

    tracing_subscriber::registry()
        .with(fmt::layer())
        .with(EnvFilter::from_default_env())
        .init();

    match args.command {
        Commands::Debug {
            input_data,
            gene_column,
            model,
        } => {
            info!("Input data: {:?}", input_data);
            let data = read_expression_data(&input_data, gene_column);
            let model = load_model(model);
            al_gore_rhythm(data, model, PaddingMode::Zero);
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
    let data = data.group_by([col(gene_col.clone())]).agg([col("*").sum()]);
    ExpressionData {
        data,
        gene_column: gene_col,
        data_columns,
    }
}

pub fn load_model(model: cli::MetabolicModel) -> Model {
    // TODO: remove hardcoded paths
    let path = match model {
        cli::MetabolicModel::Recon1Mat => "../compass/Resources/Metabolic Models/RECON1_mat",
        cli::MetabolicModel::Recon2Mat => "../compass/Resources/Metabolic Models/RECON2_mat",
    };
    let path = PathBuf::from(path);
    info!(
        "Loading model from: {}",
        path.file_name().unwrap().to_str().unwrap()
    );
    // TODO: stop hardcoding mouse species
    parse_mat_model(&path, Species::MusMusculus)
}

#[derive(Debug, Clone, Copy)]
pub enum PaddingMode {
    Zero,
    NaN,
}

pub fn al_gore_rhythm(mut data: ExpressionData, model: Model, padding_mode: PaddingMode) {
    println!("{:#?}", data.data.clone().first().collect().unwrap());
    // This is a somewhat ugly way to construct things, but whatever
    // The goal here is to transform the data such that only genes from the metabolic model are present
    // And they any missing genes are padded with zeros (or NaNs, if configured that way)
    // Then the gene evaluation can be done by a simple indexing operation.

    // Create a full list of the symbols the model can process along with the associated gene index
    // Join the gene index frame with the expression data
    // Then sort and aggregate the data by gene index
    let mut symbol_vec = Vec::new();
    let mut index_vec = Vec::new();
    println!("Number of genes {}", model.genes().len());
    for (i, gene) in model.genes().iter().enumerate() {
        for symbol in gene.get_symbols() {
            symbol_vec.push(symbol.to_lowercase());
            index_vec.push(i as u32);
        }
    }

    // TODO: Alternative naming scheme for reserved keys?
    const GENE_INDEX_KEY: &str = "_gene_index";
    let gene_df = DataFrame::new(vec![
        Column::new(data.gene_column.clone().into(), symbol_vec),
        Column::new(GENE_INDEX_KEY.into(), index_vec),
    ])
    .unwrap()
    .lazy();
    println!("{:#?}", gene_df.clone().collect().unwrap());

    // Turn the gene column into a list of strings
    // Keep around the regular one for the join
    const GENE_LIST_KEY: &str = "_gene_list";
    let df = data.data.clone().select([
        col(data.gene_column.clone())
            .cast(DataType::List(Box::new(DataType::String)))
            .alias(GENE_LIST_KEY),
        col("*"),
    ]);

    // Here we
    // 1. Select the gene symbols that appear in the model
    // 2. Sort by the model's gene index
    // 3. Sum the data columns by the gene index
    let df = df
        .join(
            gene_df,
            [col(data.gene_column.clone()).str().to_lowercase()],
            [col(data.gene_column.clone()).str().to_lowercase()],
            JoinArgs::new(JoinType::Right),
        )
        .sort([GENE_INDEX_KEY], Default::default())
        .group_by([col(GENE_INDEX_KEY)])
        .agg([
            // TODO: Not propagating the NaNs here. Unsure if that's the best approach.
            cols(data.data_columns.clone())
                .fill_null(0.0)
                .fill_nan(0.0)
                .sum(),
            col(GENE_LIST_KEY).flatten(),
        ]);
    // TODO: Other metabolic and_functions?
    // For now, default to AND=mean and OR=sum. Which should be doable as a matrix multiplication, no?
    let df = df.cache();
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
    for gene in data.data_columns.clone() {
        // TODO: With polars, should I be cloning the lazy frame like this?
        // Then I can do this with only one column in memory at a time
        let col_df = df.clone().select([col(gene.clone())]).collect().unwrap();
        assert!(col_df.width() == 1);
        // The column height is limited by the number of genes in the model
        // So should be fine to just collect it all at once
        let col = col_df.column(&gene).unwrap().f64().unwrap();
        // Handling nulls here again technically.
        let data = col.iter().map(|f| f.unwrap_or(0.0)).collect::<Vec<_>>();
        let eval = penalties::GeneRuleEval {
            gene_expr: &data,
            or_op: penalties::GeneOrSum,
            and_op: penalties::GeneAndMean,
        };
        let rxn_expr: Series = model
            .reactions()
            .iter()
            .map(|rxn| {
                if let Some(rule) = rxn.rule() {
                    eval.evaluate(rule)
                } else {
                    // Hmm, is this the correct way to fill things without a rule?
                    0.0
                }
            })
            .map(|f| (f + 1.0).log2())
            .map(|f| 1.0 / (1.0 + f))
            .collect();
        rxn_expr_cols.push(rxn_expr.with_name(gene).into());
    }
    let mut rxn_expr_df = DataFrame::new(rxn_expr_cols).unwrap();
    println!("{:#?}", rxn_expr_df);
    let penalty_file = std::fs::File::create(PathBuf::from("penalties.csv")).unwrap();
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
