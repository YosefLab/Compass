mod cli;
//mod penalties;

use std::path::{Path, PathBuf};

use clap::Parser;
use cli::{Commands, CompassCli};
use gsmm::{
    mat::parse_mat_model,
    model::{self, Model, Species},
};
use polars::prelude::*;

fn main() {
    let args = CompassCli::parse();
    match args.command {
        Commands::Debug {
            input_data,
            gene_column,
            model,
        } => {
            println!("Debugging with input data: {:?}", input_data);
            let data = read_expression_data(&input_data, gene_column);
            let model = load_model(model);
            al_gore_rhythm(data, model);
        }
    }
    // Read the expression file.

    // Preprocess the data
}

pub struct ExpressionData {
    data: LazyFrame,
    gene_column: String,
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
    let gene_col = match gene_column {
        Some(col) => col,
        None => {
            let schema = data.clone().first().collect_schema().unwrap();
            let (col, dtype) = schema.get_at_index(0).expect("No columns in the schema");
            if let DataType::String = dtype {
                col.clone().to_string()
            } else {
                panic!("The first column {col:?} is dtype {dtype:?} rather than String. Expected a list of gene symbols.");
            }
        }
    };
    let data = data.group_by([col(&gene_col)]).agg([col("*").sum()]);
    ExpressionData {
        data,
        gene_column: gene_col,
    }
}

pub fn load_model(model: cli::MetabolicModel) -> Model {
    // TODO: remove hardcoded paths
    let path = match model {
        cli::MetabolicModel::Recon1Mat => "../compass/Resources/Metabolic Models/RECON1_mat",
        cli::MetabolicModel::Recon2Mat => "../compass/Resources/Metabolic Models/RECON2_mat",
    };
    let path = PathBuf::from(path);
    println!(
        "Loading model from: {}",
        path.file_name().unwrap().to_str().unwrap()
    );
    // TODO: stop hardcoding mouse species
    parse_mat_model(&path, Species::MusMusculus)
}

pub fn al_gore_rhythm(mut data: ExpressionData, model: Model) {
    println!("{:#?}", data.data.clone().first().collect().unwrap());
    // This is a somewhat ugly way to construct things, but whatever
    // Create a full list of the symbols along with the associated gene index
    // Join the gene index frame with the expression data
    // Then sort and aggregate the data by gene index
    let mut symbol_vec = Vec::new();
    let mut index_vec = Vec::new();
    for (i, gene) in model.genes().iter().enumerate() {
        for symbol in gene.get_symbols() {
            symbol_vec.push(symbol.to_lowercase());
            index_vec.push(i as u32);
        }
    }
    const GENE_INDEX_KEY: &str = "gene_index";
    let gene_df = DataFrame::new(vec![
        Column::new(data.gene_column.clone().into(), symbol_vec),
        Column::new(GENE_INDEX_KEY.into(), index_vec),
    ])
    .unwrap()
    .lazy();
    println!("{:#?}", gene_df.clone().collect().unwrap());
    let tmp_df = data.data.clone().join(
        gene_df,
        [col(&data.gene_column).str().to_lowercase()],
        [col(&data.gene_column)],
        JoinArgs::new(JoinType::Inner),
    );
    println!("{:#?}", tmp_df.collect().unwrap());
}
