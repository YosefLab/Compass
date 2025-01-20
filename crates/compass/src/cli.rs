use std::path::PathBuf;

use clap::{Parser, Subcommand, ValueEnum};

#[derive(Parser)]
pub struct CompassCli {
    #[command(subcommand)]
    pub command: Commands,
}

#[derive(Subcommand)]
pub enum Commands {
    /// Developer thingy to go test stuff
    Debug {
        input_data: PathBuf,
        /// Optional column name for the gene symbols.
        #[arg(short, long)]
        gene_column: Option<String>,
        /// Metabolic model to use
        #[arg(short, long)]
        model: MetabolicModel,
    },
}

#[derive(Copy, Clone, PartialEq, Eq, ValueEnum)]
pub enum MetabolicModel {
    Recon2Mat,
    Recon1Mat,
}
