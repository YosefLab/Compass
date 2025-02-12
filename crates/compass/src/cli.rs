use std::path::PathBuf;

use clap::{Args, Parser, Subcommand, ValueEnum};
use gsmm::model::Species;

#[derive(Parser)]
pub struct CompassCli {
    #[command(subcommand)]
    pub command: Commands,
}

#[derive(Subcommand)]
pub enum Commands {
    /// Developer thingy to go test stuff
    Demo {
        input_data: PathBuf,
        /// Optional column name for the gene symbols.
        #[arg(short, long)]
        gene_column: Option<String>,
        #[command(flatten)]
        model: MetabolicModelConfig,
        #[arg(short, long)]
        output: PathBuf,
    },
    ModelDebug {
        #[command(flatten)]
        model: MetabolicModelConfig,
        #[command(subcommand)]
        mode: ModelDebugMode,
    },
}

#[derive(Args, Debug, Clone, PartialEq, Eq)]
pub struct MetabolicModelConfig {
    #[arg(value_enum, short, long)]
    pub model: MetabolicModel,
    #[arg(value_enum, short, long, default_value_t = Species::MusMusculus)]
    pub species: Species,
    #[arg(long, default_value_t = false)]
    pub isoform_summing: bool,
}

#[derive(Debug, Copy, Clone, PartialEq, Eq, ValueEnum)]
pub enum MetabolicModel {
    Recon2Mat,
    Recon1Mat,
}

#[derive(Clone, PartialEq, Eq, Subcommand)]
pub enum ModelDebugMode {
    /// List all reactions in the model.
    ListReactions,
    /// List info about a reaction in the model.
    ReactionInfo { reaction: String },
}

impl MetabolicModel {
    pub fn name(&self) -> &'static str {
        match self {
            MetabolicModel::Recon2Mat => "Recon2Mat",
            MetabolicModel::Recon1Mat => "Recon1Mat",
        }
    }
}
