use std::path::PathBuf;

use clap::{Args, Parser, Subcommand};
use compass::model::{
    mat::parse_mat_model,
    model::{Species},
};

// TODO: Remove hardcoded paths lol
#[allow(dead_code)]
const RECON1_XML_PATH: &str = "../compass/Resources/Metabolic Models/RECON1_xml/RECON1.xml";
#[allow(dead_code)]
const RECON1_MAT_PATH: &str = "../compass/Resources/Metabolic Models/RECON1_mat";
// rust_sbml fails to parse this one due to dc:creator.
// const RECON2_2_PATH: &str = "../compass/Resources/Metabolic Models/RECON2.2/MODEL1603150001.xml";
#[allow(dead_code)]
const RECON2_MAT_PATH: &str = "../compass/Resources/Metabolic Models/RECON2_mat";

#[derive(Parser)]
#[command(name = "compass-gsmm")]
#[command(about = "COMPASS - Characterizing metabolism from single-cell RNA-seq", long_about = None)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand)]
enum Commands {
    /// Debug GSMM model parsing and structure
    DebugGsmm(DebugGsmmArgs),
}

#[derive(Args)]
struct DebugGsmmArgs {
    /// Path to the model directory (for MAT format) or file (for XML format)
    #[arg(short, long)]
    model_path: PathBuf,

    /// Species for the model
    #[arg(short, long, default_value = "homo-sapiens")]
    species: String,

    /// Model format (mat or xml)
    #[arg(short, long, default_value = "mat")]
    format: String,
}

pub fn main() {
    let cli = Cli::parse();

    match &cli.command {
        Commands::DebugGsmm(args) => debug_gsmm(args),
    }
}

fn debug_gsmm(args: &DebugGsmmArgs) {
    println!("Debugging GSMM model...");
    println!("Model path: {}", args.model_path.display());
    println!("Species: {}", args.species);
    println!("Format: {}", args.format);

    let species = match args.species.to_lowercase().as_str() {
        "homo-sapiens" | "human" => Species::HomoSapiens,
        "mus-musculus" | "mouse" => Species::MusMusculus,
        _ => {
            eprintln!("Unknown species: {}. Using HomoSapiens as default.", args.species);
            Species::HomoSapiens
        }
    };

    match args.format.to_lowercase().as_str() {
        "mat" => {
            let model = parse_mat_model(&args.model_path, species);
            println!(
                "\nParsed model with:\n  {} genes\n  {} reactions\n  {} metabolites\n  {} non-zero S matrix entries",
                model.genes.len(),
                model.reactions.len(),
                model.metabolites.len(),
                model.s_mat.nnz(),
            );
        }
        "xml" => {
            println!("XML parsing not yet implemented");
        }
        _ => {
            eprintln!("Unknown format: {}. Supported formats: mat, xml", args.format);
        }
    }
}

pub fn sbml_parse() {
    let model_text = std::fs::read_to_string(RECON1_XML_PATH).unwrap();
    let now = std::time::Instant::now();
    let recon1 = rust_sbml::Model::parse(&model_text).unwrap();
    println!("Parsed model in {:?}", now.elapsed());

    println!("Number of reactions: {}", recon1.reactions.len());
    println!("Number of metabolites: {}", recon1.species.len());
    recon1.reactions.iter().take(5).for_each(|r| {
        println!("Reaction: {r:#?}");
    });

    println!("Parsed models");
}
