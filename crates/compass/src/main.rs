use compass::model::{
    mat::parse_mat_model,
    model::{Species},
};

// TODO: Remove hardcoded paths lol
const RECON1_XML_PATH: &str = "../compass/Resources/Metabolic Models/RECON1_xml/RECON1.xml";
const RECON1_MAT_PATH: &str = "../compass/Resources/Metabolic Models/RECON1_mat";
// rust_sbml fails to parse this one due to dc:creator.
// const RECON2_2_PATH: &str = "../compass/Resources/Metabolic Models/RECON2.2/MODEL1603150001.xml";
const RECON2_MAT_PATH: &str = "../compass/Resources/Metabolic Models/RECON2_mat";

pub fn main() {
    println!("Hello from GSMM!");

    let species = Species::MusMusculus;
    let top_dir = std::path::PathBuf::from(RECON1_MAT_PATH);
    let model = parse_mat_model(&top_dir, species);

    println!(
        "Parsed model with {} genes, {} reactions, {} metabolites, and {} non-zero S matrix entries",
        model.genes.len(),
        model.reactions.len(),
        model.metabolites.len(),
        model.s_mat.nnz(),
    );
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
