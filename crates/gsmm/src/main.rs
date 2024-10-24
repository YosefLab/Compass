use std::fs::read_to_string;

use itertools::{izip, Group};

mod mat;

// TODO: Remove hardcoded paths lol
const RECON1_XML_PATH: &str = "../compass/Resources/Metabolic Models/RECON1_xml/RECON1.xml";
const RECON1_MAT_PATH: &str = "../compass/Resources/Metabolic Models/RECON1_mat";
// rust_sbml fails to parse this one due to dc:creator.
// const RECON2_2_PATH: &str = "../compass/Resources/Metabolic Models/RECON2.2/MODEL1603150001.xml";
const RECON2_MAT_PATH: &str = "../compass/Resources/Metabolic Models/RECON2_mat";

pub enum Species {
    HomoSapiens,
    MusMusculus,
}

#[derive(Debug)]
pub struct Gene {
    id: String,
    non_i: u32,
    name: String,
    alt_symbols: Vec<String>,
}

pub fn main() {
    println!("Hello from GSMM!");

    let species = Species::HomoSapiens;

    let top_dir = std::path::PathBuf::from(RECON1_MAT_PATH);
    let model_dir = top_dir.join("model");

    // Technically are numerical, but because they're fixed point we'll treat them as strings for now.
    let genes = serde_json::from_str::<Vec<String>>(
        &read_to_string(model_dir.join("model.genes.json")).unwrap(),
    )
    .unwrap();

    let gtx = serde_json::from_str::<Vec<u32>>(
        &read_to_string(top_dir.join("non2uniqueEntrez.json")).unwrap(),
    )
    .unwrap();
    assert_eq!(genes.len(), gtx.len());

    let (gene_symbols, gene_alt_symbols) = match species {
        Species::HomoSapiens => {
            let gene_symbols: Vec<String> = serde_json::from_str(
                &read_to_string(top_dir.join("uniqueHumanGeneSymbol.json")).unwrap(),
            )
            .unwrap();
            let gene_alt_symbols = Vec::new();
            (gene_symbols, gene_alt_symbols)
        }
        Species::MusMusculus => {
            let gene_symbols: Vec<String> = serde_json::from_str(
                &read_to_string(top_dir.join("uniqueMouseGeneSymbol.json")).unwrap(),
            )
            .unwrap();
            let gene_alt_symbols: Vec<Vec<String>> = serde_json::from_str(
                &read_to_string(top_dir.join("uniqueMouseGeneSymbol_all.json")).unwrap(),
            )
            .unwrap();
            (gene_symbols, gene_alt_symbols)
        }
    };

    let genes = genes
        .iter()
        .zip(gtx)
        .map(|(id, idx)| {
            let non_i = idx - 1;
            Gene {
                id: id.to_string(),
                non_i, // TODO: This index might be unused after this point
                name: gene_symbols[non_i as usize].clone(),
                alt_symbols: gene_alt_symbols
                    .get(non_i as usize)
                    .unwrap_or(&Vec::new())
                    .to_vec(),
            }
        })
        .collect::<Vec<Gene>>();

    println!(
        "Number of genes: {}. First {:?}.",
        genes.len(),
        genes.first()
    );

    let rxns = serde_json::from_str::<Vec<String>>(
        &read_to_string(model_dir.join("model.rxns.json")).unwrap(),
    )
    .unwrap();
    println!(
        "Number of reactions: {}. First {:?}",
        rxns.len(),
        rxns.first()
    );

    let rxn_names = serde_json::from_str::<Vec<String>>(
        &read_to_string(model_dir.join("model.rxnNames.json")).unwrap(),
    )
    .unwrap();
    assert_eq!(rxns.len(), rxn_names.len());

    let lb =
        serde_json::from_str::<Vec<f64>>(&read_to_string(model_dir.join("model.lb.json")).unwrap())
            .unwrap();
    assert_eq!(rxns.len(), lb.len());
    let ub =
        serde_json::from_str::<Vec<f64>>(&read_to_string(model_dir.join("model.ub.json")).unwrap())
            .unwrap();
    assert_eq!(rxns.len(), ub.len());

    let subsystems = serde_json::from_str::<Vec<String>>(
        &read_to_string(model_dir.join("model.subSystems.json")).unwrap(),
    )
    .unwrap();
    assert_eq!(rxns.len(), subsystems.len());

    let rules = serde_json::from_str::<Vec<String>>(
        &read_to_string(model_dir.join("model.rules.json")).unwrap(),
    )
    .unwrap();
    assert_eq!(rxns.len(), rules.len());

    let rules = rules
        .iter()
        .map(|rule| parse_rule(&rule))
        .collect::<Vec<_>>();
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct GeneRuleNode {
    id: u32,
}

// TODO: add reference to the text source for debugging?
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Token {
    LeftParen,
    RightParen,
    Or,
    And,
    Gene(GeneRuleNode),
}

#[derive(Debug, Clone)]
pub enum GeneAssociation {
    Gene(GeneRuleNode),
    Or {
        left: Box<GeneAssociation>,
        right: Box<GeneAssociation>,
    },
    And {
        left: Box<GeneAssociation>,
        right: Box<GeneAssociation>,
    },
}

pub fn tokenize(rule: &str) -> Vec<Token> {
    let mut tokens = Vec::new();

    let mut citer = rule.char_indices();
    while let Some((ind, c)) = citer.next() {
        match c {
            '(' => {
                tokens.push(Token::LeftParen);
            }
            ')' => {
                tokens.push(Token::RightParen);
            }
            ' ' => {} // Ignore extra whitespace.
            'x' => {
                // Look for x(i) where i is a number
                let (_paren_ind, paren_expected) = citer.next().expect("Expected '(' after x");
                assert_eq!(paren_expected, '(');
                let (num_start, _num) = citer.next().expect("Expected number after x(");
                let num_end = 'digits: loop {
                    let (ind, c) = citer.next().expect("Looking for terminating ) after x(");
                    if c == ')' {
                        break 'digits ind;
                    }
                };
                let num_text = &rule[num_start..num_end];
                let num = num_text.parse::<u32>().unwrap_or_else(|e| {
                    panic!("Failed to parse number from {num_text}: {e}\n{rule}");
                });
                tokens.push(Token::Gene(GeneRuleNode { id: num }));
            }
            '|' => {
                tokens.push(Token::Or);
            }
            '&' => {
                tokens.push(Token::And);
            }
            c => {
                panic!("Unexpected character: {c} at index {ind}.\n{rule}");
            }
        }
    }
    tokens
}

pub fn parse_rule(rule: &str) -> Option<GeneAssociation> {
    let tokens = tokenize(rule);
    if tokens.is_empty() {
        return None;
    }

    let parsed = parse_tokens(&tokens);

    //println!("{:?}", parsed);
    parsed
}

enum Operations {
    Or,
    And,
    LeftParen,
}

fn parse_tokens(tokens: &[Token]) -> Option<GeneAssociation> {
    let mut operations = Vec::new();
    let mut output = Vec::new();
    let mut token_iter = tokens.iter();
    loop {
        match token_iter.next() {
            Some(Token::LeftParen) => {
                operations.push(Operations::LeftParen);
                //(res, tokens) = parse_tokens(&tokens[1..]);
                //assert!(tokens.first() == Some(&Token::RightParen));
            }
            Some(Token::RightParen) => {
                // Pop until we find the left paren
                loop {
                    match operations.pop() {
                        Some(Operations::LeftParen) => {
                            break;
                        }
                        Some(op @ Operations::Or) | Some(op @ Operations::And) => {
                            let right = Box::new(output.pop().expect("Expected right operand"));
                            let left = Box::new(output.pop().expect("Expected left operand"));
                            let op = match op {
                                Operations::Or => GeneAssociation::Or { left, right },
                                Operations::And => GeneAssociation::And { left, right },
                                _ => unreachable!(),
                            };
                            output.push(op);
                        }
                        None => {
                            panic!("Unmatched right paren");
                        }
                    }
                }
            }
            Some(Token::Or) => {
                operations.push(Operations::Or);
            }
            Some(Token::And) => {
                operations.push(Operations::And);
            }
            Some(Token::Gene(gene_rule_node)) => {
                output.push(GeneAssociation::Gene(*gene_rule_node));
            }
            None => {
                break;
            }
        }
    }
    loop {
        match operations.pop() {
            Some(op @ Operations::Or) | Some(op @ Operations::And) => {
                let right = Box::new(output.pop().expect("Expected right operand"));
                let left = Box::new(output.pop().expect("Expected left operand"));
                let op = match op {
                    Operations::Or => GeneAssociation::Or { left, right },
                    Operations::And => GeneAssociation::And { left, right },
                    _ => unreachable!(),
                };
                output.push(op);
            }
            None => {
                break;
            }
            Some(Operations::LeftParen) => {
                panic!("Unmatched left paren");
            }
        }
    }
    assert!(output.len() <= 1);
    output.first().cloned()
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
