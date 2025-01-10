use std::{collections::BTreeMap, fmt::Display, fs::read_to_string, mem};

use itertools::{izip, Group};

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

#[derive(Debug, Clone)]
pub struct Gene {
    id: String,
    non_i: u32,
    name: String,
    alt_symbols: Vec<String>,
}

pub fn main() {
    println!("Hello from GSMM!");

    let species = Species::HomoSapiens;

    let top_dir = std::path::PathBuf::from(RECON2_MAT_PATH);
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

    let index = rxns.iter().position(|r| r == "2OXOADOXm").unwrap();
    //println!("{index}: {}", &rules[index]);

    let rules_0 = rules
        .iter()
        .map(|rule| parse_rule(&rule).as_ref().map(GeneAssociation::collapse))
        .collect::<Vec<_>>();
    let rule_0 = rules_0[index].as_ref().unwrap();

    let rules_1 = rules
        .iter()
        .map(|rule| {
            if rule.len() > 0 {
                Some(tree_to_association(&tokens_to_tree(&tokenize(rule))))
            } else {
                None
            }
        })
        .collect::<Vec<_>>();
    let rule_1 = rules_1[index].as_ref().unwrap();

    let gene_info = genes
        .iter()
        .enumerate()
        .map(|(id, g)| (GeneId { id: id as u32 }, g.clone()))
        .collect();

    println!(
        "{}",
        GeneAssociationWithInfo {
            association: rule_0,
            info: &gene_info,
        }
    );
    println!(
        "{}",
        GeneAssociationWithInfo {
            association: rule_1,
            info: &gene_info,
        }
    );
    /*let mut count = 0;
    for (i, (a, b)) in rules_0.iter().zip(rules_1.iter()).enumerate() {
        if a != b {
            count += 1;
            if let (Some(x), Some(y)) = (a,b) {
                println!("{}", rules[i]);
                println!("{}", GeneAssociationWithInfo {
                    association: x,
                    info: &gene_info,
                });
                println!("{}", GeneAssociationWithInfo {
                    association: y,
                    info: &gene_info,
                });
                return;
            }

        }
    }
    println!("{count} out of {}", rules_0.len());*/
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub struct GeneId {
    id: u32,
}

// TODO: add reference to the text source for debugging?
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Token {
    LeftParen,
    RightParen,
    Or,
    And,
    Gene(GeneId),
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum GeneAssociationBinary {
    Gene(GeneId),
    Or {
        left: Box<GeneAssociationBinary>,
        right: Box<GeneAssociationBinary>,
    },
    And {
        left: Box<GeneAssociationBinary>,
        right: Box<GeneAssociationBinary>,
    },
}

pub struct GeneAssociationWithInfo<'a> {
    association: &'a GeneAssociation,
    info: &'a BTreeMap<GeneId, Gene>,
}

impl<'a> Display for GeneAssociationWithInfo<'a> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        self.association
            .display_with_gene_info(f, 0, Some(self.info))
    }
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum GeneAssociation {
    Gene(GeneId),
    Or(Vec<GeneAssociation>),
    And(Vec<GeneAssociation>),
}

impl GeneAssociation {
    pub fn collapse(other: &GeneAssociationBinary) -> Self {
        match other {
            GeneAssociationBinary::Gene(gene_id) => Self::Gene(*gene_id),
            GeneAssociationBinary::Or { left, right } => {
                let mut operands = Vec::new();
                let left = Self::collapse(left);
                if let Self::Or(mut operands_left) = left {
                    operands.append(&mut operands_left);
                } else {
                    operands.push(left);
                }
                let right = Self::collapse(right);
                if let Self::Or(mut operands_right) = right {
                    operands.append(&mut operands_right);
                } else {
                    operands.push(right);
                }
                Self::Or(operands)
            }
            GeneAssociationBinary::And { left, right } => {
                let mut operands = Vec::new();
                let left = Self::collapse(left);
                if let Self::And(mut operands_left) = left {
                    operands.append(&mut operands_left);
                } else {
                    operands.push(left);
                }
                let right = Self::collapse(right);
                if let Self::And(mut operands_right) = right {
                    operands.append(&mut operands_right);
                } else {
                    operands.push(right);
                }
                Self::And(operands)
            }
        }
    }

    fn display_with_gene_info(
        &self,
        f: &mut std::fmt::Formatter<'_>,
        depth: usize,
        info: Option<&BTreeMap<GeneId, Gene>>,
    ) -> std::fmt::Result {
        write!(f, "{}", " ".repeat(depth * 4))?;
        match self {
            GeneAssociation::Gene(gene_id) => match info.map(|m| m.get(gene_id)) {
                Some(Some(g)) => writeln!(f, "{gene_id:?}: {g:?}")?,
                Some(None) => writeln!(f, "{gene_id:?}: Missing info")?,
                None => writeln!(f, "{gene_id:?}")?,
            },
            GeneAssociation::And(vec) | GeneAssociation::Or(vec) => {
                if matches!(self, GeneAssociation::And(..)) {
                    writeln!(f, "and")?
                } else {
                    writeln!(f, "or")?
                }
                for expr in vec {
                    expr.display_with_gene_info(f, depth + 1, info)?;
                }
            }
        }
        Ok(())
    }
}

impl Display for GeneAssociation {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        self.display_with_gene_info(f, 0, None)
    }
}

#[derive(Debug)]
pub enum OrOp {
    Sum,
}

impl OrOp {
    pub fn apply(&self, left: Option<f64>, right: Option<f64>) -> f64 {
        let nan_to_zero = |x: f64| if x.is_nan() { 0.0 } else { x };
        match self {
            OrOp::Sum => {
                let left = left.map(nan_to_zero).unwrap_or(0.0);
                let right = right.map(nan_to_zero).unwrap_or(0.0);
                left + right
            }
        }
    }

    pub fn apply_many(&self, operands: impl Iterator<Item = Option<f64>>) -> f64 {
        let nan_to_zero = |x: f64| if x.is_nan() { 0.0 } else { x };
        match self {
            OrOp::Sum => operands.map(|x| x.map(nan_to_zero).unwrap_or(0.0)).sum(),
        }
    }
}

#[derive(Debug)]
pub enum AndOp {
    /// Propagates nans
    Min,
    /// Treats nans as 0
    Mean,
    // Because I am representating operations as a binary tree, I do not support median.
}

impl AndOp {
    pub fn apply(&self, left: Option<f64>, right: Option<f64>) -> f64 {
        let nan_to_zero = |x: f64| if x.is_nan() { 0.0 } else { x };
        match self {
            AndOp::Min => {
                let left = left.unwrap_or(0.0);
                let right = right.unwrap_or(0.0);
                left.min(right)
            }
            AndOp::Mean => {
                let left = left.map(nan_to_zero).unwrap_or(0.0);
                let right = right.map(nan_to_zero).unwrap_or(0.0);
                (left + right) / 2.0
            }
        }
    }

    pub fn apply_many(&self, operands: impl Iterator<Item = Option<f64>>) -> f64 {
        let nan_to_zero = |x: f64| if x.is_nan() { 0.0 } else { x };
        match self {
            AndOp::Min => operands
                .map(|x| x.unwrap_or(0.0))
                .fold(f64::INFINITY, |acc, val| acc.min(val)),
            AndOp::Mean => {
                let (mut sum, mut len) = (0.0, 0);
                for operand in operands {
                    sum += operand.map(nan_to_zero).unwrap_or(0.0);
                    len += 1;
                }
                if len == 0 {
                    0.0
                } else {
                    sum / len as f64
                }
            }
        }
    }
}

pub struct GeneAssociationEvaluator<'a> {
    gene_expr: &'a BTreeMap<GeneId, f64>,
    or_op: OrOp,
    and_op: AndOp,
}

impl GeneAssociationEvaluator<'_> {
    pub fn evaluate_binary(&self, node: &GeneAssociationBinary) -> Option<f64> {
        match node {
            GeneAssociationBinary::Gene(gene) => self.gene_expr.get(gene).copied(),
            GeneAssociationBinary::Or { left, right } => {
                let left = self.evaluate_binary(left);
                let right = self.evaluate_binary(right);
                Some(self.or_op.apply(left, right))
            }
            GeneAssociationBinary::And { left, right } => {
                let left = self.evaluate_binary(left);
                let right = self.evaluate_binary(right);
                Some(self.and_op.apply(left, right))
            }
        }
    }

    pub fn evaluate(&self, node: &GeneAssociation) -> Option<f64> {
        match node {
            GeneAssociation::Gene(gene_id) => self.gene_expr.get(gene_id).copied(),
            GeneAssociation::Or(vec) => {
                Some(self.or_op.apply_many(vec.iter().map(|x| self.evaluate(x))))
            }
            GeneAssociation::And(vec) => {
                Some(self.and_op.apply_many(vec.iter().map(|x| self.evaluate(x))))
            }
        }
    }
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
                // The numbers in the json are treated as 1-indexed, so subtract that here
                tokens.push(Token::Gene(GeneId { id: num - 1 }));
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

pub fn parse_rule(rule: &str) -> Option<GeneAssociationBinary> {
    let tokens = tokenize(rule);
    if tokens.is_empty() {
        return None;
    }

    let parsed = parse_tokens(&tokens);
    parsed
}

enum Operations {
    Or,
    And,
    LeftParen,
}

fn parse_tokens(tokens: &[Token]) -> Option<GeneAssociationBinary> {
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
                                Operations::Or => GeneAssociationBinary::Or { left, right },
                                Operations::And => GeneAssociationBinary::And { left, right },
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
            Some(Token::Gene(gene_id)) => {
                output.push(GeneAssociationBinary::Gene(*gene_id));
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
                    Operations::Or => GeneAssociationBinary::Or { left, right },
                    Operations::And => GeneAssociationBinary::And { left, right },
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

#[derive(Debug)]
enum TokenTree {
    Gene(GeneId),
    Or,
    And,
    Tree(Vec<TokenTree>),
}

/// Parse in the same way compass does currently
/// Prioritize OR over AND and use the middle OR to keep trees balanced.
fn tokens_to_tree(tokens: &[Token]) -> Vec<TokenTree> {
    let mut stack: Vec<Vec<TokenTree>> = Vec::new();
    let mut curr: Vec<TokenTree> = Vec::new();
    for token in tokens {
        match token {
            Token::LeftParen => {
                // New group, put curr on the stack of groups. Hopefully we find a right paren later.
                stack.push(mem::take(&mut curr));
            }
            Token::RightParen => {
                // End of group. We need there to be some previous left paren.
                let mut prev = match stack.pop() {
                    Some(prev) => prev,
                    v => panic!("Expected a group, not {v:?}"),
                };
                prev.push(TokenTree::Tree(mem::take(&mut curr)));
                curr = prev;
            }
            Token::Or => curr.push(TokenTree::Or),
            Token::And => curr.push(TokenTree::And),
            Token::Gene(gene_id) => curr.push(TokenTree::Gene(gene_id.clone())),
        }
    }

    assert!(stack.len() == 0, "Unclosed left paren");
    curr
}

fn tree_to_association(tree: &[TokenTree]) -> GeneAssociation {
    if tree.len() == 1 {
        match &tree[0] {
            TokenTree::Gene(gene_id) => return GeneAssociation::Gene(gene_id.clone()),
            TokenTree::Tree(nodes) => return tree_to_association(&nodes),
            TokenTree::Or | TokenTree::And => panic!("Invalid Token Tree"),
        }
    }
    // Iterate over all of the tokens to determine which operators appear in a group.
    // Also verify the operators and values alternate.
    let (mut or_count, mut and_count) = (0, 0);
    for (i, node) in tree.iter().enumerate() {
        match node {
            TokenTree::Or | TokenTree::And => {
                if matches!(node, TokenTree::Or) {
                    or_count += 1;
                } else {
                    and_count += 1
                }
                assert_eq!(i % 2, 1, "{i} % 2 != 1");
            }
            TokenTree::Gene(_) | TokenTree::Tree(_) => {
                assert_eq!(i % 2, 0, "{i} % 2 != 0");
            }
        }
    }

    if or_count > 0 && and_count > 0 {
        // Partition on a middle or operation and recurse.
        let (index, _) = tree
            .iter()
            .enumerate()
            .filter(|(_i, node)| matches!(node, TokenTree::Or))
            .nth(or_count / 2)
            .expect("There must exist such an or");
        let left = tree_to_association(&tree[..index]);
        let right = tree_to_association(&tree[index + 1..]);

        GeneAssociation::Or(vec![left, right])
    } else {
        assert!(
            or_count + and_count > 0,
            "Node without operation {or_count} {and_count}. {tree:?}"
        );
        let nodes = tree
            .iter()
            .step_by(2)
            .map(|node| match node {
                TokenTree::Gene(gene_id) => GeneAssociation::Gene(gene_id.clone()),
                TokenTree::Tree(nodes) => tree_to_association(nodes),
                TokenTree::Or | TokenTree::And => unreachable!(),
            })
            .collect();
        if or_count > 0 {
            GeneAssociation::Or(nodes)
        } else {
            GeneAssociation::And(nodes)
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

#[cfg(test)]
mod test {
    use super::*;

    const RULE_OR: &str = "(x(20)) | (x(17)) | (x(19)) | (x(18))";
    const TOKENS_OR: &[Token] = &[
        Token::LeftParen,
        Token::Gene(GeneId { id: 20 }),
        Token::RightParen,
        Token::Or,
        Token::LeftParen,
        Token::Gene(GeneId { id: 17 }),
        Token::RightParen,
        Token::Or,
        Token::LeftParen,
        Token::Gene(GeneId { id: 19 }),
        Token::RightParen,
        Token::Or,
        Token::LeftParen,
        Token::Gene(GeneId { id: 18 }),
        Token::RightParen,
    ];

    const RULE_AND: &str = "(x(21)) & (x(18)) & (x(22)) & (x(16))";
    const TOKENS_AND: &[Token] = &[
        Token::LeftParen,
        Token::Gene(GeneId { id: 21 }),
        Token::RightParen,
        Token::And,
        Token::LeftParen,
        Token::Gene(GeneId { id: 18 }),
        Token::RightParen,
        Token::And,
        Token::LeftParen,
        Token::Gene(GeneId { id: 22 }),
        Token::RightParen,
        Token::And,
        Token::LeftParen,
        Token::Gene(GeneId { id: 16 }),
        Token::RightParen,
    ];

    const RULE_BOTH: &str = "(x(1)) & (x(3)) | (x(4)) & (x(7))";
    const TOKENS_BOTH: &[Token] = &[
        Token::LeftParen,
        Token::Gene(GeneId { id: 1 }),
        Token::RightParen,
        Token::And,
        Token::LeftParen,
        Token::Gene(GeneId { id: 3 }),
        Token::RightParen,
        Token::Or,
        Token::LeftParen,
        Token::Gene(GeneId { id: 4 }),
        Token::RightParen,
        Token::And,
        Token::LeftParen,
        Token::Gene(GeneId { id: 7 }),
        Token::RightParen,
    ];

    #[test]
    pub fn test_tokenizer() {
        let tokens = super::tokenize(RULE_OR);
        assert_eq!(tokens, TOKENS_OR);
        let tokens = super::tokenize(RULE_AND);
        assert_eq!(tokens, TOKENS_AND);
        let tokens = super::tokenize(RULE_BOTH);
        assert_eq!(tokens, TOKENS_BOTH);
    }

    #[test]
    pub fn test_parser() {
        let or_expected = GeneAssociationBinary::Or {
            left: Box::new(GeneAssociationBinary::Gene(GeneId { id: 20 })),
            right: Box::new(GeneAssociationBinary::Or {
                left: Box::new(GeneAssociationBinary::Gene(GeneId { id: 17 })),
                right: Box::new(GeneAssociationBinary::Or {
                    left: Box::new(GeneAssociationBinary::Gene(GeneId { id: 19 })),
                    right: Box::new(GeneAssociationBinary::Gene(GeneId { id: 18 })),
                }),
            }),
        };
        let or_parsed = parse_tokens(&TOKENS_OR).unwrap();
        assert_eq!(or_parsed, or_expected);

        let and_expected = GeneAssociationBinary::And {
            left: Box::new(GeneAssociationBinary::Gene(GeneId { id: 21 })),
            right: Box::new(GeneAssociationBinary::And {
                left: Box::new(GeneAssociationBinary::Gene(GeneId { id: 18 })),
                right: Box::new(GeneAssociationBinary::And {
                    left: Box::new(GeneAssociationBinary::Gene(GeneId { id: 22 })),
                    right: Box::new(GeneAssociationBinary::Gene(GeneId { id: 16 })),
                }),
            }),
        };
        let and_parsed = parse_tokens(&TOKENS_AND).unwrap();
        assert_eq!(and_parsed, and_expected);

        let both_expected = GeneAssociationBinary::And {
            left: Box::new(GeneAssociationBinary::Gene(GeneId { id: 1 })),
            right: Box::new(GeneAssociationBinary::Or {
                left: Box::new(GeneAssociationBinary::Gene(GeneId { id: 3 })),
                right: Box::new(GeneAssociationBinary::And {
                    left: Box::new(GeneAssociationBinary::Gene(GeneId { id: 4 })),
                    right: Box::new(GeneAssociationBinary::Gene(GeneId { id: 7 })),
                }),
            }),
        };
        let both_parsed = parse_tokens(&TOKENS_BOTH).unwrap();
        assert_eq!(both_parsed, both_expected);
    }

    #[test]
    pub fn test_eval() {
        let expr = [(1, 5.0), (3, 2.0), (4, 3.0), (7, 6.0)]
            .iter()
            .copied()
            .map(|(id, val)| (GeneId { id }, val))
            .collect::<BTreeMap<_, _>>();
        let evaluator = GeneAssociationEvaluator {
            gene_expr: &expr,
            or_op: OrOp::Sum,
            and_op: AndOp::Mean,
        };
        let rule = parse_rule(RULE_BOTH).unwrap();
        // Hmm, this seems a bit off. I suppose my AST is constructed assuming associativity,
        // but the OR function is not neccesarily associative.
        // Ie ((x + y) / 2 + z) / 2 is not the same as (x + (y + z) / 2) / 2
        assert_eq!(
            evaluator.evaluate_binary(&rule),
            Some((5.0 + (2.0 + (3.0 + 6.0) / 2.0)) / 2.0)
        );

        let rule_flat = GeneAssociation::collapse(&rule);
        assert_eq!(
            evaluator.evaluate(&rule_flat),
            Some((5.0 + (2.0 + (3.0 + 6.0) / 2.0)) / 2.0)
        );
    }

    #[test]
    pub fn test_eval_associative() {
        let expr = [(1, 5.0), (2, 2.0), (3, 6.0)]
            .iter()
            .copied()
            .map(|(id, val)| (GeneId { id }, val))
            .collect::<BTreeMap<_, _>>();
        let evaluator = GeneAssociationEvaluator {
            gene_expr: &expr,
            or_op: OrOp::Sum,
            and_op: AndOp::Mean,
        };
        let rule = parse_rule("(x(1)) & (x(2)) & (x(3))").unwrap();
        const VALUE: f64 = (((6.0 + 2.0) / 2.0) + 5.0) / 2.0;
        assert_eq!(evaluator.evaluate_binary(&rule), Some(VALUE));

        let rule_flat = GeneAssociation::collapse(&rule);
        const FLAT_VALUE: f64 = (5.0 + 2.0 + 6.0) / 3.0;
        assert_eq!(evaluator.evaluate(&rule_flat), Some(FLAT_VALUE));

        assert!(VALUE != FLAT_VALUE);
    }
}
