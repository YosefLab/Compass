//! Parses the _mat folder with GSMM info.
//! Contains two ways of parsing the GPR
//!     1. Prefers OR over AND and appears to be correct. This is how the python code does it.
//!     2. Parses the operations, applying operations left to right.

use std::{fs::read_to_string, mem, path::Path};

use crate::model::{Gene, GeneAssociation, GeneId, Model, Species};

pub fn parse_mat_model(top_dir: &Path, species: Species) -> Model {
    let model_dir = top_dir.join("model");

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
            // No alt symbols for human in the mat models
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

    let rules_text = serde_json::from_str::<Vec<String>>(
        &read_to_string(model_dir.join("model.rules.json")).unwrap(),
    )
    .unwrap();
    assert_eq!(rxns.len(), rules_text.len());

    let rules = rules_text
        .iter()
        .map(|ja| {
            if ja.len() > 0 {
                Some(TokenTree::tree_to_association(&TokenTree::from_tokens(
                    Token::tokenize(ja).into_iter(),
                )))
            } else {
                None
            }
        })
        .collect::<Vec<_>>();
    todo!()
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

#[derive(Debug)]
enum TokenTree {
    Gene(GeneId),
    Or,
    And,
    Tree(Vec<TokenTree>),
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

impl Token {
    fn tokenize(rule: &str) -> Vec<Self> {
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
                    let num = num_text.parse::<usize>().unwrap_or_else(|e| {
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
}

impl TokenTree {
    /// Parse in the same way compass does currently
    /// Prioritize OR over AND and use the middle OR to keep trees balanced.
    fn from_tokens(tokens: impl Iterator<Item = Token>) -> Vec<Self> {
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

    fn tree_to_association(tree: &[Self]) -> GeneAssociation {
        if tree.len() == 1 {
            match &tree[0] {
                TokenTree::Gene(gene_id) => return GeneAssociation::Gene(gene_id.clone()),
                TokenTree::Tree(nodes) => return Self::tree_to_association(&nodes),
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
            let left = Self::tree_to_association(&tree[..index]);
            let right = Self::tree_to_association(&tree[index + 1..]);

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
                    TokenTree::Tree(nodes) => Self::tree_to_association(nodes),
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
}

impl GeneAssociationBinary {
    fn parse_tokens(mut tokens: impl Iterator<Item = Token>) -> Option<Self> {
        // Internal type only for the operations stack.
        enum Operations {
            Or,
            And,
            LeftParen,
        }
        let mut operations = Vec::new();
        let mut output = Vec::new();
        loop {
            match tokens.next() {
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
                    output.push(GeneAssociationBinary::Gene(gene_id.clone()));
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

    fn flatten(&self) -> GeneAssociation {
        match self {
            GeneAssociationBinary::Gene(gene_id) => GeneAssociation::Gene(*gene_id),
            GeneAssociationBinary::Or { left, right } => {
                let mut operands = Vec::new();
                let left = left.flatten();
                if let GeneAssociation::Or(mut operands_left) = left {
                    operands.append(&mut operands_left);
                } else {
                    operands.push(left);
                }
                let right = right.flatten();
                if let GeneAssociation::Or(mut operands_right) = right {
                    operands.append(&mut operands_right);
                } else {
                    operands.push(right);
                }
                GeneAssociation::Or(operands)
            }
            GeneAssociationBinary::And { left, right } => {
                let mut operands = Vec::new();
                let left = left.flatten();
                if let GeneAssociation::And(mut operands_left) = left {
                    operands.append(&mut operands_left);
                } else {
                    operands.push(left);
                }
                let right = right.flatten();
                if let GeneAssociation::And(mut operands_right) = right {
                    operands.append(&mut operands_right);
                } else {
                    operands.push(right);
                }
                GeneAssociation::And(operands)
            }
        }
    }
}
