//! Parses the _mat folder with GSMM info.
//! Contains two ways of parsing the GPR
//!     1. Prefers OR over AND and appears to be correct. This is how the python code does it.
//!     2. Parses the operations, applying operations left to right.

use std::{fs::read_to_string, mem, path::Path};

use itertools::izip;

use crate::model::model::StoichiometricMatrix;

use super::model::{Gene, GeneAssociation, GeneIndex, Metabolite, Model, Reaction, Species};

fn parse_json_file<T: serde::de::DeserializeOwned>(path: &Path) -> T {
    match read_to_string(path) {
        Err(e) => panic!("Failed to read {}: {}", path.display(), e),
        Ok(s) => {
            return serde_json::from_str::<T>(&s).unwrap_or_else(|e| {
                panic!("Failed to parse {}: {}", path.display(), e);
            });
        }
    }
}

pub fn parse_mat_model(top_dir: &Path, species: Species) -> Model {
    let model_dir = top_dir.join("model");

    let genes = parse_json_file::<Vec<String>>(&model_dir.join("model.genes.json"));

    let gtx = parse_json_file::<Vec<u32>>(&top_dir.join("non2uniqueEntrez.json"));
    assert_eq!(genes.len(), gtx.len());

    let (gene_symbols, gene_alt_symbols) = match species {
        Species::HomoSapiens => {
            let gene_symbols = parse_json_file::<Vec<String>>(&top_dir.join("uniqueHumanGeneSymbol.json"));
            // No alt symbols for human in the mat models
            let gene_alt_symbols = Vec::new();
            (gene_symbols, gene_alt_symbols)
        }
        Species::MusMusculus => {
            let gene_symbols = parse_json_file::<Vec<String>>(&top_dir.join("uniqueMouseGeneSymbol.json"));
            let gene_alt_symbols = parse_json_file::<Vec<Vec<String>>>(&top_dir.join("uniqueMouseGeneSymbol_all.json"));
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
                // Note that gene_alt_symbols may be empty, therefore cannot zip
                alt_symbols: gene_alt_symbols
                    .get(non_i as usize)
                    .unwrap_or(&Vec::new())
                    .to_vec(),
            }
        })
        .collect::<Vec<Gene>>();

    let rxns = parse_json_file::<Vec<String>>(&model_dir.join("model.rxns.json"));

    let rxn_names = parse_json_file::<Vec<String>>(&model_dir.join("model.rxnNames.json"));
    assert_eq!(rxns.len(), rxn_names.len());

    let lb = parse_json_file::<Vec<f64>>(&model_dir.join("model.lb.json"));
    assert_eq!(rxns.len(), lb.len());
    let ub = parse_json_file::<Vec<f64>>(&model_dir.join("model.ub.json"));
    assert_eq!(rxns.len(), ub.len());

    let subsystems = parse_json_file::<Vec<String>>(&model_dir.join("model.subSystems.json"));
    assert_eq!(rxns.len(), subsystems.len());

    let rules_text = parse_json_file::<Vec<String>>(&model_dir.join("model.rules.json"));
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
    assert_eq!(rxns.len(), rules.len());

    let reactions: Vec<Reaction> = izip!(
        rxns,
        rxn_names,
        lb,
        ub,
        subsystems,
        rules
    ).map(|(id, name, lb, ub, subsystem, rule)| {
        Reaction {
            id,
            name,
            lb,
            ub,
            subsystem,
            rule,
        }
    }).collect();

    let mets = parse_json_file::<Vec<String>>(&model_dir.join("model.mets.json"));

    let met_formulas = parse_json_file::<Vec<String>>(&model_dir.join("model.metFormulas.json"));
    assert_eq!(mets.len(), met_formulas.len());

    // Not in RECON1_mat, is in RECON2_mat
    /*let kegg_ids = serde_json::from_str::<Vec<String>>(
        &read_to_string(model_dir.join("model.metKeggID.json")).unwrap(),
    )
    .unwrap();
    assert_eq!(mets.len(), kegg_ids.len());*/

    let met_names = parse_json_file::<Vec<String>>(&model_dir.join("model.metNames.json"));
    assert_eq!(mets.len(), met_names.len());

    let metabolites = izip!(mets, met_names, met_formulas)
        .map(|(id, name, formula)| Metabolite {
            id,
            name,
            formula,
        })
        .collect::<Vec<_>>();

    let mut s_mat_coords = parse_json_file::<Vec<(usize, usize, f64)>>(&model_dir.join("model.S.json"));
    // Adjust for 1-based indexing
    for (row, col, _v) in s_mat_coords.iter_mut() {
        assert!(*row > 0 && *row <= metabolites.len(), "Row index {row} out of bounds");
        *row -= 1;
        assert!(*col > 0 && *col <= reactions.len(), "Col index {col} out of bounds");
        *col -= 1;
    }

    Model {
        genes,
        reactions,
        metabolites,
        s_mat: StoichiometricMatrix { coordinates: s_mat_coords },
        species,
    }
}

// TODO: add reference to the text source for debugging?
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Token {
    LeftParen,
    RightParen,
    Or,
    And,
    Gene(GeneIndex),
}

// The token tree applies parentheses, but no other semantics.
#[derive(Debug)]
enum TokenTree {
    Gene(GeneIndex),
    Or,
    And,
    TreeVec(Vec<TokenTree>),
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
                    tokens.push(Token::Gene(GeneIndex { id: num - 1 }));
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
                    prev.push(TokenTree::TreeVec(mem::take(&mut curr)));
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
                TokenTree::TreeVec(nodes) => return Self::tree_to_association(&nodes),
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
                TokenTree::Gene(_) | TokenTree::TreeVec(_) => {
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
                    TokenTree::TreeVec(nodes) => Self::tree_to_association(nodes),
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

#[cfg(test)]
pub mod tests {
    use super::*;

    #[test]
    pub fn test_tokenizer_or() {
        const RULE_OR: &str = "(x(20)) | (x(17)) | (x(19)) | (x(18))";
        const TOKENS_OR: &[Token] = &[
            Token::LeftParen,
            Token::Gene(GeneIndex { id: 19 }),
            Token::RightParen,
            Token::Or,
            Token::LeftParen,
            Token::Gene(GeneIndex { id: 16 }),
            Token::RightParen,
            Token::Or,
            Token::LeftParen,
            Token::Gene(GeneIndex { id: 18 }),
            Token::RightParen,
            Token::Or,
            Token::LeftParen,
            Token::Gene(GeneIndex { id: 17 }),
            Token::RightParen,
        ];

        let tokens = Token::tokenize(RULE_OR);
        assert_eq!(tokens, TOKENS_OR);
    }

    #[test]
    pub fn test_tokenizer_and() {
        const RULE_AND: &str = "(x(21)) & (x(18)) & (x(22)) & (x(16))";
        const TOKENS_AND: &[Token] = &[
            Token::LeftParen,
            Token::Gene(GeneIndex { id: 20 }),
            Token::RightParen,
            Token::And,
            Token::LeftParen,
            Token::Gene(GeneIndex { id: 17 }),
            Token::RightParen,
            Token::And,
            Token::LeftParen,
            Token::Gene(GeneIndex { id: 21 }),
            Token::RightParen,
            Token::And,
            Token::LeftParen,
            Token::Gene(GeneIndex { id: 15 }),
            Token::RightParen,
        ];

        let tokens = Token::tokenize(RULE_AND);
        assert_eq!(tokens, TOKENS_AND);
    }

    #[test]
    pub fn test_tokenizer_both() {
        const RULE_BOTH: &str = "(x(1)) & (x(3)) | (x(4)) & (x(7))";
        const TOKENS_BOTH: &[Token] = &[
            Token::LeftParen,
            Token::Gene(GeneIndex { id: 0 }),
            Token::RightParen,
            Token::And,
            Token::LeftParen,
            Token::Gene(GeneIndex { id: 2 }),
            Token::RightParen,
            Token::Or,
            Token::LeftParen,
            Token::Gene(GeneIndex { id: 3 }),
            Token::RightParen,
            Token::And,
            Token::LeftParen,
            Token::Gene(GeneIndex { id: 6 }),
            Token::RightParen,
        ];

        let tokens = Token::tokenize(RULE_BOTH);
        assert_eq!(tokens, TOKENS_BOTH);
    }

    #[test]
    #[expect(non_snake_case, reason = "using name of chemical reaction")]
    pub fn test_recon1_2OXOADOXm() {
        const EXAMPLE_RULE: &str =
            "(x(972)) & (x(318) & x(1662)) & (x(319)) | (x(973)) & (x(318) & x(1662)) & (x(319))";

        /* Here is the python code result for comparison:
        or: 
            and: 
                gene: 4967.1 {'OGDH'}
                and: 
                    gene: 1738.1 {'DLD'}
                    gene: 8050.1 {'PDHX'}
                gene: 1743.1 {'DLST'}
            and: 
                gene: 4967.2 {'OGDH'}
                and: 
                    gene: 1738.1 {'DLD'}
                    gene: 8050.1 {'PDHX'}
                gene: 1743.1 {'DLST'}
        */
        #[rustfmt::skip]
        let expected : GeneAssociation = 
        GeneAssociation::Or(vec![
            GeneAssociation::And(vec![
                GeneAssociation::Gene(GeneIndex { id: 971 }),
                GeneAssociation::And(vec![
                    GeneAssociation::Gene(GeneIndex { id: 317 }),
                    GeneAssociation::Gene(GeneIndex { id: 1661 }),
                ]),
                GeneAssociation::Gene(GeneIndex { id: 318 }),
            ]),
            GeneAssociation::And(vec![
                GeneAssociation::Gene(GeneIndex { id: 972 }),
                GeneAssociation::And(vec![
                    GeneAssociation::Gene(GeneIndex { id: 317 }),
                    GeneAssociation::Gene(GeneIndex { id: 1661 }),
                ]),
                GeneAssociation::Gene(GeneIndex { id: 318 }),
            ]),
        ]);
        let tokens = Token::tokenize(EXAMPLE_RULE);
        let tree = TokenTree::from_tokens(tokens.into_iter());
        let association = TokenTree::tree_to_association(&tree);
        assert_eq!(association, expected);
    }

    #[test]
    #[expect(non_snake_case, reason = "using name of chemical reaction")]
    pub fn test_recon1_NaKt() {
        const EXAMPLE_RULE: &str =
            "(x(941) & x(938)) | (x(937) & x(451)) | (x(943) & x(936)) | (x(941) & x(940)) | (x(942) & x(938)) | (x(942) & x(936)) | (x(942) & x(937)) | (x(941) & x(936)) | (x(451) & x(936)) | (x(451) & x(940)) | (x(941) & x(937))";

        /* Here is the python code result for comparison:
        or: 
            and: 
                gene: 481.1 {'ATP1B1'}
                gene: 478.1 {'ATP1A3'}
            and: 
                gene: 477.1 {'ATP1A2'}
                gene: 23439.1 {'ATP1B4'}
            and: 
                gene: 483.1 {'ATP1B3'}
                gene: 476.1 {'ATP1A1'}
            and: 
                gene: 481.1 {'ATP1B1'}
                gene: 480.1 {'ATP1A4'}
            and: 
                gene: 482.1 {'ATP1B2'}
                gene: 478.1 {'ATP1A3'}
            and: 
                gene: 482.1 {'ATP1B2'}
                gene: 476.1 {'ATP1A1'}
            and: 
                gene: 482.1 {'ATP1B2'}
                gene: 477.1 {'ATP1A2'}
            and: 
                gene: 481.1 {'ATP1B1'}
                gene: 476.1 {'ATP1A1'}
            and: 
                gene: 23439.1 {'ATP1B4'}
                gene: 476.1 {'ATP1A1'}
            and: 
                gene: 23439.1 {'ATP1B4'}
                gene: 480.1 {'ATP1A4'}
            and: 
                gene: 481.1 {'ATP1B1'}
                gene: 477.1 {'ATP1A2'}
        */
        #[rustfmt::skip]
        let expected : GeneAssociation = 
        GeneAssociation::Or(vec![
            GeneAssociation::And(vec![
                GeneAssociation::Gene(GeneIndex { id: 940 }),
                GeneAssociation::Gene(GeneIndex { id: 937 }),
            ]),
            GeneAssociation::And(vec![
                GeneAssociation::Gene(GeneIndex { id: 936 }),
                GeneAssociation::Gene(GeneIndex { id: 450 }),
            ]),
            GeneAssociation::And(vec![
                GeneAssociation::Gene(GeneIndex { id: 942 }),
                GeneAssociation::Gene(GeneIndex { id: 935 }),
            ]),
            GeneAssociation::And(vec![
                GeneAssociation::Gene(GeneIndex { id: 940 }),
                GeneAssociation::Gene(GeneIndex { id: 939 }),
            ]),
            GeneAssociation::And(vec![
                GeneAssociation::Gene(GeneIndex { id: 941 }),
                GeneAssociation::Gene(GeneIndex { id: 937 }),
            ]),
            GeneAssociation::And(vec![
                GeneAssociation::Gene(GeneIndex { id: 941 }),
                GeneAssociation::Gene(GeneIndex { id: 935 }),
            ]),
            GeneAssociation::And(vec![
                GeneAssociation::Gene(GeneIndex { id: 941 }),
                GeneAssociation::Gene(GeneIndex { id: 936 }),
            ]),
            GeneAssociation::And(vec![
                GeneAssociation::Gene(GeneIndex { id: 940 }),
                GeneAssociation::Gene(GeneIndex { id: 935 }),
            ]),
            GeneAssociation::And(vec![
                GeneAssociation::Gene(GeneIndex { id: 450 }),
                GeneAssociation::Gene(GeneIndex { id: 935 }),
            ]),
            GeneAssociation::And(vec![
                GeneAssociation::Gene(GeneIndex { id: 450 }),
                GeneAssociation::Gene(GeneIndex { id: 939 }),
            ]),
            GeneAssociation::And(vec![
                GeneAssociation::Gene(GeneIndex { id: 940 }),
                GeneAssociation::Gene(GeneIndex { id: 936 }),
            ]),
        ]);
        let tokens = Token::tokenize(EXAMPLE_RULE);
        let tree = TokenTree::from_tokens(tokens.into_iter());
        let association = TokenTree::tree_to_association(&tree);
        assert_eq!(association, expected);
    }

    #[test]
    #[expect(non_snake_case, reason = "using name of chemical reaction")]
    pub fn test_recon1_PNTEH() {
        const EXAMPLE_RULE: &str =
    "(x(1675)) & (x(1676)) & (x(1677)) | (x(1675)) & (x(1676)) & (x(1678)) | (x(1675)) & (x(1679)) & (x(1680)) | (x(1675)) & (x(1678)) & (x(1679)) | (x(1676)) & (x(1675)) & (x(1680)) | (x(1675)) & (x(1677)) & (x(1679))";

        /*  Here is the python code result for comparison:
        or: 
            or: 
                or: 
                    and: 
                        gene: 8876.1 {'VNN1'}
                        gene: 8875.2 {''}
                        gene: 55350.3 {''}
                    and: 
                        gene: 8876.1 {'VNN1'}
                        gene: 8875.2 {''}
                        gene: 55350.2 {''}
                and: 
                    gene: 8876.1 {'VNN1'}
                    gene: 8875.1 {''}
                    gene: 55350.1 {''}
            or: 
                or: 
                    and: 
                        gene: 8876.1 {'VNN1'}
                        gene: 55350.2 {''}
                        gene: 8875.1 {''}
                    and: 
                        gene: 8875.2 {''}
                        gene: 8876.1 {'VNN1'}
                        gene: 55350.1 {''}
                and: 
                    gene: 8876.1 {'VNN1'}
                    gene: 55350.3 {''}
                    gene: 8875.1 {''}
        */
        #[rustfmt::skip]
        let expected : GeneAssociation =
        GeneAssociation::Or(vec![
            GeneAssociation::Or(vec![
                GeneAssociation::Or(vec![
                    GeneAssociation::And(vec![
                        GeneAssociation::Gene(GeneIndex { id: 1674 }),
                        GeneAssociation::Gene(GeneIndex { id: 1675 }),
                        GeneAssociation::Gene(GeneIndex { id: 1676 }),
                    ]),
                    GeneAssociation::And(vec![
                        GeneAssociation::Gene(GeneIndex { id: 1674 }),
                        GeneAssociation::Gene(GeneIndex { id: 1675 }),
                        GeneAssociation::Gene(GeneIndex { id: 1677 }),
                    ]),
                ]),
                GeneAssociation::And(vec![
                    GeneAssociation::Gene(GeneIndex { id: 1674 }),
                    GeneAssociation::Gene(GeneIndex { id: 1678 }),
                    GeneAssociation::Gene(GeneIndex { id: 1679 }),
                ]),
            ]),
            GeneAssociation::Or(vec![
                GeneAssociation::Or(vec![
                    GeneAssociation::And(vec![
                        GeneAssociation::Gene(GeneIndex { id: 1674 }),
                        GeneAssociation::Gene(GeneIndex { id: 1677 }),
                        GeneAssociation::Gene(GeneIndex { id: 1678 }),
                    ]),
                    GeneAssociation::And(vec![
                        GeneAssociation::Gene(GeneIndex { id: 1675 }),
                        GeneAssociation::Gene(GeneIndex { id: 1674 }),
                        GeneAssociation::Gene(GeneIndex { id: 1679 }),
                    ]),
                ]),
                GeneAssociation::And(vec![
                    GeneAssociation::Gene(GeneIndex { id: 1674 }),
                    GeneAssociation::Gene(GeneIndex { id: 1676 }),
                    GeneAssociation::Gene(GeneIndex { id: 1678 }),
                ]),
            ]),
        ]);

        let tokens = Token::tokenize(EXAMPLE_RULE);
        let tree = TokenTree::from_tokens(tokens.into_iter());
        let association = TokenTree::tree_to_association(&tree);
        assert_eq!(association, expected);
    }
}
