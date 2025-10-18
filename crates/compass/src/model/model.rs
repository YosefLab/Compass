use std::{collections::BTreeMap, fmt::Display};

pub struct Model {
    pub species: Species,
    // N.B. using Vec's here for efficiency, but they should be indexed by newtypes.
    // e.g. MetaboliteId, ReactionId, GeneId
    pub genes: Vec<Gene>,
    pub reactions: Vec<Reaction>,
    pub metabolites: Vec<Metabolite>,
    pub s_mat: StoichiometricMatrix,
}

pub struct StoichiometricMatrix {
    // using COO format for now
    pub(super) coordinates: Vec<(usize, usize, f64)>, // (row, col, value)
}

#[derive(Debug)]
pub enum Species {
    HomoSapiens,
    MusMusculus,
}

#[derive(Debug)]
pub struct Metabolite {
    pub id: String,
    pub name: String,
    pub formula: String,
}

#[derive(Debug)]
pub struct Reaction {
    pub lb: f64,
    pub ub: f64,
    // TODO: replace String with a newtype?
    pub id: String,
    pub name: String,
    pub subsystem: String,
    // reactants and products are stored in S matrix, so leave them out for now.
    //pub reactants: Vec<StoichiometricValue>,
    //pub products: Vec<StoichiometricValue>,
    pub rule: Option<GeneAssociation>,
    // reverse reaction representation?
}

#[derive(Debug)]
pub struct StoichiometricValue {
    pub metabolite: MetaboliteIndex,
    // should be positive.
    pub coefficient: f64,
}

#[derive(Debug, Clone)]
pub struct Gene {
    pub(super) id: String,
    pub(super) non_i: u32,
    pub(super) name: String,
    pub(super) alt_symbols: Vec<String>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub struct GeneIndex {
    pub(super) id: usize,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub struct MetaboliteIndex {
    pub(super) id: usize,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum GeneAssociation {
    Gene(GeneIndex),
    Or(Vec<GeneAssociation>),
    And(Vec<GeneAssociation>),
}

pub struct GeneAssociationEvaluator<'a> {
    pub gene_expr: &'a BTreeMap<GeneIndex, f64>,
    pub or_op: OrOp,
    pub and_op: AndOp,
}

pub struct GeneAssociationWithInfo<'a> {
    pub association: &'a GeneAssociation,
    pub info: &'a Vec<Gene>,
}

#[derive(Debug)]
pub enum AndOp {
    /// Propagates nans
    Min,
    /// Treats nans as 0
    Mean,
}

#[derive(Debug)]
pub enum OrOp {
    Sum,
}

impl GeneAssociation {
    fn display_with_gene_info(
        &self,
        f: &mut std::fmt::Formatter<'_>,
        depth: usize,
        info: Option<&Vec<Gene>>,
    ) -> std::fmt::Result {
        write!(f, "{}", " ".repeat(depth * 4))?;
        match self {
            GeneAssociation::Gene(gene_id) => match info.map(|m| &m[gene_id.id]) {
                Some(g) => writeln!(f, "{gene_id:?}: {g:?}")?,
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

impl<'a> Display for GeneAssociationWithInfo<'a> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        self.association
            .display_with_gene_info(f, 0, Some(self.info))
    }
}

impl<'a> GeneAssociationEvaluator<'a> {
    pub fn evaluate(&self, node: &GeneAssociation) -> Option<f64> {
        match node {
            GeneAssociation::Gene(gene_id) => self.gene_expr.get(gene_id).copied(),
            GeneAssociation::Or(vec) => {
                Some(self.or_op.apply(vec.iter().map(|x| self.evaluate(x))))
            }
            GeneAssociation::And(vec) => {
                Some(self.and_op.apply(vec.iter().map(|x| self.evaluate(x))))
            }
        }
    }
}

impl OrOp {
    pub fn apply(&self, operands: impl Iterator<Item = Option<f64>>) -> f64 {
        let nan_to_zero = |x: f64| if x.is_nan() { 0.0 } else { x };
        match self {
            OrOp::Sum => operands.map(|x| x.map(nan_to_zero).unwrap_or(0.0)).sum(),
        }
    }
}

impl AndOp {
    pub fn apply(&self, operands: impl Iterator<Item = Option<f64>>) -> f64 {
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

impl StoichiometricMatrix {
    pub fn nnz(&self) -> usize {
        self.coordinates.len()
    }
}

#[cfg(test)]
mod tests {
    use rand::{prelude::*};
    use rand_xoshiro::Xoroshiro128Plus;

    use super::*;

    #[test]
    pub fn test_binary_eval() {
        let mut rng = Xoroshiro128Plus::seed_from_u64(133771331);
        let rule = GeneAssociation::Or(vec![
            GeneAssociation::Gene(GeneIndex { id: 0 }),
            GeneAssociation::And(vec![
                GeneAssociation::Gene(GeneIndex { id: 1 }),
                GeneAssociation::Gene(GeneIndex { id: 2 }),
            ]),
            GeneAssociation::Gene(GeneIndex { id: 3 }),
        ]);
        let test_inner = |vals: [f64; 4]| {
            let expr = vals
                .iter()
                .copied()
                .enumerate()
                .map(|(i, v)| (GeneIndex { id: i }, v))
                .collect();
            let evaluator = GeneAssociationEvaluator {
                gene_expr: &expr,
                or_op: OrOp::Sum,
                and_op: AndOp::Mean,
            };
            let result = evaluator.evaluate(&rule);
            let expected = vals[0] + (vals[1] + vals[2]) / 2.0 + vals[3];
            assert_eq!(result, Some(expected));
        };
        for _ in 0..10 {
            let vals = std::array::from_fn(|_| rng.random_range(0..100usize) as f64 / 10.0);
            println!("Testing with vals: {:?}", vals);
            test_inner(vals);
        }
    }

    #[test]
    pub fn test_nested_eval() {
        let mut rng = Xoroshiro128Plus::seed_from_u64(133771331);
        let rule = GeneAssociation::And(vec![
            GeneAssociation::Gene(GeneIndex { id: 0 }),
            GeneAssociation::Or(vec![
                GeneAssociation::Gene(GeneIndex { id: 0 }),
                GeneAssociation::Gene(GeneIndex { id: 1 }),
                GeneAssociation::Gene(GeneIndex { id: 2 }),
            ]),
            GeneAssociation::And(vec![
                GeneAssociation::Gene(GeneIndex { id: 3 }),
                GeneAssociation::Gene(GeneIndex { id: 4 }),
            ]),
            GeneAssociation::Or(vec![
                GeneAssociation::Gene(GeneIndex { id: 1 }),
                GeneAssociation::Gene(GeneIndex { id: 3 }),
            ]),
        ]);
        let test_inner = |vals: [f64; 5]| {
            let expr = vals
                .iter()
                .copied()
                .enumerate()
                .map(|(i, v)| (GeneIndex { id: i }, v))
                .collect();
            let evaluator = GeneAssociationEvaluator {
                gene_expr: &expr,
                or_op: OrOp::Sum,
                and_op: AndOp::Mean,
            };
            let result = evaluator.evaluate(&rule);
            let expected = (vals[0]
                + (vals[0] + vals[1] + vals[2])
                + ((vals[3] + vals[4]) / 2.0)
                + (vals[1] + vals[3]))
                / 4.0;
            assert_eq!(result, Some(expected));
        };
        for _ in 0..10 {
            let vals = std::array::from_fn(|_| rng.random_range(0..100usize) as f64 / 10.0);
            println!("Testing with vals: {:?}", vals);
            test_inner(vals);
        }
    }
}
