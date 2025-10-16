use std::{collections::BTreeMap, fmt::Display};

pub struct Model {
    pub species: Species,
    // N.B. using Vec's here for efficiency, but they should be indexed by newtypes.
    // e.g. MetaboliteId, ReactionId, GeneId
    pub genes: Vec<Gene>,
    pub reactions: Vec<Reaction>,
    pub metabolites: Vec<Metabolite>,
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

pub struct GeneAssociationWithInfo<'a> {
    pub(super) association: &'a GeneAssociation,
    pub(super) info: &'a BTreeMap<GeneIndex, Gene>,
}

#[derive(Debug)]
pub enum AndOp {
    /// Propagates nans
    Min,
    /// Treats nans as 0
    Mean,
    // Because I am representating operations as a binary tree, I do not support median.
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
        info: Option<&BTreeMap<GeneIndex, Gene>>,
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

impl<'a> Display for GeneAssociationWithInfo<'a> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        self.association
            .display_with_gene_info(f, 0, Some(self.info))
    }
}