use std::{collections::BTreeMap, fmt::Display};

pub struct Model {
    species: Species,
    // TODO: Use a polars dataframe for this instead? Or somet other data structure.
    genes: BTreeMap<GeneId, Gene>,
}

pub enum Species {
    HomoSapiens,
    MusMusculus,
}


#[derive(Debug, Clone)]
pub struct Gene {
    pub(super) id: String,
    pub(super) non_i: u32,
    pub(super) name: String,
    pub(super) alt_symbols: Vec<String>,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub struct GeneId {
    pub(super) id: usize,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum GeneAssociation {
    Gene(GeneId),
    Or(Vec<GeneAssociation>),
    And(Vec<GeneAssociation>),
}

pub struct GeneAssociationWithInfo<'a> {
    pub(super) association: &'a GeneAssociation,
    pub(super) info: &'a BTreeMap<GeneId, Gene>,
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

impl<'a> Display for GeneAssociationWithInfo<'a> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        self.association
            .display_with_gene_info(f, 0, Some(self.info))
    }
}