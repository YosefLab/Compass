use std::{
    collections::{BTreeMap, BTreeSet},
    fmt::Display,
};

pub struct Model {
    pub config: ModelConfig,
    /// A list of genes.
    pub(super) genes: Vec<Gene>,
    pub(super) reactions: Vec<Reaction>,
    pub(super) metabolites: Vec<Metabolite>,
    pub(super) subsystems: Vec<Subsystem>,
    pub(super) s_matrix: StoichiometricMatrix,
}

#[derive(Debug)]
pub struct ModelConfig {
    pub model_name: String,
    pub species: Species,
    pub isoform_summing: bool,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
#[cfg_attr(feature = "clap", derive(clap::ValueEnum))]
pub enum Species {
    HomoSapiens,
    MusMusculus,
}

/// A chemical reaction.
/// For the reactants, use the S matrix.
#[derive(Debug)]
pub struct Reaction {
    pub(super) name: String,
    pub(super) rule: Option<GeneAssociation>,
    pub(super) subsystem: SubsystemId,
    /// Lower bound
    pub(super) lb: f64,
    /// Upper bound.
    pub(super) ub: f64,
}

#[derive(Debug, Clone)]
pub struct Gene {
    /// The entrez? gene id. Pretty sure its entrez.
    pub(super) entrez: String,
    pub(super) non_i: u32,
    pub(super) symbols: Vec<String>,
}

#[derive(Debug, Clone)]
pub struct Metabolite {
    pub(super) name: String,
    pub(super) display_name: Option<String>,
    pub(super) formula: String,
    // TODO: Add kegg id or something?
}

#[derive(Debug, Clone)]
pub struct Subsystem {
    pub(super) name: String,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub struct GeneId {
    pub(super) id: usize,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub struct MetabId {
    pub(super) index: usize,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord)]
pub struct SubsystemId {
    pub(super) index: usize,
}

pub struct StoichiometricMatrix {
    // TODO: struct of arrays here?
    // This is rxn
    pub(super) data: Vec<StoichiometricEntry>,
}

pub(crate) struct StoichiometricEntry {
    pub(super) rxn: usize,
    pub(super) metab: MetabId,
    pub(super) value: f64,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum GeneAssociation {
    Gene(GeneId),
    Or(Vec<GeneAssociation>),
    And(Vec<GeneAssociation>),
}

pub struct GeneAssociationWithInfo<'a> {
    pub(super) association: &'a GeneAssociation,
    pub(super) info: &'a Vec<Gene>,
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

impl GeneId {
    /// Returns the index of the gene in the list of metabolic model's genes.
    pub fn index(&self) -> usize {
        self.id
    }
}

impl GeneAssociation {
    /// Note that `info` is expected to be vector where the index is the gene id.
    fn display_with_gene_info(
        &self,
        f: &mut std::fmt::Formatter<'_>,
        depth: usize,
        info: Option<&Vec<Gene>>,
    ) -> std::fmt::Result {
        write!(f, "{}", " ".repeat(depth * 4))?;
        match self {
            GeneAssociation::Gene(gene_id) => match info.map(|m| m.get(gene_id.id)) {
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

    /// Removes instances where the same gene is OR'd or AND'd together with itself.
    /// This can occur where there are multiple isoforms with the same gene symbol.
    pub fn remove_isoform_summing(&self, genes: &Vec<Gene>) -> Self {
        match self {
            GeneAssociation::Gene(g) => GeneAssociation::Gene(g.clone()),
            GeneAssociation::And(vec) | GeneAssociation::Or(vec) => {
                let mut new = Vec::new();
                let mut seen = BTreeSet::new();
                for expr in vec {
                    match expr {
                        v @ GeneAssociation::Gene(gene_id) => {
                            if seen.insert(genes[gene_id.id].name()) {
                                new.push(v.clone());
                            }
                        }
                        v => {
                            new.push(v.remove_isoform_summing(genes));
                        }
                    }
                }
                match self {
                    GeneAssociation::And(..) => GeneAssociation::And(new),
                    GeneAssociation::Or(..) => GeneAssociation::Or(new),
                    _ => unreachable!(),
                }
            }
        }
    }
}

impl<'a> Display for GeneAssociationWithInfo<'a> {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        self.association
            .display_with_gene_info(f, 0, Some(self.info))
    }
}

impl Model {
    pub fn species(&self) -> Species {
        self.config.species
    }

    pub fn genes(&self) -> &[Gene] {
        &self.genes
    }

    pub fn reactions(&self) -> &[Reaction] {
        &self.reactions
    }

    pub fn metabolites(&self) -> &[Metabolite] {
        &self.metabolites
    }

    pub fn reaction_display<'a>(&'a self, reaction: &'a Reaction) -> String {
        use std::fmt::Write;
        // TODO: Returning a string here is less than ideal.
        // Also, this does not check whether the reaction is in the model.

        let mut s = String::new();
        writeln!(s, "{}", reaction.name).unwrap();
        let subsystem_name = self.subsystems[reaction.subsystem.index].name.clone();
        writeln!(s, "Subsystem: {}", subsystem_name).unwrap();
        writeln!(s, "Bounds: [{}, {}]", reaction.lb, reaction.ub).unwrap();
        if let Some(rule) = &reaction.rule {
            let rule = GeneAssociationWithInfo {
                association: rule,
                info: &self.genes,
            };
            writeln!(s, "Rule:\n{}", rule).unwrap();
        } else {
            writeln!(s, "Rule: None").unwrap();
        }
        s
    }
}

impl Reaction {
    pub fn name(&self) -> &str {
        &self.name
    }

    pub fn rule(&self) -> Option<&GeneAssociation> {
        self.rule.as_ref()
    }

    pub fn subsystem(&self) -> SubsystemId {
        self.subsystem
    }

    pub fn lower_bound(&self) -> f64 {
        self.lb
    }

    pub fn upper_bound(&self) -> f64 {
        self.ub
    }
}

impl Gene {
    /// Lists all symbols for this gene, including the primary symbol.
    pub fn get_symbols(&self) -> &[String] {
        &self.symbols
    }

    /// Returns the primary symbol of the gene.
    pub fn name(&self) -> &str {
        &self.symbols[0]
    }
}
