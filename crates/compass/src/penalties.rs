use tracing::debug;

use gsmm::model::GeneAssociation;

pub struct GeneRuleEval<'a, OR, AND> {
    pub gene_expr: &'a [f64],
    pub or_op: OR,
    pub and_op: AND,
}

impl<'a, OR: GeneRuleOp, AND: GeneRuleOp> GeneRuleEval<'_, OR, AND> {
    pub fn evaluate(&self, node: &GeneAssociation) -> f64 {
        self.evaluate_debug(node, false)
    }

    pub fn evaluate_debug(&self, node: &GeneAssociation, debug: bool) -> f64 {
        match node {
            GeneAssociation::Gene(gene_id) => {
                let expr = self.gene_expr[gene_id.index()];
                if debug {
                    debug!("Gene {gene_id:?} has expr {expr}");
                }
                expr
            }
            GeneAssociation::Or(vec) => self
                .or_op
                .apply(vec.iter().map(|x| self.evaluate_debug(x, debug))),
            GeneAssociation::And(vec) => self
                .and_op
                .apply(vec.iter().map(|x| self.evaluate_debug(x, debug))),
        }
    }
}

pub trait GeneRuleOp {
    fn apply(&self, operands: impl Iterator<Item = f64>) -> f64;
}

pub struct GeneOrSum;
pub struct GeneAndMean;

impl GeneRuleOp for GeneOrSum {
    fn apply(&self, operands: impl Iterator<Item = f64>) -> f64 {
        operands.sum()
    }
}

impl GeneRuleOp for GeneAndMean {
    fn apply(&self, operands: impl Iterator<Item = f64>) -> f64 {
        let (sum, count) = operands.fold((0.0, 0usize), |(sum, count), x| (sum + x, count + 1));
        if count == 0 {
            0.0
        } else {
            sum / count as f64
        }
    }
}
