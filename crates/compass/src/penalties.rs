use gsmm::model::GeneAssociation;

pub struct GeneRuleEval<'a, OR, AND> {
    pub gene_expr: &'a [f64],
    pub or_op: OR,
    pub and_op: AND,
}

impl<'a, OR: GeneRuleOp, AND: GeneRuleOp> GeneRuleEval<'_, OR, AND> {
    pub fn evaluate(&self, node: &GeneAssociation) -> f64 {
        match node {
            GeneAssociation::Gene(gene_id) => self.gene_expr[gene_id.index()],
            GeneAssociation::Or(vec) => self.or_op.apply(vec.iter().map(|x| self.evaluate(x))),
            GeneAssociation::And(vec) => self.and_op.apply(vec.iter().map(|x| self.evaluate(x))),
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
