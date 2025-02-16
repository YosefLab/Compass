from . import load_metabolic_model, init_model
from .MetabolicModel import sum_wo_nan, min_w_nan
import pandas as pd

# Mucking about with this code while on a machine without cplex.
if __name__ == '__main__':
    print("hello world")
    species = "mus_musculus"
    model = load_metabolic_model("RECON2_mat", species)
    model.remove_isoform_summing()
    expression = pd.read_csv("~/repos/Compass/compass/Resources/Test-Data/tsv_format/expression.tsv", sep="\t", index_col=0)
    expression.index = expression.index.astype('str').str.upper()
    sample = "Ob-DHA-e_S154_L007_R1_001"
    sample_expression = expression[sample]
    print(sample_expression)
    #assert len(model.reactions) == 7440
    #assert len(model.species) == 5063
    
    rxns = { k:r for (k, r) in model.reactions.items() if r.gene_associations }
    #print(len(rxns))
    #print(next(iter(rxns)))
    #print(model.reactions["34DHOXPEGOX"].gene_associations)
    #print(model.reactions["ACCOAC"].gene_associations)
    #print(model.reactions["AKGDm"].gene_associations)
    #print(model.reactions["ALCD1"].gene_associations)
    # One with OR -> AND -> OR. Potentially a concern for associativity of the parsing?
    print(model.reactions["ENO"].gene_associations)
    print(model.reactions["ENO"].eval_expression(sample_expression, min_w_nan, sum_wo_nan))
    # Big ones, probably get deeply nested
    #"ATPasel"
    #"ATPS4m"




