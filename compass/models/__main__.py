from . import load_metabolic_model, init_model

# Mucking about with this code while on a machine without cplex.
if __name__ == '__main__':
    print("hello world")
    model = load_metabolic_model("RECON1_mat")
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
    print(model.reactions["2OXOADOXm"].gene_associations)
    # Big ones, probably get deeply nested
    #"ATPasel"
    #"ATPS4m"




