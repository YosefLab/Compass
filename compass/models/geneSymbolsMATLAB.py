"""
Based on functions in geneSymbols.py.

Statically generates for MATLAB models:

non2uniqueEntrez.json
uniqueHumanGeneSymbol.json

"""

import random
import gzip
import json
import re
import os
from collections import OrderedDict

RESOURCE_DIR = os.path.join("..", "Resources")
MODEL_DIR = os.path.join(RESOURCE_DIR, "Metabolic Models")

def detect_type(gene_list):
    """
    Takes the input MetabolicModel and detects the
    type of gene symbols in its genes' name field

    returns either 'HGNC', 'Entrez', 'Ensembl', or 'Symbol'
    """

    # accumulate all model genes
    all_genes = set(gene_list)

    # Test N genes for their type
    N = 100
    test_genes = random.sample(all_genes, N)

    ens_re = re.compile("ENSG\d+")
    entrez_re = re.compile("\d+")
    hgnc_re = re.compile("HGNC:\d+")

    ens_count = 0
    entrez_count = 0
    hgnc_count = 0

    for gene in test_genes:
        if ens_re.match(gene):
            ens_count += 1
        if entrez_re.match(gene):
            entrez_count += 1
        if hgnc_re.match(gene):
            hgnc_count += 1

    if ens_count == N:
        return "Ensembl"
    elif entrez_count == N:
        return "Entrez"
    elif hgnc_count == N:
        return "HGNC"
    else:
        return "Symbol"


def resolve_genes(gene_list):
    """
    Ensures that the gene names in a model correspond with hgnc gene symbols
    and are not some other identifier.

    And generates the mappings of non2uniqueEntrez and uniqueHumanGeneSymbol
    """
    non2uniqueEntrez = []
    uniqueHumanGeneSymbol = []

    field_type = detect_type(gene_list)

    if field_type == "Symbol":  # No work to do
        return

    # Load gene symbol information

    hgnc_file = os.path.join(
        RESOURCE_DIR, "Genes", "HGNC.json.gz"
    )

    data_dict = json.loads(gzip.open(hgnc_file).read().decode('utf-8'))
    genes = data_dict['response']['docs']

    # Key dictionary based on the identifiers that we currently have
    if field_type == "HGNC":
        key_field = "hgnc_id"
    elif field_type == "Entrez":
        key_field = "entrez_id"
    elif field_type == "Ensembl":
        key_field = "ensembl_gene_id"
    else:
        raise Exception("Error: invalid field_type")

    gene_dict = {x.get(key_field): x for x in genes}

    for gid in gene_list:
        gene_info = gene_dict.get(gid)
        if gene_info:
            gene_name = gene_info['symbol'].upper()
            try:
                i = uniqueHumanGeneSymbol.index(gene_name)
            except ValueError:
                uniqueHumanGeneSymbol.append(gene_name)
                i = len(uniqueHumanGeneSymbol)                
        else:
            gene_name = ""
            uniqueHumanGeneSymbol.append(gene_name)
            i = len(uniqueHumanGeneSymbol)
        non2uniqueEntrez.append(i)
    
    return (non2uniqueEntrez, uniqueHumanGeneSymbol)

model_name = "RECON3_mat"
model_dir = os.path.join(MODEL_DIR, model_name)
geneFilePath = os.path.join(model_dir, "model", "model.genes.json")


with open(geneFilePath) as geneFile:
    gene_list = json.load(geneFile)

gene_list = [gene[:-2] for gene in gene_list]

non2uniqueEntrez, uniqueHumanGeneSymbol = resolve_genes(gene_list)

with open(os.path.join(model_dir, "non2uniqueEntrez.json"), "w") as f:
    json.dump(non2uniqueEntrez, f, indent=2)

with open(os.path.join(model_dir, "uniqueHumanGeneSymbol.json"), "w") as f:
    json.dump(uniqueHumanGeneSymbol, f, indent=2)