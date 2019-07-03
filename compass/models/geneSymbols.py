"""
Helps to resolve gene symbols by using information from HGNC
"""
from __future__ import print_function, division, absolute_import


import random
import gzip
import json
import re
import os
from ..globals import RESOURCE_DIR
from .MetabolicModel import Gene, Association
from collections import defaultdict
import pandas as pd


this_directory = os.path.dirname(os.path.abspath(__file__))


def detect_type(model):
    """
    Takes the input MetabolicModel and detects the
    type of gene symbols in its genes' name field

    returns either 'HGNC', 'Entrez', 'Ensembl', or 'Symbol'
    """

    # accumulate all model genes
    all_genes = set()
    for reaction in model.reactions.values():
        all_genes |= set(reaction.list_genes())

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


def resolve_genes(model):
    """
    Ensures that the gene names in a model correspond with hgnc gene symbols
    and are not some other identifier.

    Also fills in the alt_symbols on genes if conversion occurs
    """

    field_type = detect_type(model)

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

    gene_dict = {x[key_field]: x for x in genes}

    # Now, recursively crawl the gene association structure
    #   and update all gene entries

    def update_association(assoc):
        if assoc.type == "gene":
            gene = assoc.gene
            if gene.name in gene_dict:
                gene_info = gene_dict[gene.name]
                gene.name = gene_info['symbol'].upper()
                if 'alias_symbol' in gene_info:
                    gene.alt_symbols = [
                        x.upper() for x in gene_info['alias_symbol']
                    ]
        else:
            for child in assoc.children:
                update_association(child)

    for reaction in model.reactions.values():
        if reaction.gene_associations is not None:
            update_association(reaction.gene_associations)


def load_mgi():
    """
    Loads the ortho2human, ortho2mouse dictionaries
    from the MGI exported file
    """

    ortho2mouse = defaultdict(set)
    ortho2human = defaultdict(set)

    mgi_file = os.path.join(
        RESOURCE_DIR, "Genes", "HOM_MouseHumanSequence.rpt.gz")

    data = pd.read_csv(mgi_file, sep="\t")

    hgs = data.loc[data['Common Organism Name'] == 'human']
    mgs = data.loc[data['Common Organism Name'] == 'mouse, laboratory']

    for ortho_id, hg in zip(hgs['HomoloGene ID'], hgs['Symbol']):
        ortho2human[ortho_id].add(hg.upper())

    for ortho_id, mg in zip(mgs['HomoloGene ID'], mgs['Symbol']):
        ortho2mouse[ortho_id].add(mg.upper())

    return ortho2human, ortho2mouse


def convert_species(model, target_species):

    if target_species == "homo_sapiens":
        return  # Nothing to do, assume input is human

    if target_species != "mus_musculus":
        raise NotImplementedError("Can only convert to mus_musculus")

    ortho2human, ortho2mouse = load_mgi()

    # Invert the human dictionary
    human2ortho = defaultdict(set)
    for ortho_id in ortho2human:
        for hg in ortho2human[ortho_id]:
            human2ortho[hg].add(ortho_id)

    # Now, recursively crawl the gene association structure
    #   and update all gene entries

    def update_association(assoc):
        if assoc.type == "gene":
            gene = assoc.gene
            if gene.name in human2ortho:
                ortho_id = human2ortho[gene.name]

                mouse_genes = set()

                for oid in ortho_id:
                    if oid in ortho2mouse:
                        mouse_genes |= ortho2mouse[oid]

                if len(mouse_genes) == 0:
                    gene.name = ""
                    gene.alt_symbols = []
                    return

                if len(mouse_genes) == 1:
                    gene.name = list(mouse_genes)[0]
                    gene.alt_symbols = []
                    return

                # About 30 human genes map to multiple mouse genes
                # Replace them with an OR association

                assoc.type = 'or'
                assoc.gene = None
                assoc.children = []
                for mg in mouse_genes:
                    new_gene = Gene()
                    new_gene.id = mg
                    new_gene.name = mg

                    new_assoc = Association()
                    new_assoc.type = 'gene'
                    new_assoc.gene = new_gene

                    assoc.children.append(new_assoc)

        else:
            for child in assoc.children:
                update_association(child)

    for reaction in model.reactions.values():
        if reaction.gene_associations is not None:
            update_association(reaction.gene_associations)
