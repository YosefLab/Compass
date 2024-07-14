import os
import json
import pandas as pd
import csv
from functools import cmp_to_key
import shutil

from compass.globals import EXCHANGE_LIMIT
from compass.models import load_metabolic_model, init_model

# NOTE: this script generates models for user-defined meta-subsystems based on the RECON2 format
# For each defined meta-subsystem, exchange reactions are added for all (non-currency) metabolites associated
# with the current meta-subsystem. For each metabolite, associated reactions that are not part of the current subsystem
# are also included

PATH_2_RECON_2_MAT = '/home/eecs/charleschien101/Compass/compass/Resources/Metabolic Models/RECON2_mat'
model_dir = os.path.join(PATH_2_RECON_2_MAT, 'model')

def smat_cmp(a, b):
    if a[1] > b[1]:
        return 1
    elif a[1] == b[1]:
        if a[0] > b[0]:
            return 1
        else:
            return -1
    else:
        return -1
    
def custom_json_dump(data):
    json_str = '[\n'
    for sublist in data:
        json_str += '    ' + json.dumps(sublist) + ',\n'
    json_str = json_str.rstrip(',\n') + '\n]'
    return json_str

def partition_model(args):

    # Parse user input
    with open(args['select_meta_subsystems']) as f:
        text = f.readlines()
    text = [line.strip() for line in text]
    entries = []
    for line in text:
        line = line.split(':')
        newline = [line[0]] + line[1].split(';')
        newline = [item.strip() for item in newline]
        entries.append(newline)

    selected_meta_subsystems = []
    selected_meta_subsystem_names = []
    meta_subsystems = {}

    for entry in entries:
        ms = entry[0]
        meta_subsystem_rxn = ms[ms.find("(") + 1 : ms.find(")")]
        meta_subsystem_name = ms[:ms.find("(") - 1]
        selected_meta_subsystems.append(meta_subsystem_rxn)
        selected_meta_subsystem_names.append(meta_subsystem_name)
        meta_subsystems[meta_subsystem_rxn] = entry[1:]

    # Create directories for meta subsystem models
    meta_subsystem_models_dir = os.path.join(args['output_dir'], 'meta_subsystem_models')
    if os.path.exists(meta_subsystem_models_dir) == False:
        os.mkdir(meta_subsystem_models_dir)

    for meta_subsystem in selected_meta_subsystems:
        output_dir = os.path.join(meta_subsystem_models_dir, f'{meta_subsystem}_mat')
        if os.path.exists(output_dir) == False:
            os.mkdir(output_dir)
        if os.path.exists(os.path.join(output_dir, 'model')) == False:
            os.mkdir(os.path.join(output_dir, 'model'))
        if os.path.exists(os.path.join(output_dir, 'media')) == False:
            os.mkdir(os.path.join(output_dir, 'media'))

    media = args['media']

    # Load reaction files
    with open(os.path.join(model_dir, f'model.lb.json')) as fin:
        recon2_lbs = json.load(fin)
    with open(os.path.join(model_dir, f'model.rev.json')) as fin:
        recon2_revs = json.load(fin)
    with open(os.path.join(model_dir, f'model.rules.json')) as fin:
        recon2_rules = json.load(fin)
    with open(os.path.join(model_dir, f'model.rxnECNumbers.json')) as fin:
        recon2_rxnECNumbers = json.load(fin)
    with open(os.path.join(model_dir, f'model.rxnKeggID.json')) as fin:
        recon2_rxnKeggIDs = json.load(fin)
    with open(os.path.join(model_dir, f'model.rxnNames.json')) as fin:
        recon2_rxnNames = json.load(fin)
    with open(os.path.join(model_dir, f'model.rxns.json')) as fin:
        recon2_rxns = json.load(fin)
    with open(os.path.join(model_dir, f'model.subSystems.json')) as fin:
        recon2_subSystems = json.load(fin)
    with open(os.path.join(model_dir, f'model.ub.json')) as fin:
        recon2_ubs = json.load(fin)

    # Load metabolite files
    with open(os.path.join(model_dir, f'model.metFormulas.json')) as fin:
        recon2_metFormulas = json.load(fin)
    with open(os.path.join(model_dir, f'model.metKeggID.json')) as fin:
        recon2_metKeggIDs = json.load(fin)
    with open(os.path.join(model_dir, f'model.metNames.json')) as fin:
        recon2_metNames = json.load(fin)
    with open(os.path.join(model_dir, f'model.mets.json')) as fin:
        recon2_mets = json.load(fin)

    recon2_model = load_metabolic_model(args['model'], species=args['species'])
    # Read in RECON2 SMAT
    recon2_smat = recon2_model.getSMAT()
    recon2_smat_transposed = recon2_model.getSMAT_transposed()

    # Read in core reactions
    with open('/home/eecs/charleschien101/Compass/compass/Resources/Metabolic Models/RECON2_mat/model/core_reactions.txt') as f:
        core_reactions = f.readlines()
    core_reactions = [r.strip() for r in core_reactions]

    model_dirs = []

    # Reaction-level information for each meta-subsystem
    meta_subsystem_lbs = {}
    meta_subsystem_ubs = {}
    meta_subsystem_revs = {}
    meta_subsystem_rules = {}
    meta_subsystem_rxnECNumbers = {}
    meta_subsystem_rxnKeggIDs = {}
    meta_subsystem_rxnNames = {}
    meta_subsystem_rxns = {}
    meta_subsystem_subSystems = {}
    meta_subsystem_smat = {}

    # Metabolite-level information for each meta-subsystem
    meta_subsystem_metFormulas = {}
    meta_subsystem_metKeggIDs = {}
    meta_subsystem_metNames = {}
    meta_subsystem_mets = {}

    # Metadata
    meta_subsystem_rxnMetas = {}
    meta_subsystem_metMetas = {}

    # Generate model for each meta subsystem
    for i in range(len(selected_meta_subsystems)):

        # Select current meta subsystem
        cur_meta_subsystem = selected_meta_subsystems[i]

        # Specify output directory for current meta subsystem
        output_dir = os.path.join(meta_subsystem_models_dir, f'{cur_meta_subsystem}_mat')
        model_dirs.append(output_dir)

        # Compute rxns
        cur_meta_subsystem_rxns = []
        for rxn, subSystem in zip(recon2_rxns, recon2_subSystems):
            if rxn not in core_reactions:
                continue

            if subSystem in meta_subsystems[cur_meta_subsystem]:
                cur_meta_subsystem_rxns.append(rxn)
        cur_meta_subsystem_rxns.sort()

        # Compute remaining reaction information
        cur_meta_subsystem_lbs = []
        cur_meta_subsystem_revs = []
        cur_meta_subsystem_rules = []
        cur_meta_subsystem_rxnECNumbers = []
        cur_meta_subsystem_rxnKeggIDs = []
        cur_meta_subsystem_rxnNames = []
        cur_meta_subsystem_subSystems = []
        cur_meta_subsystem_ubs = []

        for rxn in cur_meta_subsystem_rxns:
            idx = recon2_rxns.index(rxn)
            cur_meta_subsystem_lbs.append(recon2_lbs[idx])
            cur_meta_subsystem_revs.append(recon2_revs[idx])
            cur_meta_subsystem_rules.append(recon2_rules[idx])
            cur_meta_subsystem_rxnECNumbers.append(recon2_rxnECNumbers[idx])
            cur_meta_subsystem_rxnKeggIDs.append(recon2_rxnKeggIDs[idx])
            cur_meta_subsystem_rxnNames.append(recon2_rxnNames[idx])
            cur_meta_subsystem_subSystems.append(recon2_subSystems[idx])
            cur_meta_subsystem_ubs.append(recon2_ubs[idx])

        # Compute metabolites that are associated with meta-subsystem
        cur_meta_subsystem_mets = []
        for rxn in cur_meta_subsystem_rxns:
            smat_row = recon2_smat_transposed[rxn]
            for met, coef in smat_row:
                assert met[-3:] == '[m]' or met[-3:] == '[c]' or met[-3:] == '[e]'
                cur_meta_subsystem_mets.append(met)
        cur_meta_subsystem_mets = list(set(cur_meta_subsystem_mets))
        cur_meta_subsystem_mets.sort()

        # Compute metFormulas, metKeggIDs, and metNames
        cur_meta_subsystem_metFormulas = []
        cur_meta_subsystem_metKeggIDs = []
        cur_meta_subsystem_metNames = []
        for met in cur_meta_subsystem_mets:
            idx = recon2_mets.index(met)
            cur_meta_subsystem_metFormulas.append(recon2_metFormulas[idx])
            cur_meta_subsystem_metKeggIDs.append(recon2_metKeggIDs[idx])
            cur_meta_subsystem_metNames.append(recon2_metNames[idx])

        # Compute SMAT
        # SMAT indices correspond to reactions and metabolites of current meta-subsystem
        cur_meta_subsystem_smat = []
        for rxn in cur_meta_subsystem_rxns:
            smat_row = recon2_smat_transposed[rxn]
            rxn_idx = cur_meta_subsystem_rxns.index(rxn) + 1
            for met, coef in smat_row:
                met_idx = cur_meta_subsystem_mets.index(met) + 1
                cur_meta_subsystem_smat.append([met_idx, rxn_idx, coef])

        smat_cmp_key = cmp_to_key(smat_cmp)
        cur_meta_subsystem_smat.sort(key=smat_cmp_key)

        meta_subsystem_lbs[cur_meta_subsystem] = cur_meta_subsystem_lbs
        meta_subsystem_ubs[cur_meta_subsystem] = cur_meta_subsystem_ubs
        meta_subsystem_revs[cur_meta_subsystem] = cur_meta_subsystem_revs
        meta_subsystem_rules[cur_meta_subsystem] = cur_meta_subsystem_rules
        meta_subsystem_rxnECNumbers[cur_meta_subsystem] = cur_meta_subsystem_rxnECNumbers
        meta_subsystem_rxnKeggIDs[cur_meta_subsystem] = cur_meta_subsystem_rxnKeggIDs
        meta_subsystem_rxnNames[cur_meta_subsystem] = cur_meta_subsystem_rxnNames
        meta_subsystem_rxns[cur_meta_subsystem] = cur_meta_subsystem_rxns
        meta_subsystem_subSystems[cur_meta_subsystem] = cur_meta_subsystem_subSystems
        meta_subsystem_smat[cur_meta_subsystem] = cur_meta_subsystem_smat

        meta_subsystem_metFormulas[cur_meta_subsystem] = cur_meta_subsystem_metFormulas
        meta_subsystem_metKeggIDs[cur_meta_subsystem] = cur_meta_subsystem_metKeggIDs
        meta_subsystem_metNames[cur_meta_subsystem] = cur_meta_subsystem_metNames
        meta_subsystem_mets[cur_meta_subsystem] = cur_meta_subsystem_mets

        rxnMeta = pd.read_csv(os.path.join(model_dir, 'rxnMeta.txt'), delimiter='\t', quoting=csv.QUOTE_NONE)
        cur_meta_subsystem_rxnMeta = rxnMeta[rxnMeta['id'].isin(cur_meta_subsystem_rxns)].reset_index(drop=True).sort_values('id')
        meta_subsystem_rxnMetas[cur_meta_subsystem] = cur_meta_subsystem_rxnMeta
        assert len(cur_meta_subsystem_rxnMeta) == len(cur_meta_subsystem_rxns)

        metMeta = pd.read_csv(os.path.join(model_dir, 'metMeta.txt'), delimiter='\t', quoting=csv.QUOTE_NONE)
        cur_meta_subsystem_metMeta = metMeta[metMeta['id'].isin(cur_meta_subsystem_mets)].reset_index(drop=True).sort_values('id')
        meta_subsystem_rxnMetas[cur_meta_subsystem] = cur_meta_subsystem_rxnMeta
        assert len(cur_meta_subsystem_metMeta) == len(cur_meta_subsystem_mets)

        # Temporarily save .json files
        # Intermediary outputs are used to construct model to determine uptake/secretion metabolites
        # Files will be overwritten at the end
        with open(os.path.join(output_dir, 'model', f'model.lb.json'), 'w') as file:
            json.dump(cur_meta_subsystem_lbs, file, indent=4, separators=(',', ':'))

        with open(os.path.join(output_dir, 'model', f'model.rev.json'), 'w') as file:
            json.dump(cur_meta_subsystem_revs, file, indent=4, separators=(',', ':'))

        with open(os.path.join(output_dir, 'model', f'model.rules.json'), 'w') as file:
            json.dump(cur_meta_subsystem_rules, file, indent=4, separators=(',', ':'))

        with open(os.path.join(output_dir, 'model', f'model.rxnECNumbers.json'), 'w') as file:
            json.dump(cur_meta_subsystem_rxnECNumbers, file, indent=4, separators=(',', ':'))
            
        with open(os.path.join(output_dir, 'model', f'model.rxnKeggID.json'), 'w') as file:
            json.dump(cur_meta_subsystem_rxnKeggIDs, file, indent=4, separators=(',', ':'))

        with open(os.path.join(output_dir, 'model', f'model.rxnNames.json'), 'w') as file:
            json.dump(cur_meta_subsystem_rxnNames, file, indent=4, separators=(',', ':'))

        with open(os.path.join(output_dir, 'model', f'model.rxns.json'), 'w') as file:
            json.dump(cur_meta_subsystem_rxns, file, indent=4, separators=(',', ':'))

        with open(os.path.join(output_dir, 'model', f'model.subSystems.json'), 'w') as file:
            json.dump(cur_meta_subsystem_subSystems, file, indent=4, separators=(',', ':'))

        with open(os.path.join(output_dir, 'model', f'model.ub.json'), 'w') as file:
            json.dump(cur_meta_subsystem_ubs, file, indent=4, separators=(',', ':'))

        with open(os.path.join(output_dir, 'model', f'model.metFormulas.json'), 'w') as file:
            json.dump(cur_meta_subsystem_metFormulas, file, indent=4, separators=(',', ':'))

        with open(os.path.join(output_dir, 'model', f'model.metKeggID.json'), 'w') as file:
            json.dump(cur_meta_subsystem_metKeggIDs, file, indent=4, separators=(',', ':'))

        with open(os.path.join(output_dir, 'model', f'model.metNames.json'), 'w') as file:
            json.dump(cur_meta_subsystem_metNames, file, indent=4, separators=(',', ':'))

        with open(os.path.join(output_dir, 'model', f'model.mets.json'), 'w') as file:
            json.dump(cur_meta_subsystem_mets, file, indent=4, separators=(',', ':'))

        formatted_json_str = custom_json_dump(cur_meta_subsystem_smat)
        with open(os.path.join(output_dir, 'model', f'model.S.json'), 'w') as file:
            file.write(formatted_json_str)

        cur_meta_subsystem_rxnMeta.to_csv(os.path.join(output_dir, 'model', 'rxnMeta.txt'), sep='\t', index=False)
        cur_meta_subsystem_metMeta.to_csv(os.path.join(output_dir, 'model', 'metMeta.txt'), sep='\t', index=False)

        shutil.copy(os.path.join(model_dir, 'model.genes.json'), os.path.join(output_dir, 'model', 'model.genes.json'))
        shutil.copy(os.path.join(PATH_2_RECON_2_MAT, 'media', f'{media}.json'), os.path.join(output_dir, 'media', f'{media}.json'))

        gene_filenames = ['mouseSynonyms', 'non2uniqueEntrez', 'uniqueEntrez2Non', 'uniqueHumanEntrez', 'uniqueHumanGeneSymbol', 'uniqueMouseGeneSymbol_all', 'uniqueMouseGeneSymbol']
        for filename in gene_filenames:
            shutil.copy(os.path.join(PATH_2_RECON_2_MAT, f'{filename}.json'), os.path.join(output_dir, f'{filename}.json'))

    # **********************************************************************

    # Add exchange reaction for all metabolites associated with subsystem

    for meta_subsystem in selected_meta_subsystems:
        meta_subsystem_model = init_model(f'{meta_subsystem}_mat', species=args['species'], exchange_limit=EXCHANGE_LIMIT, media=args['media'],
                                        metabolic_model_dir=meta_subsystem_models_dir)
        
        meta_subsystem_model_smat = meta_subsystem_model.getSMAT()
        
        rxn_idx = len(meta_subsystem_rxns[meta_subsystem]) + 1

        for met in meta_subsystem_mets[meta_subsystem]:

            # First check if there is already an exchange reaction associated with the current metabolite
            has_exchange = False
            smat_row = meta_subsystem_model_smat[met]

            # Check all associated reactions
            for rxn_id, coef in smat_row:
                rxn = meta_subsystem_model.reactions[rxn_id]
                if rxn.is_exchange:
                    has_exchange = True
                    print(f'{met} has exchange reaction {rxn_id}')
                    break

            if has_exchange == True:
                continue

            # If metabolite does not have associated exchange reaction, add uptake/secretion reaction pair
            exchange_lb = -1000
            exchange_ub = 1000
            exchange_rev = 0
            exchange_rxnECNumber = ""
            exchange_rxnKeggID = ""
            exchange_rxnName = f'{met}_EXCHANGE_{meta_subsystem}'
            exchange_rxn = f'{met}_EXCHANGE_{meta_subsystem}'
            exchange_subSystem = meta_subsystem
            exchange_rule = ""

            # smat
            met_idx = meta_subsystem_mets[meta_subsystem].index(met) + 1
            exchange_smat = [met_idx, rxn_idx, -1.0]
            rxn_idx += 1

            meta_subsystem_lbs[meta_subsystem].append(exchange_lb)
            meta_subsystem_revs[meta_subsystem].append(exchange_rev)
            meta_subsystem_rules[meta_subsystem].append(exchange_rule)
            meta_subsystem_rxnECNumbers[meta_subsystem].append(exchange_rxnECNumber)
            meta_subsystem_rxnKeggIDs[meta_subsystem].append(exchange_rxnKeggID)
            meta_subsystem_rxnNames[meta_subsystem].append(exchange_rxnName)
            meta_subsystem_rxns[meta_subsystem].append(exchange_rxn)
            meta_subsystem_subSystems[meta_subsystem].append(exchange_subSystem)
            meta_subsystem_ubs[meta_subsystem].append(exchange_ub)
            meta_subsystem_smat[meta_subsystem].append(exchange_smat)

            idx = len(meta_subsystem_rxnMetas[meta_subsystem])
            meta_subsystem_rxnMetas[meta_subsystem].loc[idx] = [
                meta_subsystem_rxns[meta_subsystem][idx],
                meta_subsystem_subSystems[meta_subsystem][idx],
                meta_subsystem_rxnNames[meta_subsystem][idx],
                meta_subsystem_rxnKeggIDs[meta_subsystem][idx],
                '',
                0.0,
                0.0,
                'NA',
                meta_subsystem_rxnECNumbers[meta_subsystem][idx],
                meta_subsystem_rxnNames[meta_subsystem][idx]
            ]

    # **********************************************************************

    # For each non-currency metabolite in the meta-subsystem, add all associated reactions
    # Also add exchange reactions for newly added metabolites

    with open('/home/eecs/charleschien101/Compass/compass/Resources/Metabolic Models/currency_mets.txt') as file:
        currency_mets = file.readlines()
    currency_mets = [m.strip() for m in currency_mets]

    for meta_subsystem in selected_meta_subsystems:

        output_dir = os.path.join(meta_subsystem_models_dir, f'{meta_subsystem}_mat')

        # Reactions and metabolites that are added to the meta-subsystem after adding all associated reactions
        new_rxns = set()
        new_mets = set()

        for met in meta_subsystem_mets[meta_subsystem]:

            if met in currency_mets:
                continue

            associated_smat_row = recon2_smat[met]

            for rxn, coef in associated_smat_row:
                if rxn not in core_reactions:
                    continue
                if rxn in meta_subsystem_rxns[meta_subsystem]:
                    continue

                new_rxns.add(rxn)
            
                associated_smat_transpose_row = recon2_smat_transposed[rxn]
                for met, coef in associated_smat_transpose_row:
                    assert met[-3:] == '[m]' or met[-3:] == '[c]' or met[-3:] == '[e]'
                    if met in meta_subsystem_mets[meta_subsystem]:
                        continue
                    new_mets.add(met)

        new_rxns = list(new_rxns)
        new_rxns.sort()
        new_mets = list(new_mets)
        new_mets.sort()

        # Add new metabolites to meta-subsystem
        for new_met in new_mets:
            meta_subsystem_mets[meta_subsystem].append(new_met)

            idx = recon2_mets.index(new_met)
            meta_subsystem_metFormulas[meta_subsystem].append(recon2_metFormulas[idx])
            meta_subsystem_metKeggIDs[meta_subsystem].append(recon2_metKeggIDs[idx])
            meta_subsystem_metNames[meta_subsystem].append(recon2_metNames[idx])

        # Add associated reactions to meta-subsystem
        for new_rxn in new_rxns:
            meta_subsystem_rxns[meta_subsystem].append(new_rxn)

            idx = recon2_rxns.index(new_rxn)
            meta_subsystem_lbs[meta_subsystem].append(recon2_lbs[idx])
            meta_subsystem_revs[meta_subsystem].append(recon2_revs[idx])
            meta_subsystem_rules[meta_subsystem].append(recon2_rules[idx])
            meta_subsystem_rxnECNumbers[meta_subsystem].append(recon2_rxnECNumbers[idx])
            meta_subsystem_rxnKeggIDs[meta_subsystem].append(recon2_rxnKeggIDs[idx])
            meta_subsystem_rxnNames[meta_subsystem].append(recon2_rxnNames[idx])
            meta_subsystem_subSystems[meta_subsystem].append(recon2_subSystems[idx])
            meta_subsystem_ubs[meta_subsystem].append(recon2_ubs[idx])

            # Compute SMAT
            # SMAT indices correspond to reactions and metabolites of updated meta-subsystem
            smat_row = recon2_smat_transposed[new_rxn]
            rxn_idx = len(meta_subsystem_rxns[meta_subsystem])
            assert rxn_idx == meta_subsystem_rxns[meta_subsystem].index(new_rxn) + 1
            for met, coef in smat_row:
                met_idx = meta_subsystem_mets[meta_subsystem].index(met) + 1
                meta_subsystem_smat[meta_subsystem].append([met_idx, rxn_idx, coef])

            smat_cmp_key = cmp_to_key(smat_cmp)
            meta_subsystem_smat[meta_subsystem].sort(key=smat_cmp_key)

            meta_idx = len(meta_subsystem_rxnMetas[meta_subsystem])
            meta_subsystem_rxnMetas[meta_subsystem].loc[meta_idx] = [
                meta_subsystem_rxns[meta_subsystem][meta_idx],
                meta_subsystem_subSystems[meta_subsystem][meta_idx],
                meta_subsystem_rxnNames[meta_subsystem][meta_idx],
                meta_subsystem_rxnKeggIDs[meta_subsystem][meta_idx],
                '',
                0.0,
                0.0,
                'NA',
                meta_subsystem_rxnECNumbers[meta_subsystem][meta_idx],
                meta_subsystem_rxnNames[meta_subsystem][meta_idx]
            ]

        # Add exchange reactions for new metabolites
        # No need to consider existing exchange reactions since it is impossible for
        # newly added metabolites to be associated with exchange reactions
        rxn_idx = len(meta_subsystem_rxns[meta_subsystem]) + 1

        for new_met in new_mets:

            exchange_lb = -1000
            exchange_ub = 1000
            exchange_rev = 0
            exchange_rxnECNumber = ""
            exchange_rxnKeggID = ""
            exchange_rxnName = f'{new_met}_EXCHANGE_{meta_subsystem}'
            exchange_rxn = f'{new_met}_EXCHANGE_{meta_subsystem}'
            exchange_subSystem = meta_subsystem
            exchange_rule = ""

            # smat
            met_idx = meta_subsystem_mets[meta_subsystem].index(new_met) + 1
            exchange_smat = [met_idx, rxn_idx, -1.0]
            rxn_idx += 1

            meta_subsystem_lbs[meta_subsystem].append(exchange_lb)
            meta_subsystem_revs[meta_subsystem].append(exchange_rev)
            meta_subsystem_rules[meta_subsystem].append(exchange_rule)
            meta_subsystem_rxnECNumbers[meta_subsystem].append(exchange_rxnECNumber)
            meta_subsystem_rxnKeggIDs[meta_subsystem].append(exchange_rxnKeggID)
            meta_subsystem_rxnNames[meta_subsystem].append(exchange_rxnName)
            meta_subsystem_rxns[meta_subsystem].append(exchange_rxn)
            meta_subsystem_subSystems[meta_subsystem].append(exchange_subSystem)
            meta_subsystem_ubs[meta_subsystem].append(exchange_ub)
            meta_subsystem_smat[meta_subsystem].append(exchange_smat)

            idx = len(meta_subsystem_rxnMetas[meta_subsystem])
            meta_subsystem_rxnMetas[meta_subsystem].loc[idx] = [
                meta_subsystem_rxns[meta_subsystem][idx],
                meta_subsystem_subSystems[meta_subsystem][idx],
                meta_subsystem_rxnNames[meta_subsystem][idx],
                meta_subsystem_rxnKeggIDs[meta_subsystem][idx],
                '',
                0.0,
                0.0,
                'NA',
                meta_subsystem_rxnECNumbers[meta_subsystem][idx],
                meta_subsystem_rxnNames[meta_subsystem][idx]
            ]

        assert len(meta_subsystem_rxnMetas[meta_subsystem]) == len(meta_subsystem_rxns[meta_subsystem])

        metMeta = pd.read_csv(os.path.join(model_dir, 'metMeta.txt'), delimiter='\t', quoting=csv.QUOTE_NONE)
        meta_subsystem_metMetas[meta_subsystem] = metMeta[metMeta['id'].isin(meta_subsystem_mets[meta_subsystem])].reset_index(drop=True)
        meta_subsystem_metMetas[meta_subsystem].sort_values(by='id', key=lambda column: column.map(lambda e: meta_subsystem_mets[meta_subsystem].index(e)), inplace=True)
        assert len(meta_subsystem_metMetas[meta_subsystem]) == len(meta_subsystem_mets[meta_subsystem])

        # Save final .json files
        with open(os.path.join(output_dir, 'model', f'model.lb.json'), 'w') as file:
            json.dump(meta_subsystem_lbs[meta_subsystem], file, indent=4, separators=(',', ':'))

        with open(os.path.join(output_dir, 'model', f'model.rev.json'), 'w') as file:
            json.dump(meta_subsystem_revs[meta_subsystem], file, indent=4, separators=(',', ':'))

        with open(os.path.join(output_dir, 'model', f'model.rules.json'), 'w') as file:
            json.dump(meta_subsystem_rules[meta_subsystem], file, indent=4, separators=(',', ':'))

        with open(os.path.join(output_dir, 'model', f'model.rxnECNumbers.json'), 'w') as file:
            json.dump(meta_subsystem_rxnECNumbers[meta_subsystem], file, indent=4, separators=(',', ':'))
            
        with open(os.path.join(output_dir, 'model', f'model.rxnKeggID.json'), 'w') as file:
            json.dump(meta_subsystem_rxnKeggIDs[meta_subsystem], file, indent=4, separators=(',', ':'))

        with open(os.path.join(output_dir, 'model', f'model.rxnNames.json'), 'w') as file:
            json.dump(meta_subsystem_rxnNames[meta_subsystem], file, indent=4, separators=(',', ':'))

        with open(os.path.join(output_dir, 'model', f'model.rxns.json'), 'w') as file:
            json.dump(meta_subsystem_rxns[meta_subsystem], file, indent=4, separators=(',', ':'))

        with open(os.path.join(output_dir, 'model', f'model.subSystems.json'), 'w') as file:
            json.dump(meta_subsystem_subSystems[meta_subsystem], file, indent=4, separators=(',', ':'))

        with open(os.path.join(output_dir, 'model', f'model.ub.json'), 'w') as file:
            json.dump(meta_subsystem_ubs[meta_subsystem], file, indent=4, separators=(',', ':'))

        with open(os.path.join(output_dir, 'model', f'model.metFormulas.json'), 'w') as file:
            json.dump(meta_subsystem_metFormulas[meta_subsystem], file, indent=4, separators=(',', ':'))

        with open(os.path.join(output_dir, 'model', f'model.metKeggID.json'), 'w') as file:
            json.dump(meta_subsystem_metKeggIDs[meta_subsystem], file, indent=4, separators=(',', ':'))

        with open(os.path.join(output_dir, 'model', f'model.metNames.json'), 'w') as file:
            json.dump(meta_subsystem_metNames[meta_subsystem], file, indent=4, separators=(',', ':'))

        with open(os.path.join(output_dir, 'model', f'model.mets.json'), 'w') as file:
            json.dump(meta_subsystem_mets[meta_subsystem], file, indent=4, separators=(',', ':'))

        formatted_json_str = custom_json_dump(meta_subsystem_smat[meta_subsystem])
        with open(os.path.join(output_dir, 'model', f'model.S.json'), 'w') as file:
            file.write(formatted_json_str)

        meta_subsystem_rxnMetas[meta_subsystem].to_csv(os.path.join(output_dir, 'model', 'rxnMeta.txt'), sep='\t', index=False)
        meta_subsystem_metMetas[meta_subsystem].to_csv(os.path.join(output_dir, 'model', 'metMeta.txt'), sep='\t', index=False)

        shutil.copy(os.path.join(model_dir, 'model.genes.json'), os.path.join(output_dir, 'model', 'model.genes.json'))
        shutil.copy(os.path.join(PATH_2_RECON_2_MAT, 'media', f'{media}.json'), os.path.join(output_dir, 'media', f'{media}.json'))

        gene_filenames = ['mouseSynonyms', 'non2uniqueEntrez', 'uniqueEntrez2Non', 'uniqueHumanEntrez', 'uniqueHumanGeneSymbol', 'uniqueMouseGeneSymbol_all', 'uniqueMouseGeneSymbol']
        for filename in gene_filenames:
            shutil.copy(os.path.join(PATH_2_RECON_2_MAT, f'{filename}.json'), os.path.join(output_dir, f'{filename}.json'))

    # **********************************************************************

    model_names = [f'{meta_subsystem}_mat' for meta_subsystem in selected_meta_subsystems]

    return meta_subsystem_models_dir, model_names