import os
import json
import pandas as pd
import csv
from functools import cmp_to_key
import shutil

from ..globals import EXCHANGE_LIMIT
from compass.models import load_metabolic_model, init_model

# NOTE: this script generates models for user-defined meta-subsystems based on the RECON2 format
# For each defined meta-subsystem, exchange reactions are added for all metabolites associated with the current meta-subsystem

PATH_2_RECON2_MAT = '/home/eecs/charleschien101/Compass/compass/Resources/Metabolic Models/RECON2_mat'
model_dir = os.path.join(PATH_2_RECON2_MAT, 'model')

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


    # All .json files
    all_json_filenames = ['genes', 'lb', 'metFormulas', 'metKeggID', 'metNames', 'mets', 'rev', 'rules', 'rxnECNumbers', 'rxnKeggID', 'rxnNames', 'rxns', 'subSystems', 'ub']
    # Files that need to be modified
    json_filenames = ['lb', 'rev', 'rules', 'rxnECNumbers', 'rxnKeggID', 'rxnNames', 'rxns', 'subSystems', 'ub']
    # Files that can be copied directly
    json_copy_filenames = ['genes', 'metFormulas', 'metKeggID', 'metNames', 'mets']

    subsystem_data = {}
    for meta_subsystem in selected_meta_subsystems:
        subsystem_data[meta_subsystem] = {}
        for subsystem in meta_subsystems[meta_subsystem]:
            subsystem_data[meta_subsystem][subsystem] = {}
            for filename in json_filenames:
                subsystem_data[meta_subsystem][subsystem][filename] = []
            subsystem_data[meta_subsystem][subsystem]['smat'] = {}

    with open(os.path.join(model_dir, f'model.lb.json')) as fin:
        lbs = json.load(fin)
    with open(os.path.join(model_dir, f'model.rev.json')) as fin:
        revs = json.load(fin)
    with open(os.path.join(model_dir, f'model.rules.json')) as fin:
        rules = json.load(fin)
    with open(os.path.join(model_dir, f'model.rxnECNumbers.json')) as fin:
        rxnECNumbers = json.load(fin)
    with open(os.path.join(model_dir, f'model.rxnKeggID.json')) as fin:
        rxnKeggIDs = json.load(fin)
    with open(os.path.join(model_dir, f'model.rxnNames.json')) as fin:
        rxnNames = json.load(fin)
    with open(os.path.join(model_dir, f'model.rxns.json')) as fin:
        rxns = json.load(fin)
    with open(os.path.join(model_dir, f'model.subSystems.json')) as fin:
        subsystems = json.load(fin)
    with open(os.path.join(model_dir, f'model.ub.json')) as fin:
        ubs = json.load(fin)

    with open(os.path.join(model_dir, 'core_reactions.txt')) as fin:
        core_rxns = fin.readlines()
    core_rxns = [rxn.strip() for rxn in core_rxns]

    # Read in data for core reactions
    for meta_subsystem in selected_meta_subsystems:
        for i, subsystem in enumerate(subsystems):
            '''if rxns[i] not in core_rxns:
                continue'''
            if subsystem in meta_subsystems[meta_subsystem]:
                subsystem_data[meta_subsystem][subsystem]['lb'].append(lbs[i])
                subsystem_data[meta_subsystem][subsystem]['rev'].append(revs[i])
                subsystem_data[meta_subsystem][subsystem]['rules'].append(rules[i])
                subsystem_data[meta_subsystem][subsystem]['rxnECNumbers'].append(rxnECNumbers[i])
                subsystem_data[meta_subsystem][subsystem]['rxnKeggID'].append(rxnKeggIDs[i])
                subsystem_data[meta_subsystem][subsystem]['rxnNames'].append(rxnNames[i])
                subsystem_data[meta_subsystem][subsystem]['rxns'].append(rxns[i])
                subsystem_data[meta_subsystem][subsystem]['subSystems'].append(subsystems[i])
                subsystem_data[meta_subsystem][subsystem]['ub'].append(ubs[i])

    model = load_metabolic_model(args['model'], species=args['species'])

    # Read in RECON2 SMAT
    smat_transposed = model.getSMAT_transposed()

    for meta_subsystem in selected_meta_subsystems:
        for subsystem in meta_subsystems[meta_subsystem]:
            for rxn in subsystem_data[meta_subsystem][subsystem]['rxns']:
                subsystem_data[meta_subsystem][subsystem]['smat'][rxn] = smat_transposed[rxn]

    model_dirs = []

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
    meta_subsystem_mets = {}
    meta_subsystem_metas = {}

    # Generate model for each meta subsystem
    for i in range(len(selected_meta_subsystems)):

        # Select current meta subsystem
        cur_meta_subsystem = selected_meta_subsystems[i]
        #remaining_meta_subsystems = selected_meta_subsystems.copy()
        #remaining_meta_subsystems.pop(remaining_meta_subsystems.index(cur_meta_subsystem))

        cur_meta_subsystem_name = selected_meta_subsystem_names[i]
        #remaining_meta_subsystem_names = selected_meta_subsystem_names.copy()
        #remaining_meta_subsystem_names.pop(remaining_meta_subsystem_names.index(cur_meta_subsystem_name))

        # Specify output directory for current meta subsystem
        output_dir = os.path.join(meta_subsystem_models_dir, f'{cur_meta_subsystem}_mat')
        model_dirs.append(output_dir)

        # Compute lower bounds
        cur_meta_subsystem_lbs = []
        for subsystem in meta_subsystems[cur_meta_subsystem]:
            cur_meta_subsystem_lbs += subsystem_data[cur_meta_subsystem][subsystem]['lb']

        # Compute upper bounds
        cur_meta_subsystem_ubs = []
        for subsystem in meta_subsystems[cur_meta_subsystem]:
            cur_meta_subsystem_ubs += subsystem_data[cur_meta_subsystem][subsystem]['ub']

        # Compute rev
        cur_meta_subsystem_revs = []
        for subsystem in meta_subsystems[cur_meta_subsystem]:
            cur_meta_subsystem_revs += subsystem_data[cur_meta_subsystem][subsystem]['rev']

        # Compute ECNumber
        cur_meta_subsystem_rxnECNumbers = []
        for subsystem in meta_subsystems[cur_meta_subsystem]:
            cur_meta_subsystem_rxnECNumbers += subsystem_data[cur_meta_subsystem][subsystem]['rxnECNumbers']

        # Compute KeggIDs
        cur_meta_subsystem_rxnKeggIDs = []
        for subsystem in meta_subsystems[cur_meta_subsystem]:
            cur_meta_subsystem_rxnKeggIDs += subsystem_data[cur_meta_subsystem][subsystem]['rxnKeggID']

        # Compute rxn names
        cur_meta_subsystem_rxnNames = []
        for subsystem in meta_subsystems[cur_meta_subsystem]:
            cur_meta_subsystem_rxnNames += subsystem_data[cur_meta_subsystem][subsystem]['rxnNames']

        # Compute rxns
        cur_meta_subsystem_rxns = []
        for subsystem in meta_subsystems[cur_meta_subsystem]:
            cur_meta_subsystem_rxns += subsystem_data[cur_meta_subsystem][subsystem]['rxns']

        # Compute subsystems
        cur_meta_subsystem_subSystems = []
        for subsystem in meta_subsystems[cur_meta_subsystem]:
            cur_meta_subsystem_subSystems += subsystem_data[cur_meta_subsystem][subsystem]['subSystems']

        # Compute rules
        cur_meta_subsystem_rules = []
        for subsystem in meta_subsystems[cur_meta_subsystem]:
            cur_meta_subsystem_rules += subsystem_data[cur_meta_subsystem][subsystem]['rules']

        # Compute smat
        with open(os.path.join(model_dir, f'model.mets.json')) as fin:
            mets = json.load(fin)

        cur_meta_subsystem_smat = []
        for subsystem in meta_subsystems[cur_meta_subsystem]:
            for rxn, smat_row in subsystem_data[cur_meta_subsystem][subsystem]['smat'].items():
                rxn_idx = cur_meta_subsystem_rxns.index(rxn) + 1
                for met, coef in smat_row:
                    met_idx = mets.index(met) + 1
                    cur_meta_subsystem_smat.append([met_idx, rxn_idx, coef])

        smat_cmp_key = cmp_to_key(smat_cmp)
        cur_meta_subsystem_smat.sort(key=smat_cmp_key)

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

        formatted_json_str = custom_json_dump(cur_meta_subsystem_smat)
        with open(os.path.join(output_dir, 'model', f'model.S.json'), 'w') as file:
            file.write(formatted_json_str)

        # Copy remaining files
        for filename in json_copy_filenames:
            shutil.copy(os.path.join(model_dir, f'model.{filename}.json'), os.path.join(output_dir, 'model', f'model.{filename}.json'))

        shutil.copy(os.path.join(model_dir, 'metMeta.txt'), os.path.join(output_dir, 'model', 'metMeta.txt'))
        shutil.copy(os.path.join(PATH_2_RECON2_MAT, 'media', f'{media}.json'), os.path.join(output_dir, 'media', f'{media}.json'))

        gene_filenames = ['mouseSynonyms', 'non2uniqueEntrez', 'uniqueEntrez2Non', 'uniqueHumanEntrez', 'uniqueHumanGeneSymbol', 'uniqueMouseGeneSymbol_all', 'uniqueMouseGeneSymbol']
        for filename in gene_filenames:
            shutil.copy(os.path.join(PATH_2_RECON2_MAT, f'{filename}.json'), os.path.join(output_dir, f'{filename}.json'))

        # cur_meta_subsystem_mets does not order metabolites as in model.mets.json
        cur_meta_subsystem_mets = []
        for met_idx, rxn_idx, coef in cur_meta_subsystem_smat:
            cur_meta_subsystem_mets.append(mets[met_idx - 1])
        cur_meta_subsystem_mets = list(set(cur_meta_subsystem_mets))

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
        meta_subsystem_mets[cur_meta_subsystem] = cur_meta_subsystem_mets

        rxn_meta = pd.read_csv(os.path.join(model_dir, 'rxnMeta.txt'), delimiter='\t', quoting=csv.QUOTE_NONE)
        #cur_meta_subsystem_meta = rxn_meta[(rxn_meta['subsystem'].isin(meta_subsystems[cur_meta_subsystem])) & (rxn_meta['id'].isin(core_rxns))].reset_index(drop=True)
        cur_meta_subsystem_meta = rxn_meta[rxn_meta['subsystem'].isin(meta_subsystems[cur_meta_subsystem])].reset_index(drop=True)
        meta_subsystem_metas[cur_meta_subsystem] = cur_meta_subsystem_meta


    # ****************************************************************

    with open(os.path.join(model_dir, f'model.mets.json')) as fin:
        mets = json.load(fin)

    for meta_subsystem in selected_meta_subsystems:

        rxn_idx = len(meta_subsystem_rxns[meta_subsystem]) + 1

        for met in meta_subsystem_mets[meta_subsystem]:

            exchange_lb = -1000
            exchange_ub = 1000
            exchange_rev = 0
            exchange_rxnECNumber = ""
            exchange_rxnKeggID = ""
            exchange_rxnName = f'{met}_EXCHANGE ({meta_subsystem})'
            exchange_rxn = f'{met}_EXCHANGE_{meta_subsystem}'
            exchange_subSystem = meta_subsystem
            exchange_rule = ""

            # smat
            met_idx = mets.index(met) + 1
            exchange_smat = [met_idx, rxn_idx, 1.0]
            rxn_idx += 1

            
            print(exchange_rxn, exchange_smat, exchange_rule)

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

            idx = len(meta_subsystem_metas[meta_subsystem])
            meta_subsystem_metas[meta_subsystem].loc[idx] = [
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

        output_dir = os.path.join(meta_subsystem_models_dir, f'{meta_subsystem}_mat')
        if os.path.exists(output_dir) == False:
            os.mkdir(output_dir)
        if os.path.exists(os.path.join(output_dir, 'model')) == False:
            os.mkdir(os.path.join(output_dir, 'model'))
        if os.path.exists(os.path.join(output_dir, 'media')) == False:
            os.mkdir(os.path.join(output_dir, 'media'))

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

        formatted_json_str = custom_json_dump(meta_subsystem_smat[meta_subsystem])
        with open(os.path.join(output_dir, 'model', f'model.S.json'), 'w') as file:
            file.write(formatted_json_str)

        # Update reaction metadata file
        meta_subsystem_metas[meta_subsystem].to_csv(os.path.join(output_dir, 'model', 'rxnMeta.txt'), sep='\t', index=False)

        # Copy remaining files
        for filename in json_copy_filenames:
            shutil.copy(os.path.join(model_dir, f'model.{filename}.json'), os.path.join(output_dir, 'model', f'model.{filename}.json'))

        shutil.copy(os.path.join(model_dir, 'metMeta.txt'), os.path.join(output_dir, 'model', 'metMeta.txt'))
        shutil.copy(os.path.join(PATH_2_RECON2_MAT, 'media', f'{media}.json'), os.path.join(output_dir, 'media', f'{media}.json'))

        gene_filenames = ['mouseSynonyms', 'non2uniqueEntrez', 'uniqueEntrez2Non', 'uniqueHumanEntrez', 'uniqueHumanGeneSymbol', 'uniqueMouseGeneSymbol_all', 'uniqueMouseGeneSymbol']
        for filename in gene_filenames:
            shutil.copy(os.path.join(PATH_2_RECON2_MAT, f'{filename}.json'), os.path.join(output_dir, f'{filename}.json'))

        # **********************************************************************

    model_names = [f'{meta_subsystem}_mat' for meta_subsystem in selected_meta_subsystems]

    return meta_subsystem_models_dir, model_names