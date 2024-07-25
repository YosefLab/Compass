import os
import pandas as pd
import shutil
import libsbml
from libsbml import SBMLDocument, SBMLNamespaces

from compass.globals import EXCHANGE_LIMIT, MODEL_DIR
from compass.models import load_metabolic_model, init_model

PATH_2_MOUSE_1 = os.path.join(MODEL_DIR, 'Mouse1')
PATH_2_CORE_RXNS = os.path.join(MODEL_DIR, 'Mouse1', 'core_reactions.txt')
PATH_2_CORE_RXNS_META = os.path.join(MODEL_DIR, 'Mouse1', 'core_reactions_md.csv')
PATH_2_CURRENCY_METS = os.path.join(MODEL_DIR, 'Mouse1', 'currency_mets.txt')


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
        output_dir = os.path.join(meta_subsystem_models_dir, meta_subsystem)
        if os.path.exists(output_dir) == False:
            os.mkdir(output_dir)
        if os.path.exists(os.path.join(output_dir, 'media')) == False:
            os.mkdir(os.path.join(output_dir, 'media'))

    media = args['media']

    # *****************************************************************************

    # Select reactions that belong to subsystem
    mouse1_model = load_metabolic_model('Mouse1', 'mus_musculus')
    # Read in Mouse1 SMAT
    mouse1_smat = mouse1_model.getSMAT()
    mouse1_smat_transposed = mouse1_model.getSMAT_transposed()

    sbmlDocument = libsbml.readSBMLFromFile(os.path.join(PATH_2_MOUSE_1, 'Mouse-GEM.xml'))
    xml_model = sbmlDocument.model
    
    with open(PATH_2_CORE_RXNS) as f:
        core_rxns = f.readlines()
    core_rxns = [r.strip() for r in core_rxns]

    core_rxn_meta = pd.read_csv(PATH_2_CORE_RXNS_META)

    model_dirs = []

    # Metabolites and reactions for each meta-subsystem
    meta_subsystem_rxn_ids = {}
    meta_subsystem_met_ids = {}

    # SBML model for each meta-subsystem
    meta_subsystem_sbml_doc = {}

    # Generate model for each meta subsystem
    for i in range(len(selected_meta_subsystems)):

        # Select current meta subsystem
        cur_meta_subsystem = selected_meta_subsystems[i]

        # Specify output directory for current meta subsystem
        output_dir = os.path.join(meta_subsystem_models_dir, cur_meta_subsystem)
        model_dirs.append(output_dir)

        # Compute rxns
        cur_meta_subsystem_rxn_ids = []
        for rxn, subSystem in zip(core_rxn_meta['ID'], core_rxn_meta['SUBSYSTEM']):
            if subSystem in meta_subsystems[cur_meta_subsystem]:
                cur_meta_subsystem_rxn_ids.append(rxn)

        cur_meta_subsystem_rxn_ids.sort()

        # Compute metabolites that are associated with meta-subsystem
        cur_meta_subsystem_met_ids = []

        for rxn_id in cur_meta_subsystem_rxn_ids:

            rxn = mouse1_model.reactions[rxn_id]

            # reactants
            for met_id, coefficient in rxn.reactants.items():
                metabolite = mouse1_model.species[met_id]
                assert metabolite.compartment in ['c', 'e', 'm']

            # products
            for met_id, coefficient in rxn.products.items():
                metabolite = mouse1_model.species[met_id]
                assert metabolite.compartment in ['c', 'e', 'm']

            cur_meta_subsystem_met_ids.append(met_id)

        cur_meta_subsystem_met_ids = list(set(cur_meta_subsystem_met_ids))
        cur_meta_subsystem_met_ids.sort()

        # Create SBML model for meta-subsystem
        sbmlNamespaces = SBMLNamespaces(3, 1)
        sbmlNamespaces.addPkgNamespace('fbc', 2)
        sbmlNamespaces.addPkgNamespace('groups', 1)

        sbmlDoc = SBMLDocument(sbmlNamespaces)
        sbmlDoc.setPackageRequired('fbc', False)
        sbmlDoc.setPackageRequired('groups', False)

        cur_meta_subsystem_model = sbmlDoc.createModel()

        # Import packages
        for gene_product in xml_model.getPlugin(0).getListOfGeneProducts():
            cur_meta_subsystem_model.getPlugin(0).addGeneProduct(gene_product)
        for objective in xml_model.getPlugin(0).getListOfObjectives():
            cur_meta_subsystem_model.getPlugin(0).addObjective(objective)
        for group in xml_model.getPlugin(1).getListOfGroups():
            cur_meta_subsystem_model.getPlugin(1).addGroup(group)

        # Set IDs
        cur_meta_subsystem_model.setMetaId(f'{cur_meta_subsystem}_Model')
        cur_meta_subsystem_model.setId(f'{cur_meta_subsystem}_Model')

        # Set unit definitions
        for unit_def in xml_model.getListOfUnitDefinitions():
            cur_meta_subsystem_model.addUnitDefinition(unit_def)

        # Set compartments
        for compartment in xml_model.getListOfCompartments():
            if compartment.getId() not in ['c', 'e', 'm']:
                continue
            cur_meta_subsystem_model.addCompartment(compartment)

        # Add metabolites associated with current meta-subsystem
        for species in xml_model.getListOfSpecies():
            if species.getId() not in cur_meta_subsystem_met_ids:
                continue
            cur_meta_subsystem_model.addSpecies(species)

        # Set ub/lb parameters
        for param in xml_model.getListOfParameters():
            cur_meta_subsystem_model.addParameter(param)

        # Add reactions associated with current meta-subsystem
        for reaction in xml_model.getListOfReactions():
            if reaction.getId() not in cur_meta_subsystem_rxn_ids:
                continue
            cur_meta_subsystem_model.addReaction(reaction)

        meta_subsystem_rxn_ids[cur_meta_subsystem] = cur_meta_subsystem_rxn_ids
        meta_subsystem_met_ids[cur_meta_subsystem] = cur_meta_subsystem_met_ids
        meta_subsystem_sbml_doc[cur_meta_subsystem] = sbmlDoc

        # Temporarily save model
        libsbml.writeSBMLToFile(sbmlDoc, os.path.join(output_dir, f'{cur_meta_subsystem}_model.xml'))

        shutil.copy(os.path.join(PATH_2_MOUSE_1, 'media', f'{media}.json'), os.path.join(output_dir, 'media', f'{media}.json'))

    # **********************************************************************

    # Add exchange reaction for all metabolites associated with subsystem

    for meta_subsystem in selected_meta_subsystems:
        meta_subsystem_model = init_model(meta_subsystem, species='mus_musculus', exchange_limit=EXCHANGE_LIMIT, media=args['media'],
                                    metabolic_model_dir=meta_subsystem_models_dir)
        
        meta_subsystem_xml_model = meta_subsystem_sbml_doc[meta_subsystem].model

        meta_subsystem_model_smat = meta_subsystem_model.getSMAT()

        for met in meta_subsystem_met_ids[meta_subsystem]:

            # First check if there is already an exchange reaction associated with the current metabolite
            has_exchange = False
            smat_row = meta_subsystem_model_smat[met]

            # Check all associated reactions
            for rxn_id, coef in smat_row:
                rxn = meta_subsystem_model.reactions[rxn_id]
                if rxn.is_exchange:
                    has_exchange = True
                    break

            if has_exchange == True:
                continue

            # Add exchange reaction
            exchange_rxn = meta_subsystem_xml_model.createReaction()
            exchange_rxn.setId(f'{met}_EXCHANGE_{meta_subsystem}')
            exchange_rxn.setName(f'{met}_EXCHANGE_{meta_subsystem}')
            met_obj = meta_subsystem_xml_model.getSpecies(met)
            exchange_rxn.addReactant(met_obj, 1.0)
            exchange_rxn.getPlugin('fbc').setUpperFluxBound('FB3N1000')
            exchange_rxn.getPlugin('fbc').setLowerFluxBound('FB1N1000')

            meta_subsystem_rxn_ids[meta_subsystem].append(f'{met}_EXCHANGE_{meta_subsystem}')

        # Specify output directory for current meta subsystem
        output_dir = os.path.join(meta_subsystem_models_dir, meta_subsystem)

        # Save list of reactions to run COMPASS on
        # One-hop neighbor reactions and additional exchange reactions do not need to be computed
        with open(os.path.join(output_dir, f'{meta_subsystem}_rxns.txt'), 'w') as f:
            for rxn_id in meta_subsystem_rxn_ids[meta_subsystem]:
                f.write(f'{rxn_id}\n')

    # **********************************************************************

    # For each non-currency metabolite in the meta-subsystem, add all associated reactions
    # Also add exchange reactions for newly added metabolites

    with open(PATH_2_CURRENCY_METS) as file:
        currency_mets = file.readlines()
    currency_mets = [m.strip() for m in currency_mets]

    for meta_subsystem in selected_meta_subsystems:

        meta_subsystem_xml_model = meta_subsystem_sbml_doc[meta_subsystem].model

        output_dir = os.path.join(meta_subsystem_models_dir, meta_subsystem)

        # Reactions and metabolites that are added to the meta-subsystem after adding all associated reactions
        new_rxn_ids = set()
        new_met_ids = set()

        for met_id in meta_subsystem_met_ids[meta_subsystem]:

            if met_id in currency_mets:
                continue

            associated_smat_row = mouse1_smat[met_id]

            for rxn_id, coef in associated_smat_row:
                if rxn_id not in core_rxns:
                    continue
                if rxn_id in meta_subsystem_rxn_ids[meta_subsystem]:
                    continue

                new_rxn_ids.add(rxn_id)
            
                associated_smat_transpose_row = mouse1_smat_transposed[rxn_id]
                for new_met_id, coef in associated_smat_transpose_row:
                    metabolite = mouse1_model.species[new_met_id]
                    assert metabolite.compartment in ['c', 'e', 'm']
                    if new_met_id in meta_subsystem_met_ids[meta_subsystem]:
                        continue
                    new_met_ids.add(new_met_id)

        new_rxn_ids = list(new_rxn_ids)
        new_rxn_ids.sort()
        new_met_ids = list(new_met_ids)
        new_met_ids.sort()

        # Add new metabolites to meta-subsystem
        for new_met_id in new_met_ids:
            meta_subsystem_met_ids[meta_subsystem].append(new_met_id)

            met_obj = xml_model.getSpecies(new_met_id)
            meta_subsystem_xml_model.addSpecies(met_obj)

        # Add associated reactions to meta-subsystem
        for new_rxn_id in new_rxn_ids:
            meta_subsystem_rxn_ids[meta_subsystem].append(new_rxn_id)

            rxn_obj = xml_model.getReaction(new_rxn_id)
            meta_subsystem_xml_model.addReaction(rxn_obj)

        # Add exchange reactions for new metabolites
        # No need to consider existing exchange reactions since it is impossible for
        # newly added metabolites to be associated with exchange reactions

        for new_met_id in new_met_ids:

            # Add exchange reaction
            exchange_rxn = meta_subsystem_xml_model.createReaction()
            exchange_rxn.setId(f'{new_met_id}_EXCHANGE_{meta_subsystem}')
            exchange_rxn.setName(f'{new_met_id}_EXCHANGE_{meta_subsystem}')
            met_obj = meta_subsystem_xml_model.getSpecies(new_met_id)
            exchange_rxn.addReactant(met_obj, 1.0)
            exchange_rxn.getPlugin('fbc').setUpperFluxBound('FB3N1000')
            exchange_rxn.getPlugin('fbc').setLowerFluxBound('FB1N1000')

        # Save final model
        libsbml.writeSBMLToFile(meta_subsystem_sbml_doc[meta_subsystem], os.path.join(output_dir, f'{meta_subsystem}_model.xml'))
        
        # Save media file
        shutil.copy(os.path.join(PATH_2_MOUSE_1, 'media', f'{media}.json'), os.path.join(output_dir, 'media', f'{media}.json'))

        # Save metadata file
        cur_meta_subsystem_rxn_meta = core_rxn_meta[core_rxn_meta['ID'].isin(meta_subsystem_rxn_ids[meta_subsystem])].reset_index(drop=True).sort_values('ID')
        cur_meta_subsystem_rxn_meta.to_csv(os.path.join(output_dir, f'{meta_subsystem}_rxn_meta.csv'), index=False)

    # **********************************************************************

    model_names = selected_meta_subsystems

    return meta_subsystem_models_dir, model_names