"""
Create a simple metabolic model and test

-> m1 - R1 -> m2 - R2 - m3 - R3 - m4 ->
"""

import pandas as pd
import numpy as np
from mflux.models import (MetabolicModel, Reaction,
                          Species, Association, Gene)
from tqdm import tqdm

MODEL_UB = 10000
MODEL_LB = -10000

# Define the reactions
for x in tqdm(range(1000)):

    r1 = Reaction()
    r1.id = 'R1'
    r1.upper_bound = MODEL_UB
    r1.lower_bound = 0

    r2 = Reaction()
    r2.id = 'R2'
    r2.upper_bound = np.random.randint(MODEL_UB)
    r2.lower_bound = np.random.randint(MODEL_LB, 0)

    r3 = Reaction()
    r3.id = 'R3'
    r3.upper_bound = np.random.randint(MODEL_UB)
    r3.lower_bound = np.random.randint(MODEL_LB, 0)

    # Exchage Reactions
    e1 = Reaction()
    e1.id = 'E1'
    e1.upper_bound = np.random.randint(MODEL_UB)
    e1.lower_bound = 0

    e4 = Reaction()
    e4.id = 'E4'
    e4.upper_bound = np.random.randint(MODEL_UB)
    e4.lower_bound = 0

    reactions = {x.id: x for x in [r1, r2, r3, e1, e4]}

    # Define metabolites

    m1 = Species()
    m1.id = 'M1'

    m2 = Species()
    m2.id = 'M2'

    m3 = Species()
    m3.id = 'M3'

    m4 = Species()
    m4.id = 'M4'

    species = {x.id: x for x in [m1, m2, m3, m4]}

    # Connect reactions and metabolites

    e1.products.update({m1.id: 1.0})

    r1.reactants.update({m1.id: 1.0})
    r1.products.update({m2.id: 1.0})

    r2.reactants.update({m2.id: 1.0})
    r2.products.update({m3.id: 1.0})

    r3.reactants.update({m3.id: 1.0})
    r3.products.update({m4.id: 1.0})

    e4.reactants.update({m4.id: 1.0})

    # Create some genes
    gene1 = Gene()
    gene1.id = 'a'
    gene1.name = 'a'

    gene2 = Gene()
    gene2.id = 'b'
    gene2.name = 'b'

    # Associate with some reactions
    assc1 = Association()
    assc1.type = 'gene'
    assc1.gene = gene1
    r1.gene_associations = assc1

    assc2 = Association()
    assc2.type = 'gene'
    assc2.gene = gene2
    r2.gene_associations = assc2

    # Create the model object
    model = MetabolicModel('TestNetwork1')
    model.reactions = reactions
    model.species = species

    # Create blank expression object
    expression = pd.DataFrame()
    expression['S1'] = pd.Series([1, 2], index=['a', 'b'])

    # What do we expect?
    # Reactions
    # E1, E4 only go forward
    # R2 and R3 are found to only go forward
    from mflux.compass import EXCHANGE_LIMIT
    model.limitExchangeReactions(limit=EXCHANGE_LIMIT)
    model.make_unidirectional()

    expected_max = {}
    expected_max['R1_pos'] = min(
        model.reactions['E1_pos'].upper_bound,
        model.reactions['R1_pos'].upper_bound,
        model.reactions['R2_pos'].upper_bound,
        model.reactions['R3_pos'].upper_bound,
        model.reactions['E4_pos'].upper_bound)

    expected_max['R2_pos'] = expected_max['R1_pos']
    expected_max['R3_pos'] = expected_max['R2_pos']
    expected_max['R2_neg'] = 0
    expected_max['R3_neg'] = 0

    expected_max['E1_pos'] = expected_max['R1_pos']
    expected_max['E4_pos'] = expected_max['R1_pos']

    # Metabolite Uptake

    expected_max['M2_UPTAKE'] = min(
        model.reactions['R2_pos'].upper_bound,
        model.reactions['R3_pos'].upper_bound,
        model.reactions['E4_pos'].upper_bound,
        EXCHANGE_LIMIT)
    expected_max['M3_UPTAKE'] = min(
        model.reactions['R3_pos'].upper_bound,
        model.reactions['E4_pos'].upper_bound,
        EXCHANGE_LIMIT)
    expected_max['M4_UPTAKE'] = 0

    # - Limiting factor is E4 rate
    # - Should E1 limit be honored here?
    # Metabolite secretion
    # - Limited by E1 uptake

    expected_max['M1_SECRETION'] = 0
    expected_max['M2_SECRETION'] = min(
        model.reactions['E1_pos'].upper_bound,
        model.reactions['R1_pos'].upper_bound)

    expected_max['M3_SECRETION'] = min(
        model.reactions['E1_pos'].upper_bound,
        model.reactions['R1_pos'].upper_bound,
        model.reactions['R2_pos'].upper_bound)

    from mflux import compass
    compass._clear_cache(model)
    out = compass.run_compass(model, expression)

    # Evaluate results
    preproc = compass._load_cache(model)
    errors = 0
    # print()
    for x in preproc:
        if x not in expected_max:
            errors += 1
            print('Unexpected reaction: ', x)
        elif preproc[x] != expected_max[x]:
            errors += 1
            print('Error in ' + x + ': ')
            print('   Expected:', expected_max[x])
            print('     Actual:', preproc[x])

    for x in expected_max:
        if x not in preproc:
            errors += 1
            print('Missing reaction: ', x)

    if errors > 0:
        print('Total Errors:', errors)
    else:
        # print('-------SUCCESS---------')
        pass

    # Questions
    # During reaction maximization, lets get all the exchange reactions
    # Then for metabolite secretion/uptake, only preproc added reactions
    # Only close other added reactions
