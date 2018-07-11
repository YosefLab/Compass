"""
Create a simple metabolic model and test

This model has a branch in it

m1 - R1 - m3
             \
               R3 -> m5
             /
m2 - R2 - m4

R3 compues reaction:
    2 m4 + m3 = m5
"""

import pandas as pd
import numpy as np
from compass.models import (MetabolicModel, Reaction,
                          Species, Association, Gene)
from tqdm import tqdm

MODEL_UB = 10000
MODEL_LB = -10000

# Define the reactions

for x in tqdm(range(1000)):
    r1 = Reaction()
    r1.id = 'R1'
    r1.upper_bound = np.random.randint(1, MODEL_UB)
    r1.lower_bound = np.random.randint(MODEL_LB, 0)

    r2 = Reaction()
    r2.id = 'R2'
    r2.upper_bound = MODEL_UB
    r2.lower_bound = np.random.randint(MODEL_LB, 0)

    r3 = Reaction()
    r3.id = 'R3'
    r3.upper_bound = np.random.randint(1, MODEL_UB)
    r3.lower_bound = 0

    # Exchage Reactions
    e1 = Reaction()
    e1.id = 'E1'
    e1.upper_bound = np.random.randint(1, MODEL_UB)
    e1.lower_bound = np.random.randint(MODEL_LB, 0)

    e2 = Reaction()
    e2.id = 'E2'
    e2.upper_bound = np.random.randint(1, MODEL_UB)
    e2.lower_bound = np.random.randint(MODEL_LB, 0)

    e5 = Reaction()
    e5.id = 'E5'
    e5.upper_bound = np.random.randint(1, MODEL_UB)
    e5.lower_bound = np.random.randint(MODEL_LB, 0)

    reactions = {x.id: x for x in [r1, r2, r3, e1, e2, e5]}

    # Define metabolites

    m1 = Species()
    m1.id = 'M1'

    m2 = Species()
    m2.id = 'M2'

    m3 = Species()
    m3.id = 'M3'

    m4 = Species()
    m4.id = 'M4'

    m5 = Species()
    m5.id = 'M5'

    species = {x.id: x for x in [m1, m2, m3, m4, m5]}

    # Connect reactions and metabolites

    e1.products.update({m1.id: 1.0})
    e2.products.update({m2.id: 1.0})

    r1.reactants.update({m1.id: 1.0})
    r1.products.update({m3.id: 1.0})

    r2.reactants.update({m2.id: 1.0})
    r2.products.update({m4.id: 1.0})

    r3.reactants.update({m3.id: 1.0, m4.id: 2})
    r3.products.update({m5.id: 1.0})

    e5.reactants.update({m5.id: 1.0})

    # Create some genes
    gene1 = Gene()
    gene1.id = 'a'
    gene1.name = 'a'

    gene2 = Gene()
    gene2.id = 'b'
    gene2.name = 'b'

    gene3 = Gene()
    gene3.id = 'c'
    gene3.name = 'c'

    # Associate with some reactions
    assc1 = Association()
    assc1.type = 'gene'
    assc1.gene = gene1
    r1.gene_associations = assc1

    assc2 = Association()
    assc2.type = 'gene'
    assc2.gene = gene2
    r2.gene_associations = assc2

    assc3 = Association()
    assc3.type = 'gene'
    assc3.gene = gene3
    r3.gene_associations = assc3

    # Create the model object
    model = MetabolicModel('TestNetwork2')
    model.reactions = reactions
    model.species = species

    # Create blank expression object
    expression = pd.DataFrame()
    expression['S1'] = pd.Series([1, 2, 1], index=['a', 'b', 'c'])

    from compass.compass import EXCHANGE_LIMIT
    model.limitExchangeReactions(limit=EXCHANGE_LIMIT)
    model.make_unidirectional()

    # What do we expect?
    # Reactions
    # R1 and R2 are made irreversible because of R3
    #    - _neg versions should be pruned
    # R1 max should be limted by e2 (half as much)
    # R2 max should be limited by e2
    expected_max = {}
    expected_max['R1_pos'] = min(
        model.reactions['R1_pos'].upper_bound,
        model.reactions['R2_pos'].upper_bound / 2,
        model.reactions['E2_pos'].upper_bound / 2,
        model.reactions['E1_pos'].upper_bound,
        model.reactions['R3_pos'].upper_bound,
        model.reactions['E5_pos'].upper_bound)

    expected_max['R2_pos'] = min(
        model.reactions['R2_pos'].upper_bound,
        model.reactions['R1_pos'].upper_bound * 2,
        model.reactions['R3_pos'].upper_bound * 2,
        model.reactions['E1_pos'].upper_bound * 2,
        model.reactions['E2_pos'].upper_bound,
        model.reactions['E5_pos'].upper_bound * 2)

    expected_max['R3_pos'] = min(
        expected_max['R1_pos'],
        expected_max['R2_pos'] / 2)

    expected_max['R1_neg'] = 0
    expected_max['R2_neg'] = 0

    expected_max['E1_pos'] = expected_max['R1_pos']
    expected_max['E2_pos'] = expected_max['R2_pos']
    expected_max['E5_pos'] = expected_max['R3_pos']

    expected_max['E1_neg'] = 0
    expected_max['E2_neg'] = 0
    expected_max['E5_neg'] = 0

    # Metabolite Uptake

    expected_max['M3_UPTAKE'] = min(
        min(
            model.reactions['E2_pos'].upper_bound / 2,
            model.reactions['R2_pos'].upper_bound / 2,
            model.reactions['R3_pos'].upper_bound,
            model.reactions['E5_pos'].upper_bound) +
        min(
            model.reactions['E1_neg'].upper_bound,
            model.reactions['R1_neg'].upper_bound),
        EXCHANGE_LIMIT)

    expected_max['M4_UPTAKE'] = min(
        min(
            model.reactions['E1_pos'].upper_bound * 2,
            model.reactions['R1_pos'].upper_bound * 2,
            model.reactions['R3_pos'].upper_bound * 2,
            model.reactions['E5_pos'].upper_bound * 2) +
        min(
            model.reactions['E2_neg'].upper_bound,
            model.reactions['R2_neg'].upper_bound),
        EXCHANGE_LIMIT)

    # - Limiting factor is E4 rate
    # - Should E1 limit be honored here?
    # Metabolite secretion
    # - Limited by E1 uptake
    expected_max['M3_SECRETION'] = min(
        model.reactions['R1_pos'].upper_bound,
        model.reactions['E1_pos'].upper_bound)

    expected_max['M4_SECRETION'] = min(
        model.reactions['R2_pos'].upper_bound,
        model.reactions['E2_pos'].upper_bound)

    from compass import compass
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
        pass
        # print('-------SUCCESS---------')
