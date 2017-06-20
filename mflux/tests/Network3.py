"""
Create a simple metabolic model and test

This model can choose between two routes for a metabolite

                R3 - m3 -
               /
 - m1 - R1 -> m2
               \
                R4 - m4 -
               /
           - 2*m5

"""

import pandas as pd
import numpy as np
from tqdm import tqdm
from mflux.models import (MetabolicModel, Reaction,
                          Species, Association, Gene)
from mflux import compass

MODEL_UB = 10000
MODEL_LB = -10000


def create_model():

    # Define the reactions

    r1 = Reaction()
    r1.id = 'R1'
    r1.upper_bound = np.random.randint(1, MODEL_UB)
    r1.lower_bound = 0

    r3 = Reaction()
    r3.id = 'R3'
    r3.upper_bound = np.random.randint(1, MODEL_UB)
    r3.lower_bound = np.random.randint(MODEL_LB, 0)

    r4 = Reaction()
    r4.id = 'R4'
    r4.upper_bound = np.random.randint(1, MODEL_UB)
    r4.lower_bound = np.random.randint(MODEL_LB, 0)

    # Exchage Reactions
    e1 = Reaction()
    e1.id = 'E1'
    e1.upper_bound = np.random.randint(1, MODEL_UB)
    e1.lower_bound = np.random.randint(MODEL_LB, 0)

    e3 = Reaction()
    e3.id = 'E3'
    e3.upper_bound = np.random.randint(1, MODEL_UB)
    e3.lower_bound = np.random.randint(MODEL_LB, 0)

    e4 = Reaction()
    e4.id = 'E4'
    e4.upper_bound = np.random.randint(1, MODEL_UB)
    e4.lower_bound = np.random.randint(MODEL_LB, 0)

    e5 = Reaction()
    e5.id = 'E5'
    e5.upper_bound = np.random.randint(1, MODEL_UB)
    e5.lower_bound = np.random.randint(MODEL_LB, 0)

    reactions = {x.id: x for x in [r1, r3, r4, e1, e3, e4, e5]}

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
    e5.products.update({m5.id: 1.0})
    e3.reactants.update({m3.id: 1.0})
    e4.reactants.update({m4.id: 1.0})

    r1.reactants.update({m1.id: 1.0})
    r1.products.update({m2.id: 1.0})

    r3.reactants.update({m2.id: 1.0})
    r3.products.update({m3.id: 1.0})

    r4.reactants.update({m2.id: 1.0, m5.id: 2.0})
    r4.products.update({m4.id: 1.0})

    # Create some genes
    gene1 = Gene()
    gene1.id = 'g1'
    gene1.name = 'g1'

    gene4 = Gene()
    gene4.id = 'g4'
    gene4.name = 'g4'

    gene3 = Gene()
    gene3.id = 'g3'
    gene3.name = 'g3'

    # Associate with some reactions
    assc1 = Association()
    assc1.type = 'gene'
    assc1.gene = gene1
    r1.gene_associations = assc1

    assc4 = Association()
    assc4.type = 'gene'
    assc4.gene = gene4
    r4.gene_associations = assc4

    assc3 = Association()
    assc3.type = 'gene'
    assc3.gene = gene3
    r3.gene_associations = assc3

    # Create the model object
    model = MetabolicModel('TestNetwork3')
    model.reactions = reactions
    model.species = species

    from mflux.compass import EXCHANGE_LIMIT
    model.limitExchangeReactions(limit=EXCHANGE_LIMIT)
    model.make_unidirectional()

    # What do we expect?
    # Reactions
    # Max R1 limited by e1
    # Max R3 same as R1
    # Max R4 limited by e5
    # Note: R3/R4 are actually reversible
    # Metabolite Uptake
    # m1 uptake limit same as R1 limit
    # m2 uptake is sum of e3 limit + e5 limit / 2
    # m3 limit is same as e3 limit
    # m4 limit is same as e4 limit
    # m5 limit is R4 limit x 2
    # Metabolite secretion
    # m1 limit same as e1 limit
    # m2 limit same as R1 limit
    # ...incomplete

    expected_max = {}
    expected_max['R1_pos'] = min(
        model.reactions['R1_pos'].upper_bound,
        model.reactions['E1_pos'].upper_bound,
        min(
            model.reactions['R3_pos'].upper_bound,
            model.reactions['E3_pos'].upper_bound) +
        min(
            model.reactions['R4_pos'].upper_bound,
            model.reactions['E4_pos'].upper_bound,
            model.reactions['E5_pos'].upper_bound / 2)
    )

    expected_max['R3_pos'] = min(
        min(
            model.reactions['R1_pos'].upper_bound,
            model.reactions['E1_pos'].upper_bound) +
        min(
            model.reactions['E4_neg'].upper_bound,
            model.reactions['R4_neg'].upper_bound,
            model.reactions['E5_neg'].upper_bound / 2
        ),
        model.reactions['R3_pos'].upper_bound,
        model.reactions['E3_pos'].upper_bound
    )

    expected_max['R3_neg'] = min(
        model.reactions['R3_neg'].upper_bound,
        model.reactions['E3_neg'].upper_bound,
        model.reactions['R4_pos'].upper_bound,
        model.reactions['E4_pos'].upper_bound,
        model.reactions['E5_pos'].upper_bound / 2
    )

    expected_max['R4_pos'] = min(
        model.reactions['R4_pos'].upper_bound,
        model.reactions['E4_pos'].upper_bound,
        model.reactions['E5_pos'].upper_bound / 2,
        min(
            model.reactions['R1_pos'].upper_bound,
            model.reactions['E1_pos'].upper_bound
        ) +
        min(
            model.reactions['R3_neg'].upper_bound,
            model.reactions['E3_neg'].upper_bound
        )
    )

    expected_max['R4_neg'] = min(
        model.reactions['R4_neg'].upper_bound,
        model.reactions['E4_neg'].upper_bound,
        model.reactions['E5_neg'].upper_bound / 2,
        min(
            model.reactions['R3_pos'].upper_bound,
            model.reactions['E3_pos'].upper_bound
        )
    )

    expected_max['E1_pos'] = expected_max['R1_pos']
    expected_max['E1_neg'] = 0
    expected_max['E3_pos'] = expected_max['R3_pos']
    expected_max['E3_neg'] = expected_max['R3_neg']
    expected_max['E4_pos'] = expected_max['R4_pos']
    expected_max['E4_neg'] = expected_max['R4_neg']
    expected_max['E5_pos'] = expected_max['R4_pos'] * 2
    expected_max['E5_neg'] = expected_max['R4_neg'] * 2

    # Metabolite Uptake

    expected_max['M2_UPTAKE'] = min(
        EXCHANGE_LIMIT,
        min(
            model.reactions['R3_pos'].upper_bound,
            model.reactions['E3_pos'].upper_bound
        ) +
        min(
            model.reactions['R4_pos'].upper_bound,
            model.reactions['E4_pos'].upper_bound,
            model.reactions['E5_pos'].upper_bound / 2,
        )
    )

    # Metabolite secretion

    expected_max['M2_SECRETION'] = min(
        model.reactions['R1_pos'].upper_bound,
        model.reactions['E1_pos'].upper_bound
    ) + \
        min(
        model.reactions['R4_neg'].upper_bound,
        model.reactions['E4_neg'].upper_bound,
        model.reactions['E5_neg'].upper_bound / 2
    ) + \
        min(
        model.reactions['R3_neg'].upper_bound,
        model.reactions['E3_neg'].upper_bound
    )

    return model, expected_max


def run_model(model, expression, clear_cache=True):

    if clear_cache:
        compass._clear_cache(model)

    out = compass.run_compass(model, expression)
    preproc = compass._load_cache(model)

    return out, preproc


def eval_preproc(expected_max, preproc):

    # Evaluate results
    errors = 0
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

    return errors

# Create a Default expression object
expression = pd.DataFrame()
expression['S1'] = pd.Series([1, 1, 5], index=['g1', 'g3', 'g4'])


# Run 1000 times and look for errors
for x in tqdm(range(1000)):
    model, expected_max = create_model()
    results, preproc = run_model(model, expression)
    errors = eval_preproc(expected_max, preproc)

# Compare results between runs
all_results = []
all_preproc = []

N = 10
for i in range(N):

    if i == N-1:
        clear_cache = False  # Use the cache for the last one
    else:
        clear_cache = True

    results, preproc = run_model(model, expression, clear_cache)
    all_results.append(results)
    all_preproc.append(preproc)


# Find differences
all_preproc = [pd.Series(x) for x in all_preproc]
all_rxns = [x[0] for x in all_results]
all_uptake = [x[1] for x in all_results]
all_secretion = [x[2] for x in all_results]

all_preproc = pd.concat(all_preproc, axis=1, ignore_index=True)
print('\nPreprocessing Variance:  ')
print(all_preproc.var(axis=1))

all_rxns = pd.concat(all_rxns, axis=1, ignore_index=True)
print('\nReaction Variance:  ')
print(all_rxns.var(axis=1))

all_uptake = pd.concat(all_uptake, axis=1, ignore_index=True)
print('\nUptake Variance:  ')
print(all_uptake.var(axis=1))

all_secretion = pd.concat(all_secretion, axis=1, ignore_index=True)
print('\nSecretion Variance:  ')
print(all_secretion.var(axis=1))


# Does expression work the way it should?
# Create a Mixed expression object
multi_expression = pd.DataFrame()
multi_expression['S1'] = pd.Series([10, 2, 5], index=['g1', 'g3', 'g4'])
multi_expression['S2'] = pd.Series([10, 5, 2], index=['g1', 'g3', 'g4'])



model, expected_max = create_model()
results, preproc = run_model(model, multi_expression)
