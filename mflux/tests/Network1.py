"""
Create a simple metabolic model and test
"""

from mflux.models import (MetabolicModel, Reaction,
                          Species, Association, Gene, Compartment)

MODEL_UB =  10000
MODEL_LB = -10000

# Define the reactions

r1 = Reaction()
r1.id = 'R1'
r1.upper_bound = MODEL_UB
r1.lower_bound = 0

r2 = Reaction()
r2.id = 'R2'
r2.upper_bound = MODEL_UB
r2.lower_bound = MODEL_LB

r3 = Reaction()
r3.id = 'R3'
r3.upper_bound = MODEL_UB
r3.lower_bound = MODEL_LB

# Exchage Reactions
e1 = Reaction()
e1.id = 'E1'
e1.upper_bound = MODEL_UB
e1.lower_bound = 0

e2 = Reaction()
e2.id = 'E2'
e2.upper_bound = MODEL_UB
e2.lower_bound = 0

reactions = {x.id: x for x in [r1, r2, r3, e1, e2]}

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

e2.reactants.update({m4.id: 1.0})

# Create the model object
model = MetabolicModel('TestNetwork1')
model.reactions = reactions
model.species = species

