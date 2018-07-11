

def clean_reactions(model):
    """
    Fixes reactions so that species are not listed as
    both reactants and products
    """

    for rr in model.reactions.values():
        overlap = set(rr.products.keys()) & set(rr.reactants.keys())
        for species in overlap:
            if rr.reactants[species] > rr.products[species]:
                rr.reactants[species] = (
                    rr.reactants[species] - rr.products[species]
                )
                del rr.products[species]
            elif rr.reactants[species] < rr.products[species]:
                rr.products[species] = (
                    rr.products[species] - rr.reactants[species]
                )
                del rr.reactants[species]
            else:
                del rr.products[species]
                del rr.reactants[species]


def limit_maximum_flux(model, new_limit):
    """
    Clips flux limits so that the max value
    is equal to new_limit
    """

    if new_limit < 0:
        new_limit = new_limit * -1

    old_limit = model.maximum_flux

    if old_limit > new_limit:
        for rr in model.reactions.values():

            if abs(rr.upper_bound) > new_limit:
                sign = 1 if rr.upper_bound >= 0 else -1
                rr.upper_bound = new_limit*sign

            if abs(rr.lower_bound) > new_limit:
                sign = 1 if rr.lower_bound >= 0 else -1
                rr.lower_bound = new_limit*sign
    else:
        for rr in model.reactions.values():

            if abs(rr.upper_bound) == old_limit:
                sign = 1 if rr.upper_bound >= 0 else -1
                rr.upper_bound = new_limit*sign

            if abs(rr.lower_bound) > old_limit:
                sign = 1 if rr.lower_bound >= 0 else -1
                rr.lower_bound = new_limit*sign

    model._calc_max_flux()

