import numpy as np
import pandas as pd
from .extensions import tsne_utils


def eval_reaction_penalties_shared(model, expression,
                                   sample_index, lambda_,
                                   perplexity, symmetric_kernel,
                                   and_function):
    """
    Determines reaction penalties, on the given model, for
    the given expression data

    This version allows for blending of expression data between
    the individual sample and its neighbors.

    Parameters
    ==========
    model: mflux.models.MetabolicModel
        Model to use for reaction penalty evaluation

    expression: pandas.DataFrame
        Expression data, one row per gene, one column per sample

    sample_index: int
        Which column to evaluate penalties for

    lambda_: float
        The degree of blending.  0 (no sharing) to 1 (only use group)

    perplexity: int
        Effective number of neighbors for TSNE kernel

    symmetric_kernel: bool
        Whether or not to symmetrize the TSNE kernel (takes longer)

    and_function: str
        Which 'and' function to use for aggregating GPR associations
    """

    assert lambda_ >= 0 and lambda_ <= 1

    # Compute reaction expression for each sample

    reaction_penalties = []
    for name, expression_data in expression.iteritems():
        reaction_penalties.append(
            eval_reaction_penalties(model, expression_data, and_function)
        )

    reaction_penalties = pd.concat(reaction_penalties, axis=1)

    # Compute weights between samples
    if symmetric_kernel:
        weights = sample_weights_tsne_symmetric(
            expression, perplexity, sample_index)
    else:
        weights = sample_weights_tsne_single(
            expression, perplexity, sample_index)

    # Compute weights between samples
    sample_reaction_penalties = reaction_penalties.iloc[:, sample_index]
    weighted_reaction_penalties = reaction_penalties.dot(
        weights.reshape((-1, 1))
    ).iloc[:, 0]

    result = ((1-lambda_)*sample_reaction_penalties +
              lambda_*weighted_reaction_penalties)

    return result


def eval_reaction_penalties(model, expression_data, and_function):
    # type: (mflux.models.MetabolicModel, pandas.Series)
    """
    Determines reaction penalties, on the given model, for
    the given expression data
    """

    reaction_expression = model.getReactionExpression(
        expression_data, and_function=and_function)

    reaction_expression = pd.Series(reaction_expression)
    reaction_expression[pd.isnull(reaction_expression)] = 0
    reaction_expression = np.log2(reaction_expression + 1)

    reaction_penalties = 1 / (1 + reaction_expression)

    reaction_penalties = reaction_penalties.dropna()

    return reaction_penalties


def sample_weights_tsne_single(data, perplexity, i):
    """
    Calculates p-values between samples using the procedure
    from tSNE

    data: numpy.ndarray
        data matrix with samples as columns
    perplexity: float
        binary search perplexity target
    i : int
        The index of the data point for which to calculate p-vals
    """

    if isinstance(data, pd.DataFrame):
        data = data.values

    # Calculate affinities (distance-squared) between samples
    diff = data - data[:, [i]]
    aff = (diff**2).sum(axis=0)
    aff = aff.astype('float32')

    # Run the tsne perplexity procedure
    pvals = tsne_utils.binary_search_perplexity_single(
        affinities=aff,
        neighbors=None,
        desired_perplexity=perplexity,
        verbose=0, sample=i)

    return pvals


def sample_weights_tsne_symmetric(data, perplexity, i):
    """
    Calculates p-values between samples using the procedure
    from tSNE

    As this version uses symmetric p-values, it must calculate
    the p-values for every sample vs every other sample

    data: numpy.ndarray
        data matrix with samples as columns
    perplexity: float
        binary search perplexity target
    i : int
        The index of the data point for which to calculate p-vals
    """

    if isinstance(data, pd.DataFrame):
        data = data.values

    # Calculate affinities (distance-squared) between samples

    sumData2 = np.sum(data**2, axis=0, keepdims=True)
    aff = -2*np.dot(data.T, data)
    aff += sumData2
    aff = aff.T
    aff += sumData2
    np.fill_diagonal(aff, 0)
    aff = aff.astype('float32')

    # Run the tsne perplexity procedure
    pvals = tsne_utils.binary_search_perplexity(
        affinities=aff,
        neighbors=None,
        desired_perplexity=perplexity,
        verbose=0
    )

    # Symmetrize the pvals
    pvals += pvals.T
    pvals /= 2

    return pvals[i, :]


def sample_weights_bio(bio_groups):
    """
    Computes sample weights based on biological group

    Weights are 1 if samples in the same bio_group and 0 otherwise

    bio_groups: list of str
        labels for each bio group
    """

    # Convert labels in numbers
    label_dict = {}
    i = 0
    for x in bio_groups:
        if x not in label_dict:
            label_dict[x] = i
            i = i+1

    label_vec = [label_dict[x] for x in bio_groups]
    label_vec = np.array([label_vec])  # creates a 1xN array

    # output distance matrix
    sumData2 = np.sum(label_vec**2, axis=0, keepdims=True)
    aff = -2*np.dot(label_vec.T, label_vec)
    aff += sumData2
    aff = aff.T
    aff += sumData2

    pvals = (aff == 0).astype('float')

    return pvals
