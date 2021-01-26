import numpy as np
import pandas as pd
import logging
from .extensions import tsne_utils
from .. import utils
from .. import models
from ..globals import EXCHANGE_LIMIT
from sklearn.neighbors import NearestNeighbors
from sklearn.decomposition import PCA


def eval_reaction_penalties(expression_file, model, media,
                            species, args):
    """
    Main entry-point for evaluating reaction penalties

    Determines reaction penalties, on the given model, for
    the given expression data

    This version allows for blending of expression data between
    the individual sample and its neighbors.

    Parameters
    ==========
    expression_file: str
        Path to file containing the expression data

    model: compass.models.MetabolicModel
        Model to use for reaction penalty evaluation

    media : str or None
        Name of media to use

    species: str
        Species to use

    args: dict
        Additional arguments

    Returns
    =======
    penalties: pandas.DataFrame, Reactions x Cells
        Penalty values for each reaction in each cell

    """

    # Unpack extra arguments
    lambda_ = args['lambda']
    num_neighbors = args['num_neighbors']
    symmetric_kernel = args['symmetric_kernel']
    input_weights_file = args['input_weights']
    penalty_diffusion_mode = args['penalty_diffusion']
    and_function = args['and_function']
    input_knn = args['input_knn']
    output_knn = args['output_knn']
    latent_input = args['latent_space']

    expression = utils.read_data(expression_file) #pd.read_csv(expression_file, sep='\t', index_col=0)
    expression.index = expression.index.str.upper()  # Gene names to upper

    # If genes exist with duplicate symbols
    # Need to aggregate them out
    if not expression.index.is_unique:
        expression.index.name = "GeneSymbol"
        expression = expression.reset_index()
        expression = expression.groupby("GeneSymbol").sum()

    model = models.init_model(model, species=species,
                              exchange_limit=EXCHANGE_LIMIT,
                              media=media)

    # Evaluate reaction penalties
    input_weights = None
    if input_weights_file:
        input_weights = pd.read_csv(input_weights_file, sep='\t', index_col=0)
        # ensure same cell labels
        if len(input_weights.index & expression.columns) != \
                input_weights.shape[0]:
            raise Exception("Input weights file rows must have same sample "
                            "labels as expression columns")
        if len(input_weights.columns & expression.columns) != \
                input_weights.shape[1]:
            raise Exception("Input weights file columns must have same sample "
                            "labels as expression columns")

        input_weights = input_weights.loc[expression.columns, :] \
            .loc[:, expression.columns]

    reaction_penalties = eval_reaction_penalties_shared(
        model, expression, lambda_,
        num_neighbors=num_neighbors, symmetric_kernel=symmetric_kernel,
        and_function=and_function,
        penalty_diffusion_mode=penalty_diffusion_mode,
        input_weights=input_weights, 
        input_knn=input_knn, output_knn=output_knn, latent_input=latent_input)

    return reaction_penalties


def eval_reaction_penalties_shared(model, expression,
                                   lambda_, num_neighbors,
                                   symmetric_kernel, and_function,
                                   penalty_diffusion_mode,
                                   input_weights=None, 
                                   input_knn=None, output_knn=None,
                                   latent_input=None):
    """
    Determines reaction penalties, on the given model, for
    the given expression data

    This version allows for blending of expression data between
    the individual sample and its neighbors.

    Parameters
    ==========
    model: compass.models.MetabolicModel
        Model to use for reaction penalty evaluation

    expression: pandas.DataFrame
        Expression data, one row per gene, one column per sample

    sample_index: int
        Which column to evaluate penalties for

    lambda_: float
        The degree of blending.  0 (no sharing) to 1 (only use group)

    num_neighbors: int
        Effective number of neighbors for TSNE kernel

    symmetric_kernel: bool
        Whether or not to symmetrize the TSNE kernel (takes longer)

    and_function: str
        Which 'and' function to use for aggregating GPR associations

    penalty_diffusion_mode: str
        Either 'gaussian' or 'knn'.  Which mode to use to share penalty
        values between cells.  Ignored if `input_weights` is provided.

    input_weights: pandas.DataFrame
        Cells X Cells weights matrix used instead of computing one
        in this function call.
    """

    assert lambda_ >= 0 and lambda_ <= 1

    # Compute reaction expression for each sample

    reaction_expression = []
    for name, expression_data in expression.iteritems():
        reaction_expression.append(
            eval_reaction_expression_single(
                model, expression_data, and_function
            )
        )

    reaction_expression = pd.concat(reaction_expression, axis=1)

    # Compute weights between samples
    if input_weights is not None:
        weights = input_weights.values
    elif lambda_ == 0:
        weights = np.zeros(
            (reaction_expression.shape[1], reaction_expression.shape[1])
        )
    else:
        if latent_input is not None:
            expression_data = pd.read_csv(latent_input, sep='\t', index_col=0)
        # log scale and PCA expresion
        else:
            log_expression = np.log2(expression+1)
            model = PCA(n_components=min(
                log_expression.shape[0], log_expression.shape[1], 20)
            )
            pca_expression = model.fit_transform(log_expression.T).T
            pca_expression = pd.DataFrame(pca_expression,
                                        columns=expression.columns)
            expression_data = pca_expression
        if penalty_diffusion_mode == 'gaussian':
            weights = sample_weights_tsne_symmetric(
                expression_data, num_neighbors, symmetric_kernel)
        elif penalty_diffusion_mode == 'knn':
            weights = sample_weights_knn(expression_data, num_neighbors, input_knn, output_knn)
        else:
            raise ValueError(
                'Invalid value for penalty_diffusion_mode: {}'
                .format(penalty_diffusion_mode)
            )

    # Compute weights between samples
    neighborhood_reaction_expression = reaction_expression.dot(weights.T)
    neighborhood_reaction_expression.columns = reaction_expression.columns

    reaction_penalties = 1/(1+reaction_expression)
    neighborhood_reaction_penalties = 1/(1+neighborhood_reaction_expression)

    result = ((1-lambda_)*reaction_penalties +
              lambda_*neighborhood_reaction_penalties)

    result.index.name = "Reaction"

    return result


def eval_reaction_expression_single(model, expression_data, and_function):
    # type: (compass.models.MetabolicModel, pandas.Series)
    """
    Determines reaction penalties, on the given model, for
    the given expression data.

    Operates on a single cell
    """

    reaction_expression = model.getReactionExpression(
        expression_data, and_function=and_function)

    reaction_expression = pd.Series(reaction_expression,
                                    name=expression_data.name)

    reaction_expression[pd.isnull(reaction_expression)] = 0
    reaction_expression = np.log2(reaction_expression + 1)

    return reaction_expression


def sample_weights_tsne_symmetric(data, perplexity, symmetric):
    """
    Calculates p-values between samples using the procedure
    from tSNE

    As this version uses symmetric p-values, it must calculate
    the p-values for every sample vs every other sample

    data: pandas.DataFrame
        data matrix with samples as columns
    perplexity: float
        binary search perplexity target
    symmetric: boolean
        whether or not to symmetrize the weights
    """

    columns = data.columns
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
    if symmetric:
        pvals += pvals.T
        pvals /= 2

    # Make the rows sum to 1
    pvals /= pvals.sum(axis=1, keepdims=True)

    pvals = pd.DataFrame(pvals, index=columns, columns=columns)

    return pvals

def sample_weights_knn(data, num_neighbors, input_knn=None, output_knn=None):
    """
    Calculates cell-to-cell weights using knn on data

    data: pandas.DataFrame
        data matrix with samples as columns
    num_neighbors: float
        binary search perplexity target
    """

    columns = data.columns
    data = data.values
    ind = None
    if input_knn:
        logger = logging.getLogger('compass')
        ind_df = pd.read_csv(input_knn, sep='\t', index_col=0)
        #TBD: Check if the knn graph is specified as 2d array or adjacency list
        if ind_df.shape[0] != data.shape[1]: #incompatible dimensions for data, will recompute kNN
            logger.info("Input KNN has incompatible dimensions, recomputing KNN")
        elif not np.all(ind_df.index == columns): #See if the reactions need to be reindexed
            logger.info("Input KNN does not have the same samples as data, recomputing KNN")
        else:
            logger.info("KNN successfully loaded")
            ind = ind_df.to_numpy()
    if ind is None:
        nn = NearestNeighbors(n_neighbors=num_neighbors)
        nn.fit(data.T)
        ind = nn.kneighbors(return_distance=False)

    if output_knn:
        knn_df = pd.DataFrame(ind, index=columns)
        knn_df.to_csv(output_knn, sep='\t')

    weights = np.zeros((data.shape[1], data.shape[1]))

    weights[np.arange(data.shape[1]), ind.T] = 1
    weights /= num_neighbors  # So weights sum to 1

    weights = pd.DataFrame(weights, columns=columns, index=columns)

    return weights

def compute_knn(args):
    expression = utils.read_data(args['data'])
    expression.index = expression.index.str.upper()  # Gene names to upper

    # If genes exist with duplicate symbols
    # Need to aggregate them out
    if not expression.index.is_unique:
        expression.index.name = "GeneSymbol"
        expression = expression.reset_index()
        expression = expression.groupby("GeneSymbol").sum()
        
    log_expression = np.log2(expression+1)
    model = PCA(n_components=min(
        log_expression.shape[0], log_expression.shape[1], 20)
    )
    pca_expression = model.fit_transform(log_expression.T).T
    pca_expression = pd.DataFrame(pca_expression,columns=expression.columns)

    nn = NearestNeighbors(n_neighbors=args['num_neighbors'])
    nn.fit(pca_expression.T)
    ind = nn.kneighbors(return_distance=False)

    knn_df = pd.DataFrame(ind, index=expression.columns)
    knn_df.to_csv(args['output_knn'], sep='\t')


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

    