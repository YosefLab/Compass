"""
    Reimplementation of VISION algorithm
    For more details see https://www.nature.com/articles/s41467-019-12235-0
"""
import numpy as np
import pandas as pd
from sklearn.neighbors import NearestNeighbors
from sklearn.decomposition import PCA
from sklearn.cluster import KMeans
from scipy import sparse
import igraph
import leidenalg
from ..globals import PCA_SEED, LEIDEN_SEED, KMEANS_SEED
from ..utils import read_knn_ind, read_knn_dist

def microcluster(exprData, cellsPerPartition=10,
                         filterInput = "fano",
                         filterThreshold = None,
                         filterNumMad = 2,
                         latentSpace = None, inputKnn = None, 
                         inputKnnDistances = None, K = None, 
                         n_jobs = 1):
    
    if filterThreshold is None:
        filterThreshold = int(exprData.shape[1] * 0.05)
    if K is None:
        K = int(np.sqrt(exprData.shape[1]))
    
    #The latent space is not relevant if there is already a specified knn computed
    if inputKnn is not None:
        ind = read_knn_ind(inputKnn, exprData)
        if inputKnnDistances is not None:
            dist = read_knn_dist(inputKnnDistances, exprData)
        else:
            dist = np.zeros(ind.shape)
            #Could be vectorized using something like np.fromfunction but I am worried that numpy would copy each relevant entry of the array (ie an array of n-samples * k * dimensions of samples)
            for i in range(ind.shape[0]):
                dist[i] = np.linalg.norm(exprData.values.T[i] - exprData.values.T[ind[i]], axis=1)
        #Adj is the adjacency matrix of the knn graph
        adj = sparse.csr_matrix((dist.ravel(), (np.repeat(np.arange(ind.shape[0]), ind.shape[1]), ind.ravel())))

        
        if latentSpace is not None:
            res = pd.read_csv(latentSpace, sep='\t', index_col=0).T
        else: 
            res = exprData #Assume that the clustering is done on raw gene expression if not specified
    else:
        if latentSpace is not None:
            res = pd.read_csv(latentSpace, sep='\t', index_col=0).T
        else:
            log_expression = np.log2(exprData+1)
            
            #Determine latent space genes
            log_expression_filtered = applyFilters(log_expression, filterInput, 
                                                filterThreshold, filterNumMad)

            if len(log_expression_filtered) == 0:
                raise Exception("0 genes were selected with current threshold, set a lower one and rerun")

            model = PCA(n_components=min(log_expression_filtered.shape[0], log_expression_filtered.shape[1], 20),
                        random_state = PCA_SEED)
            if pd.__version__ >= '0.24':
                pca_expression = model.fit_transform(log_expression.to_numpy().T).T
            else:
                pca_expression = model.fit_transform(log_expression.values.T).T
            pca_expression = pd.DataFrame(pca_expression,
                                        columns=exprData.columns)
            res = pca_expression

        #Compute knn on PCA with K = min(30, K)
        nn = NearestNeighbors(n_neighbors=min(K, 30), n_jobs=n_jobs)
        nn.fit(res.T)
        dist, ind = nn.kneighbors()
        
        adj = nn.kneighbors_graph(mode='distance') #Should sparse graph of csr_format

    sigma = np.square(np.median(dist, axis=1)) #sigma <- apply(d, 1, function(x) quantile(x, c(.5))[[1]])
    adj.data = np.square(adj.data)
    for i in range(adj.shape[0]):
        adj[i] /= sigma[i]
    adj.data = np.exp(-adj.data)

    cl = leidenalg.find_partition(igraph.Graph.Weighted_Adjacency(adj), leidenalg.ModularityVertexPartition, seed=LEIDEN_SEED)

    clusters = {d:[] for d in np.unique(cl.membership)}
    for i in range(len(cl.membership)):
        clusters[cl.membership[i]].append(i)

    pools = readjust_clusters(clusters, res, cellsPerPartition = cellsPerPartition)
    #Conversion here fixes downstream type errors with numpy integers
    pools = {int(x):pools[x] for x in pools}  
    return pools

def filterGenesNovar(data):
    return data[data.var(axis=1) > 0]

def filterGenesThreshold(data, threshold):
    return data[(data > 0).sum(axis=1) >= threshold]

def filterGenesFano(data, num_mad=2):
    #Could sample columns at random to improve runtime
    mu = data.mean(axis=1).ravel()
    fano = data.var(axis=1).ravel() / mu
    
    mu_argsort = np.argsort(mu)
    mu_sort = mu[mu_argsort]
    fano_sort = fano[mu_argsort]
    
    N_QUANTS = 30
    #Q = np.linspace(0, len(mu), N_QUANTS+1, dtype=int)
    #This is just so that it exactly matches VISION code.
    m = int(len(mu) / N_QUANTS)
    Q = [i*m for i in range(N_QUANTS)] + [len(mu)]
    genePassList = []
    for i in range(N_QUANTS):
        mu_quant = mu_sort[Q[i]:Q[i+1]]
        mu_quant[mu_quant == 0] = 1
        fano_quant = fano_sort[Q[i]:Q[i+1]]
        
        mad_quant = np.median(np.abs(fano_quant - np.median(fano_quant)))
        gene_passes_quant = fano_quant > (np.median(fano_quant) + num_mad * mad_quant)
        
        genePassList.append(np.nonzero(gene_passes_quant)[0] + Q[i])
        
    gene_passes = np.concatenate(genePassList).ravel()
    return data.iloc[gene_passes]
    
def applyFilters(expr, filterInput, threshold, num_mad):
    if filterInput == "novar":
        return filterGenesNovar(expr)
    elif filterInput == "threshold":
        return filterGenesThreshold(expr, threshold)
    elif filterInput == "fano":
        expr = filterGenesThreshold(expr, threshold)
        return filterGenesFano(expr, num_mad)
    else:
        raise Exception("Filter not recognized")

def readjust_clusters(clusters, data, cellsPerPartition=100):
    NUM_PARTITIONS = round(data.shape[1] / cellsPerPartition)
    EPSILON = .15

    currPart = len(clusters)
    
    while currPart < ((1 - EPSILON)*NUM_PARTITIONS):
        clusterList = {}
        cluster_offset = 0
        for i in range(len(clusters)):
            # Apply kmeans clustering to existing cluster
            currCl = clusters[i]
            subData = data.iloc[:,currCl].T
            if len(currCl) > cellsPerPartition:
                nCl = KMeans(n_clusters=round(subData.shape[0] / cellsPerPartition), random_state=KMEANS_SEED).fit(subData)
            else:
                nCl = KMeans(n_clusters=1, random_state=KMEANS_SEED).fit(subData)
            newClust = nCl.predict(subData)
            # Gather cluster vector to list of clusters
            
            for j in range(len(newClust)):
                n = newClust[j] + cluster_offset
                sample_n = currCl[j]
                if n in clusterList:
                    clusterList[n].append(sample_n)
                else:
                    clusterList[n] = [sample_n]

            # Now add to cluster offset for next re-clustering
            cluster_offset = cluster_offset + max(newClust) + 1

        currPart = len(clusterList)
        clusters = clusterList


    return(clusters)

def pool_matrix_cols(data, pools):
    groups = {}
    for g in pools:
        for i in pools[g]:
            groups[data.columns[i]] = g
    groups = pd.DataFrame.from_dict(groups, orient='index', columns=['compass_microcluster']).T
    return data.append(groups).T.groupby("compass_microcluster").mean().T.rename(mapper=lambda x: 'cluster_'+str(int(x)), axis=1)

def unpool_columns(pooled_data, pools, data):
    unpooled_cols = []
    for cluster in pools:
        for sample in pools[cluster]:
            unpooled_cols.append(pooled_data.iloc[:,cluster].rename(data.columns[sample]))
    df = pd.concat(unpooled_cols, axis=1)
    df = df[data.columns]
    return df