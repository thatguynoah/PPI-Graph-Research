# This is called a helper function file
# It contains a number of functions that we can import into the experimental file
# Files like this are nice because they're compact, easy to understand, and importable
# Please try to work on SuperImposeCorrelation to the best of your ability
# I understand this is a challenging task and you may walk out of the hour with very little to show - THAT IS OK
# as long as you are thinking and trying to debug/read documentation you're way ahead of the curve - this is very complex data science
# unfortunately, we don't have testing files right now, which makes coding tentative and difficult, just take your best swing

# this file will design the assets needed to complete the project - I will list the functions and their application
# 1. STRINGtoGraph - turns stringdb output (raw) into a networkx graph
# 2. SuperImposeCorrelation - take the correlation values and apply weights to each connection
# 3. TopBetweenness - returns n number of top ranked betweenness genes based on weighted graph

# Note:
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5244536/ - study finds spearman's better for differential correlation (slightly different)

# import the dependencies - note that you may need to install these with the pip command in cmd. If you have trouble ask Abby or google it
import numpy as np
import pandas as pd
import networkx as nx
from scipy.stats import spearmanr
# from numba import jit

# develop STRINGtoGraph
def STRINGtoGraph(S):
    # change column 1 name so that it does not include
    try:
        S = S.rename(columns = {'#node1': 'node1'})
        S = S.iloc[:, 0:2]
    except:
        pass

    # instantiate graph
    net = nx.Graph()

    # find all unique nodes in file and make into nodes
    nodeList = list(S.node1.unique())
    net.add_nodes_from(nodeList)

    # create a list of edges and incorporate them
    x, y = S.shape
    for r in range(x):
        row = tuple(S.iloc[r, :])
        net.add_edge(*row)

    # return statement
    return net

# develop SuperImposeCorrelationMatrix
# net is the network created by the first function
# expression is the normalized expression values associated with each sample
# classification is the status 0 (control) or 1 (disease)
# we may need to numba this (REALLY GOOD THING TO KNOW HOW TO DO ANYWAYS)
#@jit(nopython=True) # same as njit, may want to test speed of this, make sure it's faster
def SuperImposeCorrelation(net, expression, classification):
    # step 1 - get list of control samples
    controls = list(classification.loc[classification['status'] == 0].index)

    # step 2 - create list of node names in network
    nodeNames = list(net.nodes)  # get names of nodes

    # step 3 - from expression matrix, extract list of nodes (just unique names)
    geneNames = list(set(list(expression.index)))

    # step 4 - find all genes in nodeNames that are in geneNames
    nodeSet = set(nodeNames)
    correlationGenes = list(nodeSet.intersection(geneNames))
    popGenes = list(nodeSet.symmetric_difference(correlationGenes))

    # step 5 - remove popGenes from the network net
    net.remove_nodes_from(popGenes)

    # step 6 - construct a correlation matrix (square)
    CorrMat = pd.DataFrame(index=correlationGenes, columns=correlationGenes)

    # step 7 - fill in the correlation matrix with spearman values
    ctrlExpr = expression.loc[:, controls] # get the expression of the control samples
    for node in correlationGenes:
        adjGenes = list(net.adj[node]) # get nodes adjacent to node x
        try:
            contExpressionList1 = list(ctrlExpr.loc[node, :].iloc[0, :]) # only calculate this once per cycle; takes first index of occurrence
        except:
            contExpressionList1 = list(ctrlExpr.loc[node, :])
        # fill in CorrMat at necessary locations
        for gene in adjGenes:
            try:
                contExpressionList2 = list(ctrlExpr.loc[gene, :].iloc[0, :]) # gene; takes the first index of occurrence
            except:
                contExpressionList2 = list(ctrlExpr.loc[gene, :])
            correlationValue, _ = spearmanr(contExpressionList1, contExpressionList2)
            CorrMat.loc[node, gene] = correlationValue  # ! make sure this line works (may be depricated)

    # step 8 - populate the network with the matrix values squared
    for node1 in correlationGenes:
        adjGenes = list(net.adj[node1])
        for node2 in adjGenes:
            if CorrMat.loc[node1, node2] != np.nan:
                net[node1][node2]['weight'] = (CorrMat.loc[node1, node2])**2

    return net

# develop TopBetweenness
def TopBetweenness(net, n=10): # n = number of genes and values to return
    topGenes = [] # will contain top n genes
    centralityValues = [] # will contain corresponding centrality values

    centrality = dict(nx.betweenness_centrality(net, weight='weight')) # get dictionary of weights
    sortedCentralityList = sorted(centrality.items(), key=lambda x: x[1], reverse=True) # make list and order highest to lowest centrality (!) check

    returnData = sortedCentralityList[0:n] # return the n number of top genes in this datastruct

    for gene, value in returnData: # seperate each of the values into respective lists for easier access
        topGenes.append(gene)
        centralityValues.append(value)

    # return the list
    return topGenes, centralityValues


# eigenvector centrality: measure of node's influence in graph
# betweenness: how much of a bridge it is to other genes, eigenvector: total influence of node
# hierarchal clustering could be way to organize genetic function into phenotypic response (especially with ML)
def TopEigen(net, n=10):
    topGenes = [] # wll contain top n genes
    centralityValues = [] # will contain corresponding centrality values

    eigen = dict(nx.eigenvector_centrality(net, weight='weight')) #dictionary of eigen centralities
    sortedEigenCentralityList = sorted(eigen.items(), key=lambda x: x[1], reverse=True) # sort the list in descending value

    returnData = sortedEigenCentralityList[0:n]

    for gene, value in returnData:
        topGenes.append(gene)
        centralityValues.append(value)

    # return the list
    return topGenes, centralityValues

# find eigen and betweenness of select genes


# testing - do not run, you will get errors
if __name__=="__main__":
    pass
    # df = pd.read_csv(r"C:\Users\noahb\Downloads\string_interactions.tsv", delim_whitespace=True)
    # G = STRINGtoGraph(df)