# Statistical Significance of Genetic Functional Networks Using Edge Density 

## Goal
Compute the statistical significance of a population of subnetworks from a
set of FA loci compared to the STRING database. 

## Description
Input is two tab-delimited files named 'Input.gmt' and 'STRING.txt'. Input.gmt is a 
tab-delimited file formatted as Broad Institute’s Gene Matrix Transformed (GMT). STRING.txt 
is a version of the STRING database of known and predicted protein-protein interactions. 
Each line in the string file represents an edge in the network between the protein genes 
formatted protein'\t'protein'\t'weight'\n' where weight is the strength of their functional 
similarity. 
Random subnetworks are created from the Input.gmt.txt file where there is one random node
picked from each loci. For each node in the loci subnetwork an equivalent node in the 
STRING database is picked. Equivalent nodes from the STRING database are determined by 
organizing the nodes into quantile bins based on their edge density and picking a node 
that has a similar density to the one in the loci subnetwork.
These subnetworks are then compared for statistical significance to determine if the genes
chosen from the loci subnetwork are more functionally connected than a random set of genes.
This is done by computing the empirical p-value based on the subnetworks densities.

## Install
scipy
matplotlib.pyplot
operator
functools 
networkx
nxviz

## Usage
#### Python Usage
```python
import networkCreation, fileParsing, networkVisualization, geneScoring

inputF = 'input.gmt.txt'
stringF = 'STRING.txt'

topOverallGenes = False

# read in networks
lociLists = fileParsing.readInput(inputF)
interactions = fileParsing.makeInteractionNetwork(stringF)
network = fileParsing.makeNetwork(lociLists, interactions)

# make loci subnetworks
# make loci subnetworks
lociSubN = networkCreation.makeLociSubnetworks(5000, network, lociLists)

# calculate gene scores and sort genes by score
geneScores = geneScoring.getGeneScores(lociSubN, lociLists, network)
geneAvg = geneScoring.getGeneScoreAvg(geneScores)
networkSorted = sorted(geneAvg, key=lambda k: geneAvg[k], reverse=True)

# get top numGenes from each loci
# make network with genes
if not topOverallGenes:
    genes = geneScoring.getTopLociGenes(geneAvg, lociLists, 3)
    visualNetwork = networkVisualization.makeCrossLociNetwork(genes, network, lociLists)

# get top numGenes regardless of loci
# make network with genes
if topOverallGenes:
    genes = networkSorted[:10]
    visualNetwork = networkVisualization.makeCrossLociNetwork(genes, network, lociLists)

# write to output file the genes loci and gene score of genes in the network
networkVisualization.outputGeneScores(geneAvg, genes, 'topGeneScores.txt', lociLists)
# make graph with specified genes
# 1. gene score - size of node
# 2. loci of gene - color of node
# 3. weight of edge - darkness of edge
# make network between loci genes no edges between genes in same loci
graph = networkVisualization.makeGraph(visualNetwork, lociLists, geneAvg)
networkVisualization.visualizeGraph(graph, args.topGenes)
```

#### Command Line Usage
```commandline
$ python main.py yourInputFile.gmt.txt
$ python main.py input.gmt.txt --topGenes=True --numGenes=10
$ 
```
#### Example of Top Genes with Loci and Gene Score

#### Example Figure of Prix Fixe subnetwork

## Input
1. Input.gmt
- disjoint gene sets
- tab-delimited file Input.gmt
- First two columns describe row of gene set
- format: (https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats/)
2. STRING.txt
- STRING database of known and predicted protein-protein interactions
- tab-delimited
- each line represents an edge in network between two genes
- weighted by strength of functional similarity

## Output 

