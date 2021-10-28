# Statistical Significance of Genetic Functional Networks Using Edge Density 

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

## Goal
Compute the statistical significance of a population of subnetworks from a 
set of FA loci compared to the STRING database. 

## Install
scipy
matplotlib.pyplot
operator
functools 

## Usage
#### Python Usage
```python
import networkCreation, fileParsing, statistics

inputF = 'input.gmt.txt'
stringF = 'STRING.txt'

# read in networks
lociLists = fileParsing.readInput(inputF)
interactions = fileParsing.makeInteractionNetwork(stringF)
network = fileParsing.makeNetwork(lociLists, interactions)

# make loci subnetworks
lociSubN = networkCreation.makeLociSubnetworks(5000, network, lociLists)

numBins = 128
# make bins for coFunctional subnetwork creation
qNetworkBins = networkCreation.makeQuantileBins(interactions, numBins)
# fNetworkBins = makeFixedBins(interactions, numBins)

# make coFunctional random subnetworks
coFSubnetworks = networkCreation.makeCoFSubnetworks(interactions, qNetworkBins, lociSubN)

# calculate the pvalue
# probability edges using cof distribution is greater than avg of loci edged divided by # of random networks
pval = statistics.empiricalPVal(lociSubN, coFSubnetworks)

# make a graph showing the edge density distributions
coFDensities = []
for network in coFSubnetworks:
    coFDensities.append(statistics.calcEdgeDensity(network))

lociDensities = []
for network in lociSubN:
    lociDensities.append(statistics.calcEdgeDensity(network))

statistics.overlappingHistogram(coFDensities, lociDensities)

print('P-val : ', pval)

```

#### Command Line Usage
```commandline
$ python main.py yourInputFile.gmt.txt
$ python main.py Input.gmt.txt --interactions STRING.txt 
$ python main.py yourInputFile.gmt.txt --numBins 100 --numSubnetworks 1000
```
#### Example Figure of Subnetwork Densities
![pval edge densityies](https://user-images.githubusercontent.com/22487858/137910105-d4a6deab-ab47-49ea-b379-22f4956b8986.png)

#### Example Figure of P-Val and Random Subnetwork Density Distribution
![pval](https://user-images.githubusercontent.com/22487858/137910088-39abf6ee-49b9-40ae-99a0-4b248dac8abf.png)

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
P-Value for edge density distribution of the subnetworks. Edge density histogram  for the two networks.
Edge density histogram for the random subnetworks with a dashed line showing the average edge density
for the FA loci subnetworks and the corresponding p-value.
