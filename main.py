# Author: Katherina Cortes
# Date: October 30, 2021
# Purpose: Take a tab delimited .gmz file of gene sets, a given tab-delimited STRING file
#   create a sub network of the gene interactions from the input file using the STRING file
#   get statistical significance

import argparse, logging, random, time
import networkCreation, fileParsing, statistics, geneScoring

# arguments:
#   - input file
#   - STRING file
#   - number of bins for cof network default=128
#   - number of networks to run
#   - fixed or quantile


parser = argparse.ArgumentParser(description='Get the statistical significance of the connection of genes vs random '
                                             'genes.')
parser.add_argument('genesFile', metavar='genes', type=str, default='Input.gmt.txt',
                    help='the input file of genes')
parser.add_argument('--interactionsFile', metavar='gene_interactions', type=str, default='STRING.txt',
                    help='the input file of gene interactions')
parser.add_argument('--numBins', type=int, default=128, help='the number of bins to separate edge densities into')
parser.add_argument('--numSubnetworks', type=int, default=5000, help='the number of subnetworks to make')

args = parser.parse_args()
print(args)

def main():
    start = time.time()
    random.seed(5)
    #inputF = 'input.gmt.txt'
    #stringF = 'STRING.txt'

    # read in networks
    lociLists = fileParsing.readInput(args.genesFile)
    interactions = fileParsing.makeInteractionNetwork(args.interactionsFile)
    network = fileParsing.makeNetwork(lociLists, interactions)


    networkSorted = sorted(network, key=lambda k: len(network[k]), reverse=True)

    highest = 0
    for i in networkSorted:
        if highest < 10:
            inLocis = {}
            #print(i,' : ', len(network[i]))
            for gene in network[i]:
                numL = 0
                for o in lociLists:
                    numL +=1
                    #if i in o:
                    #    print('LOCI: ', numL)
                    if gene in o:
                        if numL not in inLocis:
                            inLocis[numL] = 1
                        else:
                            inLocis[numL] +=1
            #print(inLocis)
        highest +=1

    # make loci subnetworks
    lociSubN = networkCreation.makeLociSubnetworks(args.numSubnetworks, network, lociLists)


    numBins = args.numBins
    # make bins for coFunctional subnetwork creation
    #qNetworkBins = networkCreation.makeQuantileBins(interactions, numBins)
    # fNetworkBins = makeFixedBins(interactions, numBins)

    # calculate gene scores
    geneScores = geneScoring.getGeneScores(lociSubN, lociLists, network)

    print(len(network))
    numG = 0
    for n in network:
        if len(network[n]) > 2:
            numG +=1
    print(numG)

    c = 0
    for gene in geneScores:
        if max(geneScores[gene]) > 1:
            #print(network[gene])
            if len(network[gene]) < max(geneScores[gene]):
                print('What the actual fuck ', gene)
            c +=1

    print(c)

    geneAvg = geneScoring.getGeneScoreAvg(geneScores)

    count = 0
    for a in geneAvg:
        if geneAvg[a] > 0:
            count += 1
    print(count)


    networkSorted = sorted(geneAvg, key=lambda k: geneAvg[k], reverse=True)

    highest = 0
    for i in networkSorted:
        if highest < 10:
            print(i,' : ', geneAvg[i], ' , ', max(geneScores[i]))
        highest +=1

    # make coFunctional random subnetworks
    # need 1000 populations where populations are the 5000 networks
    # //TODO
    #coFSubnetworks = networkCreation.makeCoFSubnetworks(interactions, qNetworkBins, lociSubN)

    # calculate the pvalue
    # probability edges using cof distribution is greater than avg of loci edged divided by # of random networks
    #pval = statistics.empiricalPVal(lociSubN, coFSubnetworks)

    # make a graph showing the edge density distributions
    #coFDensities = []
    #for network in coFSubnetworks:
    #    coFDensities.append(statistics.calcEdgeDensity(network))

    #lociDensities = []
    #for network in lociSubN:
    #    lociDensities.append(statistics.calcEdgeDensity(network))

    #statistics.overlappingHistogram(coFDensities, lociDensities)

    #print('P-val : ', pval)
    print(time.time() -start)

main()