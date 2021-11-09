# Author: Katherina Cortes
# Date: October 28,2021
# Purpose: Calculate gene score

# gene score is number of connections to other genes as defined by input

# Assumption: each gene only appears once in one locus
# genes cant be connected to themselves
#

# @param lociSubNs: List of dictionaries
# @param lociL: list of the locis
# @param interactions:
#
# @returns geneScores:
def getGeneScores(lociSubNs, lociL, interactions):
    geneScores = {}

    # get subNetwork
    for subN in lociSubNs:
        # pick a gene from a loci
        for lociGene in subN:
            # figure out which loci it comes from
            for l in lociL:
                if lociGene in l:
                    loci = l
                    break
            # get all the genes in the specified loci
            for gene in loci:
                # get number of connections in the subNetwork with chosen gene
                # make sure not looking for a connection with loci gene
                geneS = 0
                # check if there are edges with other genes in subnetwork
                for diffGene in subN:
                    # make sure not checking for edge with the gene from same loci
                    if diffGene != lociGene and diffGene in interactions[gene]:
                        geneS += 1

                if gene in geneScores:
                    geneScores[gene].append(geneS)
                else:
                    geneScores[gene] = [geneS]

    return geneScores


# Average the gene Score lists
# @param geneScores: a dictionary of lists
# @returns geneSAvg:
def getGeneScoreAvg(geneScores):
    geneSAvg = {}

    for gene in geneScores:
        geneSAvg[gene] = sum(geneScores[gene])/len(geneScores[gene])

    return geneSAvg

