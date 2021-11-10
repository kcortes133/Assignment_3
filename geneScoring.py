# Author: Katherina Cortes
# Date: October 28,2021
# Purpose: Calculate gene score

# gene score is number of connections to other genes as defined by input

# Assumption: each gene only appears once in one locus
# genes cant be connected to themselves
# @param lociSubNs: List of dictionaries
# @param lociLists: list of list of genes separated by loci
# @param interactions:
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
# @param geneScores: a dictionary of lists of gene scores for each gene
# @returns geneSAvg: dictionary of average gene score for each gene
def getGeneScoreAvg(geneScores):
    geneSAvg = {}

    for gene in geneScores:
        geneSAvg[gene] = sum(geneScores[gene])/len(geneScores[gene])

    return geneSAvg

# @param geneAvg: dictionary of genes and gene scores
# @param lociLists: list of list of genes separated by loci
# @param numGenes: number of genes to get from each loci
# @returns genes: list of numGenes from each loci
def getTopLociGenes(geneAvg, lociLists, numGenes):
    genes = []
    # account for numGenes being more than length of loci list
    for loci in lociLists:
        lociScores = {}
        for gene in loci:
            lociScores[gene] = geneAvg[gene]
        lociSorted = sorted(lociScores, key=lambda k: lociScores[k], reverse=True)

        n = numGenes
        if numGenes > len(lociSorted):
            n = len(lociSorted)
        genes.extend(lociSorted[:n])

    return genes