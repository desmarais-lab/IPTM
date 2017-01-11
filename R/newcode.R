library(mvtnorm)
library(MCMCpack)
library(inline)
library(Rcpp)

setwd('/Users/bomin8319/Desktop/SU16/Rcode/Vance')
#Inputs
########################
#1. edgelist : matrix of 3 columns (col1:sender, col2:receiver, col3=time in unix.time format)
load("Vance_edge.RData") 
edge <- edge
edge[,3]<-edge[,3]/3600

#2. nodelist : matrix with the first column being the ID of nodes (ID's starting from 1), NOTE: rest columns containing node info are optional -later
load("Vance_node.RData") 
node <- as.numeric(as.vector(node[,1]))

#3. textlist: list (of length=number of edges in total) containing the words in each edge (i.e. document)
load("Vance_text.RData") 
textlist <- textlist

#4. vocabulary: all vocabularies used over the corpus
load("Vance_vocab.RData")
vocabulary <- vocabulary
########################

###########################################################

unix.time(VanceMCMC<- MCMC(edge, node, textlist, vocabulary, 2, 10, 0.5, outer=2, n1=3, n2=3, n3=3300, burn=300, thin=3, opt=TRUE, seed=1))
save(VanceMCMC, file="VanceMCMC_asym1.RData")

