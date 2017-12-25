library(IPTM)
library(FastGP)
library(MCMCpack)
library(LaplacesDemon)
load('/Users/bomin8319/Desktop/IPTM/code/Dare_full/Darenew.RData')
# 762 - 
attach(Dare)
Dare$node = 1:nrow(Dare$node)
#Dare$edge = Dare$edge[-which(sapply(Dare$text, function(d){length(d)})==0)]
#Dare$text = Dare$text[-which(sapply(Dare$text, function(d){length(d)})==0)]
Dare$edge = lapply(Dare$edge, function(x){x[1:3]})

setwd('~/Desktop/IPTM/code/')
nIP = 2
K = 25
for (i in 1:5){
    set.seed(i)
    Daretest = IPTM.inference(edge = Dare$edge, node = Dare$node, textlist = Dare$text, vocab= Dare$vocab, nIP = nIP, K = K,
    sigma.Q = c(0.0001, 0.0001, 5, 0.1), alpha = 2, mvec = rep(1/K, K), beta = 2, prior.b = list(rep(0, 24), 10*diag(24)),
    prior.delta = c(-50, 5),  prior.eta = list(c(rep(30, length(Dare$node)), rep(0,2)), 10*diag(length(Dare$node)+2)), prior.tau = 5,
    Outer = 3, Inner = c(2200,2200, 550), burn = c(200,200, 50), thin = c(2,2,1),
    netstat = c("dyadic", "degree", "triadic"), timestat = c("dayofweek", "timeofday"), optimize = TRUE, initial = NULL)
    filename = paste0("Dare_full_",nIP,"_",K,"_ver",i,".RData")
    save(Daretest, file = filename)
}

#PPC generate
setwd('~/Desktop/IPTM/code/PPC')
for (i in 1:5) {
	 filename = paste0("Dare_full_",2,"_",25,"_ver",i,".RData")
	load(filename)
	nIP = 2
	K = 25
	PPC = IPTM.PPC(Out = 20, edge = Dare$edge, node = Dare$node, textlist = Dare$text, vocab = Dare$vocab,
                   nIP = nIP, K = K, netstat = c("dyadic", "degree", "triadic"), timestat = c("dayofweek", "timeofday"),
                   inference = Daretest)
	filename = paste0("PPC", i, "_nIP2_K25_Dare.RData")	
	save(PPC, file = filename)	
}
