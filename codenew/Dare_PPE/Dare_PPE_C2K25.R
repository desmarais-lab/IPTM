library(IPTM)
load('Darenew.RData')
attach(Dare)
Dare$node = 1:nrow(Dare$node)
Dare$edge = lapply(Dare$edge, function(x){x[1:3]})

set.seed(1)
missing = matrix(0, nrow = length(Dare$edge), ncol = 3)
missing[sample(391:2247, 371, replace = FALSE), 1] = 1
missing[sample(391:2247, 371, replace = FALSE), 2] = 1
missing[sample(391:2247, 371, replace = FALSE), 3] = 1

nIP = 2
K = 25
for (i in 1:5){
    set.seed(i)
    Daretest = IPTM.inference(edge = Dare$edge, node = Dare$node, textlist = Dare$text, vocab= Dare$vocab, nIP = nIP, K = K,
    sigma.Q = c(0.0001, 0.001, 5, 1), alpha = 2, mvec = rep(1/K, K), beta = 2, prior.b = list(rep(0, 24), 10*diag(24)),
    prior.delta = c(-1000, 5),  prior.eta = list(rep(5, length(Dare$node)+2), 10*diag(length(Dare$node)+2)), prior.tau = c(1,1),
    Outer = 3, Inner = c(2200,2200, 1100), burn = c(200,200, 100), thin = c(2,2,1),
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
