library(IPTM)

load('/Users/bomin8319/Desktop/IPTM/paper/code/Darenew.RData')
# 762 - 
attach(Dare)
Dare$node = 1:nrow(Dare$node)
Dare$edge = Dare$edge[-which(sapply(Dare$text, function(d){length(d)})==0)]
Dare$text = Dare$text[-which(sapply(Dare$text, function(d){length(d)})==0)]
mintime = Dare$edge[[1]][[3]]
for (n in 1:length(Dare$edge)){
  Dare$edge[[n]][3] = (Dare$edge[[n]][[3]] - mintime) / 3600
}
Dare$edge = lapply(Dare$edge, function(x){x[1:3]})

load("/Users/bomin8319/Desktop/IPTM/topic coherence/Daretest_IPTM_K20.RData")
nDocs = length(Daretest$edge2)
node = Dare$node
vocabulary = Dare$vocab
nIP = 2
K = 20
alpha = 2
 mvec = rep(1/K, K)
betas = 2
nvec = rep(1/length(vocabulary), length(vocabulary))
	b = lapply(1:nIP, function(IP) {
        rowMeans(Daretest$B[[IP]])
    })
    delta = mean(Daretest$D)
    currentC = Daretest$C
currentZ = Daretest$Z
iJi = Daretest$iJi
netstat = c("intercept", "dyadic")
base.edge = Dare$edge[-Daretest$edge2]
base.text = Dare$text[-Daretest$edge2]

PPC = GenerateDocs.PPC(nDocs, node, vocabulary, nIP, K, alpha, mvec, betas, nvec, iJi, b, delta, currentC, netstat, base.edge, base.text, currentZ, 1)