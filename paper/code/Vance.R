library(IPTM)
load('/Users/bomin8319/Desktop/IPTM/paper/code/Vancenew.RData')
attach(Vance)
Vance$node = 1:nrow(Vance$node)
mintime = Vance$edge[[1]][[3]]
for (n in 1:length(Vance$edge)){
  Vance$edge[[n]][3] = (Vance$edge[[n]][[3]] - mintime) / 3600
}
Vance$edge = Vance$edge[-which(sapply(Vance$text, function(d){length(d)})==0)]
Vance$text = Vance$text[-which(sapply(Vance$text, function(d){length(d)})==0)]

Vance$edge = lapply(Vance$edge, function(x){x[1:3]})
Vancetest <- IPTM_inference.data(Vance$edge, Vance$node, Vance$text, Vance$vocab, nIP = 2, K = 5, sigma_Q = c(0.01, 1),
                       alpha = 2, mvec = rep(1/5, 5), betas = 2, nvec = rep(1/620, 620), prior.b.mean = c(-3, rep(0, 24)), 
                       prior.b.var = 0.1 * diag(25), prior.delta = c(0, 1), out = 1, n_B = 5500, n_d = 550, burn = c(500, 50), 
                       thinning = c(10, 1), netstat = c("intercept", "dyadic", "degree", "triadic"), optimize = TRUE)

nDocs = length(Vancetest$edge2)
node = Vance$node
vocabulary = Vance$vocab
nIP = 2
K = 5
alpha = 2
 mvec = rep(1/5, 5)
betas = 2
nvec = rep(1/620, 620)
	b = lapply(1:nIP, function(IP) {
        rowMeans(Vancetest$B[[IP]])
    })
    delta = mean(Vancetest$D)
    currentC = Vancetest$C
currentZ = Vancetest$Z
iJi = Vancetest$iJi
netstat = c("intercept", "dyadic", "degree", "triadic")
base.edge = Vance$edge[-Vancetest$edge2]
base.text = Vance$text[-Vancetest$edge2]

PPC = GenerateDocs.PPC(nDocs, node, vocabulary, nIP, K, alpha, mvec, betas, nvec, iJi, b, delta, currentC, netstat, base.edge, base.text, currentZ, 10)