library(IPTM)
load('/Users/bomin8319/Desktop/IPTM/enron/Enron.RData')
Enrontest <- IPTM_inference.data(Enron$edge, Enron$node, Enron$text, Enron$vocab, nIP = 2, K = 5, sigma_Q = c(0.01, 1),
                       alpha = 2, mvec = rep(1/5, 5), betas = 2, nvec = rep(1/7418, 7418), prior.b.mean = c(-1, rep(0, 6)), 
                       prior.b.var = 0.1 * diag(7), prior.delta = c(0, 1), out = 1, n_B = 5500, n_d = 550, burn = c(500, 50), 
                       thinning = c(10, 1), netstat = c("intercept", "dyadic"), optimize = TRUE)

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