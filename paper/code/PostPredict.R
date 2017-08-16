library(IPTM)
library(mvtnorm)
library(MCMCpack)
set.seed(100)
nDocs = 5
node = 1:4
vocabulary = c("hi", "hello", "fine", "bye", "what")
nIP = 2
K = 4
nwords = 5
alpha = 2
mvec = rep(1/4, 4)
betas = 2
nvec = rep(1/5, 5)
netstat = c("intercept", "dyadic", "triadic", "degree")
P = 1 * ("intercept" %in% netstat) + 3 * (2 * ("dyadic" %in% netstat) + 4 * ("triadic" %in% netstat) + 2 *("degree" %in% netstat))
prior.b.mean = c(-3, rep(0, P-1))
prior.b.var = 0.05 * diag(P)
prior.delta = c(2.5, 0.0001)
sigma_Q = c(0.01, 0.001)
niters = c(1, 5500, 500, 500, 5)

b = lapply(1:nIP, function(IP) {
    prior.b.mean
    })
delta = prior.delta[1]
currentC = sample(1L:nIP, K, replace = TRUE)
supportD = gibbs.measure.support(length(node) - 1)
base.data = GenerateDocs.Gibbs(400, node, vocabulary, nIP, K, nwords, alpha, mvec, betas, nvec, b, delta, currentC, netstat, base.edge = list(),  base.text = list(), base = TRUE, support = supportD) 
base.edge = base.data$edge	   
base.text = base.data$text

D = 381
O = 5
R = 1
edge = base.edge
textlist = base.text
out = 5
n_B = 500
n_d = 50
burn = c(50, 5)
thinning = c(10, 5)
try = IPTM_predict.data(D, O, R, edge, node, textlist, vocabulary, nIP, K, sigma_Q, alpha, mvec, betas, nvec, 
							 prior.b.mean, prior.b.var, prior.delta, 
                                out, n_B, n_d, burn, thinning, netstat, plot = FALSE, optimize = FALSE, niter = c(1,1))
save(try, file = "try.RData")                              