library(IPTM)
library(mvtnorm)
library(MCMCpack)
set.seed(1234)
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
netstat = c("intercept", "dyadic")
P = 1 * ("intercept" %in% netstat) + 3 * (2 * ("dyadic" %in% netstat) + 4 * ("triadic" %in% netstat) + 2 *("degree" %in% netstat))
prior.b.mean = c(-3, rep(0, P-1))
prior.b.var = 0.05 * diag(P)
prior.delta = c(2.5, 0.0001)
sigma_Q = c(0.01, 0.001)
niters = c(1, 3300, 100, 300, 5)
b = lapply(1:nIP, function(IP) {
    prior.b.mean
    })
delta = prior.delta[1]
currentC = sample(1L:nIP, K, replace = TRUE)
supportD = gibbs.measure.support(length(node) - 1)
base.data = GenerateDocs.Gibbs(100, node, vocabulary, nIP, K, nwords, alpha, mvec, betas, nvec, b, delta, currentC, netstat, base.edge = list(),  base.text = list(), base = TRUE, support = supportD) 
base.edge = base.data$edge	   
base.text = base.data$text

sigma_Q = c(0.1, 1)
niters = c(1, 2, 2, 0, 1)
nDocs = 1
iJicompare<- Comparison.Gibbs(1000, nDocs, node, vocabulary, nIP, K, nwords, alpha, mvec, betas, nvec, 
					prior.b.mean, prior.b.var, prior.delta, sigma_Q, niters, netstat, base.edge, base.text, seed = 100, generate_trace_plots = FALSE)
