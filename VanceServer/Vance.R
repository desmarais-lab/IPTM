library(mvtnorm)
library(MCMCpack)
library(entropy)
library(Rcpp)
library(RcppArmadillo)
library(lda)
library(reshape)
library(code)
library(combinat)
source('core.R')
sourceCpp('sampler.cpp')
load('Vancenew.RData')
attach(Vance)
Vance$node = 1:nrow(Vance$node)
mintime = Vance$edge[[1]][[3]]
for (n in 1:length(Vance$edge)){
  Vance$edge[[n]][3] = (Vance$edge[[n]][[3]] - mintime) / 3600
}
Vance$edge = lapply(Vance$edge, function(x){x[1:3]})
set.seed(1)
Vancetest <- IPTM_inference.data(Vance$edge, Vance$node, Vance$text, Vance$vocab, nIP = 2, K = 5, sigma_Q = c(0.01, 1),
                       alpha = 2, mvec = rep(1/5, 5), betas = 2, nvec = rep(1/620, 620), prior.b.mean = c(-3, rep(0, 24)), 
                       prior.b.var = 0.1 * diag(25), prior.delta = c(0, 1), out = 1000, n_B = 11000, n_d = 500, burn = 1000, 
                       thinning = 20, netstat = c("intercept", "dyadic", "degree", "triadic"), plot = FALSE, optimize = TRUE)
save(Vancetest, file="VanceResult.RData")
