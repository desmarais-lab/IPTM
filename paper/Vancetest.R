library(Rcpp)
library(RcppArmadillo)
sourceCpp('/Users/bomin8319/Desktop/IPTM-master/pkg/src/sampler.cpp')
library(IPTM)
load('/Users/bomin8319/Desktop/IPTMpaper/data/Vancenew.RData')
attach(Vance)
Vance$node = 1:nrow(Vance$node)
for (n in 1:length(Vance$edge)){
  Vance$edge[[n]][3] = Vance$edge[[n]][[3]] / 3600
}
Vance$edge = lapply(Vance$edge, function(x){x[1:3]})

Vancetest <- IPTM_Inference(Vance$edge, Vance$node, Vance$text, Vance$vocab, nIP = 2, K = 5, sigma_Q = 0.02,
                       alpha = 2, mvec = rep(1/5, 5), betas = 2, nvec = rep(1/620, 620), prior.b.mean = c(-3, rep(0, 24)), 
                       prior.b.var = diag(25), prior.eta = c(0, 1), out = 100, n1 = 1, n2 = 1, n3 = 550, burn = 50, 
                       thin = 5, netstat = c("degree", "dyadic","triadic"), seed = 1, plot = TRUE, optimize = FALSE)

TableBetaIP(Vancetest)

PlotBetaIP(Vancetest)

par(mfrow = c(1,2))
PlotTopicIP(Vancetest, 10)

PlotTopic(Vancetest, 10)

TableWordIP(Vancetest, 10, text, vocab)


