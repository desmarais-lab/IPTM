library(IPTM)
library(FastGP)
library(MCMCpack)
library(LaplacesDemon)
load('Darenew.RData')
attach(Dare)
Dare$node = 1:nrow(Dare$node)
Dare$edge = lapply(Dare$edge, function(x){x[1:3]})

nIP = 3
K = 10
for (i in 5:5){
    set.seed(i)
    Daretest = IPTM.inference(edge = Dare$edge, node = Dare$node, textlist = Dare$text, vocab= Dare$vocab, nIP = nIP, K = K,
    sigma.Q = c(0.0001, 0.0001, 5, 0.1), alpha = 2, mvec = rep(1/K, K), beta = 2, prior.b = list(rep(0, 24), 10*diag(24)),
    prior.delta = c(-50, 5),  prior.eta = list(c(rep(30, length(Dare$node)), rep(0,2)), 10*diag(length(Dare$node)+2)), prior.tau = 5,
    Outer = 50, Inner = c(2200,2200, 550), burn = c(200,200, 50), thin = c(2,2,1),
    netstat = c("dyadic", "degree", "triadic"), timestat = c("dayofweek", "timeofday"), optimize = TRUE, initial = NULL)
               filename = paste0("Dare_full_",nIP,"_",K,"_ver",i,".RData")
               save(Daretest, file = filename)
}
