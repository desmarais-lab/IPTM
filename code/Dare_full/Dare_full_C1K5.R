library(IPTM)
library(FastGP)
library(MCMCpack)
library(LaplacesDemon)
load('Darenew.RData')
attach(Dare)
Dare$node = 1:nrow(Dare$node)
Dare$edge = lapply(Dare$edge, function(x){x[1:3]})

nIP = 1
K = 5
for (i in 1:5){
    set.seed(i)
    Daretest = IPTM.inference.noIP(edge = Dare$edge, node = Dare$node, textlist = Dare$text, vocab= Dare$vocab,
               nIP = nIP, K = K, sigma.Q = c(0.001, 0.001, 5, 0.5), alpha = 2, mvec = rep(1/K, K), beta = 2,
               prior.b = list(rep(0, 24), 10*diag(24)), prior.delta = c(-10, 5),  prior.eta = list(rep(10, length(Dare$node)+2), 10*diag(length(Dare$node)+2)), prior.tau = 5,
               Outer = 50, Inner = c(2200,2200, 1100), burn = c(200,200, 100), thin = c(2,2,1),
               netstat = c("dyadic", "degree", "triadic"), timestat = c("dayofweek", "timeofday"), optimize = TRUE, initial = NULL)
               filename = paste0("Dare_full_",nIP,"_",K,"_ver",i,".RData")
               save(Daretest, file = filename)
}
