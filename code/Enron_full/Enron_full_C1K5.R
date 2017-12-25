library(IPTM)
library(FastGP)
library(MCMCpack)
library(LaplacesDemon)
load('Enron.RData')
attach(Enron)

nIP = 1
K = 5
for (i in 1:5){
    set.seed(i)
    Enrontest = IPTM.inference.noIP(edge = Enron$edge, node = Enron$node, textlist = Enron$text, vocab= Enron$vocab,
               nIP = nIP, K = K, sigma.Q = c(0.001, 0.001, 5, 0.5), alpha = 2, mvec = rep(1/K, K), beta = 2,
               prior.b = list(rep(0, 24), 10*diag(24)), prior.delta = c(-10, 5),  prior.eta = list(rep(10, length(Enron$node)+2), 10*diag(length(Enron$node)+2)), prior.tau = 5,
               Outer = 3, Inner = c(2200,2200, 1100), burn = c(200,200, 100), thin = c(2,2,1),
               netstat = c("dyadic", "degree", "triadic"), timestat = c("dayofweek", "timeofday"), optimize = TRUE, initial = NULL, timeunit = (3600*24), tz = "America/Los_Angeles")
               filename = paste0("Enron_full_",nIP,"_",K,"_ver",i,".RData")
               save(Enrontest, file = filename)
}
