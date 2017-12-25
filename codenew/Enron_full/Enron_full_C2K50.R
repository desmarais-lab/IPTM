library(IPTM)
library(FastGP)
library(MCMCpack)
library(LaplacesDemon)
load('Enron.RData')
attach(Enron)

nIP = 2
K = 50
for (i in 1:5){
    set.seed(i)
    Enrontest = IPTM.inference(edge = Enron$edge, node = Enron$node, textlist = Enron$text, vocab= Enron$vocab,
               nIP = nIP, K = K, sigma.Q = c(0.0001, 0.0001, 5, 0.1), alpha = 2, mvec = rep(1/K, K), beta = 2,
               prior.b = list(rep(0, 24), 10*diag(24)), prior.delta = c(-50, 5),  prior.eta = list(c(rep(30, length(Enron$node)), rep(0,2)), 10*diag(length(Enron$node)+2)), prior.tau = 5,
               Outer = 50, Inner = c(2200,2200, 550), burn = c(200,200, 50), thin = c(2,2,1),
               netstat = c("dyadic", "degree", "triadic"), timestat = c("dayofweek", "timeofday"), optimize = TRUE, initial = NULL, timeunit = (3600*24), tz = "America/Los_Angeles")
               filename = paste0("Enron_full_",nIP,"_",K,"_ver",i,".RData")
               save(Enrontest, file = filename)
}


