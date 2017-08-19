library(IPTM)
library(mvtnorm)
library(MCMCpack)

load('/Users/bomin8319/Desktop/IPTM/paper/code/Darenew.RData')
# 762 - 
attach(Dare)
Dare$node = 1:nrow(Dare$node)
Dare$text = Dare$text[762:length(Dare$edge)]
Dare$edge = Dare$edge[762:length(Dare$edge)]
Dare$edge = Dare$edge[-which(sapply(Dare$text, function(d){length(d)})==0)]
Dare$text = Dare$text[-which(sapply(Dare$text, function(d){length(d)})==0)]
mintime = Dare$edge[[1]][[3]]
for (n in 1:length(Dare$edge)){
  Dare$edge[[n]][3] = (Dare$edge[[n]][[3]] - mintime) / 3600
}
Dare$edge = lapply(Dare$edge, function(x){x[1:3]})


#D = sample(319:1456, 1, replace = FALSE)
D = 350
O = 1
R = 10
out = 1
n_B = 1500
n_d = 150
burn = c(1000,50)
thinning = c(10, 5)
#load("/Users/bomin8319/Desktop/IPTM/paper/code/Daretest.RData")
#output = Daretest
#initial = list(C = output$C, D = output$D[200], B = lapply(1:2, function(IP){output$B[[IP]][,500]}), Z = output$Z)
try = IPTM_predict.data(D, O, R, Dare$edge, Dare$node, Dare$text, Dare$vocab, nIP = 2, K = 20, sigma_Q = c(0.00005, 0.0001), 
						alpha = 2, mvec = rep(1/20, 20), betas = 2, nvec = rep(1/length(Dare$vocab), length(Dare$vocab)),
						 prior.b.mean = c(-3, rep(0, 24)), prior.b.var =  1 * diag(25), prior.delta = c(0, 1), 
                         out, n_B, n_d, burn, thinning, netstat = c("intercept", "dyadic", "triadic", "degree"), plot = FALSE, optimize = FALSE, 
                          niter = c(5,1))
save(try, file = "try.RData")                              