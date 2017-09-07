library(IPTM)
load('/Users/bomin8319/Desktop/IPTM/paper/code/Darenew.RData')
attach(Dare)
Dare$node = 1:nrow(Dare$node)
Dare$edge = Dare$edge[-which(sapply(Dare$text, function(d){length(d)})==0)]
Dare$text = Dare$text[-which(sapply(Dare$text, function(d){length(d)})==0)]
mintime = Dare$edge[[1]][[3]]
for (n in 1:length(Dare$edge)){
  Dare$edge[[n]][3] = (Dare$edge[[n]][[3]] - mintime) / 3600
}
Dare$edge = lapply(Dare$edge, function(x){x[1:3]})



Daretest <- IPTM_inference.data(Dare$edge, Dare$node, Dare$text, Dare$vocab, nIP = 2, K = 50, sigma_Q = c(0.005, 0.1),
                        alpha = 2, mvec = rep(1/50, 50), betas = 2, nvec = rep(1/length(Dare$vocab), length(Dare$vocab)), 
                        prior.b.mean = c(-3, rep(0, 6)), 
                       prior.b.var = 1 * diag(7), prior.delta = c(0, 1), out = 1000, n_B = 15000, n_d = 1500, burn = c(10000,500), 
                       thinning = c(10,5), netstat = c("intercept", "dyadic"), optimize = TRUE, initial = NULL)

save(Daretest, file = "Daretest.RData")
