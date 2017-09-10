library(IPTM)
load('/Users/bomin8319/Desktop/IPTM/paper/code/Darenew.RData')
# 762 - 
attach(Dare)
Dare$node = 1:nrow(Dare$node)
#Dare$text = Dare$text[762:length(Dare$edge)]
#Dare$edge = Dare$edge[762:length(Dare$edge)]
Dare$edge = Dare$edge[-which(sapply(Dare$text, function(d){length(d)})==0)]
Dare$text = Dare$text[-which(sapply(Dare$text, function(d){length(d)})==0)]
mintime = Dare$edge[[1]][[3]]
for (n in 1:length(Dare$edge)){
  Dare$edge[[n]][3] = (Dare$edge[[n]][[3]] - mintime) / 3600
}
Dare$edge = lapply(Dare$edge, function(x){x[1:3]})


Daretest <- IPTM_inference.LDA(Dare$edge, Dare$node, Dare$text, Dare$vocab, nIP = 2, K = 2, sigma_Q = c(0.0000005, 0.1),
                        alpha = 2, mvec = rep(1/2, 2), betas = 2, nvec = rep(1/length(Dare$vocab), length(Dare$vocab)), 
                        prior.b.mean = c(-1, rep(0, 24)),
                       prior.b.var = 1 * diag(25), prior.delta = c(0, 1), out = 100, n_B = 15000, n_d = 1500, burn = c(5000,500),
                       thinning = c(20,5), netstat = c("intercept", "dyadic", "degree", "triadic"), optimize = TRUE)
Daretest_LDA = Daretest
save(Daretest_LDA, file = "Daretest_LDA_K2.RData")

