
library(IPTM)
load('/Users/bomin8319/Desktop/IPTM/code/Darenew.RData')
# 762 - 
attach(Dare)
Dare$node = 1:nrow(Dare$node)
Dare$edge = Dare$edge[-which(sapply(Dare$text, function(d){length(d)})==0)]
Dare$text = Dare$text[-which(sapply(Dare$text, function(d){length(d)})==0)]
Dare$edge = lapply(Dare$edge, function(x){x[1:3]})


nIP = 2
K = 25

for (i in 1:5){
    set.seed(i)
    Daretest = IPTM.inference(edge = Dare$edge, node = Dare$node, textlist = Dare$text, vocab= Dare$vocab, nIP = nIP, K = K,
    sigma.Q = c(0.001, 0.01, 0.01), alpha = 2, mvec = rep(1/K, K), beta = 2, prior.b = list(rep(0, 24), diag(24)), 
    prior.delta = c(-3.5, 1),  prior.eta = list(rep(0, length(Dare$node)+2), diag(length(Dare$node)+2)), prior.tau = c(2, 1),
    Outer = 3, Inner = c(22000,22000, 1100), burn = c(2000,2000, 100), thin = c(20,20,1), 
    netstat = c("dyadic", "degree", "triadic"), timestat = c("dayofweek", "timeofday"),
    optimize = TRUE, initial = NULL)
    filename = paste0("Dare_full_",nIP,"_",K,"_ver",i,".RData")
    save(Daretest, file = filename)
}
