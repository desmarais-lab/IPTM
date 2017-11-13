library(IPTM)
load('/Users/bomin8319/Desktop/IPTM/code/Darenew.RData')
# 762 - 
attach(Dare)
Dare$node = 1:nrow(Dare$node)
Dare$edge = Dare$edge[-which(sapply(Dare$text, function(d){length(d)})==0)]
Dare$text = Dare$text[-which(sapply(Dare$text, function(d){length(d)})==0)]
mintime = Dare$edge[[1]][[3]]
for (n in 1:length(Dare$edge)){
  Dare$edge[[n]][3] = (Dare$edge[[n]][[3]] - mintime) / 3600
}
Dare$edge = lapply(Dare$edge, function(x){x[1:3]})


nIP = 2
K = 25

for (i in 1:5){
    set.seed(i)
    Daretest = IPTM.inference(edge = Dare$edge, node = Dare$node, textlist = Dare$text, vocab= Dare$vocab, nIP = nIP, K = K,
    sigma.Q = c(0.01, 0.01), alpha = 2, mvec = rep(1/K, K), beta = 2, prior.b = list(rep(0, 25), diag(25)), 
    prior.delta = c(-3.5, 1),  prior.eta = list(rep(0, 27), diag(27)), prior.tau = c(1, 1),
    Outer = 2, Inner = c(550,550), burn = c(50,50), thin = c(1,1), 
    netstat = c("intercept", "dyadic", "degree", "triadic"), timestat = c("dayofweek", "timeofday"),
    optimize = TRUE, initial = NULL)
    filename = paste0("Dare_full_",nIP,"_",K,"_ver",i,".RData")
    save(Daretest, file = filename)
}
