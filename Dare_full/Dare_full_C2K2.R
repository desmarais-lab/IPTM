library(IPTM)
load('Darenew.RData')
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
K = 2

for (i in 1:5){
    set.seed(i)
Daretest = IPTM_inference.data(edge = Dare$edge, node = Dare$node, textlist = Dare$text, vocabulary = Dare$vocab, nIP = nIP, K = K,
						sigma_Q = c(0.001, 0.01), alpha = 2, mvec = rep(1/K, K),
						betas = 2, nvec = rep(1/length(Dare$vocab), length(Dare$vocab)), prior.b.mean = rep(0,25), prior.b.var = diag(25),
						prior.delta = c(-3.5, 1), out = 100, n_B = 5500, n_d = 550, burn = c(500, 50), thinning = c(10, 1),
						netstat = c("intercept", "dyadic", "degree", "triadic"), optimize = TRUE, initial = NULL)
filename = paste0("Dare_full_",nIP,"_",K,"_ver",i,".RData")
save(Daretest, file = filename)
}
