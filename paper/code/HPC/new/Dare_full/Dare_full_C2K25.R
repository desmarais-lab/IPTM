library(IPTM2)
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
K = 25

for (i in 1:5){
    set.seed(i)
    Daretest = IPTM_inference.data(edge = Dare$edge, node = Dare$node, textlist = Dare$text, vocabulary = Dare$vocab, nIP = nIP, K = K,
    sigma_Q = c(0.001, 0.01, 0.05), alpha = 2, mvec = rep(1/K, K),
    betas = 2, nvec = rep(1/length(Dare$vocab), length(Dare$vocab)), prior.b.mean = rep(0,24), prior.b.var = diag(24),
    prior.delta = c(-3.5, 1), out = 50, n_B = 5500, n_d = 550, burn = c(500, 50), thinning = c(10, 1),
    netstat = c("dyadic", "degree", "triadic"), optimize = TRUE, initial = NULL)
    filename = paste0("Dare_full_",nIP,"_",K,"_ver",i,".RData")
    save(Daretest, file = filename)
}


Daretest =  Daretestnew
save(Daretest, file = "/Users/bomin8319/Desktop/Daretest.RData")
load("/Users/bomin8319/Desktop/Dare_full_2_25_ver5.RData")

	initial = list()
	initial$alpha =  Daretest$alpha[dim(Daretest$alpha)[1],1]
    initial$mvec = Daretest$mvec[dim(Daretest$mvec)[1],]
    initial$delta = Daretest$D[length(Daretest$D)]
    initial$eta = Daretest$E[length(Daretest$E)]
    initial$b = lapply(1:nIP, function(IP) {
        Daretest$B[[IP]][,ncol(Daretest$B[[IP]])]
    })
    initial$C = Daretest$C
    initial$Z = Daretest$Z
    initial$bmat = Daretest$B
    initial$dmat = Daretest$D
    initial$emat = Daretest$E
    initial$proposal.var = Daretest$proposal.var
    initial$sigma_Q = Daretest$sigma_Q
    initial$iJi = Daretest$iJi


 Daretestnew = IPTM_inference.data(edge = Dare$edge, node = Dare$node, textlist = Dare$text, vocabulary = Dare$vocab, nIP = nIP, K = K,
    sigma_Q = c(0.001, 0.01, 0.05), alpha = 2, mvec = rep(1/K, K),
    betas = 2, nvec = rep(1/length(Dare$vocab), length(Dare$vocab)), prior.b.mean = rep(0,24), prior.b.var = diag(24),
    prior.delta = c(-3.5, 1), out = 1, n_B = 550000, n_d = 5500, burn = c(500000, 500), thinning = c(100, 10),
    netstat = c("dyadic", "degree", "triadic"), optimize = TRUE, initial = initial)
