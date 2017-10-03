set.seed(1)
selectD = sample(1105:2210, 200, replace = FALSE)
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
K = 75

predictIPTM = function(doc) {
	D = selectD[doc]
	Daretest = IPTM_inference.data(edge = Dare$edge[1:(D-1)], node = Dare$node, textlist = Dare$text[1:(D-1)], vocabulary = Dare$vocab, nIP = nIP, K = K,
						sigma_Q = c(0.001, 0.01, 0.05), alpha = 2, mvec = rep(1/K, K),
						betas = 2, nvec = rep(1/length(Dare$vocab), length(Dare$vocab)), prior.b.mean = rep(0,24), prior.b.var = diag(24),
						prior.delta = c(-3.5, 1), out = 25, n_B = 5500, n_d = 550, burn = c(500, 50), thinning = c(10, 1),
						netstat = c("dyadic", "degree", "triadic"), optimize = TRUE, initial = NULL)
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
						
PPE = IPTM_predict.data(D, O = 10, R = 50, edge = Dare$edge, node = Dare$node, textlist = Dare$text, vocabulary = Dare$vocab, nIP = nIP, K = K,
						sigma_Q = c(0.001, 0.01, 0.05), alpha = 2, mvec = rep(1/K, K),
						betas = 2, nvec = rep(1/length(Dare$vocab), length(Dare$vocab)), prior.b.mean = rep(0,24), prior.b.var = diag(24),
						prior.delta = c(-3.5, 1), out = 2, n_B = 5500, n_d = 550, burn = c(500, 50), thinning = c(10, 1),
						netstat = c("dyadic", "degree", "triadic"), optimize = TRUE, initial = initial)
filename = paste0("PPE_",D,"_",nIP,"_",K,".RData")
save(PPE, file = filename)
}

for (doc in 151:200) {
	predictIPTM(doc)
}
