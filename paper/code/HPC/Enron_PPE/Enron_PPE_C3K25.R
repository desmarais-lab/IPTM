set.seed(1)
selectD = sample(1963:3925, 200, replace = FALSE)
library(IPTM)
load('Enron.RData')
nIP = 3
K = 25



predictIPTM = function(doc) {
	D = selectD[doc]
	Daretest = IPTM_inference.data(edge = Enron$edge[1:(D-1)], node = Enron$node, textlist = Enron$text[1:(D-1)], vocabulary = Enron$vocab, nIP = nIP, K = K,
						sigma_Q = c(0.001, 0.01), alpha = 2, mvec = rep(1/K, K),
						betas = 2, nvec = rep(1/length(Enron$vocab), length(Enron$vocab)), prior.b.mean = rep(0,7), prior.b.var = diag(7),
						prior.delta = c(-3.5, 1), out = 50, n_B = 5500, n_d = 550, burn = c(500, 50), thinning = c(10, 1),
						netstat = c("intercept", "dyadic"), optimize = TRUE, initial = NULL)
	initial = list()
	initial$alpha =  Daretest$alpha[dim(Daretest$alpha)[1],1]
    initial$mvec = Daretest$mvec[dim(Daretest$mvec)[1],]
    initial$delta = Daretest$D[length(Daretest$D)]
    initial$b = lapply(1:nIP, function(IP) {
        Daretest$B[[IP]][,ncol(Daretest$B[[IP]])]
    })
    initial$C = Daretest$C
    initial$Z = Daretest$Z
    initial$bmat = Daretest$B
    initial$dmat = Daretest$D
    intial$sigma_Q = Daretest$sigma_Q
    initial$proposal.var = Daretest$proposal.var
    initial$iJi = Daretest$iJi
						
PPE = IPTM::IPTM_predict.data(D, O = 10, R = 50, edge = Enron$edge, node = Enron$node, textlist = Enron$text, vocabulary = Enron$vocab, nIP = nIP, K = K,
						sigma_Q = c(0.001, 0.01), alpha = 2, mvec = rep(1/K, K),
						betas = 2, nvec = rep(1/length(Enron$vocab), length(Enron$vocab)), prior.b.mean = rep(0,7), prior.b.var = diag(7),
						prior.delta = c(-3.5, 1), out = 2, n_B = 5500, n_d = 550, burn = c(500, 50), thinning = c(10, 1),
						netstat = c("intercept", "dyadic"), optimize = TRUE, initial = initial)
filename = paste0("PPE_",d,"_",nIP,"_",K,".RData")
save(PPE, file = filename)
}

for (doc in 1:200) {
	predictIPTM(doc)
}
