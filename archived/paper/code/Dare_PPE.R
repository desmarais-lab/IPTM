library(IPTM)
load('/Users/bomin8319/Desktop/IPTM/paper/code/Darenew.RData')
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

load("/Users/bomin8319/Desktop/IPTM/paper/code/HPC/topic coherence/tDaretest_IPTM_K20.RData")
nIP = 2
K = 20
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
initial$proposal.var = Daretest$proposal.var
initial$sigma_Q = Daretest$sigma_Q
initial$iJi = Daretest$iJi


selectD = sample(1105:2210, 200, replace = FALSE)


predictIPTM = function(d) {
	initial$Z = Daretest$Z[1:(d-1)]
	initial$iJi = Daretest$iJi[1:(d-1)]
    PPE = IPTM_predict.data(d, O = 1, R = 10, edge = Dare$edge, node = Dare$node, textlist = Dare$text,   vocabulary = Dare$vocab, nIP = 2, K = 20,
						sigma_Q = Daretest$sigma_Q, alpha = Daretest$alpha[dim(Daretest$alpha)[1],1], mvec = Daretest$mvec[dim(Daretest$mvec)[1],],
						betas = 2, nvec = rep(1/length(Dare$vocab), length(Dare$vocab)), prior.b.mean = rep(0,7), prior.b.var = diag(7),
						prior.delta = c(-3.5, 1), out = 2, n_B = 5500, n_d = 550, burn = c(500, 50), thinning = c(10, 1), 
						netstat = c("intercept", "dyadic"), optimize = TRUE, initial = initial)
filename = paste0("PPE", d, ".RData")	
save(PPE, file = filename)
}


library(foreach)
library(doParallel)
cl = makeCluster(4)
registerDoParallel(cl)
foreach (doc = 1:200) %dopar% predictIPTM(selectD[doc])
stopCluster(cl)
