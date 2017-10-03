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

setwd('~/Desktop/IPTM/Dare_full')
for (i in 1:5) {
    filename = paste0("Dare_full_",2,"_",10,"_ver",i,".RData")
	load(filename)
	nIP = 2
	K = 10
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
	initial$iJi = Daretest$iJi
	initial$sigma_Q = Daretest$sigma_Q
	PPC = IPTM_check.data(O = 20, edge = Dare$edge, node = Dare$node, textlist = Dare$text, vocabulary = Dare$vocab, nIP = nIP, K = K,
						sigma_Q = Daretest$sigma_Q, alpha = Daretest$alpha[dim(Daretest$alpha)[1],1], mvec = Daretest$mvec[dim(Daretest$mvec)[1],],
						betas = 2, nvec = rep(1/length(Dare$vocab), length(Dare$vocab)), prior.b.mean = rep(0,25), prior.b.var = diag(25),
						prior.delta = c(0, 1), out = 1, n_B = 5500, n_d = 550, burn = c(500, 50), thinning = c(10, 1), netstat = c("intercept", "dyadic", "degree", "triadic"), 
						optimize = TRUE, initial = initial, Daretest)
	filename = paste0("PPC", i, "_nIP2_K10_Dare.RData")	
	save(PPC, file = filename)	
}


for (i in 1:5) {
  filename = paste0("Dare_full_",2,"_",25,"_ver",i,".RData")
  load(filename)
  nIP = 2
  K = 25
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
  initial$iJi = Daretest$iJi
  initial$sigma_Q = Daretest$sigma_Q
  PPC = IPTM_check.data(O = 20, edge = Dare$edge, node = Dare$node, textlist = Dare$text, vocabulary = Dare$vocab, nIP = nIP, K = K,
                        sigma_Q = Daretest$sigma_Q, alpha = Daretest$alpha[dim(Daretest$alpha)[1],1], mvec = Daretest$mvec[dim(Daretest$mvec)[1],],
                        betas = 2, nvec = rep(1/length(Dare$vocab), length(Dare$vocab)), prior.b.mean = rep(0,25), prior.b.var = diag(25),
                        prior.delta = c(0, 1), out = 1, n_B = 5500, n_d = 550, burn = c(500, 50), thinning = c(10, 1), netstat = c("intercept", "dyadic", "degree", "triadic"), 
                        optimize = TRUE, initial = initial, Daretest)
  filename = paste0("PPC", i, "_nIP2_K25_Dare.RData")	
  save(PPC, file = filename)	
}
