library(IPTM)
load('Enron.RData')

setwd('~/Desktop/Enron_full')
for (i in 1:5) {
    filename = paste0("Enron_full_",2,"_",10,"_ver",i,".RData")
	load(filename)
	nIP = 2
	K = 10
	initial = list()
	initial$alpha =  Enrontest$alpha[dim(Enrontest$alpha)[1],1]
	initial$mvec = Enrontest$mvec[dim(Enrontest$mvec)[1],]
	initial$delta = Enrontest$D[length(Enrontest$D)]	
	initial$b = lapply(1:nIP, function(IP) {
	      Enrontest$B[[IP]][,ncol(Enrontest$B[[IP]])]
    	})
	initial$C = Enrontest$C
	initial$Z = Enrontest$Z
	initial$bmat = Enrontest$B	
	initial$dmat = Enrontest$D
	initial$proposal.var = Enrontest$proposal.var
	initial$iJi = Enrontest$iJi
	initial$sigma_Q = Enrontest$sigma_Q
	PPC = IPTM_check.data(O = 20, edge = Enron$edge, node = Enron$node, textlist = Enron$text, vocabulary =Enron$vocab, nIP = nIP, K = K,
						sigma_Q = Enrontest$sigma_Q, alpha = Enrontest$alpha[dim(Enrontest$alpha)[1],1], mvec = Enrontest$mvec[dim(Enrontest$mvec)[1],],
						betas = 2, nvec = rep(1/length(Enron$vocab), length(Enron$vocab)), prior.b.mean = rep(0,7), prior.b.var = diag(7),
						prior.delta = c(-3.5, 1), out = 1, n_B = 5500, n_d = 550, burn = c(500, 50), thinning = c(10, 1), netstat = c("intercept", "dyadic"), 
						optimize = TRUE, initial = initial)
	filename = paste0("PPC", i, ".RData")	
	save(PPC, file = filename)	
}
