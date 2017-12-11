library(IPTM)
setwd('~/Desktop/IPTM/code')
load('Darenew.RData')
# 762 - 
attach(Dare)
Dare$node = 1:nrow(Dare$node)
Dare$edge = Dare$edge[-which(sapply(Dare$text, function(d){length(d)})==0)]
Dare$text = Dare$text[-which(sapply(Dare$text, function(d){length(d)})==0)]
Dare$edge = lapply(Dare$edge, function(x){x[1:3]})

for (i in 1:1) {
    filename = paste0("Dare_full_",2,"_",25,"_ver",i,".RData")
	load(filename)
	nIP = 2
	K = 25
	initial = list()
	initial$alpha =  Daretest$alpha[length(Daretest$alpha)]
	initial$mvec = Daretest$mvec[dim(Daretest$mvec)[1],]
	initial$delta = Daretest$delta[length(Daretest$delta)]
    initial$sigma2_tau = Daretest$sigma2_tau[length(Daretest$sigma2_tau)]
	initial$b = lapply(1:nIP, function(IP) {
        Daretest$b[[IP]][,ncol(Daretest$b[[IP]])]
    	})
    initial$eta = lapply(1:nIP, function(IP) {
        Daretest$eta[[IP]][,ncol(Daretest$eta[[IP]])]
    })
	initial$l = Daretest$l
	initial$z = Daretest$z
	initial$bmat = Daretest$b
    initial$emat = Daretest$eta
	initial$dmat = Daretest$delta
    initial$sigma2_tau = Daretest$sigma2_tau
	initial$proposal.var = Daretest$proposal.var
	initial$iJi = Daretest$iJi
	initial$sigma.Q = Daretest$sigma.Q
	PPC = IPTM.PPC(O = 20, edge = Dare$edge, node = Dare$node, textlist = Dare$text, vocab = Dare$vocab,
                   nIP = nIP, K = K, sigma.Q = initial$sigma.Q, initial$alpha, initial$mvec, beta=2,
                   prior.b = list(rep(0, 24), diag(24)), prior.delta = c(-3.5, 1),
                   prior.eta = list(rep(0, 26), diag(26)), prior.tau = c(2, 1), Outer, Inner, burn, thin, netstat, timestat, optimize = FALSE, initial = NULL, inference
                   
                   
						 out = 1, n_B = 5500, n_d = 550, burn = c(500, 50), thinning = c(10, 1), netstat = c("dyadic", "degree", "triadic"),
						optimize = TRUE, initial = initial, Daretest)
	filename = paste0("PPC", i, "_nIP2_K25_Dare.RData")	
	save(PPC, file = filename)	
}
