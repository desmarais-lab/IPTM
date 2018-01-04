library(IPTM)
load('Darenew.RData')
attach(Dare)
Dare$node = 1:nrow(Dare$node)
Dare$edge = lapply(Dare$edge, function(x){x[1:3]})

set.seed(1)
missing = matrix(0, nrow = length(Dare$edge), ncol = 3)
missing[sample(391:2247, 371, replace = FALSE), 1] = 1
missing[sample(391:2247, 371, replace = FALSE), 2] = 1
missing[sample(391:2247, 371, replace = FALSE), 3] = 1

for (nIP in 1:2) {
for (K in c(5, 10, 20, 30, 40, 50)) {
DarePPE = list()
for (i in 1:5){
  filename = paste0("Dare_full_",nIP,"_",K,"_ver",i,".RData")
  load(filename)
  initial = list()
  initial$alpha =  Daretest$alpha[length(Daretest$alpha)[1]]
  initial$mvec = Daretest$mvec[dim(Daretest$mvec)[1],]
  initial$delta = Daretest$delta[length(Daretest$delta)]
  initial$sigma_tau = Daretest$sigma_tau[length(Daretest$sigma_tau)]
  initial$b = t(vapply(1:nIP, function(IP) {Daretest$b[[IP]][,ncol(Daretest$b[[IP]])]}, rep(0, nrow(Daretest$b[[1]]))))
  initial$eta = t(vapply(1:nIP, function(IP) {Daretest$eta[[IP]][,ncol(Daretest$eta[[IP]])]}, rep(0, nrow(Daretest$eta[[1]]))))
  initial$l = Daretest$l
  initial$z = Daretest$z
  initial$u = Daretest$u
  DarePPE[[i]] = IPTM.PPE(missing, edge = Dare$edge, node = Dare$node, textlist = Dare$text, vocab= Dare$vocab, nIP = nIP, K = K,
                     Outer = 20, netstat = c("dyadic", "degree", "triadic"), timestat = c("dayofweek", "timeofday"), initial = initial)
}
filename = paste0("Dare_PPE_",nIP,"_",K,".RData")
save(DarePPE, file = filename)
}
}



for (i in 1:5){
  filename = paste0("Dare_full_",nIP,"_",K,"_ver",i,".RData")
  load(filename)
  initial = list()
  initial$alpha =  Daretest$alpha[length(Daretest$alpha)[1]]
  initial$mvec = Daretest$mvec[dim(Daretest$mvec)[1],]
  initial$delta = Daretest$delta[length(Daretest$delta)]
  initial$sigma_tau = Daretest$sigma_tau[length(Daretest$sigma_tau)]
  initial$b = t(vapply(1:nIP, function(IP) {Daretest$b[[IP]][,ncol(Daretest$b[[IP]])]}, rep(0, nrow(Daretest$b[[1]]))))
  initial$eta = t(vapply(1:nIP, function(IP) {Daretest$eta[[IP]][,ncol(Daretest$eta[[IP]])]}, rep(0, nrow(Daretest$eta[[1]]))))
  initial$l = Daretest$l
  initial$z = Daretest$z
  initial$proposal.var1 = diag(nrow(Daretest$b[[1]]))
  initial$proposal.var2 = diag(nrow(Daretest$eta[[1]]))
  initial$sigma.Q = Daretest$sigma.Q
  initial$u = Daretest$u
  DarePPE = IPTM.inference.PPE(missing, edge = Dare$edge, node = Dare$node, textlist = Dare$text, vocab= Dare$vocab, nIP = nIP, K = K,
    sigma.Q = c(0.0001, 0.001, 5, 0.1), alpha = 2, mvec = rep(1/K, K), beta = 2, prior.b = list(rep(0, 24), 10*diag(24)),
    prior.delta = c(-1000, 5),  prior.eta = list(rep(5, length(Dare$node)+2), 10*diag(length(Dare$node)+2)), prior.tau = 5,
    Outer = 20, Inner = c(2200,2200, 1100), burn = c(200,200, 100), thin = c(2,2,1),
    netstat = c("dyadic", "degree", "triadic"), timestat = c("dayofweek", "timeofday"), optimize = TRUE, initial = initial)
    filename = paste0("Dare_PPE_",nIP,"_",K,"_ver",i,".RData")
    save(DarePPE, file = filename)
}
