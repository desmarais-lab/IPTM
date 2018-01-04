library(IPTM)
load('Enron.RData')
attach(Enron)

set.seed(1)
missing = matrix(0, nrow = length(Enron$edge), ncol = 3)
missing[sample(556:6613, 1212, replace = FALSE), 1] = 1
missing[sample(556:6613, 1212, replace = FALSE), 2] = 1
missing[sample(556:6613, 1212, replace = FALSE), 3] = 1

for (nIP in 1:2) {
for (K in c(5, 10, 25, 50, 75, 100)) {
EnronPPE = list()
for (i in 1:5){
  filename = paste0("Enron_full_",nIP,"_",K,"_ver",i,".RData")
  load(filename)
  initial = list()
  initial$alpha =  Enrontest$alpha[length(Enrontest$alpha)[1]]
  initial$mvec = Enrontest$mvec[dim(Enrontest$mvec)[1],]
  initial$delta = Enrontest$delta[length(Enrontest$delta)]
  initial$sigma_tau = Enrontest$sigma_tau[length(Enrontest$sigma_tau)]
  initial$b = t(vapply(1:nIP, function(IP) {Enrontest$b[[IP]][,ncol(Enrontest$b[[IP]])]}, rep(0, nrow(Enrontest$b[[1]]))))
  initial$eta = t(vapply(1:nIP, function(IP) {Enrontest$eta[[IP]][,ncol(Enrontest$eta[[IP]])]}, rep(0, nrow(Enrontest$eta[[1]]))))
  initial$l = Enrontest$l
  initial$z = Enrontest$z
  initial$u = Enrontest$u
  EnronPPE[[i]] = IPTM.PPE(missing, edge = Enron$edge, node = Enron$node, textlist = Enron$text, vocab= Enron$vocab, nIP = nIP, K = K,
                     Outer = 20, netstat = c("dyadic", "degree", "triadic"), timestat = c("dayofweek", "timeofday"), initial = initial)
}
filename = paste0("Enron_PPE_",nIP,"_",K,".RData")
save(EnronPPE, file = filename)
}
}


