library(IPTM)
load('Enron.RData')
attach(Enron)
missing = list()
for (i in 1:5){
    set.seed(i)
    missing[[i]] = matrix(0, nrow = length(Enron$edge), ncol = 3)
    #missing[[i]][sample(556:6613, 606, replace = FALSE), 1] = 1
    #missing[[i]][sample(556:6613, 606, replace = FALSE), 2] = 1
    missing[[i]][sample(556:6613, 606, replace = FALSE), 3] = 1
}


nIP = 1
K = 75
EnronPPE = list()
for (i in 1:5){
  filename = paste0("Enron_full_",nIP,"_",K,"_ver",1,".RData")
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
  initial$proposal.var1 = diag(nrow(Enrontest$b[[1]]))
  initial$proposal.var2 = diag(nrow(Enrontest$eta[[1]]))
  initial$sigma.Q = Enrontest$sigma.Q
  initial$u = Enrontest$u
  EnronPPE[[i]] = IPTM.inference.PPE(missing[[i]], edge = Enron$edge, node = Enron$node, textlist = Enron$text, vocab= Enron$vocab, nIP = nIP, K = K,
  sigma.Q = c(0.0001, 0.001, 5, 0.1), alpha = 2, mvec = rep(1/K, K), beta = 2, prior.b = list(rep(0, 24), 10*diag(24)),
  prior.delta = c(-50, 5),  prior.eta = list(c(rep(30, length(Enron$node)), rep(0,2)), 10*diag(length(Enron$node)+2)), prior.tau = 5, Outer = 30, Inner = c(2200,2200, 550), burn = c(200,200, 50), thin = c(2,2,1),
  netstat = c("dyadic", "degree", "triadic"), timestat = c("dayofweek", "timeofday"), optimize = TRUE, initial = initial, timeunit = (3600*24), tz = "America/Los_Angeles")
}
filename = paste0("Enron_PPE_",nIP,"_",K,"time.RData")
save(EnronPPE, file = filename)


