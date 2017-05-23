library(IPTM)
library(mvtnorm)
library(MCMCpack)
set.seed(10)
nDocs = 10
node = 1:4
vocabulary = c("hi", "hello", "fine", "bye", "what")
nIP = 2
K = 4
nwords = 4
alpha = 2
mvec = rep(1/4, 4)
betas = 2
nvec = rep(1/5, 5)
prior.b.mean = c(-2.5, rep(0, 6))
prior.b.var =  diag(7)
prior.delta = c(4, 8)
sigma_Q = c(1, 2.5)
niters = c(1, 100, 20, 0, 5)
netstat = c("intercept", "dyadic")
P = 1 * ("intercept" %in% netstat) + 3 * (2 * ("dyadic" %in% netstat) + 4 * ("triadic" %in% netstat) + 2 *("degree" %in% netstat))
b = lapply(1:nIP, function(IP) {
    c(rmvnorm(1,  prior.b.mean, prior.b.var))
  })
delta = rgamma(1, 4, 8)
currentC = sample(1:nIP, K, replace = TRUE)	 
base.data = GenerateDocs(30, node, vocabulary, nIP, K, nwords, alpha, mvec, betas, nvec, b, delta, currentC, netstat, base.edge = list(),  base.text = list(), base = TRUE) 
base.edge = base.data$edge	   
base.text = base.data$text
  
  
TryGiR2<- GiR(5*10^4, nDocs, node, vocabulary, nIP, K, nwords, alpha, mvec, betas, nvec, 
               prior.b.mean, prior.b.var, prior.delta, sigma_Q, niters, netstat, base.edge, base.text, seed = 12345)
save(TryGiR, file = "TryGiR.RData")
par(mfrow=c(3,7))
GiR_PP_Plots(TryGiR$Forward, TryGiR$Backward)

Nsamp = nrow(TryGiR$Forward)
thin = seq(from = floor(Nsamp / 5), to = Nsamp, length.out = 500)
par(mfrow=c(2,2))
matplot(TryGiR$delta[thin, ],type = 'l', col = 1:2, lty = 1, xlab = "iter", ylab = "delta")
matplot(cbind(TryGiR$Zstat1[thin],TryGiR$Zstat2[thin]) ,type = 'l', col = 1:2, lty =1, xlab = "iter", ylab = "mean(entropy(topic-token[[d]]))")


matplot(TryGiR$delta[thin, ],type = 'l', col = 2:1, lty = 1, xlab = "iter", ylab = "delta")
matplot(cbind(TryGiR$Zstat2[thin],TryGiR$Zstat1[thin]) ,type = 'l', col = 2:1, lty =1, xlab = "iter", ylab = "mean(entropy(topic-token[[d]]))")


Nsamp = nrow(TryGiR$Forward)
thin = seq(from = floor(Nsamp / 5), to = Nsamp, length.out = 500)
par(mfrow = c(3, 7))
for (p in 1:21){
matplot(cbind(TryGiR$Forward[thin,p], TryGiR$Backward[thin,p]), type = 'l', col = 1:2, lty = 1, main = colnames(TryGiR$Forward)[p], xlab = 'iter', ylab ='')
}
par(mfrow = c(3, 7))
for (p in 1:21){
matplot(cbind(TryGiR$Backward[thin,p], TryGiR$Forward[thin,p]), type = 'l', col = 2:1, lty = 1, main = colnames(TryGiR$Forward)[p], xlab = 'iter', ylab ='')
}