library(IPTM)
library(mvtnorm)
library(MCMCpack)
set.seed(100)
nDocs = 5
node = 1:4
vocabulary = c("hi", "hello", "fine", "bye", "what")
nIP = 2
K = 4
nwords = 4
alpha = 2
mvec = rep(1/4, 4)
betas = 2
nvec = rep(1/5, 5)
netstat = c("intercept", "dyadic")
P = 1 * ("intercept" %in% netstat) + 3 * (2 * ("dyadic" %in% netstat) + 4 * ("triadic" %in% netstat) + 2 *("degree" %in% netstat))
prior.b.mean = c(-3, rep(0, P-1))
prior.b.var = 0.1 * diag(P)
prior.delta = c(2.5, 0.1)
sigma_Q = c(1.25, 2)
niters = c(3, 440, 200, 40, 2)
b = lapply(1:nIP, function(IP) {
    c(rmvnorm(1,  prior.b.mean, prior.b.var))
  })
delta = rnorm(1, prior.delta[1], sqrt(prior.delta[2]))
currentC = sample(1:nIP, K, replace = TRUE)	 
supportD = gibbs.measure.support(length(node) - 1)
base.data = GenerateDocs.Gibbs(100, node, vocabulary, nIP, K, nwords, alpha, mvec, betas, nvec, b, delta, currentC, netstat, base.edge = list(),  base.text = list(), base = TRUE, support = supportD) 
base.edge = base.data$edge	   
base.text = base.data$text
TryGiR2<- GiR.Gibbs(100, nDocs, node, vocabulary, nIP, K, nwords, alpha, mvec, betas, nvec, 
					prior.b.mean, prior.b.var, prior.delta, sigma_Q, niters, netstat, base.edge, base.text, seed = 1)

save(TryGiR, file = "/Users/bomin8319/Desktop/IPTM-master/paper/TryGiR.RData")

load("TryGiR.RData")
load("TryGiR2.RData")
par(mfrow = c(1,2))
qqplot(TryGiR$[,1], TryGiR$delta[,2])
abline(0, 1, col = 'red')
qqplot(TryGiR2$delta[,1], TryGiR2$delta[,2])
abline(0, 1, col = 'red')

par(mfrow=c(5,5), oma = c(3,3,3,3), mar = c(2,1,1,1))
GiR_PP_Plots(TryGiR$Forward, TryGiR$Backward)

TryGiR = TryGiR2
Nsamp = nrow(TryGiR$Forward)
thin = seq(from = floor(Nsamp / 10), to = Nsamp, length.out = 1000)
par(mfrow=c(2,2))
matplot(TryGiR$delta[thin, ],type = 'l', col = 1:2, lty = 1, xlab = "iter", ylab = "delta")
matplot(cbind(TryGiR$Zstat1[thin],TryGiR$Zstat2[thin]) ,type = 'l', col = 1:2, lty =1, xlab = "iter", ylab = "mean(entropy(topic-token[[d]]))")


matplot(TryGiR$delta[thin, ],type = 'l', col = 2:1, lty = 1, xlab = "iter", ylab = "delta")
matplot(cbind(TryGiR$Zstat2[thin],TryGiR$Zstat1[thin]) ,type = 'l', col = 2:1, lty =1, xlab = "iter", ylab = "mean(entropy(topic-token[[d]]))")


Nsamp = nrow(TryGiR$Forward)
thin = seq(from = floor(Nsamp / 10), to = Nsamp, length.out = 500)
par(mfrow = c(3, 7))
for (p in 1:21){
matplot(cbind(TryGiR$Forward[thin,p], TryGiR$Backward[thin,p]), type = 'l', col = 1:2, lty = 1, main = colnames(TryGiR$Forward)[p], xlab = 'iter', ylab ='')
}
par(mfrow = c(3, 7))
for (p in 1:21){
matplot(cbind(TryGiR$Backward[thin,p], TryGiR$Forward[thin,p]), type = 'l', col = 2:1, lty = 1, main = colnames(TryGiR$Forward)[p], xlab = 'iter', ylab ='')
}