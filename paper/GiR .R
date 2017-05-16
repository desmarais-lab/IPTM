library(IPTM)
TryGiR<- GiR(5000, nDocs = 10, node = 1:4, prior.delta = c(4, 8), niters = c(1, 1, 1, 50, 0, 1), prior.b.mean = c(-2, rep(0, 6)), sigma_Q = c(0.25, 2), seed = 12345)

Nsamp = nrow(TryGiR$Forward)
thin = seq(from = floor(Nsamp / 5), to = Nsamp, length.out = 500)
par(mfrow=c(1,5))
matplot(TryGiR$delta[thin, ],type = 'l', col = 1:2, lty = 1, xlab = "iter", ylab = "delta")
matplot(cbind(TryGiR$entropyC1[thin,1],TryGiR$entropyC2[thin,1]) ,type = 'l', col = 1:2, lty =1, xlab = "iter", ylab = "entropy(topic-IP)")
matplot(cbind(TryGiR$entropyC1[thin,2],TryGiR$entropyC2[thin,2]) ,type = 'l', col = 1:2, lty =1, xlab = "iter", ylab = "entropy(topic-token overall)")
matplot(cbind(TryGiR$entropyZ1[thin],TryGiR$entropyZ2[thin]) ,type = 'l', col = 1:2, lty =1, xlab = "iter", ylab = "entropy(word-token)")
matplot(cbind(TryGiR$Zstat1[thin],TryGiR$Zstat2[thin]) ,type = 'l', col = 1:2, lty =1, xlab = "iter", ylab = "mean(entropy(topic-token[[d]]))")


matplot(TryGiR$delta[thin, ],type = 'l', col = 2:1, lty = 1, xlab = "iter", ylab = "delta")
matplot(cbind(TryGiR$entropyC2[thin,1],TryGiR$entropyC1[thin,1]) ,type = 'l', col = 2:1, lty =1, xlab = "iter", ylab = "entropy(topic-IP)")
matplot(cbind(TryGiR$entropyC2[thin,2],TryGiR$entropyC1[thin,2]) ,type = 'l', col = 2:1, lty =1, xlab = "iter", ylab = "entropy(topic-token overall)")
matplot(cbind(TryGiR$entropyZ2[thin],TryGiR$entropyZ1[thin]) ,type = 'l', col =  2:1, lty =1, xlab = "iter", ylab = "entropy(word-token)")
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