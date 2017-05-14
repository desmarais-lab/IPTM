library(IPTM)
TryGiR<- GiR(5*10^3, nDocs = 10, niters = c(1, 3, 3, 150, 0, 3), sigma_Q = 0.25, seed = 1)
TryGiR2<- GiR(5*10^3, nDocs = 20, niters = c(1, 3, 3, 150, 0, 3), sigma_Q = 0.25, seed = 1)

TryGiR = TryGiR2
Nsamp = nrow(TryGiR$Forward)
thin = seq(from = floor(Nsamp / 5), to = Nsamp, length.out = 500)
par(mfrow=c(2,5))
matplot(TryGiR$delta[thin, ],type = 'l', col = 1:2, lty = 1)
matplot(cbind(TryGiR$entropyC1[thin,1],TryGiR$entropyC2[thin,1]) ,type = 'l', col = 1:2, lty =1)
matplot(cbind(TryGiR$entropyC1[thin,2],TryGiR$entropyC2[thin,2]) ,type = 'l', col = 1:2, lty =1)
matplot(cbind(TryGiR$entropyZ1[thin],TryGiR$entropyZ2[thin]) ,type = 'l', col = 1:2, lty =1)
matplot(cbind(TryGiR$Zstat1[thin],TryGiR$Zstat2[thin]) ,type = 'l', col = 1:2, lty =1)


matplot(TryGiR$delta[thin, ],type = 'l', col = 2:1, lty = 1)
matplot(cbind(TryGiR$entropyC2[thin,1],TryGiR$entropyC1[thin,1]) ,type = 'l', col = 1:2, lty =1)
matplot(cbind(TryGiR$entropyC2[thin,2],TryGiR$entropyC1[thin,2]) ,type = 'l', col = 1:2, lty =1)
matplot(cbind(TryGiR$entropyZ2[thin],TryGiR$entropyZ1[thin]) ,type = 'l', col = 1:2, lty =1)
matplot(cbind(TryGiR$Zstat2[thin],TryGiR$Zstat1[thin]) ,type = 'l', col = 1:2, lty =1)


for (p in 1:14) {
	matplot(cbind(TryGiR$b1[thin ,p], TryGiR$b2[thin,p]),type = 'l', col = 1:2, lty = 1)
}


TryGiR = TryGiR2
Nsamp = nrow(TryGiR$Forward)
thin = seq(from = floor(Nsamp / 5), to = Nsamp, length.out = 500)
par(mfrow = c(2, 7))
for (p in 1:ncol(TryGiR$Forward)){
matplot(cbind(TryGiR$Backward[thin,p], TryGiR$Forward[thin,p]), type = 'l', col = 1:2, lty = 1, main = colnames(TryGiR$Forward)[p], xlab = 'iter', ylab ='')
}