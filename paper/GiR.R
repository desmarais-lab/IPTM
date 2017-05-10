library(IPTM)
TryGiR <- GiR(10^4, nDocs = 10, seed = 12)
TryGiR2 <- GiR(10^4, nDocs = 10, sigma_Q = 0.5, seed = 12)

# Nsamp = nrow(TryGiR$Forward)
# thin = seq(from = floor(Nsamp / 5), to = Nsamp, length.out = 500)
# par(mfrow=c(4, 5))
# matplot(TryGiR$delta[thin, ],type = 'l', col = 1:2, lty = 1)
# matplot(cbind(TryGiR$entropyC1[thin,1],TryGiR$entropyC2[thin,1]) ,type = 'l', col = 1:2, lty =1)
# matplot(cbind(TryGiR$entropyC1[thin,2],TryGiR$entropyC2[thin,2]) ,type = 'l', col = 1:2, lty =1)
# matplot(cbind(TryGiR$entropyZ1[thin],TryGiR$entropyZ2[thin]) ,type = 'l', col = 1:2, lty =1)
# for (p in 1:14) {
	# matplot(cbind(TryGiR$b1[thin ,p], TryGiR$b2[thin,p]),type = 'l', col = 1:2, lty = 1)
# }

Nsamp = nrow(TryGiR$Forward)
thin = seq(from = floor(Nsamp / 5), to = Nsamp, length.out = 500)
par(mfrow = c(3, 7))
for (p in 1:ncol(TryGiR$Forward)){
matplot(cbind(TryGiR$Forward[thin,p], TryGiR$Backward[thin,p]), type = 'l', col = 2:1, lty = 1, main = colnames(TryGiR$Forward)[p], xlab = 'iter', ylab ='')
}