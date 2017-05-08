library(IPTM)
TryGiR <- GiR(10^5, niters = c(1, 1, 1, 330, 30, 3), seed = 123)

Nsamp = nrow(TryGiR$Forward)
thin = seq(from = floor(Nsamp / 5), to = Nsamp, length.out = 500)
par(mfrow=c(4, 4))
matplot(TryGiR$delta[thin, ],type = 'l', col = 1:2, lty = 1)
matplot(TryGiR$entropy[thin, ],type = 'l', col = 1:2, lty =1)
for (p in 1:14) {
	matplot(cbind(TryGiR$b1[thin ,p], TryGiR$b2[thin,p]),type = 'l', col = 1:2, lty = 1)
}