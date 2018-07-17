load("Daresurv.RData")
attach(Daresurv)	

#get initial value using normal dist.
agam0 <- lm(Y~sender+D+W)


#first try lognormal
norm <- glm(log(Y)~sender+D+W, family = gaussian(link = "identity"),control = glm.control(maxit = 100), start = as.numeric(coef(agam0)))
normpred <- simulate(norm, type = "response", se =T)
par(mfrow = c(2,2))
plot(norm) 


norm <- glm(log(Y)~sender+D+W, family = gaussian(link = "identity"),control = glm.control(maxit = 100), start = as.numeric(coef(agam0)))
normpred <- simulate(norm, type = "response", se =T)
par(mfrow = c(2,2))
plot(norm) 

norm <- glm(log(Y)~sender+D+W, family = inverse.gaussian(link = "log"),control = glm.control(maxit = 100), start = as.numeric(coef(agam0)))
normpred <- simulate(norm, type = "response", se =T)
par(mfrow = c(2,2))
plot(norm) 




#try Gamma
agam <- glm(Y~sender+D+W, family = Gamma(link = "identity"),control = glm.control(maxit = 100), start = c(rep(1, 25), coef(agam0)[-1:-25]))
summary(agam)
par(mfrow = c(2,2))
plot(agam)

library(MASS)
myshape<- gamma.shape(agam)
summary(agam, dispersion = 1 / myshape$alpha)     # different dispersion parameter changes the significance level
gampred <- simulate(agam, type = "response", se =T, dispersion = 1/myshape$alpha)


par(mfrow = c(1,2))
qqplot(exp(normpred$sim_1), Daresurv$Y)
qqplot(gampred$sim_1, Daresurv$Y)  # maximum time generated is 3.1354



#Gamma has convergence issue -> how about trying W and D only?

agam <- glm(Y~D+W, family = Gamma(link = "identity"),control = glm.control(maxit = 100), start =rep(0.1, 10))
summary(agam)
par(mfrow = c(2,2))
plot(agam)

library(MASS)
myshape<- gamma.shape(agam)
summary(agam, dispersion = 1 / myshape$alpha)     # different dispersion parameter changes the significance level
gampred <- simulate(agam, type = "response", se =T, dispersion = 1/myshape$alpha)
qqplot(gampred$sim_1, Daresurv$Y)  # maximum time generated is 2.2576





par(mar=c(2, 3, 1, 1), mfrow=c(3,2),
     oma = c(1, 1, 0.5, 0.3))
     norm <- glm(log(Y)~sender+D+W, family = gaussian(link = "identity"),control = glm.control(maxit = 100), start = as.numeric(coef(agam0)))
normpred <- simulate(norm, nsim = 1000)
sim.results = unlist(exp(normpred))

all= c(sim.results, Daresurv$Y)
 uniqueValues = quantile(all,seq(0, 1, length = 100))
    qx1 = numeric(length(uniqueValues))
  	qx2 = numeric(length(uniqueValues))
	for (j in 1:length(uniqueValues)) {
  		qx1[j] = mean(sim.results <= uniqueValues[j])
  		qx2[j] = mean(Daresurv$Y<= uniqueValues[j])
  	}
qqplot(x = qx1,
           y = qx2,
           ylim = c(0, 1),
           xlim = c(0, 1),
           ylab = "Observed",
           xlab = "Simulated",
           col = "blue",
           pch = 19,
           cex = 0.25,
           main = "Lognormal",
           cex.lab = 0.25,
           cex.axis = 0.25,
           cex.main = 0.5)
    abline(0, 1, lty = 1, col = "red", lwd = 1)
    
agam <- glm(Y~sender+D+W, family = Gamma(link = "identity"),control = glm.control(maxit = 100), start = c(rep(1, 25), coef(agam0)[-1:-25]))
gampred <- simulate(agam, nsim = 1000)
sim.results = unlist(gampred)

all= c(sim.results, Daresurv$Y)
 uniqueValues = quantile(all,seq(0, 1, length = 100))
    qx1 = numeric(length(uniqueValues))
  	qx2 = numeric(length(uniqueValues))
	for (j in 1:length(uniqueValues)) {
  		qx1[j] = mean(sim.results <= uniqueValues[j])
  		qx2[j] = mean(Daresurv$Y<= uniqueValues[j])
  	}
qqplot(x = qx1,
           y = qx2,
           ylim = c(0, 1),
           xlim = c(0, 1),
           ylab = "Observed",
           xlab = "Simulated",
           col = "blue",
           pch = 19,
           cex = 0.25,
           main = "Gamma",
           cex.lab = 0.25,
           cex.axis = 0.25,
           cex.main = 0.5)
    abline(0, 1, lty = 1, col = "red", lwd = 1)

qqplot(unlist(exp(normpred)), Daresurv$Y)
abline(0,1, col= 'red')
qqplot(unlist(gampred), Daresurv$Y,)
abline(0,1, col = 'red')    

qqplot(unlist(exp(normpred)), Daresurv$Y, xlim = c(0, 70))
abline(0,1, col= 'red')
qqplot(unlist(gampred), Daresurv$Y, xlim = c(0, 70))
abline(0,1, col = 'red')    

####
hist(unlist(exp(normpred)), breaks = 50, freq = FALSE, main = " ")
lines(density(Daresurv$Y), col = 'red')
hist(unlist(gampred), breaks = 50,freq = FALSE, main = "")
lines(density(Daresurv$Y), col = 'red')
