load("Daresurv.RData")
attach(Daresurv)	

#get initial value using normal dist.
agam0 <- lm(Y~X1+X2+X3+X4+X5+X6+X7+X8+X9+X10+X11+X12+X13+X14+X15+X16+X17+X18+X19+X20+X21+X22+X23+X24+D+W)


#first try lognormal
norm <- glm(log(Y)~X1+X2+X3+X4+X5+X6+X7+X8+X9+X10+X11+X12+X13+X14+X15+X16+X17+X18+X19+X20+X21+X22+X23+X24+D+W, family = gaussian(link = "identity"),control = glm.control(maxit = 100), start = as.numeric(coef(agam0)))
normpred <- simulate(norm, type = "response", se =T)
par(mfrow = c(2,2))
plot(norm) 

#try Gamma
agam <- glm(Y~X1+X2+X3+X4+X5+X6+X7+X8+X9+X10+X11+X12+X13+X14+X15+X16+X17+X18+X19+X20+X21+X22+X23+X24+D+W, family = Gamma(link = "identity"),control = glm.control(maxit = 100), start = c(rep(1, 25), coef(agam0)[-1:-25]))
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
     norm <- glm(log(Y)~X1+X2+X3+X4+X5+X6+X7+X8+X9+X10+X11+X12+X13+X14+X15+X16+X17+X18+X19+X20+X21+X22+X23+X24+D+W, family = gaussian(link = "identity"),control = glm.control(maxit = 100), start = as.numeric(coef(agam0)))
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
    
agam <- glm(Y~X1+X2+X3+X4+X5+X6+X7+X8+X9+X10+X11+X12+X13+X14+X15+X16+X17+X18+X19+X20+X21+X22+X23+X24+D+W, family = Gamma(link = "identity"),control = glm.control(maxit = 100), start = c(rep(1, 25), coef(agam0)[-1:-25]))
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

hist(unlist(exp(normpred)), breaks = 50, freq = FALSE, main = " ")
lines(density(Daresurv$Y), col = 'red')
hist(unlist(gampred), breaks = 50,freq = FALSE, main = "")
lines(density(Daresurv$Y), col = 'red')

#How can we parameterize alpha and beta as Aaron's paper?

Gamma = function(link = "inverse") {
	 linktemp <- substitute(link)
    if (!is.character(linktemp)) 
        linktemp <- deparse(linktemp)
    okLinks <- c("inverse", "log", "identity")
    if (linktemp %in% okLinks) 
        stats <- make.link(linktemp)
    else if (is.character(link)) 
        stats <- make.link(link)
    else {
        if (inherits(link, "link-glm")) {
            stats <- link
            if (!is.null(stats$name)) 
                linktemp <- stats$name
        }
        else {
            stop(gettextf("link \"%s\" not available for gamma family; available links are %s", 
                linktemp, paste(sQuote(okLinks), collapse = ", ")), 
                domain = NA)
        }
    }
    variance <- function(mu) mu^2
    validmu <- function(mu) all(is.finite(mu)) && all(mu > 0)
    dev.resids <- function(y, mu, wt) -2 * wt * (log(ifelse(y == 
        0, 1, y/mu)) - (y - mu)/mu)
    aic <- function(y, n, mu, wt, dev) {
        n <- sum(wt)
        disp <- dev / n                                      # probably change this part... but nothing worked
        -2 * sum(dgamma(y, 1/disp, scale = mu * disp, log = TRUE) * 
            wt) + 2
    }
    initialize <- expression({
        if (any(y <= 0)) stop("non-positive values not allowed for the 'gamma' family")
        n <- rep.int(1, nobs)
        mustart <- y
    })
    simfun <- function(object, nsim) {
        wts <- object$prior.weights
        if (any(wts != 1)) 
            message("using weights as shape parameters")
        ftd <- fitted(object)
        shape <- MASS::gamma.shape(object)$alpha * wts
        rgamma(nsim * length(ftd), shape = shape, rate = shape/ftd)
    }
    structure(list(family = "Gamma", link = linktemp, linkfun = stats$linkfun, 
        linkinv = stats$linkinv, variance = variance, dev.resids = dev.resids, 
        aic = aic, mu.eta = stats$mu.eta, initialize = initialize, 
        validmu = validmu, valideta = stats$valideta, simulate = simfun), 
        class = "family")
}
