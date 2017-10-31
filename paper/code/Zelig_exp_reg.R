library(survival)
library(Zelig)


setwd("/Users/bomin8319/Desktop/IPTM/paper/code/predictiveComparisonModels/results/")
# load Dare observed data
load("Dare.RData")
# load Dare network statistics
load("netstats_Dare.RData")


# extract the number of nodes
nodes <- dim(netstats)[2]

# extract number of features
features <- dim(netstats)[4]

# build dataset for training
# sender, receivers (nodes columns), timeinc
observed_edge_data <- matrix(NA,length(Dare$edge),nodes+3)
for(d in 1:length(Dare$edge)){
  dati <- c(Dare$edge[[d]]$sender)
  receivers <- numeric(nodes)
  receivers[Dare$edge[[d]]$receiver] <- 1
  dati <- c(dati,receivers,Dare$edge[[d]]$timeinc,Dare$edge[[d]]$unixtime)
  observed_edge_data[d,] <- dati
}
colnames(observed_edge_data) <- c("sender",paste("R",1:nodes,sep=""),"timeinc","time")

observed_edge_data <- data.frame(observed_edge_data)

# observed_edge_data <- observed_edge_data[observed_edge_data$time > 384,]

first.training.obs <- min(which(observed_edge_data$time > 384))




Y <- observed_edge_data[first.training.obs:nrow(observed_edge_data), 29]

load("/Users/bomin8319/Desktop/IPTM/paper/Daretest.RData")
B <- lapply(Daretest$B, function(x) {rowMeans(x)})
C <- Daretest$C
Z <- Daretest$Z
library(IPTM2)
setwd('/Users/bomin8319/Desktop/IPTM/pkg2')
library(devtools)
document()
timestamps = vapply(Dare$edge, function(d) {
  d[[3]]
}, c(1))

edge2 = which_int(384, timestamps) : length(Dare$edge)
p.d = pdmat(Z, C, 2)
netstat = c(0, 1,1,1)
L = 3
P = L * (2 * netstat[2] + 2 * netstat[3] + 4 * netstat[4])
lambda = list()
X = list()
for (d in edge2) {
  history.t = History(Dare$edge, p.d, Dare$node, Dare$edge[[d-1]][[3]] + exp(-745))
  X[[d]] = Netstats_cpp(history.t, Dare$node, netstat)
  XB = MultiplyXBList(X[[d]], B)     
  lambda[[d]] = lambda_cpp(p.d[d,], XB)
}

library(anytime)
library(lubridate)
load("/Users/bomin8319/Desktop/IPTM/full/Darenew.RData")
wday = function(x) {
  as.POSIXlt(x)$wday
}
weeks = sapply(edge2, function(x) wday(anytime(Dare$edge[[x]]$unixtime,  tz = "America/New_York")) )
day1 = sapply(edge2, function(x) hour(ymd_hms(anytime(Dare$edge[[x]]$unixtime, tz = "America/New_York"))))
day = cut(day1, c(0,6,12,18,24), c("0-6", "6-12", "12-18", "18-24"))
day[day1==0] = "0-6"

sender = matrix(0, length(edge2), 1)
lambdaij = matrix(0, length(edge2), 1)
Xij = matrix(0, length(edge2), 24)
iter = 1

for (d in edge2) {
	sender[iter, ] = Dare$edge[[d]]$sender
 lambdaij[iter, ] = mean(lambda[[d]][Dare$edge[[d]]$sender, Dare$edge[[d]]$receiver])
 Xs = p.d[d, 1] * X[[d]][[Dare$edge[[d]]$sender]][[1]][Dare$edge[[d]]$receiver,] + p.d[d, 2] * X[[d]][[Dare$edge[[d]]$sender]][[2]][Dare$edge[[d]]$receiver,]
 if (is.matrix(Xs)) {Xs = colMeans(Xs)}
 Xij[iter, ] = Xs
   iter = iter + 1
}
Xij = Xij[-which(Y == 0),]
sender = sender[-which(Y == 0),]
Daresurv = cbind(Daresurv, Xij, sender)
names(Daresurv)[5:28] = sapply(1:24, function(x) paste0("X",x))
#save(Daresurv, file = "Daresurv.RData")
#Daresurv = data.frame(Y = Y, L = lambdaij, W = as.factor(weeks), D = as.factor(day))
#Daresurv = Daresurv[-which(Daresurv$Y == 0),]
Daresurv$logL = log(Daresurv$L)

par(mfrow = c(1,2))
z.out <- survreg(Surv(Y)~ as.factor(sender)+logL+W + D,dist= "weibull", data = Daresurv)

quantiles <- seq(0.0001,.9999,length=1000)
quantileMat <- predict(z.out,type="quantile",p=quantiles)
sim.quantiles <- sample(1:1000,nrow(Daresurv),rep=T)
simulated.times <- quantileMat[cbind(1:nrow(Daresurv),sim.quantiles)]
qqplot(simulated.times, Daresurv$Y, main = "logL+W+D")
qqline(simulated.times)

z.out <- survreg(Surv(Y)~ logL,dist= "weibull", data = Daresurv)

quantiles <- seq(0.0001,.9999,length=1000)
quantileMat <- predict(z.out,type="quantile",p=quantiles)
sim.quantiles <- sample(1:1000,nrow(Daresurv),rep=T)
simulated.times <- quantileMat[cbind(1:nrow(Daresurv),sim.quantiles)]
qqplot(simulated.times, Daresurv$Y, main = "offset(logL)+W+D")
qqline(simulated.times)



par(mfrow = c(2,2))
z.out <- survreg(Surv(Y)~ -1 +sender+ W+ D,dist= "exponential", data = Daresurv)

quantiles <- seq(0.0001,.9999,length=1000)
quantileMat <- predict(z.out,type="quantile",p=quantiles)

intervals = c()
intervals[1] = mean(quantileMat[,450] < Daresurv$Y & quantileMat[,550] > Daresurv$Y)
intervals[2] = mean(quantileMat[,400] < Daresurv$Y & quantileMat[,600] > Daresurv$Y)
intervals[3] = mean(quantileMat[,350] < Daresurv$Y & quantileMat[,650] > Daresurv$Y)
intervals[4] = mean(quantileMat[,300] < Daresurv$Y & quantileMat[,700] > Daresurv$Y)
intervals[5] = mean(quantileMat[,250] < Daresurv$Y & quantileMat[,750] > Daresurv$Y)
intervals[6] = mean(quantileMat[,200] < Daresurv$Y & quantileMat[,800] > Daresurv$Y)
intervals[7] = mean(quantileMat[,150] < Daresurv$Y & quantileMat[,850] > Daresurv$Y)
intervals[8] = mean(quantileMat[,100] < Daresurv$Y & quantileMat[,900] > Daresurv$Y)
intervals[9] = mean(quantileMat[,50] < Daresurv$Y & quantileMat[,950] > Daresurv$Y)
intervals[10] = mean(quantileMat[,1] < Daresurv$Y & quantileMat[,1000] > Daresurv$Y)
plot(c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.00), intervals, xlab = "interval", ylab ="% within", main = "Exponential")
abline(0,1)
sim.results = matrix(0, nrow(Daresurv), 1000)
for (k in 1:1000) {
sim.quantiles <- sample(1:1000,nrow(Daresurv),rep=T)
simulated.times <- quantileMat[cbind(1:nrow(Daresurv),sim.quantiles)]
sim.results[,k ] = simulated.times
}

all= c(rowMeans(sim.results), Daresurv$Y)
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
           main = "sender+X+W+D (Exponential)",
           cex.lab = 0.25,
           cex.axis = 0.25,
           cex.main = 0.5)
    abline(0, 1, lty = 1, col = "red", lwd = 1)
    
z.out <- survreg(Surv(Y)~ X1+X2+X3+X4+X5+X6+X7+X8+X9+X10+X11+X12+X13+X14+X15+X16+X17+X18+X19+X20+X21+X22+X23+X24+W + D,dist= "lognormal", data = Daresurv)

quantiles <- seq(0.0001,.9999,length=1000)
quantileMat <- predict(z.out,type="quantile",p=quantiles)
sim.results = matrix(0, nrow(Daresurv), 1000)
for (k in 1:1000) {
sim.quantiles <- sample(1:1000,nrow(Daresurv),rep=T)
simulated.times <- quantileMat[cbind(1:nrow(Daresurv),sim.quantiles)]
sim.results[,k ] = simulated.times
}


all= c(rowMeans(sim.results), Daresurv$Y)
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
           main = "sender+X+W+D",
           cex.lab = 0.25,
           cex.axis = 0.25,
           cex.main = 0.5)
    abline(0, 1, lty = 1, col = "red", lwd = 1)
    

z.out <- survreg(Surv(Y)~ -1 + log(L)+W+D,dist= "lognormal", data = Daresurv)

quantiles <- seq(0.0001,.9999,length=1000)
quantileMat <- predict(z.out,type="quantile",p=quantiles)
sim.results = matrix(0, nrow(Daresurv), 1000)
for (k in 1:1000) {
sim.quantiles <- sample(1:1000,nrow(Daresurv),rep=T)
simulated.times <- quantileMat[cbind(1:nrow(Daresurv),sim.quantiles)]
sim.results[,k ] = simulated.times
}
intervals = c()
intervals[1] = mean(quantileMat[,450] < Daresurv$Y & quantileMat[,550] > Daresurv$Y)
intervals[2] = mean(quantileMat[,400] < Daresurv$Y & quantileMat[,600] > Daresurv$Y)
intervals[3] = mean(quantileMat[,350] < Daresurv$Y & quantileMat[,650] > Daresurv$Y)
intervals[4] = mean(quantileMat[,300] < Daresurv$Y & quantileMat[,700] > Daresurv$Y)
intervals[5] = mean(quantileMat[,250] < Daresurv$Y & quantileMat[,750] > Daresurv$Y)
intervals[6] = mean(quantileMat[,200] < Daresurv$Y & quantileMat[,800] > Daresurv$Y)
intervals[7] = mean(quantileMat[,150] < Daresurv$Y & quantileMat[,850] > Daresurv$Y)
intervals[8] = mean(quantileMat[,100] < Daresurv$Y & quantileMat[,900] > Daresurv$Y)
intervals[9] = mean(quantileMat[,50] < Daresurv$Y & quantileMat[,950] > Daresurv$Y)
intervals[10] = mean(quantileMat[,1] < Daresurv$Y & quantileMat[,1000] > Daresurv$Y)
plot(c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.00), intervals, xlab = "interval", ylab ="% within", main = "Lognormal")
abline(0,1)

all= c(rowMeans(sim.results), Daresurv$Y)
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
           main = "sender+X+W+D (Lognormal)",
           cex.lab = 0.25,
           cex.axis = 0.25,
           cex.main = 0.5)
    abline(0, 1, lty = 1, col = "red", lwd = 1)
    
z.out <- survreg(Surv(Y)~ -1+X1+X2+X3+X4+X5+X6+X7+X8+X9+X10+X11+X12+X13+X14+X15+X16+X17+X18+X19+X20+X21+X22+X23+X24+sender+W + D,dist= "extreme", data = Daresurv)

quantiles <- seq(0.0001,.9999,length=1000)
quantileMat <- predict(z.out,type="quantile",p=quantiles)
sim.results = matrix(0, nrow(Daresurv), 100)
for (k in 1:100) {
sim.quantiles <- sample(1:1000,nrow(Daresurv),rep=T)
simulated.times <- quantileMat[cbind(1:nrow(Daresurv),sim.quantiles)]
sim.results[,k ] = simulated.times
}

all= c(rowMeans(sim.results), Daresurv$Y)
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
           main = "sender+X+W+D",
           cex.lab = 0.25,
           cex.axis = 0.25,
           cex.main = 0.5)
    abline(0, 1, lty = 1, col = "red", lwd = 1)
    


myloggamma = list(
name = "Loggamma",
init = function (x, weights, ...) 
{
    mean <- sum(x * weights)/sum(weights)
    var <- sum(weights * (x - mean)^2)/sum(weights)
    c(mean, var)
},
density = function (x, parms) 
{	expx = exp(x)
    k <- 2.5
    cbind(pgamma(expx, shape = k, scale = 1), 1 - pgamma(expx, 
        shape = k, scale =1), expx * dgamma(expx, shape = k, 
        scale = 1), 1 + (k - 1 - expx), 1 + 3 * (k - 1 - expx) + 
        (k - 1) * (k - 2 - expx) - expx * (k - 1 - expx))
},
quantile = function (p, parms) 
log(qgamma(p, shape = k, scale = 1)),
deviance = function (...) 
stop("deviance residuals not defined")
)
    
z.out <- survreg(Surv(log(Y))~ -1+X1+X2+X3+X4+X5+X6+X7+X8+X9+X10+X11+X12+X13+X14+X15+X16+X17+X18+X19+X20+X21+X22+X23+X24+as.factor(sender)+W + D,dist= myloggamma, data = Daresurv)

quantiles <- seq(0.0001,.9999,length=1000)
quantileMat <- exp(predict(z.out,type="quantile",p=quantiles))
sim.results = matrix(0, nrow(Daresurv), 100)
for (k in 1:100) {
sim.quantiles <- sample(1:1000,nrow(Daresurv),rep=T)
simulated.times <- quantileMat[cbind(1:nrow(Daresurv),sim.quantiles)]
sim.results[,k ] = simulated.times
}
qqplot(sim.results, Daresurv$Y)
all= c(rowMeans(sim.results), Daresurv$Y)
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
           main = "sender+X+W+D",
           cex.lab = 0.25,
           cex.axis = 0.25,
           cex.main = 0.5)
    abline(0, 1, lty = 1, col = "red", lwd = 1)
    


library(flexsurv)
z.out <- flexsurvreg(Surv(Y)~ -1+X1+X2+X3+X4+X5+X6+X7+X8+X9+X10+X11+X12+X13+X14+X15+X16+X17+X18+X19+X20+X21+X22+X23+X24+sender+W + D,dist="gamma", data = Daresurv)
mu = -1.281755
sigma = 1.706303
Q = 0.292506

predicted = model.matrix(z.out) %*% z.out$res[-(1:2),"est"]

mean.gengamma <- function(mu, sigma, Q, horizon = 100, ...) {
	surv <- function(t, ...) { 1-pgengamma(q = t, mu = mu, sigma = sigma, Q = Q, ...)}
	integrate(surv, 0, horizon, ...)$value}
what = summary(z.out, t = Daresurv$Y, fn = mean.gengamma)	
	
	
	