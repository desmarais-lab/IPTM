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


lambdaij = matrix(0, length(edge2), 1)
Xij = matrix(0, length(edge2), 24)
iter = 1
for (d in edge2) {
 lambdaij[iter, ] = mean(lambda[[d]][Dare$edge[[d]]$sender, Dare$edge[[d]]$receiver])
 Xs = p.d[d, 1] * X[[d]][[Dare$edge[[d]]$sender]][[1]][Dare$edge[[d]]$receiver,] + p.d[d, 2] * X[[d]][[Dare$edge[[d]]$sender]][[2]][Dare$edge[[d]]$receiver,]
 if (is.matrix(Xs)) {Xs = colMeans(Xs)}
 Xij[iter, ] = Xs
   iter = iter + 1
}
Xij = Xij[-which(Y == 0),]
Daresurv = cbind(Daresurv, Xij)
names(Daresurv)[5:28] = sapply(1:24, function(x) paste0("X",x))
#save(Daresurv, file = "Daresurv.RData")
#Daresurv = data.frame(Y = Y, L = lambdaij, W = as.factor(weeks), D = as.factor(day))
#Daresurv = Daresurv[-which(Daresurv$Y == 0),]
Daresurv$logL = log(Daresurv$L)

par(mfrow = c(1,2))
z.out <- survreg(Surv(Y)~ logL+W + D,dist= "exp", data = Daresurv)

quantiles <- seq(0.0001,.9999,length=1000)
quantileMat <- predict(z.out,type="quantile",p=quantiles)
sim.quantiles <- sample(1:1000,nrow(Daresurv),rep=T)
simulated.times <- quantileMat[cbind(1:nrow(Daresurv),sim.quantiles)]
qqplot(simulated.times, Daresurv$Y, main = "logL+W+D")
qqline(simulated.times)

z.out <- survreg(Surv(Y)~ offset(logL)+W + D,dist= "exp", data = Daresurv)

quantiles <- seq(0.0001,.9999,length=1000)
quantileMat <- predict(z.out,type="quantile",p=quantiles)
sim.quantiles <- sample(1:1000,nrow(Daresurv),rep=T)
simulated.times <- quantileMat[cbind(1:nrow(Daresurv),sim.quantiles)]
qqplot(simulated.times, Daresurv$Y, main = "offset(logL)+W+D")
qqline(simulated.times)



par(mfrow = c(2,2))
z.out <- survreg(Surv(Y)~ X1+X2+X3+X4+X5+X6+X7+X8+X9+X10+X11+X12+X13+X14+X15+X16+X17+X18+X19+X20+X21+X22+X23+X24+W + D,dist= "exp", data = Daresurv)

quantiles <- seq(0.0001,.9999,length=1000)
quantileMat <- predict(z.out,type="quantile",p=quantiles)
sim.quantiles <- sample(1:1000,nrow(Daresurv),rep=T)
simulated.times <- quantileMat[cbind(1:nrow(Daresurv),sim.quantiles)]
qqplot(simulated.times, Daresurv$Y, main = "X+W+D")
qqline(simulated.times)


z.out <- survreg(Surv(Y)~ X1+X2+X3+X4+X5+X6+X7+X8+X9+X10+X11+X12+X13+X14+X15+X16+X17+X18+X19+X20+X21+X22+X23+X24+W,dist= "exp", data = Daresurv)
quantiles <- seq(0.0001,.9999,length=1000)
quantileMat <- predict(z.out,type="quantile",p=quantiles)
sim.quantiles <- sample(1:1000,nrow(Daresurv),rep=T)
simulated.times <- quantileMat[cbind(1:nrow(Daresurv),sim.quantiles)]
qqplot(simulated.times, Daresurv$Y, main = "X+W")
qqline(simulated.times)


z.out <- survreg(Surv(Y)~ X1+X2+X3+X4+X5+X6+X7+X8+X9+X10+X11+X12+X13+X14+X15+X16+X17+X18+X19+X20+X21+X22+X23+X24+D,dist= "exp", data = Daresurv)
quantiles <- seq(0.0001,.9999,length=1000)
quantileMat <- predict(z.out,type="quantile",p=quantiles)
sim.quantiles <- sample(1:1000,nrow(Daresurv),rep=T)
simulated.times <- quantileMat[cbind(1:nrow(Daresurv),sim.quantiles)]
qqplot(simulated.times, Daresurv$Y, main = "X+D")
qqline(simulated.times)

z.out <- survreg(Surv(Y)~ X1+X2+X3+X4+X5+X6+X7+X8+X9+X10+X11+X12+X13+X14+X15+X16+X17+X18+X19+X20+X21+X22+X23+X24,dist= "exp", data = Daresurv)
quantiles <- seq(0.0001,.9999,length=1000)
quantileMat <- predict(z.out,type="quantile",p=quantiles)
sim.quantiles <- sample(1:1000,nrow(Daresurv),rep=T)
simulated.times <- quantileMat[cbind(1:nrow(Daresurv),sim.quantiles)]
qqplot(simulated.times, Daresurv$Y, main = "X")
qqline(simulated.times)



lambdaij = matrix(0, length(edge2), 1)
Xij = matrix(0, length(edge2), 24)
Zij =  matrix(0, length(edge2), 24)
iter = 1
for (d in edge2) {
 lambdaij[iter, ] = mean(lambda[[d]][Dare$edge[[d]]$sender, Dare$edge[[d]]$receiver])
 Xs = X[[d]][[Dare$edge[[d]]$sender]][[1]][Dare$edge[[d]]$receiver,]
 if (is.matrix(Xs)) Xs = colMeans(Xs)
 Xij[iter, ] = Xs
 Zs = X[[d]][[Dare$edge[[d]]$sender]][[2]][Dare$edge[[d]]$receiver,]
 if (is.matrix(Zs)) Zs = colMeans(Zs)
 Zij[iter,] =  Zs
   iter = iter + 1
}
Xij = Xij[-which(Y == 0),]
Zij = Zij[-which(Y == 0),]
Daresurv = cbind(Daresurv, Xij)
Daresurv = cbind(Daresurv, Zij)
names(Daresurv)[30:53] = sapply(1:24, function(x) paste0("XIP1",x))
names(Daresurv)[54:77] = sapply(1:24, function(x) paste0("XIP2",x))

z.out <- survreg(Surv(Y)~ XIP11+XIP12+XIP13+XIP14+XIP15+XIP16+XIP17+XIP18+XIP19+XIP110+XIP111+XIP112+XIP113+XIP114+XIP115+XIP116+XIP117+XIP118+XIP119+XIP120+XIP121+XIP122+XIP123+XIP124+
XIP21+XIP22+XIP23+XIP24+XIP25+XIP26+XIP27+XIP28+XIP29+XIP210+XIP211+XIP212+XIP213+XIP214+XIP215+XIP216+XIP217+XIP218+XIP219+XIP220+XIP221+XIP222+XIP223+XIP224+W + D,dist= "exp", data = Daresurv)

quantiles <- seq(0.0001,.9999,length=1000)
quantileMat <- predict(z.out,type="quantile",p=quantiles)
sim.quantiles <- sample(1:1000,nrow(Daresurv),rep=T)
simulated.times <- quantileMat[cbind(1:nrow(Daresurv),sim.quantiles)]
qqplot(simulated.times, Daresurv$Y, main = "X+W+D")
qqline(simulated.times)

