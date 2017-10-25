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
timestamps = vapply(Dare$edge, function(d) {
  d[[3]]
}, c(1))

edge2 = which_int(384, timestamps) : length(Dare$edge)
p.d = pdmat(Z, C, 2)
netstat = c(0, 1,1,1)
L = 3
P = L * (2 * netstat[2] + 2 * netstat[3] + 4 * netstat[4])
lambda = list()
for (d in edge2) {
  history.t = History(Dare$edge, p.d, Dare$node, Dare$edge[[d-1]][[3]] + exp(-745))
  X = Netstats_cpp(history.t, Dare$node, netstat)
  XB = MultiplyXBList(X, B)     
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
iter = 1
for (d in edge2) {
 lambdaij[iter, ] = mean(lambda[[d]][Dare$edge[[d]]$sender, Dare$edge[[d]]$receiver])
  iter = iter + 1
}

save(Daresurv, file = "Daresurv.RData")
Daresurv = data.frame(Y = Y, L = lambdaij, W = as.factor(weeks), D = as.factor(day))
Daresurv = Daresurv[-which(Daresurv$Y == 0),]
z.out <- zelig(Surv(Y)~ offset(L) + W + D, model = "exp", data = Daresurv)

x.out = setx(z.out, W = Daresurv$W[1], D = Daresurv$D[1] )

s.out <- sim(z.out, x = x.out)

