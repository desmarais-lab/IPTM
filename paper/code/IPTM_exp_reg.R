library(survival)
library(Zelig)
loglik_expreg <- function(logEta,y,baseLambda){
  # logEta is the log of the scaling parameter
  # y is the dependent variable
  # baseLambda is the rate from the logit model
  rate <- exp(logEta)*baseLambda
  ll <- sum(log(dexp(y,rate)))
  ll
}

mcmc_exponential <-function(y,lambda,burnin = 5000, mcmc = 15000, thin = 10){
  require(MCMCpack)
  lambda.est <- 1/mean(y)
  average.lambda <- mean(lambda)
  init <- log(lambda.est/average.lambda)
  mcmcRes <- MCMCmetrop1R(loglik_expreg,init,burnin=burnin,mcmc=mcmc,thin=thin,baseLambda=lambda,y=y)
  mcmcRes
}

set.seed(1)
testDocuments = sample(1105:2210, 200, replace = FALSE)
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

# reshape into dyadic dataset, ignore initial history
dyadic_sendrec <- matrix(NA,nrow(observed_edge_data)*(nodes-1),features+5)
rowind <- 1
for(i in 1:nrow(observed_edge_data)){
  for(j in (1:nodes)[-observed_edge_data$sender[i]]){
    datij <- c(i,observed_edge_data$sender[i],j,observed_edge_data[i,1+j],netstats[i,observed_edge_data$sender[i],j,],observed_edge_data$timeinc[i])
    dyadic_sendrec[rowind,] <- datij
    rowind <- rowind + 1
  }
}
colnames(dyadic_sendrec) <- c("document","sender","alter","recieved",paste("f",1:features),"timeinc")

# number of simulated documents
nsim = 100
times = matrix(0, length(testDocuments), nsim)
rownames(times) = testDocuments
iter = 1
for (d in testDocuments){
  
  result_file <- paste("DarelogitResultsD",d,".RData",sep="")
  load(result_file)
  
  for (s in 1:nsim) {
  beta <- training_results[sample(1:nrow(training_results),1,rep=T),]
  meanLambda <- rep(NA,nrow(observed_edge_data))
  for(i in first.training.obs:(d-1)){
    recipient_features <- rbind(dyadic_sendrec[dyadic_sendrec[,1]==i & dyadic_sendrec[,4]==1,5:(ncol(dyadic_sendrec)-1)])
    meanX <- apply(recipient_features,2,mean)
    meanLambda[i] <- exp(rbind(meanX)%*%cbind(beta))
  }

  expreg <- mcmc_exponential(observed_edge_data$timeinc[first.training.obs:(d-1)],meanLambda[first.training.obs:(d-1)],burnin=10000,mcmc=1000)
  eta <- exp(expreg[100,1])
  
  # generate times for senders of document d
    dyadic_sendrec_d = dyadic_sendrec[dyadic_sendrec[,1] == d,]
    netstatsi <- rbind(netstats[d, dyadic_sendrec_d[1,2], dyadic_sendrec_d[which(dyadic_sendrec_d[,4]==1),3],])
    lambdaid <- exp(rbind(apply(netstatsi,2,mean))%*%cbind(beta))
    times_d <- rexp(1,exp(expreg[100,1])*lambdaid)
    times[iter, s] <- times_d
  }
  iter = iter + 1
}