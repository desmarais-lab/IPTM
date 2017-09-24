loglik_logit <- function(beta,y,x){
  # beta is a vector of regression coefficients
  # y is the dependent variable
  # x is a matrix of covariates
  linear_predictor <- cbind(x)%*%cbind(beta)
  proby <- 1/(1+exp(-linear_predictor))
  ll <- sum(log(dbinom(y,1,proby)))
  ll
}

loglik_expreg <- function(logEta,y,baseLambda){
  # logEta is the log of the scaling parameter
  # y is the dependent variable
  # baseLambda is the rate from the logit model
  rate <- exp(logEta)*baseLambda
  ll <- sum(log(dexp(y,rate)+1/10000^2))
  ll
}

library(MCMCpack)

mcmc_logit <-function(y,X,burnin = 5000, mcmc = 15000, thin = 10){
  require(MCMCpack)
  init <- coef(glm(y~X-1,family="binomial"))
  mcmcRes <- MCMCmetrop1R(loglik_logit,init,burnin=burnin,mcmc=mcmc,thin=thin,x=X,y=y)
  mcmcRes
}


# test out logit
set.seed(5)
n <- 1000
X <- cbind(rnorm(n),rnorm(n))
y <- rbinom(n,1,prob=1/(1+exp(-(X%*%cbind(c(1,-1))))))
test.logit <- mcmc_logit(y,X)
glm(y~X,family="binomial")
summary(test.logit)


mcmc_exponential <-function(y,lambda,burnin = 5000, mcmc = 15000, thin = 10){
  require(MCMCpack)
  lambda.est <- 1/mean(y)
  average.lambda <- min(lambda)
  init <- log(lambda.est/average.lambda)
  mcmcRes <- MCMCmetrop1R(loglik_expreg,init,burnin=burnin,mcmc=mcmc,thin=thin,baseLambda=lambda,y=y)
  mcmcRes
}

# test out exponential model
set.seed(5)
n <- 1000
lambda <- exp(rnorm(n))
logEta <- 1.25
y <- rexp(n,lambda*exp(logEta))
test.expreg <- mcmc_exponential(y,lambda)
summary(test.expreg)


## Read in data
# load Dare observed data
load("./enron/Enron.RData")
# load Dare network statistics
load("./enron/X_Enron2.RData")

# extract the number of nodes
nodes <- dim(X)[2]

# extract number of features
features <- dim(X)[4]

# build dataset for training
# sender, receivers (nodes columns), timeinc
observed_edge_data <- matrix(NA,length(Enron$edge),nodes+3)
for(d in 1:length(Enron$edge)){
  dati <- c(Enron$edge[[d]]$sender)
  print(length(dati))
  receivers <- numeric(nodes)
  receivers[Enron$edge[[d]]$receiver] <- 1
  dati <- c(dati,receivers,Enron$edge[[d]][[4]],Enron$edge[[d]]$timestamp)
  print(length(dati))
  observed_edge_data[d,] <- dati
}
colnames(observed_edge_data) <- c("sender",paste("R",1:nodes,sep=""),"timeinc","time")

observed_edge_data <- data.frame(observed_edge_data)

#observed_edge_data <- observed_edge_data[observed_edge_data$time > 384,]

first.training.obs <- min(which(observed_edge_data$time > 384))

# reshape into dyadic dataset, ignore initial history
dyadic_sendrec <- matrix(NA,nrow(observed_edge_data)*(nodes-1),features+4)
rowind <- 1
for(i in 1:nrow(observed_edge_data)){
  for(j in (1:nodes)[-observed_edge_data$sender[i]]){
    datij <- c(i,observed_edge_data$sender[i],j,observed_edge_data[i,1+j],X[i,observed_edge_data$sender[i],j,])
    dyadic_sendrec[rowind,] <- datij
    rowind <- rowind + 1
  }
}
colnames(dyadic_sendrec) <- c("document","sender","alter","recieved",paste("f",1:features))

set.seed(1)
testDocuments = sample(1963:3925, 200, replace = FALSE)

for(d in testDocuments){

result_file <- paste("./results/EnronlogitResultsD",d,".RData",sep="")
load(result_file)
  
# number of simulated documents
nsim = 100

simulated_document_matrix <- matrix(NA,nsim,nodes+2)

colnames(simulated_document_matrix) <- c("sender_id",paste("recipient_indicator",1:nodes,sep=""),"time_increment")

for(s in 1:nsim){

beta <- training_results[sample(1:nrow(training_results),1,rep=T),]
meanLambda <- rep(NA,nrow(observed_edge_data))
for(i in first.training.obs:(d-1)){
  recipient_features <- rbind(dyadic_sendrec[dyadic_sendrec[,1]==i & dyadic_sendrec[,4]==1,5:(ncol(dyadic_sendrec))])
  meanX <- apply(recipient_features,2,mean)
  meanLambda[i] <- exp(rbind(meanX)%*%cbind(beta))
}

expreg <- mcmc_exponential(observed_edge_data$timeinc[first.training.obs:(d-1)],meanLambda[first.training.obs:(d-1)],burnin=10000,mcmc=1000)

# generate latent recipients for all senders of document d
adjacency_matrix <- matrix(0,nodes,nodes)
for(i in 1:nodes){
  recipients <- rep(0,nodes)
  while(sum(recipients)==0){
  for(j in (1:nodes)[-i]){
    pij <- 1/(1+exp(-rbind(X[d,i,j,])%*%cbind(beta)))
    if(pij > runif(1)) recipients[j] <- 1
  }
  }
  adjacency_matrix[i,] <- recipients
}


# generate times for all senders of document d
times_d <- numeric(nodes)
for(i in 1:nodes){
  netstatsi <- rbind(X[d,i,which(adjacency_matrix[i,]==1),])
  lambdaid <- exp(rbind(apply(netstatsi,2,mean))%*%cbind(beta))
  times_d[i] <- rexp(1,exp(expreg[100,1])*lambdaid)
}

simulated_sender_d <- which.min(times_d)

simulated_time_d <- times_d[simulated_sender_d]

simulated_recipients_d <- adjacency_matrix[simulated_sender_d,]

simulated_document_matrix[s,] <- c(simulated_sender_d,simulated_recipients_d,simulated_time_d)

}

simulation_result_file <- paste("./results/EnronSimulatedDocumentsD",d,".RData",sep="")

save(list="simulated_document_matrix",file=simulation_result_file)

print(d)

}







