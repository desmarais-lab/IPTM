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
  ll <- sum(log(dexp(y,rate)))
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
  average.lambda <- mean(lambda)
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
  training_data <- dyadic_sendrec[which(dyadic_sendrec[,1]<d & dyadic_sendrec[,1] >= first.training.obs),]
  training_results <- mcmc_logit(training_data[,"recieved"],training_data[,5:ncol(training_data)],burnin = 15000, mcmc = 25000,thin=50)
  result_file <- paste("./results/EnronlogitResultsD",d,".RData",sep="")
  save(list="training_results",file=result_file)
  print(d)
}



