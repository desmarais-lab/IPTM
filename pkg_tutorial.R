#install the package from Github
library(devtools)
install_github("desmarais-lab/MulticastNetwork/pkg")
library(MulticastNetwork)
#load our Montgomery county email data
load("/Users/bomin8319/Desktop/MulticastNetwork/Lab/Montgomery.RData")
names(Montgomery)

#data including sender a_d, reciever vector r_d, timestamp t_d
edge = Montgomery$edge
head(edge)

#covariates affecting "who sends to whom"
X = Montgomery$X
head(X[100, 1, , ])
X = X[,,,c(1:5)]    #use few covariates for this application

#covariates affecting "when to send"
Y = Montgomery$Y
head(Y[100, , ])
Y = Y[,,c(1:5)]    #use few covariates for this application



##########
library(IPTM)
library(combinat)
library(mvtnorm)
library(MCMCpack)
library(Rcpp)
C = 2
K = 4
V = 5
D = 50
A = 5
P = 4
Q = 3
X = array(rnorm(D*A*A*P), dim = c(D,A,A,P))
Y = array(rnorm(D*A*Q), dim = c(D,A,Q))
prior.epsilon = 1
psi = rgamma(C, prior.epsilon, prior.epsilon)
theta = matrix(rgamma(C*K, prior.epsilon, prior.epsilon), C, K)
phi = matrix(rgamma(K*V, prior.epsilon, prior.epsilon), K, V)
support = gibbs.measure.support(A-1)
prior.beta = list(mean = rep(0, P), var = diag(P))
prior.eta = list(mean = rep(0, Q), var = diag(Q))
prior.sigma2 = list(a = 4, b = 1)
Nsamp = 2500
outer = 10
inner = c(5, 5, 5)
burn = 0
#Schein test
result = matrix(NA, Nsamp, 2*(10+P+Q))
for (n in 1:Nsamp) {
	beta = rmvnorm(1, prior.beta$mean, prior.beta$var)
	eta = rmvnorm(1, prior.eta$mean, prior.eta$var)
	sigma2 = 1/rgamma(1, prior.sigma2$a, prior.sigma2$b)
	initial = Generate(D, A, psi, theta, phi, beta, eta, sigma2, X, Y, support, prior.epsilon=1)
	infer =  Inference(initial$data, X, Y, C, K, V, outer, inner, burn, prior.epsilon, prior.beta, prior.eta, prior.sigma2, 
                    initial = initial, proposal.var = c(0.01, 0.1, 0.01, 0.5), lasttime = 0)
    initial2 = Generate(D, A, infer$psi, infer$theta, infer$phi, infer$beta[outer,], infer$eta[outer,], infer$sigma2[outer,], X, Y, support, prior.epsilon=1)                  
                
	result[n, ] = c(mean(vapply(initial$data, function(x) sum(x$r_d), c(1))),
                  var(vapply(initial$data, function(x) sum(x$r_d), c(1))),
 				 mean(vapply(2:D, function(d) initial$data[[d]]$t_d-initial$data[[d-1]]$t_d, c(1))),
 				var(vapply(2:D, function(d) initial$data[[d]]$t_d-initial$data[[d-1]]$t_d, c(1))),
                  initial$beta, initial$eta, initial$sigma2, mean(initial$pi), mean(initial$psi), mean(initial$theta), mean(initial$phi), mean(initial$c_d),
                 mean(vapply(initial2$data, function(x) sum(x$r_d), c(1))),
 			var(vapply(initial2$data, function(x) sum(x$r_d), c(1))),
                  mean(vapply(2:D, function(d) initial2$data[[d]]$t_d-initial2$data[[d-1]]$t_d, c(1))),
 				var(vapply(2:D, function(d) initial2$data[[d]]$t_d-initial2$data[[d-1]]$t_d, c(1))),
                  initial2$beta, initial2$eta, initial2$sigma2, mean(initial2$pi), mean(initial2$psi), mean(initial2$theta), mean(initial2$phi), mean(initial2$c_d))		
}
par(mfrow=c(4,5))
GiR_PP_Plots(result[1:(n-1),c(1:(10+P+Q))], result[1:(n-1),c((11+P+Q):(2*(10+P+Q)))])





#run inference to estimate beta, eta, u, and sigma2
prior.beta = list(mean = rep(0, P), var = 2*diag(P))
prior.eta = list(mean = rep(0, Q), var = 2*diag(Q))
prior.sigma2 = list(a = 2, b = 1)

outer = 500
inner = c(1, 1, 1)
burn = 0

#this initial values are my results after convergence
initialval = list()
initialval$beta = c(-3.548488159, -0.108898033,  0.086928281,  0.274102964,  0.029273815)
initialval$eta = c(7.38197540,  0.12883606, -1.07065227, -0.20698218, -0.05894448)
initialval$sigma2 = 14.09288
initialval$u = lapply(1:dim(X)[1], function(d) matrix(0, dim(X)[2], dim(X)[2]))

#run infernece
Montgomery_infer = Inference(edge, X, Y, outer, inner, burn, prior.beta, prior.eta, prior.sigma2, initialval = initialval,
		  proposal.var = c(0.0001, 0.001, 0.1), timeunit = 3600, lasttime = Montgomery$lasttime, timedist = "lognormal")
names(Montgomery_infer)

# check convergence
plot(Montgomery_infer$loglike, type = 'l')

# generate data from the model estimates
Montgomery_PPC = PPC(length(edge), beta = colMeans(Montgomery_infer$beta), eta = colMeans(Montgomery_infer$eta), 
                     sigma2 = mean(Montgomery_infer$sigma2), X, Y, timeunit = 3600, u = Montgomery_infer$u, timedist = "lognormal")

# prediction experiment
# hide 10% of senders, receivers, and timestamps
set.seed(1)
missing = list()
#missingness of senders
missing[[1]] = matrix(0, nrow = dim(Y)[1], 1)    
missing[[1]][sample(1:dim(Y)[1], 62, replace = FALSE), ] = 1
#missingness of receivers
missing[[2]] = matrix(0, nrow = dim(Y)[1], A)    
missing[[2]][sample(1:(dim(Y)[1]*A), 1118, replace = FALSE)] = 1
#missingness of timestamps
missing[[3]] = matrix(0, nrow = dim(Y)[1], 1)
missing[[3]][sample(1:dim(Y)[1], 62, replace = FALSE), ] = 1

# use parameter estimates as initial values
initial = list()
initial$beta = colMeans(Montgomery_infer$beta)
initial$eta =  colMeans(Montgomery_infer$eta)
initial$u = Montgomery_infer$u
initial$sigma2 = mean(Montgomery_infer$sigma2)

#will generate 10 predictions (iterate two steps: imputation -> inference)
Montgomery_PPE = PPE(edge, missing, X, Y, 10, c(5,5,1), 0, prior.beta, prior.eta, prior.sigma2, 
                     initial = initial, proposal.var = c(0.0001, 0.001, 0.1), timeunit = 3600, 
                     lasttime = Montgomery$lasttime, MHprop.var = 0.1, timedist = "lognormal")


