library(IPTMnew)
library(FastGP)
library(MCMCpack)
library(LaplacesDemon)
set.seed(526113325)
nDocs = 20
node = 1:4
vocab = c("hi", "hello", "fine", "bye", "what")

nIP = 2
K = 4
n.d = 10
# alphas = c(25, 50, 50)
# beta = 50
# zeta = 10
alphas = c(5,5,25)
beta = 50
zeta = 5
netstat = c("dyadic")
timestat = c("timeofday", "dayofweek")
L = 3
P = 6
prior.b = list(rep(0.5, P), 0.5* diag(P))
prior.delta = c(-2, 0.1)
prior.eta = list(rep(2.5, length(node) + length(timestat)), 0.5*diag(length(node) +length(timestat)))
prior.tau = 5
sigma.Q = c(0.01, 0.0001, 0.01, 0.5)

b = matrix(c(prior.b[[1]],prior.b[[1]]), nrow = nIP, byrow = TRUE)
eta =  matrix(c(prior.eta[[1]],prior.eta[[1]]), nrow = nIP, byrow = TRUE)
delta = prior.delta[1]
sigma_tau = prior.tau
psi = rdirichlet(1, rep(zeta/nIP, nIP))
cd = vapply(1:500, function(i) which(rmultinom(1,1,psi)==1), c(1))
support = gibbs.measure.support(length(node)-1)
base.data = GenerateDocs(500, node, vocab, nIP, K, n.d, alphas, beta, zeta, b, eta, delta, sigma_tau, cd, support, netstat, timestat,
                        base.data= NULL, topic_token_assignments = NULL, 
                        backward = FALSE, base = TRUE) 


Outer = 1
#Inner = c(1,1,1)
Outer = 10
Inner = c(5,5,5)
Schein <- Schein(50, nDocs, node, vocab, nIP, K, n.d, alphas, beta, zeta,
              prior.b, prior.delta, prior.eta, prior.tau, sigma.Q, Outer, Inner,
               netstat = c("dyadic"), timestat = c("timeofday", "dayofweek"),
              base.data = base.data, generate_PP_plots = TRUE)

#######################

library(IPTMnew)
library(FastGP)
library(MCMCpack)
library(LaplacesDemon)
set.seed(526113325)
nDocs = 10
node = 1:4
vocab = c("hi", "hello", "fine", "bye", "what")

nIP = 2
K = 4
n.d = 5
alphas = c(5,5,5)
beta = 5
zeta = 5
netstat = c("dyadic")
timestat = c("timeofday", "dayofweek")
L = 3
P = 6
prior.b = list(rep(0.5, P), 0.5* diag(P))
prior.delta = c(-2, 0.1)
prior.eta = list(rep(2.5, length(node) + length(timestat)), 0.5*diag(length(node) +length(timestat)))
prior.tau = 5
sigma.Q = c(0.01, 0.0001, 0.00001, 0.5)

b = matrix(c(prior.b[[1]],prior.b[[1]]), nrow = nIP, byrow = TRUE)
eta =  matrix(c(prior.eta[[1]],prior.eta[[1]]), nrow = nIP, byrow = TRUE)
delta = prior.delta[1]
sigma_tau = prior.tau
psi = rdirichlet(1, rep(zeta/nIP, nIP))
cd = vapply(1:500, function(i) which(rmultinom(1,1,psi)==1), c(1))
support = gibbs.measure.support(length(node)-1)
base.data = GenerateDocs(500, node, vocab, nIP, K, n.d, alphas, beta, zeta, b, eta, delta, sigma_tau, cd, support, netstat, timestat,
                         base.data= NULL, topic_token_assignments = NULL, 
                         backward = FALSE, base = TRUE) 


Outer = 1
Inner = c(1,1,1)
Outer = 10
Inner = c(5,5,5)
Schein <- Schein(1000, nDocs, node, vocab, nIP, K, n.d, alphas, beta, zeta,
                 prior.b, prior.delta, prior.eta, prior.tau, sigma.Q, Outer, Inner,
                 netstat = c("dyadic"), timestat = c("timeofday", "dayofweek"),
                 base.data = base.data, generate_PP_plots = TRUE)

############################

Outer = 3
Inner = c(3300,3300,550)
burn = c(300,300, 50)
thin = c(3,3, 1)
GettingItRight <- GiR(1000, nDocs, node, vocab, nIP, K, n.d, alpha, mvec, beta,
              prior.b, prior.delta, prior.eta, prior.tau, sigma.Q, Outer, Inner, burn, thin,
              netstat = c("dyadic"), timestat = c("timeofday", "dayofweek"),
              base.edge, base.text, generate_PP_plots = TRUE)
