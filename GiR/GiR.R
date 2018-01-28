library(IPTM)
library(FastGP)
library(MCMCpack)
library(LaplacesDemon)
set.seed(526113325)
nDocs = 20
node = 1:4
vocab = c("hi", "hello", "fine", "bye", "what")

nIP = 2
K = 4
n.d = 6
alpha = 2
mvec = rep(1/4, 4)
beta = 2
netstat = c("dyadic")
timestat = c("timeofday", "dayofweek")
L = 3
P = 6
prior.b = list(rep(0, P), 0.5* diag(P))
prior.delta = c(-2, 0.1)
prior.eta = list(rep(3, length(node) + length(timestat)), 0.5*diag(length(node) +length(timestat)))
prior.tau = 5
sigma.Q = c(0.15, 0.025, 0.05, 25)

b = matrix(c(prior.b[[1]],prior.b[[1]]), nrow = nIP, byrow = TRUE)
eta =  matrix(c(prior.eta[[1]],prior.eta[[1]]), nrow = nIP, byrow = TRUE)
delta = prior.delta[1]
sigma_tau = prior.tau

l = sample(1:nIP, K, replace = TRUE)
support = gibbs.measure.support(length(node)-1)
base.data = GenerateDocs(500, node, vocab, nIP, K, n.d, alpha, mvec, beta, b, eta, delta, sigma_tau, l, support, netstat, timestat,
                        base.edge = NULL, base.text = NULL, topic_token_assignments = NULL, 
                        backward = FALSE, base = TRUE) 
base.edge = base.data$edge	   
base.text = base.data$text


Outer = 1
Inner = c(1,1,1)
burn = c(0,0,0)
thin = c(1,1,1)
Outer = 100
Inner = c(10,10,10)
burn = c(0,0,0)
thin = c(1, 1, 1)
Schein <- Schein(5000, nDocs, node, vocab, nIP, K, n.d, alpha, mvec, beta, 
              prior.b, prior.delta, prior.eta, prior.tau, sigma.Q, Outer, Inner, burn, thin,
              netstat = c("dyadic"), timestat = c("timeofday", "dayofweek"),
              base.edge, base.text, generate_PP_plots = TRUE)

Outer = 3
Inner = c(3300,3300,550)
burn = c(300,300, 50)
thin = c(3,3, 1)
GettingItRight <- GiR(1000, nDocs, node, vocab, nIP, K, n.d, alpha, mvec, beta,
              prior.b, prior.delta, prior.eta, prior.tau, sigma.Q, Outer, Inner, burn, thin,
              netstat = c("dyadic"), timestat = c("timeofday", "dayofweek"),
              base.edge, base.text, generate_PP_plots = TRUE)
