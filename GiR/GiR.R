library(IPTM)
library(FastGP)
library(MCMCpack)
library(LaplacesDemon)
set.seed(123)
nDocs = 10
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
prior.eta = list(rep(2, length(node) + 2), 0.5*diag(length(node) +2))
prior.tau = 5
sigma.Q = c(0.01, 0.01, 0.05, 5)

b = rcpp_rmvnorm(nIP, prior.b[[2]], prior.b[[1]])
eta =  rcpp_rmvnorm(nIP, prior.eta[[2]], prior.eta[[1]])
delta = rnorm(1, prior.delta[1], sqrt(prior.delta[2]))
sigma_tau = rhalfcauchy(1, prior.tau)

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
Outer = 3
Inner = c(3300,5500,1100)
burn = c(300,500, 100)
thin = c(3, 5, 2)
Schein <- Schein(1000, nDocs, node, vocab, nIP, K, n.d, alpha, mvec, beta, 
              prior.b, prior.delta, prior.eta, prior.tau, sigma.Q, Outer, Inner, burn, thin,
              netstat = c("dyadic"), timestat = c("timeofday", "dayofweek"),
              base.edge, base.text, generate_PP_plots = TRUE)



Outer = 3
Inner = c(3300,3300,550)
burn = c(300,300, 50)
thin = c(3,3, 1)
GettingItRight <- GiR(10000, nDocs, node, vocab, nIP, K, n.d, alpha, mvec, beta,
              prior.b, prior.delta, prior.eta, prior.tau, sigma.Q, Outer, Inner, burn, thin,
              netstat = c("dyadic"), timestat = c("timeofday", "dayofweek"),
              base.edge, base.text, generate_PP_plots = TRUE)
