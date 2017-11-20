library(IPTM)
library(FastGP)
library(MCMCpack)
set.seed(12345)
D = 5
node = 1:4
vocab = c("hi", "hello", "fine", "bye", "what")

nIP = 2
K = 4
n.d = 6
alpha = 2
mvec = rep(1/4, 4)
beta = 2
netstat = c("intercept","dyadic")
timestat = c("timeofday", "dayofweek")
#netstat = as.numeric(c("intercept", "degree", "dyadic", "triadic") %in% netstat)
#timestat = as.numeric(c("dayofweek","timeofday") %in% timestat)

L = 3
P = 7
prior.b = list(rep(0, P), diag(P))
prior.delta = c(-2.5, 0.1)
prior.eta = list(rep(1, P+2), diag(P+2))
prior.tau = c(2,1)
sigma.Q = c(0.01, 0.005)

b = lapply(1:nIP, function(IP) {c(rcpp_rmvnorm(1, prior.b[[2]], prior.b[[1]]))}) 
eta = lapply(1:nIP, function(IP) {c(rcpp_rmvnorm(1, prior.eta[[2]], prior.eta[[1]]))})
delta = rnorm(1, prior.delta[1], sqrt(prior.delta[2]))
sigma2_tau = 1/rgamma(1, prior.tau[1], prior.tau[2])

 l = sample(1:nIP, K, replace = TRUE)
 support = gibbs.measure.support(length(node) - 1)
base.data = GenerateDocs(1000, node, vocab, nIP, K, n.d, alpha, mvec, beta, b, eta, delta, sigma2_tau, l, support, netstat, timestat,
                        base.edge = NULL, base.text = NULL, topic_token_assignments = NULL, 
                        backward = FALSE, base = TRUE) 
base.edge = base.data$edge	   
base.text = base.data$text


Outer = 1
Inner = c(1,1)
burn = c(0,0)
thin = c(1,1)
Outer = 5
Inner = c(3300,3300)
burn = c(300,300)
thin = c(3,3)

Schein <- Schein(10000, D, node, vocab, nIP, K, n.d, alpha, mvec, beta, 
              prior.b, prior.delta, prior.eta, prior.tau, sigma.Q, Outer, Inner, burn, thin,
              netstat = c("intercept", "dyadic"), timestat = c("timeofday", "dayofweek"),
              base.edge, base.text, generate_PP_plots = TRUE)


Outer = 3
Inner = c(5500,5500)
burn = c(500,500)
thin = c(5,5)

GettingItRight <- GiR(5000, D, node, vocab, nIP, K, n.d, alpha, mvec, beta, 
              prior.b, prior.delta, prior.eta, prior.tau, sigma.Q, Outer, Inner, burn, thin,
              netstat = c("intercept", "dyadic"), timestat = c("timeofday", "dayofweek"),
              base.edge, base.text, generate_PP_plots = TRUE)
