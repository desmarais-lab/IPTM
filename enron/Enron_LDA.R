library(IPTM)
load('Enron.RData')
attach(Dare)
Enrontest <- IPTM_inference.LDA(Enron$edge, Enron$node, Enron$text, Enron$vocab, nIP = 2, K = 2, sigma_Q = c(0.0005, 0.01),
                        alpha = 2, mvec = rep(1/2, 2), betas = 2, nvec = rep(1/length(Dare$vocab), length(Dare$vocab)), 
                        prior.b.mean = c(-3, rep(0, 24)), 
                       prior.b.var = 1 * diag(25), prior.delta = c(0, 1), out = 100, n_B = 15000, n_d = 1500, burn = c(10000,500), 
                       thinning = c(10,5), netstat = c("intercept", "dyadic", "degree", "triadic"), optimize = TRUE)
Enrontest_LDA = Daretest
save(Enrontest_LDA, file = "Enrontest_LDA_K2.RData")
