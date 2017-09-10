library(IPTM)
load('/Users/bomin8319/Desktop/IPTM/enron/Enron.RData')
Enrontest <- IPTM_inference.data(Enron$edge, Enron$node, Enron$text, Enron$vocab, nIP = 2, K = 10, sigma_Q = c(0.01, 1),
                       alpha = 2, mvec = rep(1/10, 10), betas = 2, nvec = rep(1/length(Enron$vocab), length(Enron$vocab)), prior.b.mean = c(-1, rep(0, 6)),
                       prior.b.var = 0.1 * diag(7), prior.delta = c(0, 1), out = 100, n_B = 5500, n_d = 550, burn = c(500, 50),
                       thinning = c(10, 1), netstat = c("intercept", "dyadic"), optimize = TRUE)

save(Enrontest, file = "Enrontest_IPTM_K10.RData")
