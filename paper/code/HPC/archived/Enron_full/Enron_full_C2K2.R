library(IPTM)
load('Enron.RData')
nIP = 2
K = 2

for (i in 1:5){
    set.seed(i)
    Enrontest = IPTM_inference.data(edge = Enron$edge, node =Enron$node, textlist = Enron$text, vocabulary = Enron$vocab, nIP = nIP, K = K,
    sigma_Q = c(0.001, 0.01), alpha = 2, mvec = rep(1/K, K),
    betas = 2, nvec = rep(1/length(Enron$vocab), length(Enron$vocab)), prior.b.mean = rep(0,7), prior.b.var = diag(7),
    prior.delta = c(-3.5, 1), out = 50, n_B = 5500, n_d = 550, burn = c(500, 50), thinning = c(10, 1),
    netstat = c("intercept", "dyadic"), optimize = TRUE, initial = NULL)
    filename = paste0("Enron_full_",nIP,"_",K,"_ver",i,".RData")
    save(Enrontest, file = filename)
}
