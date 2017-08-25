library(IPTM)
load('/Users/bomin8319/Desktop/IPTM/paper/code/Vancenew.RData')
attach(Vance)
Vance$node = 1:nrow(Vance$node)
mintime = Vance$edge[[1]][[3]]
for (n in 1:length(Vance$edge)){
  Vance$edge[[n]][3] = (Vance$edge[[n]][[3]] - mintime) / 3600
}
Vance$edge = Vance$edge[-which(sapply(Vance$text, function(d){length(d)})==0)]
Vance$text = Vance$text[-which(sapply(Vance$text, function(d){length(d)})==0)]

Vance$edge = lapply(Vance$edge, function(x){x[1:3]})
set.seed(1)

#Vance$node = Vance$node[-c(5, 8, 16, 17)]
#delete = sapply(1:length(Vance$edge), function(d){(unlist(Vance$edge[[d]][1:2]) %in% Vance$node)})
#Vance$edge = Vance$edge[- which(sapply(1:length(delete), function(d){"FALSE" %in% delete[[d]]})>0)]
Vancetest <- IPTM_inference.data(Vance$edge, Vance$node, Vance$text, Vance$vocab, nIP = 2, K = 10, sigma_Q = c(0.01, 1),
                       alpha = 2, mvec = rep(1/10, 10), betas = 2, nvec = rep(1/620, 620), prior.b.mean = c(-5, rep(0, 24)), 
                       prior.b.var = 0.1 * diag(25), prior.delta = c(0, 1), out = 200, n_B = 5500, n_d = 550, burn = c(500, 50), 
                       thinning = c(10, 1), netstat = c("intercept", "dyadic", "degree", "triadic"), optimize = TRUE)
Vancetest_new = Vancetest
save(Vancetest_new, file = "/Users/bomin8319/Desktop/IPTM/paper/code/Vancetest_new.RData")

Vancetest <- IPTM_inference.LDA(Vance$edge, Vance$node, Vance$text, Vance$vocab, nIP = 2, K = 10, sigma_Q = c(0.01, 1),
                                 alpha = 2, mvec = rep(1/10, 10), betas = 2, nvec = rep(1/620, 620), prior.b.mean = c(-5, rep(0, 24)), 
                                 prior.b.var = 0.1 * diag(25), prior.delta = c(0, 1), out = 200, n_B = 5500, n_d = 550, burn = c(500, 50), 
                                 thinning = c(10, 1), netstat = c("intercept", "dyadic", "degree", "triadic"), optimize = TRUE)
Vancetest_LDA = Vancetest
save(Vancetest_LDA, file = "/Users/bomin8319/Desktop/IPTM/paper/code/Vancetest_LDA.RData")




TablebetaIP(Vancetest)

PlotbetaIP(Vancetest$B)


PlotTopic(Vancetest, 10)


set.seed(100)
load('/Users/bomin8319/Desktop/IPTM/paper/code/Darenew.RData')
# 762 - 
attach(Dare)
Dare$node = 1:nrow(Dare$node)
Dare$text = Dare$text[762:length(Dare$edge)]
Dare$edge = Dare$edge[762:length(Dare$edge)]
Dare$edge = Dare$edge[-which(sapply(Dare$text, function(d){length(d)})==0)]
Dare$text = Dare$text[-which(sapply(Dare$text, function(d){length(d)})==0)]
mintime = Dare$edge[[1]][[3]]
for (n in 1:length(Dare$edge)){
  Dare$edge[[n]][3] = (Dare$edge[[n]][[3]] - mintime) / 3600
}
Dare$edge = lapply(Dare$edge, function(x){x[1:3]})


Daretest <- IPTM_inference.data(Dare$edge, Dare$node, Dare$text, Dare$vocab, nIP = 3, K = 20, sigma_Q = c(0.0005, 0.01),
                        alpha = 2, mvec = rep(1/20, 20), betas = 2, nvec = rep(1/length(Dare$vocab), length(Dare$vocab)), 
                        prior.b.mean = c(-3, rep(0, 24)), 
                       prior.b.var = 1 * diag(25), prior.delta = c(0, 1), out = 100, n_B = 15000, n_d = 1500, burn = c(10000,500), 
                       thinning = c(10,5), netstat = c("intercept", "dyadic", "degree", "triadic"), plot = FALSE, optimize = TRUE)
load("/Users/bomin8319/Desktop/IPTM/paper/code/Daretest.RData")
output = Daretest
initial = list(C = output$C, D = output$D[200], B = lapply(1:2, function(IP){output$B[[IP]][,500]}), Z =output$Z)

Daretest1 <- IPTM_inference.data2(Dare$edge, Dare$node, Dare$text, Dare$vocab, nIP = 3, K = 20, sigma_Q = c(0.00005, 0.0075),
                        alpha = 2, mvec = rep(1/20, 20), betas = 2, nvec = rep(1/length(Dare$vocab), length(Dare$vocab)), 
                        prior.b.mean = c(-3, rep(0, 24)), 
                       prior.b.var = 1 * diag(25), prior.delta = c(0, 1), out = 1, n_B = 500000, n_d = 50000, burn = c(0,0), 
                       thinning =  c(1,1), netstat = c("intercept", "dyadic", "degree", "triadic"), plot = FALSE, optimize = TRUE, 
                       initial = initial)
save(Daretest, file = "/Users/bomin8319/Desktop/Daretest.RData")
save(Daretest1, file = "Daretest1.RData")

PlotbetaIP(Daretest$B)
