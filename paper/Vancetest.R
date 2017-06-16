library(IPTM)
load('/Users/bomin8319/Desktop/IPTM/VanceServer/Vancenew.RData')
attach(Vance)
Vance$node = 1:nrow(Vance$node)
mintime = Vance$edge[[1]][[3]]
for (n in 1:length(Vance$edge)){
  Vance$edge[[n]][3] = (Vance$edge[[n]][[3]] - mintime) / 3600
}
Vance$edge = lapply(Vance$edge, function(x){x[1:3]})
set.seed(1)

#Vance$node = Vance$node[-c(5, 8, 16, 17)]
#delete = sapply(1:length(Vance$edge), function(d){(unlist(Vance$edge[[d]][1:2]) %in% Vance$node)})
#Vance$edge = Vance$edge[- which(sapply(1:length(delete), function(d){"FALSE" %in% delete[[d]]})>0)]
Vancetest <- IPTM_inference.data(Vance$edge, Vance$node, Vance$text, Vance$vocab, nIP = 2, K = 5, sigma_Q = c(0.01, 1),
                       alpha = 2, mvec = rep(1/5, 5), betas = 2, nvec = rep(1/620, 620), prior.b.mean = c(-3, rep(0, 24)), 
                       prior.b.var = 0.1 * diag(25), prior.delta = c(0, 1), out = 800, n_B = 5500, n_d = 500, burn = 500, 
                       thinning = 10, netstat = c("intercept", "dyadic", "degree", "triadic"), plot = TRUE, optimize = TRUE)

TablebetaIP(Vancetest)

PlotbetaIP(Vancetest$B)

par(mfrow = c(1,2))
PlotTopicIP(Vancetest, 10)

PlotTopic(Vancetest, 5)

TableWordIP(Vancetest, 5, text, vocab)

set.seed(100)
load('/Users/bomin8319/Desktop/IPTM/paper/Darenew.RData')
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

initial = list(C = Daretest$C, D = Daretest$D[500], B = lapply(1:2, function(IP){rowMeans(Daretest$B[[IP]])}), Z = Daretest$Z)

Daretest1 <- IPTM_inference.Schein(Dare$edge, Dare$node, Dare$text, Dare$vocab, nIP = 2, K = 20, sigma_Q = c(0.0001, 0.01),
                        alpha = 2, mvec = rep(1/20, 20), betas = 2, nvec = rep(1/length(Dare$vocab), length(Dare$vocab)), 
                        prior.b.mean = c(-3, rep(0, 24)), 
                       prior.b.var = 1 * diag(25), prior.delta = c(0, 1), out = 1, n_B = 110, n_d = 1, burn = 10, 
                       thinning = 1, netstat = c("intercept", "dyadic", "degree", "triadic"), plot = FALSE, optimize = TRUE, 
                       initial = initial)

Daretest1 <- IPTM_inference.Schein(Dare$edge, Dare$node, Dare$text, Dare$vocab, nIP = 2, K = 20, sigma_Q = c(0.01, 0.01),
                        alpha = 2, mvec = rep(1/20, 20), betas = 2, nvec = rep(1/length(Dare$vocab), length(Dare$vocab)), 
                        prior.b.mean = c(-3, rep(0, 24)), 
                       prior.b.var = 1 * diag(25), prior.delta = c(0, 1), out = 1, n_B = 55000, n_d = 1, burn = 5000, 
                       thinning = 50, netstat = c("intercept", "dyadic", "degree", "triadic"), plot = FALSE, optimize = TRUE, 
                       initial = initial)
save(Daretest, file = "Daretest.RData")
save(Daretest1, file = "Daretest1.RData")

PlotbetaIP(Daretest$B)
