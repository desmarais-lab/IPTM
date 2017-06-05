library(IPTM)
load('/Users/bomin8319/Desktop/IPTMpaper/data/Vancenew.RData')
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
Vancetest2 <- IPTM_inference.data(Vance$edge, Vance$node, Vance$text, Vance$vocab, nIP = 2, K = 5, sigma_Q = c(0.01, 1),
                       alpha = 2, mvec = rep(1/5, 5), betas = 2, nvec = rep(1/620, 620), prior.b.mean = rep(0, 25), 
                       prior.b.var = diag(25), prior.delta = c(0, 1), out = 500, n_B = 550, n_d = 100, burn = 50, 
                       thinning = 5, netstat = c("intercept", "dyadic", "degree", "triadic"), plot = TRUE, optimize = TRUE)

TablebetaIP(Vancetest)

PlotbetaIP(Vancetest2$B)

par(mfrow = c(1,2))
PlotTopicIP(Vancetest, 10)

PlotTopic(Vancetest, 5)

TableWordIP(Vancetest, 5, text, vocab)


load('/Users/bomin8319/Desktop/IPTMpaper/data/Darenew.RData')
# attach(Dare)
# Sandy = c(1038:1967)
# mintime = Dare$edge[[1038]][[3]]
# for (n in 1:length(Sandy)){
  # Dare$edge[[n]][1] = Dare$edge[[Sandy[n]]][1]
  # Dare$edge[[n]][2] = Dare$edge[[Sandy[n]]][2]
  # Dare$edge[[n]][3] = (Dare$edge[[Sandy[n]]][[3]] - mintime) / 3600
  # Dare$edge[[n]][4] = Dare$edge[[Sandy[n]]][4]
# }
# Dare$edge = lapply(Dare$edge[1:930], function(x){x[1:3]})
# Dare$text = Dare$text[Sandy]
# set.seed(1)
# Daretest <- IPTM_inference.data(Dare$edge, Dare$node, Dare$text, Dare$vocab, nIP = 2, K = 5, sigma_Q = c(0.001, 0.1),
                        # alpha = 2, mvec = rep(1/5, 5), betas = 2, nvec = rep(1/length(Dare$vocab), length(Dare$vocab)), 
                        # prior.b.mean = rep(0, 25), 
                       # prior.b.var = diag(25), prior.delta = c(0, 1), out = 1, n_B = 550, n_d = 100, burn = 50, 
                       # thinning = 5, netstat = c("intercept", "dyadic", "degree", "triadic"), plot = TRUE, optimize = TRUE)
                       
attach(Dare)
Dare$node = 1:nrow(Dare$node)
mintime = Dare$edge[[1]][[3]]
for (n in 1:length(Dare$edge)){
  Dare$edge[[n]][3] = (Dare$edge[[n]][[3]] - mintime) / 3600
}
Dare$edge = lapply(Dare$edge, function(x){x[1:3]})
set.seed(1231)
Daretest <- IPTM_inference.data(Dare$edge, Dare$node, Dare$text, Dare$vocab, nIP = 2, K = 5, sigma_Q = c(0.001, 0.1),
                        alpha = 2, mvec = rep(1/5, 5), betas = 2, nvec = rep(1/length(Dare$vocab), length(Dare$vocab)), 
                        prior.b.mean = c(-10, rep(0, 24)), 
                       prior.b.var = diag(25), prior.delta = c(0, 1), out = 10, n_B = 550, n_d = 100, burn = 50, 
                       thinning = 5, netstat = c("intercept", "dyadic", "degree", "triadic"), plot = TRUE, optimize = TRUE)