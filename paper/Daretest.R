library(IPTM)
attach(Dare)
Daretest <- MCMC(edge, node, text, vocab, nIP = 3, K = 10, delta_B = 0.01, lambda = 0.05, outer = 10)

table_betaIP(Daretest)

plot_betaIP(Daretest)

plot_topicIP(Daretest, 10)

plot_topic(Daretest, 10)

table_wordIP(Daretest, 10, text, vocab)





load('/Users/bomin8319/Desktop/IPTMpaper/data/Darenew.RData')
attach(Dare)
Dare$node = 1:nrow(Dare$node)
for (n in 1:length(Dare$edge)){
  Dare$edge[[n]][3] = Dare$edge[[n]][[3]] / 3600
}
Dare$edge = lapply(Dare$edge, function(x){x[1:3]})