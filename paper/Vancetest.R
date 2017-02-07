library(IPTM)
attach(Vance)
Vancetest <- MCMC(edge, node, text, vocab, nIP = 2, K = 10, delta_B = 0.5, lambda = 0.05, outer = 10)

table_betaIP(Vancetest)

plot_betaIP(Vancetest)

plot_topicIP(Vancetest, 10)

plot_topic(Vancetest, 10)

table_wordIP(Vancetest, 10, text, vocab)


