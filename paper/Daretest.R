library(IPTM)
attach(Dare)
Daretest <- MCMC(edge, node, text, vocab, nIP = 3, K = 10, delta_B = 0.01, lambda = 0.05, outer = 10)

table_betaIP(Daretest)

plot_betaIP(Daretest)

plot_topicIP(Daretest, 10)

plot_topic(Daretest, 10)

table_wordIP(Daretest, 10, text, vocab)


