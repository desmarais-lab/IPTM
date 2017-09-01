library(IPTM)
load('/Users/bomin8319/Desktop/IPTM/paper/code/Darenew.RData')
# 762 - 
attach(Dare)
Dare$node = 1:nrow(Dare$node)
#Dare$text = Dare$text[762:length(Dare$edge)]
#Dare$edge = Dare$edge[762:length(Dare$edge)]
Dare$edge = Dare$edge[-which(sapply(Dare$text, function(d){length(d)})==0)]
Dare$text = Dare$text[-which(sapply(Dare$text, function(d){length(d)})==0)]
mintime = Dare$edge[[1]][[3]]
for (n in 1:length(Dare$edge)){
  Dare$edge[[n]][3] = (Dare$edge[[n]][[3]] - mintime) / 3600
}
Dare$edge = lapply(Dare$edge, function(x){x[1:3]})


Daretest <- IPTM_inference.LDA(Dare$edge, Dare$node, Dare$text, Dare$vocab, nIP = 2, K = 10, sigma_Q = c(0.0005, 0.01),
                        alpha = 2, mvec = rep(1/10, 10), betas = 2, nvec = rep(1/length(Dare$vocab), length(Dare$vocab)), 
                        prior.b.mean = c(-3, rep(0, 24)), 
                       prior.b.var = 1 * diag(25), prior.delta = c(0, 1), out = 1000, n_B = 15000, n_d = 1500, burn = c(10000,500), 
                       thinning = c(10,5), netstat = c("intercept", "dyadic", "degree", "triadic"), optimize = TRUE)
Daretest_LDA = Daretest
save(Daretest_LDA, file = "Daretest_LDA_K5.RData")



library(lda)
library(SpeedReader)
load('~/Desktop/IPTM/paper/code/AISTAT_version/Vancetest_LDA_K5.RData')
K = 5
top.words = TableWord(Vancetest_LDA$Z, K, Vance$text, Vance$vocab)
topic.word = matrix(0, nrow = length(Vancetest_LDA$Z), ncol = length(Vance$vocab))
rownames(topic.word) = 1:length(Vance$text)
colnames(topic.word) = Vance$vocab
for (d in seq(along = Vance$text)) {
	topic.word[d, ] = tabulate(Vance$text[[d]], length(Vance$vocab))
}
sapply(1:K, function(k) {
topic_coherence(top.words[,k], topic.word,Vance$vocab)	
})

load('~/Desktop/IPTM/paper/code/AISTAT_version/Vancetest_LDA_K5.RData')
K = 10
top.words = TableWord(Vancetest_LDA$Z, K, Vance$text, Vance$vocab)
topic.word = matrix(0, nrow = length(Vancetest_LDA$Z), ncol = length(Vance$vocab))
rownames(topic.word) = 1:length(Vance$text)
colnames(topic.word) = Vance$vocab
for (d in seq(along = Vance$text)) {
	topic.word[d, ] = tabulate(Vance$text[[d]], length(Vance$vocab))
}
sapply(1:K, function(k) {
topic_coherence(top.words[,k], topic.word,Vance$vocab)	
})
