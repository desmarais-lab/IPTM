

library(lda)
library(SpeedReader)
load('~/Desktop/IPTM/paper/code/Darenew.RData')
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

TableWord = function(Zchain, K, textlist, vocabulary, top = 20) {
  # Generate a table of token-topic assignments with high probabilities for each IP
  #
  # Args 
  #  Zchain summary of Z obtained using MCMC function
  #  K total number of topics specified by the user
  #  textlist list of text containing the words in each document
  #  vocabulary all vocabularies used over the corpus
  #
  # Returns
  #  List of table that summarize token-topic assignments for each IP
  W = length(vocabulary)
    Zsummary = list()
    topic.word = matrix(0, nrow = K, ncol = W)
    colnames(topic.word) = vocabulary
    iter = 1
    for (d in seq(along = textlist)) {
      if (length(Zchain[[d]]) > 0){
        Zsummary[[iter]] = Zchain[[d]]
        names(Zsummary[[iter]])<- vocabulary[textlist[[d]]]
        iter = iter+1
      }
    }
    topic.dist = t(tabulate(unlist(Zsummary), K)/length(unlist(Zsummary)))
    colnames(topic.dist) = c(1:K)
    top.topic = topic.dist[, order(topic.dist, decreasing = TRUE)]
    all.word = unlist(Zsummary)
    for (i in seq(along = all.word)){
      matchWZ = which(colnames(topic.word) == names(all.word[i]))
      topic.word[all.word[i], matchWZ] = topic.word[all.word[i], matchWZ] + 1
    }
    table.word = top.topic.words(topic.word, num.words = top)
    colnames(table.word) = names(top.topic)
  return(table.word)
}

Dare_topic = matrix(0, 4, 8)
colnames(Dare_topic) = as.numeric(c("1", "2", "5", "10", "20", "30", "40", "50"))

load('~/Desktop/IPTM/paper/code/HPC/topic coherence/Daretest_IPTM_K50.RData')
top.words = TableWord(Daretest$Z, 50, Dare$text, Dare$vocab)
topic.word = matrix(0, nrow = length(Daretest$Z), ncol = length(Dare$vocab))
rownames(topic.word) = 1:length(Dare$text)
colnames(topic.word) = Dare$vocab
for (d in seq(along = Dare$text)) {
	topic.word[d, ] = tabulate(Dare$text[[d]], length(Dare$vocab))
}
K =50
Dare_topic[1,8] = mean(sapply(1:K, function(k) {
topic_coherence(top.words[,k], topic.word,Dare$vocab)	
}))

load('~/Desktop/IPTM/paper/code/HPC/topic coherence/Daretest_IPTM_K40.RData')
top.words = TableWord(Daretest$Z, 40, Dare$text, Dare$vocab)
topic.word = matrix(0, nrow = length(Daretest$Z), ncol = length(Dare$vocab))
rownames(topic.word) = 1:length(Dare$text)
colnames(topic.word) = Dare$vocab
for (d in seq(along = Dare$text)) {
	topic.word[d, ] = tabulate(Dare$text[[d]], length(Dare$vocab))
}
K =40
Dare_topic[1,7] = mean(sapply(1:K, function(k) {
topic_coherence(top.words[,k], topic.word,Dare$vocab)	
}))

load('~/Desktop/IPTM/paper/code/HPC/topic coherence/Daretest_IPTM_K30.RData')
top.words = TableWord(Daretest$Z, 30, Dare$text, Dare$vocab)
topic.word = matrix(0, nrow = length(Daretest$Z), ncol = length(Dare$vocab))
rownames(topic.word) = 1:length(Dare$text)
colnames(topic.word) = Dare$vocab
for (d in seq(along = Dare$text)) {
	topic.word[d, ] = tabulate(Dare$text[[d]], length(Dare$vocab))
}
K =30
Dare_topic[1,6] = mean(sapply(1:K, function(k) {
topic_coherence(top.words[,k], topic.word,Dare$vocab)	
}))


load('~/Desktop/IPTM/paper/code/HPC/topic coherence/Daretest_IPTM_K20.RData')
top.words = TableWord(Daretest$Z, 20, Dare$text, Dare$vocab)
topic.word = matrix(0, nrow = length(Daretest$Z), ncol = length(Dare$vocab))
rownames(topic.word) = 1:length(Dare$text)
colnames(topic.word) = Dare$vocab
for (d in seq(along = Dare$text)) {
	topic.word[d, ] = tabulate(Dare$text[[d]], length(Dare$vocab))
}
K =20
Dare_topic[1,5] = mean(sapply(1:K, function(k) {
topic_coherence(top.words[,k], topic.word,Dare$vocab)	
}))


load('~/Desktop/IPTM/paper/code/HPC/topic coherence/Daretest_IPTM_K10.RData')
top.words = TableWord(Daretest$Z, 10, Dare$text, Dare$vocab)
topic.word = matrix(0, nrow = length(Daretest$Z), ncol = length(Dare$vocab))
rownames(topic.word) = 1:length(Dare$text)
colnames(topic.word) = Dare$vocab
for (d in seq(along = Dare$text)) {
	topic.word[d, ] = tabulate(Dare$text[[d]], length(Dare$vocab))
}
K = 10
Dare_topic[1,4] = mean(sapply(1:K, function(k) {
topic_coherence(top.words[,k], topic.word,Dare$vocab)	
}))



load('~/Desktop/IPTM/paper/code/HPC/topic coherence/Daretest_IPTM_K5.RData')
top.words = TableWord(Daretest$Z, 5, Dare$text, Dare$vocab)
topic.word = matrix(0, nrow = length(Daretest$Z), ncol = length(Dare$vocab))
rownames(topic.word) = 1:length(Dare$text)
colnames(topic.word) = Dare$vocab
for (d in seq(along = Dare$text)) {
	topic.word[d, ] = tabulate(Dare$text[[d]], length(Dare$vocab))
}
K = 5
Dare_topic[1,3] = mean(sapply(1:K, function(k) {
topic_coherence(top.words[,k], topic.word,Dare$vocab)	
}))


load('~/Desktop/IPTM/paper/code/HPC/topic coherence/Daretest_IPTM_K2.RData')
top.words = TableWord(Daretest$Z, 2, Dare$text, Dare$vocab)
topic.word = matrix(0, nrow = length(Daretest$Z), ncol = length(Dare$vocab))
rownames(topic.word) = 1:length(Dare$text)
colnames(topic.word) = Dare$vocab
for (d in seq(along = Dare$text)) {
	topic.word[d, ] = tabulate(Dare$text[[d]], length(Dare$vocab))
}
K = 2
Dare_topic[1,2] = mean(sapply(1:K, function(k) {
topic_coherence(top.words[,k], topic.word,Dare$vocab)	
}))

for (i in 1:length(Daretest$Z)) {
    Daretest$Z[[i]] = rep(1, length(Daretest$Z[[i]]))
}
top.words = TableWord(Daretest$Z, 1, Dare$text, Dare$vocab)
topic.word = matrix(0, nrow = length(Daretest$Z), ncol = length(Dare$vocab))
rownames(topic.word) = 1:length(Dare$text)
colnames(topic.word) = Dare$vocab
for (d in seq(along = Dare$text)) {
    topic.word[d, ] = tabulate(Dare$text[[d]], length(Dare$vocab))
}
K = 1
Dare_topic[1,1] = mean(sapply(1:K, function(k) {
    topic_coherence(top.words[,k], topic.word,Dare$vocab)
}))




load('~/Desktop/IPTM/paper/code/HPC/topic coherence/sDaretest_IPTM_K50.RData')
top.words = TableWord(Daretest$Z, 50, Dare$text, Dare$vocab)
topic.word = matrix(0, nrow = length(Daretest$Z), ncol = length(Dare$vocab))
rownames(topic.word) = 1:length(Dare$text)
colnames(topic.word) = Dare$vocab
for (d in seq(along = Dare$text)) {
    topic.word[d, ] = tabulate(Dare$text[[d]], length(Dare$vocab))
}
K =50
Dare_topic[2,8] = mean(sapply(1:K, function(k) {
    topic_coherence(top.words[,k], topic.word,Dare$vocab)
}))

load('~/Desktop/IPTM/paper/code/HPC/topic coherence/sDaretest_IPTM_K40.RData')
top.words = TableWord(Daretest$Z, 40, Dare$text, Dare$vocab)
topic.word = matrix(0, nrow = length(Daretest$Z), ncol = length(Dare$vocab))
rownames(topic.word) = 1:length(Dare$text)
colnames(topic.word) = Dare$vocab
for (d in seq(along = Dare$text)) {
    topic.word[d, ] = tabulate(Dare$text[[d]], length(Dare$vocab))
}
K =40
Dare_topic[2,7] = mean(sapply(1:K, function(k) {
    topic_coherence(top.words[,k], topic.word,Dare$vocab)
}))

load('~/Desktop/IPTM/paper/code/HPC/topic coherence/sDaretest_IPTM_K30.RData')
top.words = TableWord(Daretest$Z, 30, Dare$text, Dare$vocab)
topic.word = matrix(0, nrow = length(Daretest$Z), ncol = length(Dare$vocab))
rownames(topic.word) = 1:length(Dare$text)
colnames(topic.word) = Dare$vocab
for (d in seq(along = Dare$text)) {
    topic.word[d, ] = tabulate(Dare$text[[d]], length(Dare$vocab))
}
K =30
Dare_topic[2,6] = mean(sapply(1:K, function(k) {
    topic_coherence(top.words[,k], topic.word,Dare$vocab)
}))


load('~/Desktop/IPTM/paper/code/HPC/topic coherence/sDaretest_IPTM_K20.RData')
top.words = TableWord(Daretest$Z, 20, Dare$text, Dare$vocab)
topic.word = matrix(0, nrow = length(Daretest$Z), ncol = length(Dare$vocab))
rownames(topic.word) = 1:length(Dare$text)
colnames(topic.word) = Dare$vocab
for (d in seq(along = Dare$text)) {
    topic.word[d, ] = tabulate(Dare$text[[d]], length(Dare$vocab))
}
K =20
Dare_topic[2,5] = mean(sapply(1:K, function(k) {
    topic_coherence(top.words[,k], topic.word,Dare$vocab)
}))


load('~/Desktop/IPTM/paper/code/HPC/topic coherence/sDaretest_IPTM_K10.RData')
top.words = TableWord(Daretest$Z, 10, Dare$text, Dare$vocab)
topic.word = matrix(0, nrow = length(Daretest$Z), ncol = length(Dare$vocab))
rownames(topic.word) = 1:length(Dare$text)
colnames(topic.word) = Dare$vocab
for (d in seq(along = Dare$text)) {
    topic.word[d, ] = tabulate(Dare$text[[d]], length(Dare$vocab))
}
K = 10
Dare_topic[2,4] = mean(sapply(1:K, function(k) {
    topic_coherence(top.words[,k], topic.word,Dare$vocab)
}))



load('~/Desktop/IPTM/paper/code/HPC/topic coherence/sDaretest_IPTM_K5.RData')
top.words = TableWord(Daretest$Z, 5, Dare$text, Dare$vocab)
topic.word = matrix(0, nrow = length(Daretest$Z), ncol = length(Dare$vocab))
rownames(topic.word) = 1:length(Dare$text)
colnames(topic.word) = Dare$vocab
for (d in seq(along = Dare$text)) {
    topic.word[d, ] = tabulate(Dare$text[[d]], length(Dare$vocab))
}
K = 5
Dare_topic[2,3] = mean(sapply(1:K, function(k) {
    topic_coherence(top.words[,k], topic.word,Dare$vocab)
}))


load('~/Desktop/IPTM/paper/code/HPC/topic coherence/sDaretest_IPTM_K2.RData')
top.words = TableWord(Daretest$Z, 2, Dare$text, Dare$vocab)
topic.word = matrix(0, nrow = length(Daretest$Z), ncol = length(Dare$vocab))
rownames(topic.word) = 1:length(Dare$text)
colnames(topic.word) = Dare$vocab
for (d in seq(along = Dare$text)) {
    topic.word[d, ] = tabulate(Dare$text[[d]], length(Dare$vocab))
}
K = 2
Dare_topic[2,2] = mean(sapply(1:K, function(k) {
    topic_coherence(top.words[,k], topic.word,Dare$vocab)
}))

for (i in 1:length(Daretest$Z)) {
    Daretest$Z[[i]] = rep(1, length(Daretest$Z[[i]]))
}
top.words = TableWord(Daretest$Z, 1, Dare$text, Dare$vocab)
topic.word = matrix(0, nrow = length(Daretest$Z), ncol = length(Dare$vocab))
rownames(topic.word) = 1:length(Dare$text)
colnames(topic.word) = Dare$vocab
for (d in seq(along = Dare$text)) {
    topic.word[d, ] = tabulate(Dare$text[[d]], length(Dare$vocab))
}
K = 1
Dare_topic[2,1] = mean(sapply(1:K, function(k) {
    topic_coherence(top.words[,k], topic.word,Dare$vocab)
}))





load('~/Desktop/IPTM/paper/code/HPC/topic coherence/Daretest_LDA_K50.RData')
Daretest = Daretest_LDA
top.words = TableWord(Daretest$Z, 50, Dare$text, Dare$vocab)
topic.word = matrix(0, nrow = length(Daretest$Z), ncol = length(Dare$vocab))
rownames(topic.word) = 1:length(Dare$text)
colnames(topic.word) = Dare$vocab
for (d in seq(along = Dare$text)) {
	topic.word[d, ] = tabulate(Dare$text[[d]], length(Dare$vocab))
}
K =50
Dare_topic[3,8] = mean(sapply(1:K, function(k) {
topic_coherence(top.words[,k], topic.word,Dare$vocab)	
}))

load('~/Desktop/IPTM/paper/code/HPC/topic coherence/Daretest_LDA_K40.RData')
Daretest = Daretest_LDA
top.words = TableWord(Daretest$Z, 40, Dare$text, Dare$vocab)
topic.word = matrix(0, nrow = length(Daretest$Z), ncol = length(Dare$vocab))
rownames(topic.word) = 1:length(Dare$text)
colnames(topic.word) = Dare$vocab
for (d in seq(along = Dare$text)) {
	topic.word[d, ] = tabulate(Dare$text[[d]], length(Dare$vocab))
}
K =40
Dare_topic[3,7] = mean(sapply(1:K, function(k) {
topic_coherence(top.words[,k], topic.word,Dare$vocab)	
}))

load('~/Desktop/IPTM/paper/code/HPC/topic coherence/Daretest_LDA_K30.RData')
Daretest = Daretest_LDA
top.words = TableWord(Daretest$Z, 30, Dare$text, Dare$vocab)
topic.word = matrix(0, nrow = length(Daretest$Z), ncol = length(Dare$vocab))
rownames(topic.word) = 1:length(Dare$text)
colnames(topic.word) = Dare$vocab
for (d in seq(along = Dare$text)) {
	topic.word[d, ] = tabulate(Dare$text[[d]], length(Dare$vocab))
}
K =30
Dare_topic[3,6] = mean(sapply(1:K, function(k) {
topic_coherence(top.words[,k], topic.word,Dare$vocab)	
}))


load('~/Desktop/IPTM/paper/code/HPC/topic coherence/Daretest_LDA_K20.RData')
Daretest = Daretest_LDA
top.words = TableWord(Daretest$Z, 20, Dare$text, Dare$vocab)
topic.word = matrix(0, nrow = length(Daretest$Z), ncol = length(Dare$vocab))
rownames(topic.word) = 1:length(Dare$text)
colnames(topic.word) = Dare$vocab
for (d in seq(along = Dare$text)) {
	topic.word[d, ] = tabulate(Dare$text[[d]], length(Dare$vocab))
}
K =20
Dare_topic[3,5] = mean(sapply(1:K, function(k) {
topic_coherence(top.words[,k], topic.word,Dare$vocab)	
}))


load('~/Desktop/IPTM/paper/code/HPC/topic coherence/Daretest_LDA_K10.RData')
Daretest = Daretest_LDA
top.words = TableWord(Daretest$Z, 10, Dare$text, Dare$vocab)
topic.word = matrix(0, nrow = length(Daretest$Z), ncol = length(Dare$vocab))
rownames(topic.word) = 1:length(Dare$text)
colnames(topic.word) = Dare$vocab
for (d in seq(along = Dare$text)) {
	topic.word[d, ] = tabulate(Dare$text[[d]], length(Dare$vocab))
}
K = 10
Dare_topic[3,4] = mean(sapply(1:K, function(k) {
topic_coherence(top.words[,k], topic.word,Dare$vocab)	
}))



load('~/Desktop/IPTM/paper/code/HPC/topic coherence/Daretest_LDA_K5.RData')
Daretest = Daretest_LDA
top.words = TableWord(Daretest$Z, 5, Dare$text, Dare$vocab)
topic.word = matrix(0, nrow = length(Daretest$Z), ncol = length(Dare$vocab))
rownames(topic.word) = 1:length(Dare$text)
colnames(topic.word) = Dare$vocab
for (d in seq(along = Dare$text)) {
	topic.word[d, ] = tabulate(Dare$text[[d]], length(Dare$vocab))
}
K = 5
Dare_topic[3,3] = mean(sapply(1:K, function(k) {
topic_coherence(top.words[,k], topic.word,Dare$vocab)	
}))


load('~/Desktop/IPTM/paper/code/HPC/topic coherence/Daretest_LDA_K2.RData')
Daretest = Daretest_LDA
top.words = TableWord(Daretest$Z, 2, Dare$text, Dare$vocab)
topic.word = matrix(0, nrow = length(Daretest$Z), ncol = length(Dare$vocab))
rownames(topic.word) = 1:length(Dare$text)
colnames(topic.word) = Dare$vocab
for (d in seq(along = Dare$text)) {
	topic.word[d, ] = tabulate(Dare$text[[d]], length(Dare$vocab))
}
K = 2
Dare_topic[3,2] = mean(sapply(1:K, function(k) {
topic_coherence(top.words[,k], topic.word,Dare$vocab)	
}))

Dare_topic[3,1] = Dare_topic[1,1]



load('~/Desktop/IPTM/paper/code/HPC/topic coherence/sDaretest_LDA_K50.RData')
Daretest = Daretest_LDA
top.words = TableWord(Daretest$Z, 50, Dare$text, Dare$vocab)
topic.word = matrix(0, nrow = length(Daretest$Z), ncol = length(Dare$vocab))
rownames(topic.word) = 1:length(Dare$text)
colnames(topic.word) = Dare$vocab
for (d in seq(along = Dare$text)) {
    topic.word[d, ] = tabulate(Dare$text[[d]], length(Dare$vocab))
}
K =50
Dare_topic[4,8] = mean(sapply(1:K, function(k) {
    topic_coherence(top.words[,k], topic.word,Dare$vocab)
}))

load('~/Desktop/IPTM/paper/code/HPC/topic coherence/sDaretest_LDA_K40.RData')
Daretest = Daretest_LDA
top.words = TableWord(Daretest$Z, 40, Dare$text, Dare$vocab)
topic.word = matrix(0, nrow = length(Daretest$Z), ncol = length(Dare$vocab))
rownames(topic.word) = 1:length(Dare$text)
colnames(topic.word) = Dare$vocab
for (d in seq(along = Dare$text)) {
    topic.word[d, ] = tabulate(Dare$text[[d]], length(Dare$vocab))
}
K =40
Dare_topic[4,7] = mean(sapply(1:K, function(k) {
    topic_coherence(top.words[,k], topic.word,Dare$vocab)
}))

load('~/Desktop/IPTM/paper/code/HPC/topic coherence/sDaretest_LDA_K30.RData')
Daretest = Daretest_LDA
top.words = TableWord(Daretest$Z, 30, Dare$text, Dare$vocab)
topic.word = matrix(0, nrow = length(Daretest$Z), ncol = length(Dare$vocab))
rownames(topic.word) = 1:length(Dare$text)
colnames(topic.word) = Dare$vocab
for (d in seq(along = Dare$text)) {
    topic.word[d, ] = tabulate(Dare$text[[d]], length(Dare$vocab))
}
K =30
Dare_topic[4,6] = mean(sapply(1:K, function(k) {
    topic_coherence(top.words[,k], topic.word,Dare$vocab)
}))


load('~/Desktop/IPTM/paper/code/HPC/topic coherence/sDaretest_LDA_K20.RData')
Daretest = Daretest_LDA
top.words = TableWord(Daretest$Z, 20, Dare$text, Dare$vocab)
topic.word = matrix(0, nrow = length(Daretest$Z), ncol = length(Dare$vocab))
rownames(topic.word) = 1:length(Dare$text)
colnames(topic.word) = Dare$vocab
for (d in seq(along = Dare$text)) {
    topic.word[d, ] = tabulate(Dare$text[[d]], length(Dare$vocab))
}
K =20
Dare_topic[4,5] = mean(sapply(1:K, function(k) {
    topic_coherence(top.words[,k], topic.word,Dare$vocab)
}))


load('~/Desktop/IPTM/paper/code/HPC/topic coherence/sDaretest_LDA_K10.RData')
Daretest = Daretest_LDA
top.words = TableWord(Daretest$Z, 10, Dare$text, Dare$vocab)
topic.word = matrix(0, nrow = length(Daretest$Z), ncol = length(Dare$vocab))
rownames(topic.word) = 1:length(Dare$text)
colnames(topic.word) = Dare$vocab
for (d in seq(along = Dare$text)) {
    topic.word[d, ] = tabulate(Dare$text[[d]], length(Dare$vocab))
}
K = 10
Dare_topic[4,4] = mean(sapply(1:K, function(k) {
    topic_coherence(top.words[,k], topic.word,Dare$vocab)
}))



load('~/Desktop/IPTM/paper/code/HPC/topic coherence/sDaretest_LDA_K5.RData')
Daretest = Daretest_LDA
top.words = TableWord(Daretest$Z, 5, Dare$text, Dare$vocab)
topic.word = matrix(0, nrow = length(Daretest$Z), ncol = length(Dare$vocab))
rownames(topic.word) = 1:length(Dare$text)
colnames(topic.word) = Dare$vocab
for (d in seq(along = Dare$text)) {
    topic.word[d, ] = tabulate(Dare$text[[d]], length(Dare$vocab))
}
K = 5
Dare_topic[4,3] = mean(sapply(1:K, function(k) {
    topic_coherence(top.words[,k], topic.word,Dare$vocab)
}))


load('~/Desktop/IPTM/paper/code/HPC/topic coherence/sDaretest_LDA_K2.RData')
Daretest = Daretest_LDA
top.words = TableWord(Daretest$Z, 2, Dare$text, Dare$vocab)
topic.word = matrix(0, nrow = length(Daretest$Z), ncol = length(Dare$vocab))
rownames(topic.word) = 1:length(Dare$text)
colnames(topic.word) = Dare$vocab
for (d in seq(along = Dare$text)) {
    topic.word[d, ] = tabulate(Dare$text[[d]], length(Dare$vocab))
}
K = 2
Dare_topic[4,2] = mean(sapply(1:K, function(k) {
    topic_coherence(top.words[,k], topic.word,Dare$vocab)
}))

Dare_topic[4,1] = Dare_topic[1,1]





#############

plot(colnames(Dare_topic), colMeans(Dare_topic[3:4,]), type = 'l', col = 'red', xlab = "Number of Topics", lty = 2, ylab = "Topic Coherence", ylim = c(min(Dare_topic), max(Dare_topic)))
lines(colnames(Dare_topic), colMeans(Dare_topic[1:2,]), lty = 3,col = "blue")
legend(40, -410, pch = 21, legend = c("IPTM","LDA"), col = c("blue","red"))

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
