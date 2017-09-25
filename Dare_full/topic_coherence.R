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
    table.word = top.topic.words(topic.word, num.words = top, by.score =TRUE)
    colnames(table.word) = names(top.topic)
  return(table.word)
}

Dare_topic = matrix(0, 3, 8)
colnames(Dare_topic) = as.numeric(c("1", "2", "5", "10", "25", "50", "75", "100"))

setwd("~/Desktop/Dare_full")
for (nIP in 1:3) {
	iter = 1
	for (K in c(2, 5, 10, 25, 50, 75, 100)){
		filename = paste0("Dare_full_",nIP,"_",K,"_ver",1,".RData")
		load(filename)
		top.words = TableWord(Daretest$Z, K, Dare$text, Dare$vocab)
        topic.word = matrix(0, nrow = length(Daretest$Z), ncol = length(Dare$vocab))
		rownames(topic.word) = 1:length(Dare$text)
		colnames(topic.word) = Dare$vocab
		for (d in seq(along = Dare$text)) {
		topic.word[d, ] = tabulate(Dare$text[[d]], length(Dare$vocab))
		}
		iter = iter + 1
		Dare_topic[nIP, iter] = mean(sapply(1:K, function(k) {
		topic_coherence(top.words[,k], topic.word,Dare$vocab)	
		}))
	}
}

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
Dare_topic[,1] = mean(sapply(1:K, function(k) {
    topic_coherence(top.words[,k], topic.word,Dare$vocab)
}))


plot(colnames(Dare_topic), Dare_topic[1,], type = 'l', col = 'red', xlab = "Number of Topics", lty = 2, ylab = "Topic Coherence", ylim = c(min(Dare_topic), max(Dare_topic)), main = "Dare")
lines(colnames(Dare_topic), Dare_topic[2,], lty = 3,col = "blue")
lines(colnames(Dare_topic), Dare_topic[3,], lty = 4, col = "green")
legend(65, -405, pch = 21, legend = c("IPTM with C = 3","IPTM with C = 2","LDA"), col = c("green","blue","red"))
