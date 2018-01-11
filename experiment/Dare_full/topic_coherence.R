library(lda)
library(SpeedReader)
setwd("~/Desktop/IPTM/experiment/Dare_full")
load('Darenew.RData')
attach(Dare)
Dare$node = 1:nrow(Dare$node)
Dare$edge = lapply(Dare$edge, function(x){x[1:3]})

TableWord = function(Zchain, K, textlist, vocabulary, top = 20) {
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
    top.topic = topic.dist[, order(topic.dist[1,], decreasing = TRUE)]
    all.word = unlist(Zsummary)
    for (i in seq(along = all.word)){
      matchWZ = which(colnames(topic.word) == names(all.word[i]))
      topic.word[all.word[i], matchWZ] = topic.word[all.word[i], matchWZ] + 1
    }
    table.word = top.topic.words(topic.word, num.words = top, by.score =TRUE)
    colnames(table.word) = names(top.topic)
  return(table.word)
}


TableWord2 = function(Zchain, K, textlist, vocabulary, top = 20) {
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
    top.topic = topic.dist[, order(topic.dist[1,], decreasing = TRUE)]
    all.word = unlist(Zsummary)
    for (i in seq(along = all.word)){
        matchWZ = which(colnames(topic.word) == names(all.word[i]))
        topic.word[all.word[i], matchWZ] = topic.word[all.word[i], matchWZ] + 1
    }
    table.word = top.topic.words(topic.word, num.words = top, by.score =FALSE)
    colnames(table.word) = names(top.topic)
    return(table.word)
}



Dare_topic = matrix(0, 3, 7)
colnames(Dare_topic) = as.numeric(c( "1", "5", "10", "20", "30", "40", "50"))

for (i in 1:5) {
for (nIP in 1:3) {
	iter = 1
	for (K in c(5, 10, 20, 30, 40, 50)){
		filename = paste0("Dare_full_",nIP,"_",K,"_ver",i,".RData")
		load(filename)
		top.words = TableWord(Daretest$z, K, Dare$text, Dare$vocab)
        topic.word = matrix(0, nrow = length(Daretest$z), ncol = length(Dare$vocab))
		rownames(topic.word) = 1:length(Dare$text)
		colnames(topic.word) = Dare$vocab
		for (d in seq(along = Dare$text)) {
			if (length(Dare$text[[d]])>0) {
		topic.word[d, ] = tabulate(Dare$text[[d]], length(Dare$vocab))
		}
		}
		iter = iter + 1
		Dare_topic[nIP, iter] = Dare_topic[nIP, iter] + mean(sapply(1:K, function(k) {
		topic_coherence(top.words[,k], topic.word,Dare$vocab)	
		}))
	}
}
}


Dare_topic = Dare_topic / 5

for (i in 1:length(Daretest$z)) {
    Daretest$z[[i]] = rep(1, length(Daretest$z[[i]]))
}
top.words = TableWord2(Daretest$z, 1, Dare$text, Dare$vocab)
topic.word = matrix(0, nrow = length(Daretest$z), ncol = length(Dare$vocab))
rownames(topic.word) = 1:length(Dare$text)
colnames(topic.word) = Dare$vocab
for (d in seq(along = Dare$text)) {
	if (length(Dare$text[[d]]) > 0){
    topic.word[d, ] = tabulate(Dare$text[[d]], length(Dare$vocab))
}
}
K = 1
Dare_topic[,1] = mean(sapply(1:K, function(k) {
    topic_coherence(top.words[,k], topic.word,Dare$vocab)
}))
library(ggplot2)
library(gridExtra)


Dare_topic_new = data.frame(Coherence = c(t(Dare_topic)), nIP = as.factor(c(sapply(1:3, function(x) rep(x, 7)))), K =as.factor(rep(c(1,5, 10, 20, 30, 40, 50), 3)))
ggplot(Dare_topic_new, aes(K, Coherence, col = nIP))+geom_line(aes(group=nIP)) + geom_point(size = 1)+
ggtitle("Dare")

###################
setwd("~/Desktop/IPTM/experiment/Enron_full")
load('Enron.RData')
Enron_topic = matrix(0, 3, 7)
colnames(Enron_topic) = as.numeric(c("1", "5", "10", "25", "50", "75", "100"))

for (i in 1:5) {
for (nIP in 1:2) {
	iter = 1
	for (K in c(5, 10, 25, 50, 75, 100)){
		filename = paste0("Enron_full_",nIP,"_",K,"_ver",i,".RData")
		load(filename)
		top.words = TableWord(Enrontest$z, K, Enron$text, Enron$vocab)
        topic.word = matrix(0, nrow = length(Enrontest$z), ncol = length(Enron$vocab))
		rownames(topic.word) = 1:length(Enron$text)
		colnames(topic.word) = Enron$vocab
		for (d in seq(along = Enron$text)) {
			if (length(Enron$text[[d]] > 0)) {
			topic.word[d, ] = tabulate(Enron$text[[d]], length(Enron$vocab))
		}
		}
		iter = iter + 1
		Enron_topic[nIP, iter] = Enron_topic[nIP, iter] + mean(sapply(1:K, function(k) {
		topic_coherence(top.words[,k], topic.word,Enron$vocab)	
		}))
	}
}
}
Enron_topic = Enron_topic / 5

for (i in 1:length(Enrontest$z)) {
   Enrontest$z[[i]] = rep(1, length(Enrontest$z[[i]]))
}
top.words = TableWord2(Enrontest$z, 1, Enron$text, Enron$vocab)
topic.word = matrix(0, nrow = length(Enrontest$z), ncol = length(Enron$vocab))
rownames(topic.word) = 1:length(Enron$text)
colnames(topic.word) = Enron$vocab
for (d in seq(along =Enron$text)) {
	if (length(Enron$text[[d]] > 0)) {
    topic.word[d, ] = tabulate(Enron$text[[d]], length(Enron$vocab))
}
}
K = 1
Enron_topic[,1] = mean(sapply(1:K, function(k) {
    topic_coherence(top.words[,k], topic.word,Enron$vocab)
}))

################



Enron_topic_new = data.frame(Coherence = c(t(Enron_topic)), nIP = as.factor(c(sapply(1:3, function(x) rep(x, 7)))), K =as.factor(rep(c(1,5, 10, 25, 50, 75, 100), 3)))
ggplot(Enron_topic_new, aes(K, Coherence, col = nIP))+geom_line(aes(group=nIP)) + geom_point(size = 1)+
ggtitle("Enron")

