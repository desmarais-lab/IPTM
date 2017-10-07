library(lda)
library(SpeedReader)
setwd("~/Desktop/IPTM/full")
load('Darenew.RData')
attach(Dare)
Dare$node = 1:nrow(Dare$node)
Dare$edge = Dare$edge[-which(sapply(Dare$text, function(d){length(d)})==0)]
Dare$text = Dare$text[-which(sapply(Dare$text, function(d){length(d)})==0)]
mintime = Dare$edge[[1]][[3]]
for (n in 1:length(Dare$edge)){
  Dare$edge[[n]][3] = (Dare$edge[[n]][[3]] - mintime) / 3600
}
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



Dare_topic = matrix(0, 3, 8)
colnames(Dare_topic) = as.numeric(c("1", "2", "5", "10", "25", "50", "75", "100"))

for (i in 1:5) {
for (nIP in 1:3) {
	iter = 1
	for (K in c(2, 5, 10, 25, 50, 75, 100)){
		filename = paste0("Dare_full_",nIP,"_",K,"_ver",i,".RData")
		load(filename)
		top.words = TableWord(Daretest$Z, K, Dare$text, Dare$vocab)
        topic.word = matrix(0, nrow = length(Daretest$Z), ncol = length(Dare$vocab))
		rownames(topic.word) = 1:length(Dare$text)
		colnames(topic.word) = Dare$vocab
		for (d in seq(along = Dare$text)) {
		topic.word[d, ] = tabulate(Dare$text[[d]], length(Dare$vocab))
		}
		iter = iter + 1
		Dare_topic[nIP, iter] = Dare_topic[nIP, iter] + mean(sapply(1:K, function(k) {
		topic_coherence(top.words[,k], topic.word,Dare$vocab)	
		}))
	}
}
}

Dare_topic = Dare_topic / 5

for (i in 1:length(Daretest$Z)) {
    Daretest$Z[[i]] = rep(1, length(Daretest$Z[[i]]))
}
top.words = TableWord2(Daretest$Z, 1, Dare$text, Dare$vocab)
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

###################
load('Enron.RData')
Enron_topic = matrix(0, 3, 8)
colnames(Enron_topic) = as.numeric(c("1", "2", "5", "10", "25", "50", "75", "100"))

for (i in c(1, 5)) {
for (nIP in 1:3) {
	iter = 1
	for (K in c(2, 5, 10, 25, 50, 75)){
		filename = paste0("Enron_full_",nIP,"_",K,"_ver",i,".RData")
		load(filename)
		top.words = TableWord(Enrontest$Z, K, Enron$text, Enron$vocab)
        topic.word = matrix(0, nrow = length(Enrontest$Z), ncol = length(Enron$vocab))
		rownames(topic.word) = 1:length(Enron$text)
		colnames(topic.word) = Enron$vocab
		for (d in seq(along = Enron$text)) {
		topic.word[d, ] = tabulate(Enron$text[[d]], length(Enron$vocab))
		}
		iter = iter + 1
		Enron_topic[nIP, iter] = Enron_topic[nIP, iter] + mean(sapply(1:K, function(k) {
		topic_coherence(top.words[,k], topic.word,Enron$vocab)	
		}))
	}
}
}

Enron_topic = Enron_topic /2

for (i in 1:length(Enrontest$Z)) {
   Enrontest$Z[[i]] = rep(1, length(Enrontest$Z[[i]]))
}
top.words = TableWord2(Enrontest$Z, 1, Enron$text, Enron$vocab)
topic.word = matrix(0, nrow = length(Enrontest$Z), ncol = length(Enron$vocab))
rownames(topic.word) = 1:length(Enron$text)
colnames(topic.word) = Enron$vocab
for (d in seq(along =Enron$text)) {
    topic.word[d, ] = tabulate(Enron$text[[d]], length(Enron$vocab))
}
K = 1
Enron_topic[,1] = mean(sapply(1:K, function(k) {
    topic_coherence(top.words[,k], topic.word,Enron$vocab)
}))

################

par(mfrow = c(1,2))
plot(colnames(Dare_topic), Dare_topic[1,], type = 'b', pch = 19, col = 'blue', xlab = "Number of Topics", lty = 2, ylab = "Topic Coherence", ylim = c(min(Dare_topic), max(Dare_topic)), main = "Dare")
lines(colnames(Dare_topic), Dare_topic[2,], lty = 3, pch = 19, type = 'b', col = "green")
lines(colnames(Dare_topic), Dare_topic[3,], lty = 4, pch = 19, type = 'b',col = "purple")
legend(65, -405, pch = 19, legend = c("IPTM with C = 3","IPTM with C = 2","LDA"), col = c("purple","green","blue"))

Enron_topic = Enron_topic[, -8]
plot(colnames(Enron_topic), Enron_topic[1,], type = 'b', pch = 19, col = 'blue', xlab = "Number of Topics", lty = 2, ylab = "Topic Coherence", ylim = c(min(Enron_topic), max(Enron_topic)), main = "Enron")
lines(colnames(Enron_topic), Enron_topic[2,], lty = 3, pch = 19, type = 'b', col = "green")
lines(colnames(Enron_topic), Enron_topic[3,], lty = 4, pch = 19, type = 'b',col = "purple")
legend(47, -550, pch = 19, legend = c("IPTM with C = 3","IPTM with C = 2","LDA"), col = c("purple","green","blue"))
