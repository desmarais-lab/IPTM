library(stringr)

Clean_String <- function(string){
  # Lowercase
  temp <- tolower(string)
  #' Remove everything that is not a number or letter (may want to keep more 
  #' stuff in your actual analyses). 
  temp <- stringr::str_replace_all(temp,"[^a-zA-Z\\s]", " ")
  # Shrink down to just one white space
  temp <- stringr::str_replace_all(temp,"[\\s]+", " ")
  # Split it
  temp <- stringr::str_split(temp, " ")[[1]]
  # Get rid of trailing "" if necessary
  indexes <- which(temp == "")
  if(length(indexes) > 0){
    temp <- temp[-indexes]
  } 
  return(temp)
}

library(tm)
load("/Users/bomin8319/Desktop/IPTM/enron/enron_raw.RData")
Enron = list()
edge = list()
text = list()
reservedM<-c("re", "a", "b", "c", "d", "e", "g", "h", "q", "s", "t", "v","w","x","z",
             "l", "p", "r", "f", "n", "y", "u", "k", "j", "i", "m", "o")
for (d in 1:dim(enron)[1]) {
  receivers =as.numeric(strsplit(Corpus(VectorSource(as.character(enron[d,2])))$content, ",")[[1]])
  edge[[d]] = list(sender = enron[d, 1], receiver = receivers, timestamp = enron[d, 3] / 3600)
  words = Clean_String(enron[d, 4])
  words = words[nchar(words) >= 3]
  words = tm_map(Corpus(VectorSource(words)), removeWords, c(stopwords('english'), reservedM))$content
  words = words[!words %in% ""]
  text[[d]] = words
}
vocab = unique(unlist(text))
for (d in 1:dim(enron)[1]) {
  text[[d]] = which(vocab %in% text[[d]])
}  
sender = lapply(edge, function(i){i[[1]]})

#top_ten_writers
newsender = which(tabulate(unlist(sender)) > 100)
enron = enron[enron[,1] %in% newsender,]
edge = list()
text = list()
for (d in 1:dim(enron)[1]) {
  receivers =as.numeric(strsplit(Corpus(VectorSource(as.character(enron[d,2])))$content, ",")[[1]])
  edge[[d]] = list(sender = enron[d, 1], receiver = receivers, timestamp = enron[d, 3] / 3600)
  words = Clean_String(enron[d, 4])
  words = words[nchar(words) >= 3]
  words = tm_map(Corpus(VectorSource(words)), removeWords, c(stopwords('english'), reservedM))$content
  words = words[!words %in% ""]
  text[[d]] = words
}
vocab = unique(unlist(text))
for (d in 1:dim(enron)[1]) {
  text[[d]] = which(vocab %in% text[[d]])
} 

#only include receivers > 100 received emails
receiver = lapply(edge, function(i){i[[2]]})

newreceiver = which(tabulate(unlist(receiver)) > 100)
oldedge = edge
oldtext = text
edge = list()
text = list()
it = 1
for (d in 1:dim(enron)[1]) {
  if (sum(oldedge[[d]]$receiver %in% newreceiver)==length(oldedge[[d]]$receiver)) {
  receivers =as.numeric(strsplit(Corpus(VectorSource(as.character(enron[d,2])))$content, ",")[[1]])
  edge[[it]] = list(sender = enron[d, 1], receiver = receivers, timestamp = enron[d, 3] / 3600)
  words = Clean_String(enron[d, 4])
  words = words[nchar(words) >= 3]
  words = tm_map(Corpus(VectorSource(words)), removeWords, c(stopwords('english'), reservedM))$content
  words = words[!words %in% ""]
  text[[it]] = words
  it = it + 1
  }
}
vocab = unique(unlist(text))
delete = which(sapply(text, function(x){length(x)}) ==0 |sapply(text, function(x){length(x)}) ==6519 )
edge2 = edge[-delete]
text2 = text[-delete]
edge = edge2
text = text2
sender = lapply(edge, function(i){i[[1]]})
receiver = lapply(edge, function(i){i[[2]]})
node = union(unique(unlist(sender)), unique(unlist(receiver)))
vocab = unique(unlist(text))
text2 = list()
for (d in 1:length(edge)) {
  edge[[d]]$sender = which(node == edge[[d]]$sender)
  edge[[d]]$receiver = which(node %in% edge[[d]]$receiver)
  if (d > 1 && edge[[d]]$timestamp == edge[[d-1]]$timestamp) {
  edge[[d]]$timestamp = edge[[d]]$timestamp + exp(-20)
  }
   text2[[d]] = which(vocab %in% text[[d]])
  }

Enron = list(edge = edge, node = 1:length(node), text = text2, vocab = vocab)
save(Enron, file = "Enron.RData")
