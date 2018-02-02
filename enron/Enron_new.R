library(tm)
library(stringr)
load("enron_raw.RData")
email = enron_raw
sender = sapply(1:length(email), function(d) {email[[d]]$from})

#writers_sent_over_100
newsender = which(table(sender) >= 200)
email = email[which(sender %in% names(newsender))]

#only_include_receivers_over_100_received_emails
receiver = unlist(lapply(1:length(email), function(d) {email[[d]]$to}))
newreceiver = which(table(receiver) >= 200)
receiverlist = lapply(1:length(email), function(d) {email[[d]]$to})
include = sapply(1:length(email), function(d) {sum(receiverlist[[d]] %in% names(newreceiver))==length(receiverlist[[d]])})
email = email[include]

#good to be used for other purpose
#save(email, file = "email.RData")
#can also separate edge info and text
#edge = list()
#text = list()
#for (d in 1:length(email)) {
#    edge[[d]] = list()
#    edge[[d]]$sender = email[[d]][[1]]
#    edge[[d]]$receiver = email[[d]][[2]]
#    edge[[d]]$date = as.POSIXct(email[[d]][[3]])
#    text[[d]] = email[[d]][[4]]
#}
#save(edge, file = "edge.RData")
#save(text, file = "text.RData")

##############################################
#finalize the dataset for IPTM (e.g. use node id's and unixtime)
sender = lapply(email, function(i){i[[1]]})
receiver = lapply(email, function(i){i[[2]]})
node = union(unique(unlist(sender)), unique(unlist(receiver)))  
text = lapply(email, function(i){i[[4]]})
vocab = unique(unlist(text))
#dominated by these words (ridiculously too many)---address, email, forwarded, message---so delete those
vocab = vocab[!vocab %in% c("address", "email", "forwarded", "message")]

edge = list()
text2 = list()
for (d in 1:length(email)) {
  edge[[d]] = list()
  edge[[d]]$sender = which(node == email[[d]][[1]])           #make node id's starting from 1 to A 
  edge[[d]]$receiver = which(node %in% email[[d]][[2]])     #make node id's starting from 1 to A
  edge[[d]]$timestamp = as.numeric(as.POSIXct(email[[d]][[3]]))
  text2[[d]] = which(vocab %in% email[[d]][[4]])
}
Enron = list(edge = edge, node = 1:length(node), text = text2, vocab = vocab)
save(Enron, file = "Enron.RData")
