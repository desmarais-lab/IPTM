enron = data.frame(sender = enron_emails$from, receiver = enron_emails$to)

enron_dtm <- Corpus(VectorSource(enron_emails$body)) %>%
  DocumentTermMatrix( control = dtm.control) %>%
  removeSparseTerms(.999) %>%
  .[rowSums(as.matrix(.))>0,]
 
docterm = as.data.frame(tidy(enron_dtm)) 
docterm$document = as.numeric(docterm$document)
docterm$count = as.numeric(docterm$count)


enron_dtm2 <- Corpus(VectorSource(enron_emails$subject)) %>%
  DocumentTermMatrix( control = dtm.control) %>%
  removeSparseTerms(.999) %>%
  .[rowSums(as.matrix(.))>0,]
 
docterm2 = as.data.frame(tidy(enron_dtm2)) 
docterm2$document = as.numeric(docterm2$document)
docterm2$count = as.numeric(docterm2$count)

subject = matrix(NA, nrow = 30237, ncol = 1)
body = matrix(NA, nrow = 30237, ncol = 1)

for (i in 1:30237) {
	rawsub = docterm2[which(docterm2[,1] ==i),-1]
	rawtext = docterm[which(docterm[,1] ==i),-1]
	if (nrow(rawsub) >0) {
	subject[i, ] = paste(unlist(lapply(1:nrow(rawsub), function(x) {rep(rawsub[x,1], rawsub[x,2])})), collapse = " ")
	} 
	if (nrow(rawtext) >0) {
	body[i, ] = paste(unlist(lapply(1:nrow(rawtext), function(x) {rep(rawtext[x,1], rawtext[x,2])})), collapse = " ")
	} 
}

reservedM<-c("re", "a", "b", "c", "d", "e", "g", "h", "q", "s", "t", "v","w","x","z",
             "l", "p", "r", "f", "n", "y", "u", "k", "j", "i", "m", "o")

Clean_String <- function(string){
  # Lowercase
  temp <- tolower(string)
  #' Remove everything that is not a number or letter (may want to keep more 
  #' stuff in your actual analyses). 
  temp <- stringr::str_replace_all(temp,"[^a-zA-Z\\s]", " ")
  # Shrink down to just one white space
  temp <- stringr::str_replace_all(temp,"[\\s]+", " ")
  # Split it
  temp <- unlist(stringr::str_split(temp, " "))
  # Get rid of trailing "" if necessary
  indexes <- which(temp == "")
  if(length(indexes) > 0){
    temp <- temp[-indexes]
  } 
  return(temp)
}

enron = data.frame(sender = enron_emails$from, receiver = enron_emails$to, subject = subject, body = body)

write.csv(enron, file = "enron.csv")