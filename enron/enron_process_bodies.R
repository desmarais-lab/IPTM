library(tidyverse)
library(tidytext)
library(tm)
library(stringr)
enron_emails <- read_csv('emails.csv')

enron_emails <- enron_emails  %>%
  filter(str_detect(file, '/_sent_mail/'))
  
enron_emails <- enron_emails %>%
  mutate(message = message %>% str_replace_all('\r', '')) %>%
  mutate(header = str_sub(message, end = str_locate(message, '\n\n')[,1] -1)) %>%
  mutate(body = str_sub(message, start = str_locate(message, '\n\n')[,2] + 1) %>% 
           str_replace_all('\n|\t', ' ') %>% 
           str_replace_all('---Original Message .*', 'FORWARDED_MESSAGE') %>% 
           str_replace_all('--- Forwarded by .*', 'FORWARDED_MESSAGE') %>%
           str_replace_all('From: .*', 'FORWARDED_MESSAGE') %>%
           str_replace('To:.*', 'FORWARDED_MESSAGE') %>%
           str_replace_all('\\S*@\\S*', 'EMAIL_ADDRESS')) 
           
enron_emails <- enron_emails %>%
  mutate(date = str_extract(header, 'Date:.*') %>% 
           str_replace('Date: ', '') %>% 
           str_replace('.+, ', '') %>% 
           strptime(format = '%d %b %Y %H:%M:%S %z') %>%
           as.POSIXct()) %>%
  mutate(from = str_extract(header, 'From:.*') %>% 
           str_replace('From: ', '')) %>%
  mutate(to = header %>% str_replace_all('\n|\t', ' ') %>%
           str_extract('To:.*Subject:') %>%
           str_replace_all('To: |Subject:', '')) %>%
  mutate(subject = str_extract(header, 'Subject:.*') %>% 
           str_replace('Subject: ', '')) %>%
  mutate(xfrom = str_extract(header, 'X-From:.*') %>% 
           str_replace('X-From: ', '')) %>%
  mutate(xto = str_extract(header, 'X-To:.*') %>% 
           str_replace('X-To: ', '')) %>%
  mutate(xcc = str_extract(header, 'X-cc:.*') %>% 
           str_replace('X-cc: ', '')) %>%
  mutate(xbcc = str_extract(message, 'X-bcc:.*') %>% 
           str_replace('X-bcc: ', '')) %>%
  arrange(date)

dtm.control = list(
  tolower = T,
  removePunctuation = T,
  removeNumbers = T,
  stopwords = stopwords('english'),
  weighting = weightTf,
  seed = 0
)

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


email = list()
for (i in 1:30237) {
	email[[i]] = list()
	email[[i]]$from = enron_emails$from[i]
	email[[i]]$to = gsub(" ", "", strsplit(enron_emails$to[i], ", ")[[1]])
	email[[i]]$date = enron_emails$date[i]
	rawtext = rbind(docterm[which(docterm[,1] ==i),-1], docterm2[which(docterm2[,1] ==i),-1])
	if (nrow(rawtext) >0) {
	 email[[i]]$text = unlist(lapply(1:nrow(rawtext), function(x) {rep(rawtext[x,1], rawtext[x,2])}))
	} else {
	email[[i]]$text = integer(0)
	}
}
enron_raw = email
save(enron_raw, file = "enron_raw.RData")

#clean text once more
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

for (i in 1:30237) {
  if (length(email[[i]]$text) >0){
  words = Clean_String(email[[i]]$text)
  words = words[nchar(words) >= 3]
  words = tm_map(Corpus(VectorSource(words)), removeWords, c(stopwords('english'), reservedM))$content
  words = words[!words %in% ""]
 email[[i]]$text = words
}
}
enron_raw = email
save(enron_raw, file = "enron_raw.RData")





