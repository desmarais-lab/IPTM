library(ggplot2)
library(lubridate)
library(anytime)
library(gridExtra)
load('~/Desktop/IPTM/experiment/Dare_PPE/Darenew.RData')
senders = sapply(1:2247, function(x) Dare$edge[[x]]$sender)
receivers = sapply(1:2247, function(x) length(Dare$edge[[x]]$receiver))
time = c(Dare$edge[[1]]$unixtime, sapply(2:2247, function(x) Dare$edge[[x]]$unixtime-Dare$edge[[x-1]]$unixtime))/3600

day = sapply(1:2247, function(x) format(anytime(Dare$edge[[x]]$unixtime, tz = "America/New_York"), format = "%m/%d"))
p = list()
#emails sent by County Manager (node 1) --- total 350
byday = rep(0, length(unique(day)))
names(byday)  = unique(day)
i = 1
for (d in unique(day)) {
	docs = which(day == d)
	byday[i] = sum(senders[docs] == 1)
	i = i + 1
}

byday2 = rep(0, length(unique(day)))
names(byday2)  = unique(day)
i = 1
for (d in unique(day)) {
	docs = which(day == d)
	byday2[i] = sum(receivers[docs][senders[docs] == 1])
	i = i + 1
}

words = unlist(lapply(which(senders==1), function(x) Dare$text[[x]]))
Dare$vocab[which(tabulate(words) %in% sort(tabulate(words), decreasing = TRUE)[1:15])]

SummaryTable <- data.frame(Words1 = Dare$vocab[which(tabulate(words) %in% sort(tabulate(words), decreasing = TRUE)[1:15])][1:15])

data = data.frame(day = as.factor(unique(day)), send = byday, num_receiver = byday2)
p[[1]]<-ggplot(data=data,aes(x=day, y=send)) + geom_bar(stat="identity",fill = "steelblue") + geom_line(aes(x=day, y=num_receiver, group = 1, colour="#FF9999"), show.legend = F)+labs(x="", y="")+theme_minimal()+ scale_x_discrete(breaks = function(n) n[seq(1, length(n), by = length(n)/10)])+annotate(geom="text", size = 5, x=70, y=25, label="County Manager",color="steelblue")

summary(time[which(senders==1)])
summary(time[which(senders==4)])
#emails sent by Emergency Services official (node 4) --- total 325
byday = rep(0, length(unique(day)))
names(byday)  = unique(day)
i = 1
for (d in unique(day)) {
	docs = which(day == d)
	byday[i] = sum(senders[docs] == 4)
	i = i + 1
}
byday2 = rep(0, length(unique(day)))
names(byday2)  = unique(day)
i = 1
for (d in unique(day)) {
	docs = which(day == d)
	byday2[i] = sum(receivers[docs][senders[docs] == 4])
	i = i + 1
}

words = unlist(lapply(which(senders==4), function(x) Dare$text[[x]]))
Dare$vocab[which(tabulate(words) %in% sort(tabulate(words), decreasing = TRUE)[1:15])]

SummaryTable$Words2 <- Dare$vocab[which(tabulate(words) %in% sort(tabulate(words), decreasing = TRUE)[1:15])]

data = data.frame(day = as.factor(unique(day)), send = byday, num_receiver = byday2)
p[[2]]<-ggplot(data=data, aes(x=day, y=send)) + geom_bar(stat="identity",fill = "#009E73") + geom_line(aes(x=day, y=num_receiver, group = 1, colour="#FF9999"), show.legend = F)+labs(x="", y="")+theme_minimal()+ scale_x_discrete(breaks = function(n) n[seq(1, length(n), by = length(n)/10)])+annotate(geom="text", size = 5, x=66, y=188, label="Emergency Department",color="#009E73")


marrangeGrob(p[1:2], nrow = 2, ncol = 1, top = NULL)
             
             

names(SummaryTable) <- c(

expression("County Manager"),

              expression("Emergency Department"))


tt <- ttheme_default(colhead=list(fg_params = list(parse=TRUE)))

tbl <- tableGrob(SummaryTable[,1], rows=NULL, theme=ttheme_default())

pl = replicate(2, ggplot(), FALSE)
pl[[1]] = p[[1]]
pl[[2]] = tbl
grid.arrange(grobs= lapply(pl, "+", theme(plot.margin=margin(1,1,1,1))))
