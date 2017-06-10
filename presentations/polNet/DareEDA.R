#Dare EDA
library(anytime)
library(ggplot2)
library(MCMCpack)
library(reshape2)
library(gridExtra)
library(ggrepel)
library(RColorBrewer)

load('/Users/bomin8319/Desktop/IPTM/paper/Darenew.RData')
# 762 - 
attach(Dare)
Dare$text = Dare$text[762:length(Dare$edge)]
Dare$edge = Dare$edge[762:length(Dare$edge)]
Dare$edge = Dare$edge[-which(sapply(Dare$text, function(d){length(d)})==0)]
Dare$text = Dare$text[-which(sapply(Dare$text, function(d){length(d)})==0)]
#mintime = Dare$edge[[1]][[3]]
#for (n in 1:length(Dare$edge)){
#  Dare$edge[[n]][3] = (Dare$edge[[n]][[3]] - mintime) / 3600
#}

sender = sapply(Dare$edge, function(d){d[[1]]})
sender = data.frame(sender=sender, dept = Dare$node[sender, 3], time = sapply(Dare$edge, function(d){d[[3]]}), date = anydate(sapply(Dare$edge, function(d){d[[3]]})))
sender$date = format(sender$date, format = "%m/%d")
Dept = data.frame(Date = c(sapply(unique(sender$date), function(d){rep(as.character(d), 22)})), 
				  Department = c(sapply(unique(sender$date), function(d){sort(unique(Dare$node[,3]))
})), Send = NA)
i = 1
for (date in unique(sender$date)) {
Dept[(22*(i-1)+1) : (22*i),3]= tabulate(sender[which(sender[,4] == date),2], 22)
i = i + 1
}

	f <- ggplot(Dept, aes(Date, Send, colour = Department), show.legend=FALSE)
	f + geom_line(aes(group = Department)) + guides(col = guide_legend(nrow=22)) + scale_x_discrete(breaks = function(n) n[seq(0, length(n), by = length(n)/10)]) + theme(legend.text = element_text(size = 10)) +geom_vline(xintercept = 23, colour = "red", size = 1) +geom_vline(xintercept = 27, colour = "red", size = 1) + 
	geom_vline(xintercept = 20, colour = "red", size = 0.5, linetype = "dashed")+ annotate("text", x =25, y = 21, label = "Sandy", colour= "red" ) + 
	annotate("segment", x = 18, xend = 20, y = 15, yend = 14, colour = "red", size = 0.1, arrow = arrow()) + annotate("text", x =18, y = 15.5, label = "First Sandy", colour= "red", size = 3)  + theme(legend.position="none")

	which(Dare$vocab == "sandy")
which(Dare$vocab == "hurricane")
which(sapply(Dare$text, function(d) {49 %in% d})==TRUE)
which(sapply(Dare$text, function(d) {81 %in% d})==TRUE)

Dare$edge[[394]]
Dare$edge[[176]]


receiver = unlist(sapply(Dare$edge, function(d){d[[2]]}))
time = unlist(sapply(Dare$edge, function(d){rep(d[[3]], length(d[[2]]))}))
date = anydate(unlist(sapply(Dare$edge, function(d){rep(d[[3]], length(d[[2]]))})))
receiver  = data.frame(receiver =receiver, dept = Dare$node[receiver , 3], time = time, date = date)
receiver$date = format(receiver$date, format = "%m/%d")
Dept = data.frame(Date = c(sapply(unique(receiver$date), function(d){rep(as.character(d), 22)})), 
				  Department = c(sapply(unique(receiver$date), function(d){sort(unique(Dare$node[,3]))
})), Receive = NA)
i = 1
for (date in unique(receiver$date)) {
Dept[(22*(i-1)+1) : (22*i),3]= tabulate(receiver[which(receiver[,4] == date),2], 22)
i = i + 1
}

	f <- ggplot(Dept, aes(Date, Receive, colour = Department))
	f + geom_line(aes(group = Department)) + guides(col = guide_legend(nrow=22)) + scale_x_discrete(breaks = function(n) n[seq(0, length(n), by = length(n)/10)]) + theme(legend.text = element_text(size = 5)) +geom_vline(xintercept = 23, colour = "red", size = 1) +geom_vline(xintercept = 27, colour = "red", size = 1) + 
	geom_vline(xintercept = 20, colour = "red", size = 0.5, linetype = "dashed")+ annotate("text", x =25, y = 27, label = "Sandy", colour= "red" ) + 
	annotate("segment", x = 18, xend = 20, y = 20, yend = 19, colour = "red", size = 0.1, arrow = arrow()) + annotate("text", x =18, y = 20.5, label = "First Sandy", colour= "red", size = 3)  + theme(legend.position="none")

	


Sandy = sapply(Dare$text, function(d){sum(49 == d)})
Sandy = data.frame(Sandy=Sandy, time = sapply(Dare$edge, function(d){d[[3]]}), date = anydate(sapply(Dare$edge, function(d){d[[3]]})))
Sandy$date = format(Sandy$date, format = "%m/%d")
Hurr = sapply(Dare$text, function(d){sum(81 == d)})
Hurr = data.frame(Hurr =Hurr , time = sapply(Dare$edge, function(d){d[[3]]}), date = anydate(sapply(Dare$edge, function(d){d[[3]]})))
Hurr $date = format(Hurr$date, format = "%m/%d")

Dept = data.frame(Date = c(sapply(unique(Sandy$date), function(d){rep(as.character(d), 2)})), Word = c(sapply(unique(Sandy$date), function(d){c("Hurricane", "Sandy")})), Count = NA)
i = 1
for (date in unique(Sandy$date)) {
Dept[(2*(i-1)+1) : (2*i),3]= c(sum(Hurr[which(Hurr[,3] == date),1]), sum(Sandy[which(Sandy[,3] == date),1]))
i = i + 1
}


	f <- ggplot(Dept, aes(Date, Count, colour = Word))
	f + geom_line(aes(group = Word))+ guides(col = guide_legend(nrow=2))+ scale_x_discrete(breaks = function(n) n[seq(0, length(n), by = length(n)/10)]) +geom_vline(xintercept = 23, colour = "red", size = 1) +geom_vline(xintercept = 27, colour = "red", size = 1) + 
	geom_vline(xintercept = 20, colour = "red", size = 0.5, linetype = "dashed")+ annotate("text", x =25, y = 27, label = "Sandy", colour= "red" ) + 
	annotate("segment", x = 18, xend = 20, y = 20, yend = 19, colour = "red", size = 0.1, arrow = arrow()) + annotate("text", x =18, y = 20.5, label = "First Sandy", colour= "red", size = 3) 
	
which(Dare$vocab == "sandy")
which(Dare$vocab == "hurricane")
Sandy = which(sapply(Dare$text, function(d) {49 %in% d})==TRUE)
Hurr = which(sapply(Dare$text, function(d) {81 %in% d})==TRUE)

sapply(Sandy, function(d){sum(49 == Dare$text[[d]])})
sapply(Hurr, function(d){sum(81== Dare$text[[d]])})



library(GGally)
library(network)
library(sna)

Network = list()
Network[[1]]= Dare$edge[which(Sandy$date %in% unique(Sandy$date)[1:18])] #pre-sandy
Network[[2]]= Dare$edge[which(Sandy$date %in% unique(Sandy$date)[19:39])] #sandy 
Network[[3]]= Dare$edge[which(Sandy$date %in% unique(Sandy$date)[40:55])] #post-sandy
ggplotColours <- function(n = 6, h = c(0, 360) + 15){
  if ((diff(h) %% 360) < 1) h[2] <- h[2] - 360/n
  hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
}
titles = c("Pre-Sandy", "Sandy", "Post-Sandy")
colors = ggplotColours(22)
names(colors) = (unique(sort(Dare$node[,3])))
g= list()
for (i in 1:3) {
	edge= matrix(unlist(sapply(Network[[i]], function(d){cbind(rep(d[[1]], length(d[[2]])), d[[2]])})), ncol = 2)
	net = as.network(edge, matrix.type = "edgelist", directed = TRUE)
	network.vertex.names(net) = 1:27
	net %v% "Dept" <- as.character(Dare$node[,3])
	g[[i]] = ggnet2(net, size = "degree", arrow.size =8, color = colors[net %v% "Dept"],  color.legend = "Dept",label.color = colors[net %v% "Dept"], mode = "kamadakawai",legend.position = "none", size.min=1) + ggtitle(titles[i]) +  theme(panel.border = element_rect(color = "grey50", fill = NA),
          aspect.ratio = 1)
}
grid.arrange <- getFromNamespace("grid.arrange", asNamespace("gridExtra"))
grid.arrange(grobs = g, nrow = 1)

