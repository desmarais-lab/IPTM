setwd('/Users/bomin8319/Desktop/IPTM/Dare_full')

load('Darenew.RData')
# 762 - 
attach(Dare)
Dare$node = 1:nrow(Dare$node)
Dare$edge = Dare$edge[-which(sapply(Dare$text, function(d){length(d)})==0)]
Dare$text = Dare$text[-which(sapply(Dare$text, function(d){length(d)})==0)]
mintime = Dare$edge[[1]][[3]]
for (n in 1:length(Dare$edge)){
  Dare$edge[[n]][3] = (Dare$edge[[n]][[3]] - mintime) / 3600
}
Dare$edge = lapply(Dare$edge, function(x){x[1:3]})

o_outdegree = tabulate(sapply(Dare$edge[385:2210], function(x) {x$sender}), 27)
o_indegree = tabulate(unlist(lapply(Dare$edge[385:2210], function(x) {x$receiver})), 27)
o_multicast = tabulate(sapply(Dare$edge[385:2210], function(x) {length(x$receiver)}), 27)
o_timediff = sapply(385:2210, function(x) {Dare$edge[[x]]$unixtime - Dare$edge[[x-1]]$unixtime })
names(o_outdegree) = c(1:27)
names(o_indegree) = c(1:27)

outdegree = matrix(0, 100, 27)
indegree = matrix(0, 100, 27)
colnames(outdegree) = c(1:27)
colnames(indegree) = c(1:27)
multicast = matrix(0, 100, 27)
timediff = matrix(0, 100, length(o_timediff))
d = 1
for (i in 1:5) {
	filename = paste0("PPC", i, ".RData")
	load(filename)
	for (j in 1:20) {
		PPC_j = PPC[[j]]
		outdegree[d,] = tabulate(sapply(PPC_j$edge[385:2210], function(x) {x$sender}), 27)
		indegree[d,] = tabulate(unlist(lapply(PPC_j$edge[385:2210], function(x) {x$receiver})), 27)
		multicast[d,] =  tabulate(sapply(PPC_j$edge[385:2210], function(x) {length(x$receiver)}), 27)
		timediff[d, ] = unlist(sapply(385:2210, function(x) {PPC_j$edge[[x]][[3]]- PPC_j$edge[[x-1]][[3]]}))
	d = d + 1
	}
}

outdegree = outdegree[1:80, ]
indegree = indegree[1:80, ]
multicast = multicast[1:80,]
timediff = timediff[1:80,]


sortedoutdegree = sort(o_outdegree, decreasing = TRUE)
outdegree = outdegree[, as.numeric(names(sortedoutdegree))]
bplot(outdegree, xlab = "Actor", ylab = "Outdegree")
lines(sortedoutdegree, col = 'red')


sortedindegree = sort(o_indegree, decreasing = TRUE)
indegree = indegree[, as.numeric(names(sortedindegree))]
bplot(indegree, xlab = "Actor", ylab = "Indegree")
lines(sortedindegree, col = 'red')

bplot(multicast[,1:24],xlab = "Number of Receipient", ylab = "")
lines(o_multicast[1:24], col = 'red')

qqplot(o_timediff, timediff, xlab = "Observed Quantile", ylab = "Simulated Quantile")
abline(0, 1)

