setwd('/Users/bomin8319/Desktop/IPTM/full')
library(fields)
library(ergm)
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

o_outdegree = tabulate(sapply(Dare$edge[385:2210], function(x) {x[[1]]}), 27)
o_indegree = tabulate(unlist(lapply(Dare$edge[385:2210], function(x) {x[[2]]})), 27)
o_multicast = tabulate(sapply(Dare$edge[385:2210], function(x) {length(x[[2]])}), 27)
o_timediff = sapply(385:2210, function(x) {Dare$edge[[x]][[3]] - Dare$edge[[x-1]][[3]] })
names(o_outdegree) = c(1:27)
names(o_indegree) = c(1:27)

K = 25
o_mi = matrix(0, 5, K)
for (i in c(1,3,4,5)) {
    filename = paste0("Dare_full_",2,"_",25,"_ver",i,".RData")
    load(filename)
    o_mi[i, ] = mi_calc(Dare$text[385:2210], topics = Daretest$Z[385:2210], K, length(Dare$vocab))
}
networkmat = matrix(0, nrow = 0, ncol = 2)
        for (d in 385:2210) {
            for (l in 1:length(Dare$edge[[d]]$receiver)) {
                networkmat = rbind(networkmat, c(Dare$edge[[d]][[1]], Dare$edge[[d]][[2]][l]))
            }
        }
        net=network(networkmat,matrix.type="edgelist",directed=TRUE)
        ergmfit = ergm(net ~ edges)
        gofs = gof(ergmfit)
o_geodist = gofs$obs.dist / sum(gofs$obs.dist)
o_espart = gofs$obs.espart / sum(gofs$obs.espart)


outdegree = matrix(0, 100, 27)
indegree = matrix(0, 100, 27)
colnames(outdegree) = c(1:27)
colnames(indegree) = c(1:27)
multicast = matrix(0, 100, 27)
timediff = matrix(0, 100, length(o_timediff))
MI = matrix(0, 100, K)
geodist = matrix(0, 100, 27)
espart = matrix(0, 100, 26)


d = 1
for (i in c(1,3, 4, 5)) {
	filename = paste0("PPC", i, "_nIP2_K25_Dare.RData")
	load(filename)
    filename2 = paste0("Dare_full_",2,"_",K,"_ver",i,".RData")
    load(filename2)
   	for (j in 1:20) {
		PPC_j = PPC[[j]]
		outdegree[d,] = tabulate(sapply(PPC_j$edge[385:2210], function(x) {x[[1]]}), 27)
		indegree[d,] = tabulate(unlist(lapply(PPC_j$edge[385:2210], function(x) {x[[2]]})), 27)
		#multicast[d,] =  tabulate(sapply(PPC_j$edge[385:2210], function(x) {length(x[[2]])}), 27)
		#timediff[d, ] = unlist(sapply(385:2210, function(x) {PPC_j$edge[[x]][[3]]- PPC_j$edge[[x-1]][[3]]}))
        #MI[d, ] = mi_calc(PPC_j$text[385:2210], topics = Daretest$Z[385:2210], K, length(Dare$vocab))
        ##Z's from different inferences may have label switching problem
        #networkmat = matrix(0, nrow = 0, ncol = 2)
        #for (doc in 385:2210) {
        #    for (l in 1:length(PPC_j$edge[[doc]]$receiver)) {
        #        networkmat = rbind(networkmat, c(PPC_j$edge[[doc]][[1]], PPC_j$edge[[doc]][[2]][l]))
        #    }
        #}
        #net=network(networkmat,matrix.type="edgelist",directed=TRUE)
        #ergmfit = ergm(net ~ edges, control = control.ergm(MCMC.interval = 1, MCMC.burnin = 0, MCMC.samplesize = 1, MCMLE.maxit = 1))
        #gofs = gof(ergmfit)
		#geodist[d, ] =gofs$obs.dist / sum(gofs$obs.dist)
		#espart[d,] = gofs$obs.espart / sum(gofs$obs.espart)
	d = d + 1
	}
}

outdegree = outdegree[1:80,]
indegree = indegree[1:80,]
multicast = multicast[1:80,]
timediff = timediff[1:80,]
MI= MI[1:80,]
geodist = geodist[1:80,]
espart= espart[1:80,]

par(mfrow = c(2,3))
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

# # o_geodist2 = unlist(sapply(1:4, function(x){rep(x, o_geodist[x])}))
# geodist2 = list()
# for (d in 1:100) {
	# geodist2[[d]] =  unlist(sapply(1:4, function(x){rep(x, geodist[d,x])}))
# }
# hi = sapply(geodist2, function(x) { mean(x)})
# hist(hi)
# abline(h = mean(o_geodist2), col ='red')
# qqplot(o_geodist, geodist,  xlab = "Observed Quantile", ylab = "Simulated Quantile")
# abline(0,1)

# qqplot(o_geodist2, geodist2,  xlab = "Observed Quantile", ylab = "Simulated Quantile")
# abline(0,1)

bplot(geodist,xlab = "Geodesic distance", ylab = "Proportion")
lines(o_geodist, col = 'red')

bplot(espart,xlab = "Edgewise shared partner", ylab = "Proportion")
lines(o_espart, col = 'red')

par(mfrow = c(1,5))
for (k in 1:5) {
	hist(MI[,k], main = paste("Topic", k), xlab = "Mutual Information")
	abline(v = mean(o_mi[c(1,3,4,5), k]), col = 'red')
}

mi_calc =function(text, topics = lapply(text, function(x){as.numeric(names(x))}), K =K, W) {
	D = length(text)
	MI = rep(0, K)
	N_k = tabulate(unlist(topics), K)
	N_dk = t(sapply(topics, function(x){ tabulate(x, K)}))
	N_wk = matrix(0, W, K)
	for (w in 1:W) {
		N_wk[w, ] = tabulate(unlist(topics)[which(unlist(text) == w)], K)
	}
    for (d in 1:D) {
        text_d = text[[d]]
        for (w in unique(text_d)) {
            N_wdk = tabulate(topics[[d]][which(text_d == w)], K)
            for (k in which(N_wdk > 0)) {
                MI[k] = MI[k] + N_wdk[k] /N_k[k] * (log(N_k[k]) + log(N_wdk[k]) - log(N_dk[d,k]) - log(N_wk[w,k]))
            }
        }
    }
	return(MI)
}
