#minimal#

library(MCMCpack)
alpha = 2
alpha0 = 2
alpha1 =2
K = 4
D = 10
samp = 500
results = matrix(NA, nrow = samp, ncol = 2*K)
for (s in 1:samp) {
#assume IP is fixed
nIP = 2
cd = sapply(1:D, function(d) which(rmultinom(1, 1, rep(1/nIP, nIP))==1))

#generate topics according to the generative process
#m = rdirichlet(1, alpha0*rep(1/K, K))
#mc = rdirichlet(nIP, alpha1*m)
#theta = t(sapply(1:D, function(d){rdirichlet(1, alpha*mc[cd[d],])}))
n.d = 10
z = matrix(NA, D, n.d)
table.k = rep(0, K)
table.dk = matrix(0, K, D)
table.cd = matrix(0, K, nIP)
for (d in 1:D) {
	for (n in 1:n.d) {
	 ratio2 = (table.k+alpha0/K)/sum(table.k+alpha0/K)
     ratio1 = (table.cd[,cd[d]]+alpha1*ratio2) / sum(table.cd[,cd[d]]+alpha1*ratio2)
      #ratio1 = (table.cd[,cd[d]]+alpha1/K) / sum(table.cd[,cd[d]]+alpha1/K)
	  prob = table.dk[,d]+alpha*ratio1
      z[d,n] = which(rmultinom(1, 1, prob)==1)
      #re-fill the tables 
      table.dk[z[d,n],d] = table.dk[z[d,n],d]+1
      if (table.dk[z[d,n],d]==1) {
      table.cd[z[d,n], cd[d]] = table.cd[z[d,n], cd[d]]+1		
	  }
	  table.k = rowSums(table.cd>0)
	}
}
#this table.k.null is the true topic distribution across the corpus
table.k.null = tabulate(z, K)

#initialize N_k, N_dk, N_ck 
#z = matrix(NA, D, n.d)
#for (d in 1:D) {
#  theta = rdirichlet(1, alpha0*rep(1/K, K))
#  hi = rmultinom(n.d, 1, theta)
#  z[d,] = sapply(1:n.d, function(x) which(hi[,x]==1))
#}

#table.k = tabulate(z, K)
#table.dk = sapply(1:D, function(d) tabulate(z[d,], K))
#table.cd = sapply(1:nIP, function(IP) {tabulate(z[which(cd==IP),],K)})

#inference on z using equation (15)---without word part that we ignored for simplified model
for (i in 1:100){
  for (d in 1:D) {
    for (n in 1:n.d) {
      #excluding \dn
      table.dk[z[d,n],d] = table.dk[z[d,n],d]-1
      if (table.dk[z[d,n], d] == 0) {
      	 table.cd[z[d,n], cd[d]] = table.cd[z[d,n], cd[d]]-1
      }
      table.k = rowSums(table.cd>0)
      
      #calculate Equation (15) and update z_dn
      ratio2 = (table.k+alpha0/K)/sum(table.k+alpha0/K)
      ratio1 = (table.cd[,cd[d]]+alpha1*ratio2) / sum(table.cd[,cd[d]]+alpha1*ratio2)
      #ratio1 = (table.cd[,cd[d]]+alpha1/K) / sum(table.cd[,cd[d]]+alpha1/K)
      prob = table.dk[,d]+alpha*ratio1
      z[d,n] = which(rmultinom(1, 1, prob)==1)
      #re-fill the tables 
      table.dk[z[d,n],d] = table.dk[z[d,n],d]+1
      if (table.dk[z[d,n],d]==1) {
      table.cd[z[d,n], cd[d]] = table.cd[z[d,n], cd[d]]+1		
	  }
	  table.k = rowSums(table.cd>0)
    }
  }
}

#see the results
results[s,] = cbind(table.k.null,  tabulate(z, K))
#even in this simplified model (only topic being random variable) and started with the true topic assignments, 
#we see that the topic assignments do not converges to true distribution
}

plot(sort(results[,1]), sort(results[,5]))
abline(0,1)



##########################maximal
library(MCMCpack)
alpha = 5
alpha0 = 100
alpha1 = 50
K = 4
D = 10
samp = 500
results = matrix(NA, nrow = samp, ncol = 2*K)
for (s in 1:samp) {
#assume IP is fixed
nIP =2
cd = sapply(1:D, function(d) which(rmultinom(1, 1, rep(1/nIP, nIP))==1))

#generate topics according to the generative process
#m = rdirichlet(1, alpha0*rep(1/K, K))
#mc = rdirichlet(nIP, alpha1*m)
#theta = t(sapply(1:D, function(d){rdirichlet(1, alpha*mc[cd[d],])}))
n.d = 10
z = matrix(NA, D, n.d)
table.k = rep(0, K)
table.dk = matrix(0, K, D)
table.cd = matrix(0, K, nIP)
for (d in 1:D) {
	for (n in 1:n.d) {
	  ratio2 = (table.k+alpha0/K)/sum(table.k+alpha0/K)
      ratio1 = (table.cd[,cd[d]]+alpha1*ratio2) / sum(table.cd[,cd[d]]+alpha1*ratio2)
      #ratio1 = (table.cd[,cd[d]]+alpha1/K) / sum(table.cd[,cd[d]]+alpha1/K)
	  #ratio1 = (table.cd[,cd[d]]+alpha1*m) / sum(table.cd[,cd[d]]+alpha1*m)
	  prob = table.dk[,d]+alpha*ratio1
      z[d,n] = which(rmultinom(1, 1, prob)==1)
      #re-fill the tables 
      table.dk[z[d,n],d] = table.dk[z[d,n],d]+1
      table.cd[z[d,n], cd[d]] = table.cd[z[d,n], cd[d]]+1		
	  table.k[z[d,n]] = table.k[z[d,n]] +1
	}
}
#this table.k.null is the true topic distribution across the corpus
table.k.null = tabulate(z, K)

# #initialize N_k, N_dk, N_ck 
 # z = matrix(NA, D, n.d)
 # for (d in 1:D) {
   # theta = rdirichlet(1, alpha0*rep(1/K, K))
   # hi = rmultinom(n.d, 1, theta)
   # z[d,] = sapply(1:n.d, function(x) which(hi[,x]==1))
 # }

 # table.k = tabulate(z, K)
 # table.dk = sapply(1:D, function(d) tabulate(z[d,], K))
 # table.cd = sapply(1:nIP, function(IP) {tabulate(z[which(cd==IP),],K)})
#inference on z using equation (15)---without word part that we ignored for simplified model
for (i in 1:100){
  for (d in 1:D) {
    for (n in 1:n.d) {
      #excluding \dn
      table.dk[z[d,n],d] = table.dk[z[d,n],d]-1
      table.cd[z[d,n], cd[d]] = table.cd[z[d,n], cd[d]]-1
      table.k[z[d,n]] = table.k[z[d,n]] - 1

      #calculate Equation (15) and update z_dn
      ratio2 = (table.k+alpha0/K)/sum(table.k+alpha0/K)
      ratio1 = (table.cd[,cd[d]]+alpha1*ratio2) / sum(table.cd[,cd[d]]+alpha1*ratio2)
      #ratio1 = (table.cd[,cd[d]]+alpha1/K) / sum(table.cd[,cd[d]]+alpha1/K)
      #ratio1 = (table.cd[,cd[d]]+alpha1*mnew) / sum(table.cd[,cd[d]]+alpha1*m)
      prob = table.dk[,d]+alpha*ratio1
      z[d,n] = which(rmultinom(1, 1, prob)==1)
      #re-fill the tables 
      table.dk[z[d,n],d] = table.dk[z[d,n],d]+1
      table.cd[z[d,n], cd[d]] = table.cd[z[d,n], cd[d]]+1		
	  table.k[z[d,n]] = table.k[z[d,n]] +1
    }
  }
}

#see the results
results[s,] = cbind(table.k.null,  tabulate(z, K))
#even in this simplified model (only topic being random variable) and started with the true topic assignments, 
#we see that the topic assignments do not converges to true distribution
}

plot(sort(results[,1]), sort(results[,5]))
abline(0,1)


#######Two hierarchy#########

##################
library(MCMCpack)
alpha = 5
alpha1 = 50
K = 4
D = 10
samp = 500
results = matrix(NA, nrow = samp, ncol = 2*K)
for (s in 1:samp) {
#assume IP is fixed
nIP = 4

#generate topics according to the generative process
m = rdirichlet(1, alpha1*rep(1/K, K))
#mc = rdirichlet(nIP, alpha1*m)
theta = t(sapply(1:D, function(d){rdirichlet(1, alpha*m)}))
n.d =10
z = matrix(NA, D, n.d)
for (d in 1:D) {
  hi = rmultinom(n.d, 1, theta[d,])
  z[d,] = sapply(1:n.d, function(x) which(hi[,x]==1))
}
#this table.k.null is the true topic distribution across the corpus
table.k.null = tabulate(z, K)
table.dk.null =  sapply(1:D, function(d) tabulate(z[d,], K))
#initialize N_k, N_dk, N_ck using the true values
table.k = tabulate(z, K)
table.dk = sapply(1:D, function(d) tabulate(z[d,], K))

#inference on z using equation (15)---without word part that we ignored for simplified model
for (i in 1:100){
  for (d in 1:D) {
    for (n in 1:n.d) {
      #excluding \dn
      table.k[z[d,n]] = table.k[z[d,n]]-1
      table.dk[z[d,n],d] = table.dk[z[d,n],d]-1
      prob = table.dk[,d]+alpha * (table.k+alpha1/K)/(sum(table.k)+alpha1)
      ratio2 = (table.k+alpha1/K)/(sum(table.k)+alpha1)
      prob = table.dk[,d]+alpha *ratio2
      z[d,n] = which(rmultinom(1, 1, prob)==1)
      #re-fill the tables 
      table.k[z[d,n]] = table.k[z[d,n]]+1
      table.dk[z[d,n],d] = table.dk[z[d,n],d]+1
    }
  }
}

#see the results
results[s,] = cbind(table.k.null, table.k)
#even in this simplified model (only topic being random variable) and started with the true topic assignments, 
#we see that the topic assignments do not converges to true distribution
}

plot(sort(results[,1]), sort(results[,5]))
abline(0,1)

#################
#collapsed
library(MCMCpack)
alpha = 10
alpha1 = 10
K = 4
D = 10
samp = 500
results = matrix(NA, nrow = samp, ncol = 2*K)
for (s in 1:samp) {
#assume IP is fixed
nIP = 2

#generate topics according to the generative process
#m = rdirichlet(1, alpha1*rep(1/K, K))
#mc = rdirichlet(nIP, alpha1*m)
#theta = t(sapply(1:D, function(d){rdirichlet(1, alpha*m)}))
n.d =10
z = matrix(NA, D, n.d)
table.k = rep(0, K)
table.dk = matrix(0, K, D)
for (d in 1:D) {
	for (n in 1:n.d) {
	  prob = table.dk[,d]+alpha* (table.k+alpha1/K)/(sum(table.k)+alpha1)
      z[d,n] = which(rmultinom(1, 1, prob)==1)
	table.k[z[d,n]] = table.k[z[d,n]]+1
	table.dk[z[d,n], d] = table.dk[z[d,n], d] + 1
	}      
 }

#this table.k.null is the true topic distribution across the corpus
table.k.null = tabulate(z, K)
table.dk.null =  sapply(1:D, function(d) tabulate(z[d,], K))
#initialize N_k, N_dk, N_ck using the true values
table.k = tabulate(z, K)
table.dk = sapply(1:D, function(d) tabulate(z[d,], K))

#inference on z using equation (15)---without word part that we ignored for simplified model
for (i in 1:100){
  for (d in 1:D) {
    for (n in 1:n.d) {
      #excluding \dn
      table.k[z[d,n]] = table.k[z[d,n]]-1
      table.dk[z[d,n],d] = table.dk[z[d,n],d]-1
      ratio2 = (table.k+alpha1/K)/(sum(table.k)+alpha1)
      prob = table.dk[,d]+alpha *ratio2
      z[d,n] = which(rmultinom(1, 1, prob)==1)
      #re-fill the tables 
      table.k[z[d,n]] = table.k[z[d,n]]+1
      table.dk[z[d,n],d] = table.dk[z[d,n],d]+1
    }
  }
}

#see the results
results[s,] = cbind(table.k.null, table.k)
#even in this simplified model (only topic being random variable) and started with the true topic assignments, 
#we see that the topic assignments do not converges to true distribution
}

plot(sort(results[,1]), sort(results[,5]))
abline(0,1)



##############standard LDA########


#################
library(MCMCpack)
alpha = 5
K = 4
D = 10
samp = 1000
results = matrix(NA, nrow = samp, ncol = 2*K)
for (s in 1:samp) {
#assume IP is fixed
nIP = 4

#generate topics according to the generative process
#m = rdirichlet(1, alpha0*rep(1/K, K))
#mc = rdirichlet(nIP, alpha1*m)
theta = t(sapply(1:D, function(d){rdirichlet(1, alpha*rep(1/K, K))}))
n.d =10
z = matrix(NA, D, n.d)
for (d in 1:D) {
  hi = rmultinom(n.d, 1, theta[d,])
  z[d,] = sapply(1:n.d, function(x) which(hi[,x]==1))
}
#this table.k.null is the true topic distribution across the corpus
table.k.null = tabulate(z, K)
table.dk.null =  sapply(1:D, function(d) tabulate(z[d,], K))
#initialize N_k, N_dk, N_ck using the true values
table.k = tabulate(z, K)
table.dk = sapply(1:D, function(d) tabulate(z[d,], K))
#table.cd = sapply(1:nIP, function(IP) {tabulate(z[which(cd==IP),],K)})

#inference on z using equation (15)---without word part that we ignored for simplified model
for (i in 1:100){
  for (d in 1:D) {
    for (n in 1:n.d) {
      #excluding \dn
      table.k[z[d,n]] = table.k[z[d,n]]-1
      table.dk[z[d,n],d] = table.dk[z[d,n],d]-1
      #table.cd[z[d,n], cd[d]] = table.cd[z[d,n], cd[d]]-1
      #calculate Equation (15) and update z_dn
      #ratio2 = (table.k+alpha0/K)/(sum(table.k)+alpha0)
      #ratio1 = (table.cd[,cd[d]]+alpha1*ratio2) / sum(table.cd[,cd[d]]+alpha1*ratio2)
      prob = log(table.dk[,d]+alpha/K)
      z[d,n] = which(rmultinom(1, 1, exp(prob))==1)
      #re-fill the tables 
      table.k[z[d,n]] = table.k[z[d,n]]+1
      table.dk[z[d,n],d] = table.dk[z[d,n],d]+1
      #table.cd[z[d,n], cd[d]] = table.cd[z[d,n], cd[d]]+1
    }
  }
}

#see the results
results[s,] = cbind(table.k.null, table.k)
#even in this simplified model (only topic being random variable) and started with the true topic assignments, 
#we see that the topic assignments do not converges to true distribution
}

plot(sort(results[,1]), sort(results[,5]))
abline(0,1)


##########
#collapsed
library(MCMCpack)
alpha = 5
K = 4
D = 10
samp = 1000
results = matrix(NA, nrow = samp, ncol = 2*K)
for (s in 1:samp) {
#assume IP is fixed
nIP = 4

#generate topics according to the generative process
#m = rdirichlet(1, alpha0*rep(1/K, K))
#mc = rdirichlet(nIP, alpha1*m)
#theta = t(sapply(1:D, function(d){rdirichlet(1, alpha*rep(1/K, K))}))
n.d =10
z = matrix(NA, D, n.d)
table.k = rep(0, K)
table.dk = matrix(0, K, D)
for (d in 1:D) {
	for (n in 1:n.d) {
	  prob = table.dk[,d]+alpha*rep(1/K, K)
      z[d,n] = which(rmultinom(1, 1, prob)==1)
	table.k[z[d,n]] = table.k[z[d,n]]+1
	table.dk[z[d,n], d] = table.dk[z[d,n], d] + 1
	}      
 }
#this table.k.null is the true topic distribution across the corpus
table.k.null = tabulate(z, K)
table.dk.null =  sapply(1:D, function(d) tabulate(z[d,], K))
#initialize N_k, N_dk, N_ck using the true values
table.k = tabulate(z, K)
table.dk = sapply(1:D, function(d) tabulate(z[d,], K))
#table.cd = sapply(1:nIP, function(IP) {tabulate(z[which(cd==IP),],K)})

#inference on z using equation (15)---without word part that we ignored for simplified model
for (i in 1:100){
  for (d in 1:D) {
    for (n in 1:n.d) {
      #excluding \dn
      table.k[z[d,n]] = table.k[z[d,n]]-1
      table.dk[z[d,n],d] = table.dk[z[d,n],d]-1
      #table.cd[z[d,n], cd[d]] = table.cd[z[d,n], cd[d]]-1
      #calculate Equation (15) and update z_dn
      #ratio2 = (table.k+alpha0/K)/(sum(table.k)+alpha0)
      #ratio1 = (table.cd[,cd[d]]+alpha1*ratio2) / sum(table.cd[,cd[d]]+alpha1*ratio2)
      prob = log(table.dk[,d]+alpha/K)
      z[d,n] = which(rmultinom(1, 1, exp(prob))==1)
      #re-fill the tables 
      table.k[z[d,n]] = table.k[z[d,n]]+1
      table.dk[z[d,n],d] = table.dk[z[d,n],d]+1
      #table.cd[z[d,n], cd[d]] = table.cd[z[d,n], cd[d]]+1
    }
  }
}

#see the results
results[s,] = cbind(table.k.null, table.k)
#even in this simplified model (only topic being random variable) and started with the true topic assignments, 
#we see that the topic assignments do not converges to true distribution
}

plot(sort(results[,1]), sort(results[,5]))
abline(0,1)




#' @title AlphamvecOpt
#' @description Optimize the hyperparmeter vector (= alpha*mvec) given the current assignment of interaction patterns and topics
#'
#' @param K total number of topics specified by the user
#' @param z current state of the assignment of topics 
#' @param alpha Dirichlet concentration prior for topic distribution
#' @param mvec Dirichlet base prior for topic distribution
#' @param niter number of iterations to perfom
#'
#' @return The optimized value of the vector (= alpha*mvec)
#'
#' @export
AlphamvecOpt =  function(K, z, alpha, mvec, niter) {
	current.vec = alpha*mvec
	n.word = mapply(length, z)
	n.word.table = tabulate(n.word)
	nK.word.list = matrix(NA, nrow = length(z), ncol = K)
	for (d in seq(along = z)) {
		nK.word.list[d, ] = tabulate(z[[d]], K)
	}
	nK.word.table = lapply(1:K, function(k){
		tabulate(nK.word.list[,k])
	})
	for (i in 1:niter) {
	alpha = sum(current.vec)
	S = UpdateDenom(alpha, n.word.table)		
	s = UpdateNum(current.vec, nK.word.table)
	current.vec = current.vec*s/S
	}
	return(current.vec)	
}

znew = list()
for (d in 1:D) {
	znew[[d]] = z[d,]
}
alphamvec = AlphamvecOpt(K, znew, alpha, m, 50)
mnew = alphamvec / sum(alphamvec)
alphanew = sum(alphamvec)