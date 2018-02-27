#three level cluster LDA
library(MCMCpack)

###########################################
#maximal assumption which uses full counts#
###########################################
alpha = 5               #document-level parameter
alpha1 = 50			#cluster-level parameter
alpha0 = 100			#corpus-level parameter

K = 4
D = 10

#samp is the number of GiR samples 
samp = 5000     #takes about a minute with 1000
results = matrix(NA, nrow = samp, ncol = 2*K)
for (s in 1:samp) {
# alpha = 5              #document-level parameter
# alpha1 = 5			#cluster-level parameter
# alpha0 = 5			#corpus-level parameter
alpha = 10               #document-level parameter
alpha1 = 10			#cluster-level parameter
alpha0 = 10			#corpus-level parameter

alphas = c(alpha0, alpha1, alpha)

#assume IP is fixed
nIP =2
cd = sapply(1:D, function(d) which(rmultinom(1, 1, rep(1/nIP, nIP))==1))

#generate topics according to the generative process
#let's not use the non-collapsed version
#m = rdirichlet(1, alpha0*rep(1/K, K))
#mc = rdirichlet(nIP, alpha1*m)
#theta = t(sapply(1:D, function(d){rdirichlet(1, alpha*mc[cd[d],])}))

#collapsed cluster LDA generating process
n.d = 10
z = matrix(NA, D, n.d)
table.k = rep(0, K)
table.dk = matrix(0, K, D)
table.cd = matrix(0, K, nIP)
for (d in 1:D) {
	for (n in 1:n.d) {
	  ratio2 = (table.k+alpha0/K) / sum(table.k+alpha0/K)
      ratio1 = (table.cd[,cd[d]]+alpha1*ratio2) / sum(table.cd[,cd[d]]+alpha1*ratio2)
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

# initialize N_k, N_dk, N_ck using the true distribution
# for real data analysis, need to start from Z that is initialized using standard LDA 
#table.k = tabulate(z, K)
#table.dk = sapply(1:D, function(d) tabulate(z[d,], K))
#table.cd = sapply(1:nIP, function(IP) {tabulate(z[which(cd==IP),],K)})

#inference on z using equation (15)---without word part that we ignored for simplified model
for (i in 1:10){
	alphas = slice(alphas, c(5,1,0.1), 5)
	alpha0 = alphas[1]
	alpha1 = alphas[2]
	alpha = alphas[3]

  for (d in 1:D) {
    for (n in 1:n.d) {
      #excluding \dn
      table.dk[z[d,n],d] = table.dk[z[d,n],d]-1
      table.cd[z[d,n], cd[d]] = table.cd[z[d,n], cd[d]]-1
      table.k[z[d,n]] = table.k[z[d,n]] - 1

      #calculate Equation (15) and update z_dn
      ratio2 = (table.k+alpha0/K)/sum(table.k+alpha0/K)
      ratio1 = (table.cd[,cd[d]]+alpha1*ratio2) / sum(table.cd[,cd[d]]+alpha1*ratio2)
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
}

# similar to PPplot
par(mfrow = c(1,4))
for (j in 1:4) {
	plot(sort(results[,j]), sort(results[,j+4]))
	abline(0,1)
}
head(results)
#currently it passes the test, but adjusting alphas (e.g. (10, 10, 10)) makes it fail

###########################################
#minimal assumption which uses unique draw#
###########################################

library(MCMCpack)
alpha = 5               #document-level parameter
alpha1 = 50				#cluster-level parameter
alpha0 = 100			#corpus-level parameter

K = 4
D = 10
#samp is the number of GiR samples
samp = 1000
results = matrix(NA, nrow = samp, ncol = 2*K)
for (s in 1:samp) {

#assume IP is fixed -> need to use larger nIP to avoid uniform base measure for m
nIP = 4
cd = sapply(1:D, function(d) which(rmultinom(1, 1, rep(1/nIP, nIP))==1))

#generate topics according to the generative process
#only difference is how we count table.cd and table.k
#N_{kc} is the number of different documents belonging to c that use topic k
#N_k is the number of different interaction patterns in which k has been used
n.d = 10
z = matrix(NA, D, n.d)
table.k = rep(0, K)
table.dk = matrix(0, K, D)
table.cd = matrix(0, K, nIP)
for (d in 1:D) {
	for (n in 1:n.d) {
	 ratio2 = (table.k+alpha0/K)/sum(table.k+alpha0/K)
     ratio1 = (table.cd[,cd[d]]+alpha1*ratio2) / sum(table.cd[,cd[d]]+alpha1*ratio2)
	  prob = table.dk[,d]+alpha*ratio1
      z[d,n] = which(rmultinom(1, 1, prob)==1)
      #re-fill the tables 
      table.dk[z[d,n],d] = table.dk[z[d,n],d]+1
      #increase table.cd only if k is used in this document for the first time
      if (table.dk[z[d,n],d]==1) {
      table.cd[z[d,n], cd[d]] = table.cd[z[d,n], cd[d]]+1		
	  }
	  # rowSums gives the number of different IPs
	  table.k = rowSums(table.cd>0)
	}
}
#this table.k.null is the true topic distribution across the corpus
table.k.null = tabulate(z, K)

#this table.k.null is the true topic distribution across the corpus
table.k.null = tabulate(z, K)

# initialize N_k, N_dk, N_ck using the true distribution
# need to do nothing
# for real data analysis, need to start from Z that is initialized using standard LDA 

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
      prob = table.dk[,d]+alpha*ratio1
      z[d,n] = which(rmultinom(1, 1, prob)==1)
      #re-fill the tables in the same manner
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
}

# similar to PPplot
par(mfrow = c(1,4))
for (j in 1:4) {
	plot(sort(results[,j]), sort(results[,j+4]))
	abline(0,1)
}
head(results)
#currently it passes the test, but adjusting alphas (e.g. (10, 10, 10)) makes it fail


slice = function(alphas, sigmas, S) {
	for (s in 1:S) {
		unew = runif(1, 0, f(alphas))
		for (i in 1:length(alphas)) {
			alphas.left[i] = max(0, alphas[i] - runif(1) * sigmas[i])
			alphas.right[i] = max(0, alphas.left[i] + sigmas[i])
		}	
		while (TRUE){
		alphas.new = sapply(1:length(alphas), function(i) runif(1, alphas.left[i], alphas.right[i]))
		if (f(alphas.new) > unew) break
		for (i in 1:length(alphas)) {
			if (alphas.new[i] < alphas[i]) {
				alphas.left[i] = max(0, alphas.new[i])
			} else {
				alphas.right[i] = max(0, alphas.new[i])
			}
		}
		}
		alphas = alphas.new
	}
	return(alphas)
}


f = function(alphas) {
	alpha0 = alphas[1]
	alpha1 = alphas[2]
	alpha = alphas[3]
	#prior = sum(dgamma(alphas, 1,1))
	prob = 0
	for (d in 1:D) {
      ratio2 = (table.k+alpha0/K)/(sum(table.k)+alpha0)
      ratio1 = (table.cd[,cd[d]]+alpha1*ratio2) / (sum(table.cd[,cd[d]])+alpha1)
      prob = prob + sum(log((table.dk[,d]+alpha*ratio1) /(sum(table.dk[,d])+alpha)))
  }
	return(exp(prob))
}