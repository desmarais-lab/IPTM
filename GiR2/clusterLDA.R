library(MCMCpack)
alpha = 10
alpha0 = 10
alpha1 = 10
K = 4
D = 20
#assume IP is fixed
cd = rbinom(D, 1, c(0.5,0.5))+1

#generate topics according to the generative process
m = rdirichlet(1, alpha0*rep(1/K, K))
mc = rdirichlet(2, alpha1*m)
theta = t(sapply(1:D, function(d){rdirichlet(1, alpha*mc[cd[d],])}))
n.d = 10
z = matrix(NA, D, n.d)
for (d in 1:D) {
  hi = rmultinom(n.d, 1, theta[d,])
  z[d,] = sapply(1:n.d, function(x) which(hi[,x]==1))
}
#this table.k.null is the true topic distribution across the corpus
table.k.null = tabulate(z, K)

#initialize N_k, N_dk, N_ck using the true values
table.k = tabulate(z, K)
table.dk = sapply(1:D, function(d) tabulate(z[d,], K))
table.cd = sapply(1:2, function(IP) {tabulate(z[which(cd==IP),],K)})

#inference on z using equation (15)---without word part that we ignored for simplified model
for (i in 1:100){
  for (d in 1:D) {
    for (n in 1:n.d) {
      #excluding \dn
      table.k[z[d,n]] = table.k[z[d,n]]-1
      table.dk[z[d,n],d] = table.dk[z[d,n],d]-1
      table.cd[z[d,n], cd[d]] = table.cd[z[d,n], cd[d]]-1
      #calculate Equation (15) and update z_dn
      ratio2 = (table.k+alpha0/K)/sum(table.k+alpha0/K)
      ratio1 = (table.cd[,cd[d]]+alpha1*ratio2) / sum(table.cd[,cd[d]]+alpha1*ratio2)
      prob = log(table.dk[,d]+alpha*ratio1)
      z[d,n] = which(rmultinom(1, 1, exp(prob))==1)
      #re-fill the tables 
      table.k[z[d,n]] = table.k[z[d,n]]+1
      table.dk[z[d,n],d] = table.dk[z[d,n],d]+1
      table.cd[z[d,n], cd[d]] = table.cd[z[d,n], cd[d]]+1
    }
  }
}

#see the results
rbind(table.k.null, table.k)
#even in this simplified model (only topic being random variable) and started with the true topic assignments, 
#we see that the topic assignments do not converges to true distribution
