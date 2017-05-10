prior_eta = function(eta, s2) {
  dnorm(eta, 0, s2, log = TRUE)
}
ll_delta = function(delta, lambda = 100, N = 10) {
  N * log(delta * lambda) - log(delta * lambda + 1)
}

eta = sort(rnorm(10^4, 0, 1))

D = 200
s2 = 0.25
par(mfrow= c(1,2))
plot(pnorm(eta), prior_eta(eta, s2) + D * ll_delta(pnorm(eta), lambda = 100, N = D))


#first derivative of posterior
dprior_eta = function(eta, s2) {
 -eta/ (s2)
}
dll_delta = function(delta, lambda = 100, N = 10) {
  N / delta - lambda / (delta * lambda + 1)
}

plot(pnorm(eta), dprior_eta(eta, s2) + dll_delta(pnorm(eta), lambda = 100, N = D))


par(mfrow=c(1,3))

eta = sort(rnorm(5*10^3, 0, 1))
post = function(eta) {
  dnorm(eta, 0, 1, log = TRUE) + sum(vapply(edge2, function(d) {
    EdgeInEqZ(iJi[[d]], lambda[[d]], pnorm(eta))
  }, c(1))) / length(edge2)
}
posteta = vapply(eta, function(i){post(i)}, c(1))
plot(pnorm(eta), posteta)

# this time, 
# delta = 0.5600206
# > summary(sapply(lambda, function(d){mean(d)}))
#Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
#0.000e+00 2.000e+00 7.427e+03 4.192e+32 1.203e+13 2.767e+34 
#> summary(sapply(iJi, function(d){sum(d)}))
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#1.00    3.00   10.00   19.17   28.00  127.00 


eta = sort(rnorm(5*10^3, 0, 1))
post = function(eta) {
  dnorm(eta, 0, 1, log = TRUE) + sum(vapply(edge2, function(d) {
    EdgeInEqZ(iJi[[d]], lambda[[d]], pnorm(eta))
  }, c(1))) / length(edge2)
}
posteta = vapply(eta, function(i){post(i)}, c(1))
plot(pnorm(eta), posteta)

#this time, 
# delta = 0.003549501
#> summary(sapply(lambda, function(d){mean(d)}))
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#0.01253 0.02075 0.03899 0.03891 0.05018 0.28284 
#> summary(sapply(iJi, function(d){sum(d)}))
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#1.000   1.000   1.000   1.514   1.000  12.000       2 


eta = sort(rnorm(5*10^3, 0, 1))
post = function(eta) {
  dnorm(eta, 0, 1, log = TRUE) + sum(vapply(edge2, function(d) {
    EdgeInEqZ(iJi[[d]], lambda[[d]], pnorm(eta))
  }, c(1))) / length(edge2)
}
posteta = vapply(eta, function(i){post(i)}, c(1))
plot(pnorm(eta), posteta)

#this time,
#delta = 2.545531e-05
#> summary(sapply(lambda, function(d){mean(d)}))
#Min.   1st Qu.    Median      Mean   3rd Qu.      Max. 
#0.000e+00 6.000e+00 2.168e+03 3.445e+15 9.863e+06 2.050e+17 
#> summary(sapply(iJi, function(d){sum(d)}))
#Min. 1st Qu.  Median    Mean 3rd Qu.    Max.    NA's 
#1.000   1.000   1.000   1.473   1.000  12.000       3 
