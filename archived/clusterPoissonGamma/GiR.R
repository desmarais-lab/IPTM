library(MCMCpack)

generate = function(D, V, C, K, a_C, b, e, xi) {
	 #psi = sapply(1:C, function(c) rgamma(1, gamma0/C, xi))
	 #cd = sapply(1:D, function(d) which(rmultinom(1, 1, psi)==1))
	 psi = rdirichlet(1, xi * rep(1/C, C))
	 cd = sapply(1:D, function(d) which(rmultinom(1, 1, psi)==1))
	 pi.dc = matrix(0, D, C)
	 for (d in 1:D) {
	 	pi.dc[d,cd[d]] = rgamma(1, a_C[cd[d]], b)
	 }
	theta.ck = matrix(NA, C, K)
	for (c in 1:C) {
		theta.ck[c,] = rgamma(K, e, e)
	}
	phi.kv = matrix(NA, K, V)
	for (k in 1:K) {
		phi.kv[k,] = rgamma(V, e, e)
	}
	W.dv = matrix(NA, D, V)
	for (d in 1:D) {
		for (v in 1:V) {
			W.dv[d,v] = rpois(1, pi.dc[d,]%*%theta.ck%*%phi.kv[,v])
		}
	}
	return(list(cd = cd, W.dv = W.dv, pi.dc = pi.dc, theta.ck = theta.ck, phi.kv = phi.kv))
}	

update = function(cd, W.dv, pi.dc, theta.ck, phi.kv, a_C, b, e, xi, niter) {
	D = nrow(W.dv)
	V = ncol(W.dv)
	C = nrow(theta.ck)
	K = ncol(theta.ck)
	for (n in 1:niter) {
		alloc = allocate_naive(W.dv, pi.dc, theta.ck, phi.kv)
		W.dc = alloc$W.dc
		W.ck = alloc$W.ck
		W.kv = alloc$W.kv
		for (d in 1:D) {
			pi.dc[d,cd[d]] = rgamma(1, a_C[cd[d]] + W.dc[d,cd[d]], b + theta.ck[cd[d],] %*% rowSums(phi.kv))
			# if (sum(pi.dc[d,]) == 0) browser()
			# prob = log(tabulate(cd[-d], C) + xi/C) + sapply(1:C, function(c) {
				# dgamma(pi.dc[d,cd[d]], a_C[c] + W.dc[d,cd[d]], b + theta.ck[c,] %*% rowSums(phi.kv), log = TRUE)
				# })
			# cdnew = which(rmultinom(1, 1, exp(prob - max(prob)))==1)
			# if (cdnew != cd[d]) {
				# pi.dc[d, cdnew] = pi.dc[d, cd[d]]
				# pi.dc[d, -cdnew] = 0
				# cd[d] = cdnew
			# }
		}
		theta.ck = matrix(rgamma(C*K, e + W.ck, e + outer(colSums(pi.dc),rowSums(phi.kv))), C, K)
		for (v in 1:V) {
			phi.kv[,v] = rgamma(K, e + W.kv[,v], e + colSums(pi.dc)%*%theta.ck)
		}
	}
	return(list(cd = cd, W.dv = W.dv, pi.dc = pi.dc, theta.ck = theta.ck, phi.kv = phi.kv))
}


allocate_naive = function(W.dv, pi.dc, theta.ck, phi.kv) {
	D = nrow(W.dv)
	V = ncol(W.dv)
	C = nrow(theta.ck)
	K = ncol(theta.ck)
	W.dc = matrix(0, D, C)
	W.ck = matrix(0, C, K)
	W.kv = matrix(0, K, V)
	for (d in 1:D) {
		for (v in 1:V) {
			if (W.dv[d,v] > 0) {
				p.ck = outer(pi.dc[d,], phi.kv[,v]) * theta.ck		
				W.dvck = matrix(rmultinom(1, W.dv[d,v], p.ck), C, K)
				W.dc[d, ] = W.dc[d, ] + rowSums(W.dvck)
				W.ck = W.ck + W.dvck
				W.kv[,v] = W.kv[,v] + colSums(W.dvck)
			}
		}
	}	
	return(list(W.dc = W.dc, W.ck = W.ck, W.kv = W.kv))
}

D = 5
V = 6
C = 2
K = 4
e = 0.1
xi = 5
niter = 5
alpha = 1
beta = 0.01
Schein = function(D, V, C, K, alpha, beta, e, xi, niter, nsamp) {
	results = matrix(NA, nrow = nsamp, ncol = 14)
	for (n in 1:nsamp) {
		if (n %% 100 == 0) print(n)
		a_C = rgamma(C, alpha, beta)
		b = 1
		forward = generate(D, V, C, K, a_C, b, e, xi)
		backward = update(forward$cd, forward$W.dv, forward$pi.dc, forward$theta.ck, forward$phi.kv, a_C, b, e, xi, niter)
		results[n, 1:7] = c(mean(forward$cd), mean(c(forward$pi.dc)), var(c(forward$pi.dc)), mean(c(forward$theta.ck)), var(c(forward$theta.ck)), mean(c(forward$phi.kv)), var(c(forward$phi.kv)))
		results[n, 8:14] = c(mean(backward$cd), mean(c(backward$pi.dc)), var(c(backward$pi.dc)), mean(c(backward$theta.ck)), var(c(backward$theta.ck)), mean(c(backward$phi.kv)), var(c(backward$phi.kv)))
	}
	Forward_stats = results[,1:7]
	colnames(Forward_stats) = c("cd","mean(pi)", "var(pi)","mean(theta)", "var(theta)", "mean(phi)", "var(phi)")
	Backward_stats = results[,8:14]
	colnames(Backward_stats) = colnames(Forward_stats)
	par(mfrow = c(2,4))
	GiR_PP_Plots(Forward_stats, Backward_stats)
	return(results)
}
#' @title GiR_PP_Plots
#' @description Generate PP-plots (probability-probability plot) for Gettig it Right test
#'
#' @param Forward_stats statistics obtained from forward sampling
#' @param Backward_stats statistics obtained from backward sampling
#'
#' @return PP-plots for different GiR statistics of interest
#'
#' @export
GiR_PP_Plots = function(Forward_stats, Backward_stats) {
  nms = colnames(Forward_stats)
  
  for (i in 1:ncol(Forward_stats)) {
    all = c(Backward_stats[, i], Forward_stats[, i])
    
    quantiles = 100
    uniqueValues = quantile(all,seq(0, 1, length = quantiles))
    qx1 = numeric(length(uniqueValues))
  	qx2 = numeric(length(uniqueValues))
  		
  	for (j in 1:length(uniqueValues)) {
  		qx1[j] = mean(Forward_stats[, i] <= uniqueValues[j])
  		qx2[j] = mean(Backward_stats[, i] <= uniqueValues[j])
  	}
    
    qqplot(x = qx1,
           y = qx2,
           ylim = c(0, 1),
           xlim = c(0, 1),
           ylab = "Backward",
           xlab = "Forward",
           col = "blue",
           pch = 19,
           cex = 0.25,
           main = nms[i],
           cex.lab = 0.25,
           cex.axis = 0.25,
           cex.main = 0.5)
    abline(0, 1, lty = 1, col = "red", lwd = 1)
    
      if (nrow(Forward_stats) > 10000) {
       thinning2 = seq(from = floor(nrow(Forward_stats) / 10), to = nrow(Forward_stats), length.out = 10000)
       Forward_test2 = Forward_stats[thinning2, i]
       Backward_test2 = Backward_stats[thinning2, i]
       } else {
        Forward_test2 = Forward_stats[, i]
        Backward_test2 = Backward_stats[, i]    	
      }
    text(paste("Backward Mean:", round(mean(Backward_stats[, i]), 4),
                "\nForward Mean:", round(mean(Forward_stats[, i]), 4),
                "\nt-test p-value:", round(t.test(Backward_test2, Forward_test2)$p.value, 4),
                "\nMann-Whitney p-value:", round(wilcox.test(Backward_test2, Forward_test2)$p.value, 4)),
                x = 0.65, y = 0.15, cex = 0.4)
  }
}      

test = Schein(D, V, C, K, alpha, beta, e, xi, 5, 20000)