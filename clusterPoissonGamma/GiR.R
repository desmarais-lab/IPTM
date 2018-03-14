library(MCMCpack)

generate = function(D, V, C, K, a_C, b_C, e, xi) {
	 ksi = rdirichlet(1, xi*rep(1/C, C))
	 cd = sapply(1:D, function(d) which(rmultinom(1, 1, ksi)==1))
	#cd = sample(1:C, D, TRUE)
	theta.ck = matrix(NA, C, K)
	for (c in 1:C) {
		theta.ck[c,] = rgamma(K, a_C[c], b_C[c])
	}
	phi.kv = matrix(NA, K, V)
	for (k in 1:K) {
		phi.kv[k,] = rgamma(V, e, e)
	}
	W.dv = matrix(NA, D, V)
	for (d in 1:D) {
		for (v in 1:V) {
			W.dv[d,v] = rpois(1, sum(theta.ck[cd[d],]*phi.kv[,v]))
		}
	}
	return(list(W.dv = W.dv, cd = cd, theta.ck = theta.ck, phi.kv = phi.kv))
}	

update = function(W.dv, cd, theta.ck, phi.kv, a_C, b_C, e, xi, niter) {
	D = nrow(W.dv)
	V = ncol(W.dv)
	C = nrow(theta.ck)
	K = ncol(theta.ck)
	for (n in 1:niter) {
		alloc = allocate_naive(W.dv, cd, theta.ck, phi.kv)
		W.ck = alloc$W.ck
		W.kv = alloc$W.kv
		W.dk = alloc$W.dk
		for (c in 1:C) {
			theta.ck[c,] = rgamma(K, a_C[c] + W.ck[c,], b_C[c]+sum(cd==c)*rowSums(phi.kv))
		}
		theta.dk = t(sapply(1:D, function(d) {theta.ck[cd[d],]}))
		for (k in 1:K) {
			phi.kv[k,] = rgamma(V, e + W.kv[k,], e + colSums(theta.dk))
		}
		for (d in 1:D) {
			W.ck[cd[d],] = W.ck[cd[d],] - W.dk[d,]
			prior = log(tabulate(cd[-d], C)+ xi/C)
			thetapart = sapply(1:C, function(c) {
				cd[d] = c; 
				W.ck[cd[d],] = W.ck[cd[d],]+W.dk[d,];
				sum(sapply(1:C, function(c) {
					sum(dgamma(theta.ck[c,], a_C[c] + W.ck[c,], b_C[c]+sum(cd==c)*rowSums(phi.kv), log = TRUE))}))
				})	
			Wpart = sapply(1:C, function(c) {sum(sapply(1:V, function(v) {dpois(W.dv[d,v], sum(theta.ck[c,]*phi.kv[,v]) , log = TRUE)}))	})
			prob = prior+Wpart+thetapart
			prob = exp(prob - max(prob))
			if (sum(is.na(prob)) > 0) browser()
			cd[d] = which(rmultinom(1, 1, prob)==1)
			W.ck[cd[d],] = W.ck[cd[d],] + W.dk[d,]
		}
	}
	return(list(W.dv = W.dv, cd = cd, theta.ck = theta.ck, phi.kv = phi.kv))
}


allocate_naive = function(W.dv, cd, theta.ck, phi.kv) {
	D = nrow(W.dv)
	V = ncol(W.dv)
	C = nrow(theta.ck)
	K = ncol(theta.ck)
	W.dk = matrix(0, D, K)
	W.ck = matrix(0, C, K)
	W.kv = matrix(0, K, V)
	for (d in 1:D) {
		for (v in 1:V) {
			if (W.dv[d,v] > 0) {
				p.ck = theta.ck[cd[d],] * phi.kv[,v]
				W.dvk = tabulate(sapply(1:W.dv[d,v], function(x) which(rmultinom(1, 1, p.ck)==1)), K)
				W.dk[d, ] = W.dk[d, ] + W.dvk
				W.ck[cd[d], ] = W.ck[cd[d], ] + W.dvk
				W.kv[,v] = W.kv[,v] + W.dvk
			}
		}
	}	
	return(list(W.dk = W.dk, W.ck = W.ck, W.kv = W.kv))
}

D = 10
V = 6
C = 2
K = 2
e = 0.1
xi = 5
niter = 5
alpha = 1
beta = 0.01
Schein = function(D, V, C, K, alpha, beta, e, xi, niter, nsamp) {
	results = matrix(NA, nrow = nsamp, ncol = 10)
	for (n in 1:nsamp) {
		if (n %% 100 == 0) print(n)
		a_C = rgamma(C, alpha, beta)
		#a_C = rep(10, C)
		b_C = rep(1, C)
		forward = generate(D, V, C, K, a_C, b_C, e, xi)
		backward = update(forward$W.dv, forward$cd, forward$theta.ck, forward$phi.kv, a_C, b_C, e, xi, niter)
		results[n, 1:5] = c(mean(forward$cd), mean(c(forward$theta.ck)), var(c(forward$theta.ck)), mean(c(forward$phi.kv)), var(c(forward$phi.kv)))
		results[n, 6:10] = c(mean(backward$cd), mean(c(backward$theta.ck)), var(c(backward$theta.ck)), mean(c(backward$phi.kv)), var(c(backward$phi.kv)))
	}
	Forward_stats = results[,1:5]
	colnames(Forward_stats) = c("cd", "mean(theta)", "var(theta)", "mean(phi)", "var(phi)")
	Backward_stats = results[,6:10]
	colnames(Backward_stats) = colnames(Forward_stats)
	par(mfrow = c(2,3))
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
    if (grepl("b_", nms[i]) ) {
      quantiles = 1000
    }
    if (grepl("eta_", nms[i]) ) {
      quantiles = 1000
    }
    if (grepl("delta", nms[i]) ) {
      quantiles = 1000
    }
    if (grepl("sigma_tau", nms[i]) ) {
      quantiles = 1000
    }
    
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

test = Schein(D, V, C, K, alpha, beta, e, xi, 5, 10000)