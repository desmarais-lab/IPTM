#' @useDynLib IPTMnew
#' @import stats
#' @import grDevices
#' @import graphics
#' @import RcppArmadillo
#' @importFrom Rcpp sourceCpp
#' @importFrom reshape melt
#' @importFrom coda mcmc geweke.diag
#' @importFrom combinat permn
#' @importFrom mgcv uniquecombs
#' @importFrom lubridate wday hour
#' @importFrom LaplacesDemon dhalfcauchy rhalfcauchy
#' @importFrom truncnorm rtruncnorm dtruncnorm

tvapply = function(...) transpose(vapply(...))

#' @title gibbs.measure.support
#' @description List out the support of Gibbs measure
#'
#' @param n length of the vector to be sampled 
#'
#' @return  a 2^n x n binary matrix representing the support of the binary Gibbs measure in n elements
#'
#' @export
gibbs.measure.support = function(n) {
	gibbs.support = rbind(rep(1, n))
	for(i in 1:(n-1)){
		gibbs.mat.i = do.call('rbind',permn(c(rep(1, i), rep(0, n-i))))
		gibbs.support = rbind(gibbs.support, gibbs.mat.i)
	}
	out = as.matrix(unique(gibbs.support))
	return(out)
}

#' @title r.gibbs.measure
#' @description List out the support of Gibbs measure
#'
#' @param nsamp the number of binary vector samples to draw
#' @param lambda.i a vector of coefficients according to which each element
#' @param delta a positive real valued parameter that controls the density penalty 
#' @param support support of Gibbs measure
#'
#' @return nsamp number of samples with each row denoting binary vector
#'
#' @export
r.gibbs.measure <- function(nsamp, lambda.i, delta, support) {
	#gibbsNormalizer = prod(exp(delta+lambda.i)+1)-1
	logitNumerator = vapply(1:nrow(support), function(s) {
		sum((delta+lambda.i)*support[s,])
		}, c(1))		
	samp = lmultinom(logitNumerator)
	return(support[samp,])	
}

#' @title adaptive.MH 
#' @description adaptive Metropolis Hastings to maintain target acceptance rate
#'
#' @param sigma.Q proposal distribution variance parameter
#' @param accept.rates acceptance rate from previous iteration
#' @param target target acceptance rate
#' @param update.size size of update to be adjusted 
#' @param tol tolerance level to determine if acceptance rate is too high
#'
#' @return nsamp number of samples with each row denoting binary vector
#'
#' @export
adaptive.MH = function(sigma.Q, accept.rates, target = 0.1, update.size, tol = 0.8) {
	for (i in 1:length(sigma.Q)) {
	  if (accept.rates[i] < target) {
			sigma.Q[i] = sigma.Q[i]-update.size[i]
		}
		if (accept.rates[i] > (target+tol)) {
			sigma.Q[i] = sigma.Q[i]+update.size[i]
		}		
	}
	return(sigma.Q)
}

#' @title IPTM.inference
#' @description Iterate Markov Chain Monte Carlo (MCMC) algorithm for the interaction-partitioned topic model
#'
#' @param edge list of tie data with 3 elements (1: author, 2: recipient, 3: timestamp in unix.time format)
#' @param node vector of node id's (ID starting from 1)
#' @param textlist list of text containing the words in each document
#' @param vocab all vocabularies used over the corpus
#' @param nIP total number of interaction patterns
#' @param K total number of topics
#' @param sigma.Q proposal distribution variance parameter
#' @param alphas Dirichlet concentration prior for document-topic distribution (alpha0, alpha1, alpha)
#' @param beta Dirichlet concentration prior for topic-word distribution
#' @param zeta Dirichlet concentration prior for document-interaction-pattern distribution
#' @param prior.b prior mean and covariance of b in multivariate Normal distribution
#' @param prior.delta prior mean and variance of delta in Normal distribution
#' @param prior.eta prior mean and covariance of eta in multivariate Normal distribution
#' @param prior.tau prior shape and scale parameter of sigma_tau in inverse-Gamma distribution
#' @param Outer size of outer iterations 
#' @param Inner size of inner iteration for Metropolis-Hastings updates
#' @param netstat which type of network statistics to use ("dyadic", "triadic", "degree")
#' @param timestat additional statistics to be used for timestamps other than netstat ("sender", "receiver","timeofday", "dayofweek")
#' @param initial list of initial values user wants to assign including (delta, b, eta, cd, z, u, sigma_tau, proposal.var1, proposal.var2)
#' @param timeunit hour (= 3600) or day (=3600*24) and so on
#' @param tz timezone such as EST and PST
#'
#' @return MCMC output containing all parameter estimates
#'
#' @export

IPTM.inference = function(edge, node, textlist, vocab, nIP, K, sigma.Q, alphas, beta, zeta,
                          prior.b, prior.delta, prior.eta, prior.tau, Outer, Inner,
                          netstat, timestat, initial = NULL, timeunit = 3600, tz = "America/New_York") {
  
  # trim the edge so that we only model edges after 384 hours
  A = length(node)
  D = length(edge)
  timestamps = vapply(edge, function(d) { d[[3]] }, c(1))
  senders = vapply(edge, function(d) { d[[1]] }, c(1))
  edge.trim = which_num(384*timeunit, timestamps-timestamps[1]):D
  max.edge = max(edge.trim)
  timeinc = c(timestamps[1], timestamps[-1]-timestamps[-length(timestamps)])/timeunit
  timeinc[timeinc==0] = runif(sum(timeinc==0), 0, min(timeinc[timeinc!=0]))
  emptytext = which(sapply(textlist, function(d){length(d)})==0)
  # initialization 
  convergence = c()
  topicpart = c()
  edgepart = c()
  netstat = as.numeric(c("degree", "dyadic", "triadic") %in% netstat)
  timestat = as.numeric(c("dayofweek","timeofday") %in% timestat)
  timemat = matrix(0, nrow = D, ncol = sum(timestat))
  if (sum(timestat) > 0) {
    Sys.setenv(TZ = tz)
    time_ymd = as.POSIXct(timestamps, tz = getOption("tz"), origin = "1970-01-01")
    if (timestat[1] > 0) {
      days = vapply(time_ymd, function(d) {wday(d)}, c(1))
      days[days==1] = 8
      timemat[,1] = as.numeric(cut(days, c(1,6,8), c("weekdays","weekends")))-1
      it = 1
    }
    if (timestat[2] > 0) {
      hours = vapply(time_ymd, function(d) {hour(d)}, c(1))
      timemat[,it+1] = as.numeric(cut(hours, c(-1,12,24), c("AM", "PM")))-1
    }     
  }
  L = 3
  P = L*(2*netstat[1]+2*netstat[2]+4*netstat[3])
  Q = length(prior.eta[[1]])
  proposal.var1 = diag(P)
  proposal.var2 = diag(Q)
  V = length(vocab)
  phi = matrix(NA, K, V)
  for (k in 1:K) {
    phi[k,] = rdirichlet_cpp(1, rep(beta/V, V))
  }
  alpha0 = alphas[1]
  alpha1 = alphas[2]
  alpha = alphas[3]
  m = rdirichlet_cpp(1, rep(alpha0/K, K))
  mc = matrix(NA, nIP, K)
  for (IP in 1:nIP) {
    mc[IP,] = rdirichlet_cpp(1, alpha1*m)
  }
  psi = rdirichlet_cpp(1, rep(zeta/nIP, nIP))
  if (length(initial) == 0) {
    cd = multinom_vec(D, psi) 
    theta = tvapply(seq(along = edge), function(d) {rdirichlet_cpp(1, alpha*mc[cd[d],]) }, rep(0, K))
    delta = rnorm(1, prior.delta[1], sqrt(prior.delta[2]))
    sigma_tau = rhalfcauchy(1, prior.tau)
    b.old = rmvnorm_arma(nIP, prior.b[[1]], prior.b[[2]])
    eta.old = rmvnorm_arma(nIP, prior.eta[[1]], prior.eta[[2]])
    z = lapply(seq(along = edge), function(d) multinom_vec(max(1, length(textlist[[d]])), theta[d, ]))
    sigma.Q = sigma.Q
    u = list()
    for (d in edge.trim) {
      u[[d]] = matrix(rbinom(A^2, 1, 1/A), nrow = A, ncol = A)
      diag(u[[d]]) = 0
      u[[d]][senders[d],] = tabulateC(as.numeric(unlist(edge[[d]][2])), A)
    }
  } else {
    delta = initial$delta
    sigma_tau = initial$sigma_tau
    b.old = initial$b
    eta.old = initial$eta
    cd = initial$cd
    z = initial$z
    sigma.Q = initial$sigma.Q
    u = initial$u
  }						 
  bmat = list()
  etamat = list()
  for (IP in 1:nIP) {
    bmat[[IP]] = matrix(b.old[IP,], nrow = P, ncol = 1)
    etamat[[IP]] = matrix(eta.old[IP,], nrow = Q, ncol = 1)
  }
  deltamat = rep(delta, 1)
  sigma_taumat = rep(sigma_tau, 1)		
  mu = matrix(0, nrow = D, ncol = A)
  textlist.raw = unlist(textlist)
  accept.rates = rep(0, 4)
  hist.d = c()
  for (d in 1:D) {
  if (timestamps[d]+384*timeunit > timestamps[max.edge]) {
        hist.d[d] = max.edge
      } else {
        hist.d[d] = which_num(timestamps[d]+384*timeunit, timestamps)
      }
  }
  timeinterval = timefinder(timestamps, edge.trim, timeunit)
  X = list()
  for (d in edge.trim) {
    X[[d]] = Netstats_cpp(edge, timestamps, timeinterval[[d]], senders, cd, A, timeunit, netstat)
  }
  
  table.W = tvapply(1:K, function(k) {
      tabulateC(textlist.raw[which(unlist(z[-emptytext]) == k)], V)
        }, rep(0, V))
  totalN = sum(table.W)-1 
  table.C = tabulateC(cd, nIP)
  table.cd = vapply(1:nIP, function(IP) {
    if (sum(cd==IP) > 0) {
      tabulateC(unlist(z[which(cd == IP)]), K)
    } else {
      rep(0, K)
    }  
  }, rep(0, K))
  table.dk = vapply(1:D, function(d) {
      tabulateC(unlist(z[[d]]), K)
        }, rep(0, K))
  table.k = rowSums(table.dk)     
  #start outer iteration
  for (o in 1:Outer) {
    print(o)
    # Data augmentation
    for (d in edge.trim) {
        lambda = MultiplyXB(X[[d]], b.old[cd[d],])
        for (i in node[-senders[d]]) {
            for (j in sample(node[-i], A-1)) {
                probij = u_Gibbs(u[[d]][i, ], lambda[i,], delta, j)
                u[[d]][i, j] = lmultinom(probij)-1
            }
        }
    }
    # cd update	
     for (d in 1:D) {
     	table.C[cd[d]] = table.C[cd[d]]-1
     	table.cd[, cd[d]] = table.cd[, cd[d]] - tabulateC(z[[d]], K)
        for (IP in 1:nIP) {
         cd[d] = IP
         IPpart = log(table.C[IP] + 1 + zeta/nIP)
       	 Xnew =  Netstats_cpp(edge, timestamps, timeinterval[[d]], senders, cd, A, timeunit, netstat)
       	 Edgepart = Edgepartsum(Xnew, b.old[IP,], u[[d]], delta)
       	 munew = mu_vec(timemat[d,], A, eta.old[IP,])
         Timepart = Timepart(munew, sigma_tau, senders[d], timeinc[d])
         Topicpart = Topicpart(K, z[[d]], table.cd[,IP], table.k-tabulateC(z[[d]], K), alphas)
         const.C[IP] = Edgepart+Timepart+Topicpart+IPpart
       }
       cd[d] = lmultinom(const.C)
       table.C[cd[d]] = table.C[cd[d]] + 1
       table.cd[, cd[d]] = table.cd[, cd[d]] + tabulateC(z[[d]], K)
    }
    for (d in edge.trim) {
         X[[d]] = Netstats_cpp(edge, timestamps, timeinterval[[d]], senders, cd, A, timeunit, netstat)
	 }
      
	 #Z update
	 topicsum = 0
     for (d in 1:D) {   
       textlist.d = textlist[[d]]
 	   for (w in 1:length(z[[d]])) {
 	   		zw.old = z[[d]][w]
 	   		table.cd[zw.old, cd[d]] = table.cd[zw.old, cd[d]]-1
 	   		table.dk[zw.old, d] = table.dk[zw.old, d]-1
 	   		table.k[zw.old] = table.k[zw.old]-1
         if (length(textlist.d) > 0) {
             table.W[zw.old, textlist.d[w]] = table.W[zw.old, textlist.d[w]]-1
        	 topicword.d = TopicWord(K, table.dk[,d], table.W[,textlist.d[w]], table.cd[,cd[d]], table.k, totalN,
        	  						alphas, beta, V)
            zw.new = lmultinom(topicword.d)
            if (zw.new != zw.old) {
               z[[d]][w] = zw.new
           }
           table.W[z[[d]][w], textlist.d[w]] = table.W[z[[d]][w], textlist.d[w]]+1
           } else {
           topicword.d = TopicWord0(K, table.cd[,cd[d]], table.k, totalN, alphas, beta, V)
           zw.new = lmultinom(topicword.d)
           if (zw.new != zw.old) {
               z[[d]][w] = zw.new
           }
         }
        table.cd[z[[d]][w], cd[d]] = table.cd[z[[d]][w], cd[d]]+1
        table.dk[z[[d]][w], d] = table.dk[z[[d]][w], d]+1
 	    table.k[z[[d]][w]] = table.k[z[[d]][w]]+1
 	    topicsum = topicsum + topicword.d[z[[d]][w]]
	 	}
	 }

  # adaptive M-H   
    #if (o > 1) {
    #	accept.rates[1] = accept.rates[1]/Inner[1]
    #	accept.rates[2] = accept.rates[2]/Inner[2]
    #  accept.rates[3] = accept.rates[3]/Inner[1]
    #  accept.rates[4] = accept.rates[1]/Inner[3]
    #	sigma.Q = adaptive.MH(sigma.Q, accept.rates, update.size = 0.2*sigma.Q)
    #}
    #accept.rates = rep(0, 4)
    
    prior.old1 = priorsum(prior.b[[2]], prior.b[[1]], b.old)+
    			 dnorm(delta, prior.delta[1], sqrt(prior.delta[2]), TRUE)
    post.old1 = Edgepartsum(X[[max.edge]], b.old[cd[max.edge],], u[[max.edge]], delta)
    b.new = matrix(NA, nIP, P)
    for (inner in 1:Inner[1]) {
      for (IP in 1:nIP) {
			  b.new[IP, ] = rmvnorm_arma(1, b.old[IP,], sigma.Q[1]*proposal.var1)
	  }
      delta.new = rnorm(1, delta, sqrt(sigma.Q[4]))
      prior.new1 = priorsum(prior.b[[2]], prior.b[[1]], b.new)+
    				 dnorm(delta.new, prior.delta[1], sqrt(prior.delta[2]), TRUE)
      post.new1 = Edgepartsum(X[[max.edge]], b.new[cd[max.edge],], u[[max.edge]], delta.new)
      loglike.diff = prior.new1+post.new1-prior.old1-post.old1
      if (log(runif(1, 0, 1)) < loglike.diff) {
        b.old = b.new
        delta = delta.new
        prior.old1 = prior.new1
        post.old1 = post.new1
      #  accept.rates[1] = accept.rates[1]+1
      }
        for (IP in 1:nIP) {
          bmat[[IP]] = cbind(bmat[[IP]], b.old[IP,])
        }
        deltamat = c(deltamat, delta)
    }
     
    mu = mu_mat(timemat, A, eta.old, cd)
	  prior.old2 = priorsum(prior.eta[[2]], prior.eta[[1]], eta.old)
	  post.old2 = Timepartsum(mu, sigma_tau, senders, timeinc)
      eta.new = matrix(NA, nIP, Q)
      for (inner in 1:Inner[2]) {
          for (IP in 1:nIP) {
			  eta.new[IP, ] = rmvnorm_arma(1, eta.old[IP,], sigma.Q[2]*proposal.var2)
		  }
      mu = mu_mat(timemat, A, eta.new, cd)
      prior.new2 = priorsum(prior.eta[[2]], prior.eta[[1]], eta.old)
      post.new2 = Timepartsum(mu, sigma_tau, senders, timeinc)
      loglike.diff = prior.new2+post.new2-prior.old2-post.old2
      if (log(runif(1, 0, 1)) < loglike.diff) {
        eta.old = eta.new
        prior.old2 = prior.new2
        post.old2 = post.new2
      #  accept.rates[2] = accept.rates[2]+1
      }
      for (IP in 1:nIP) {
          etamat[[IP]] = cbind(etamat[[IP]], eta.old[IP,])
        }
    }
      
    prior.old3 = dhalfcauchy(sigma_tau, prior.tau, TRUE)
    post.old3 = post.old2
    for (inner in 1:Inner[3]) {
      sigma_tau.new = rtruncnorm(1, 0, Inf, sigma_tau, sqrt(sigma.Q[3]))
      prior.new3 = dhalfcauchy(sigma_tau.new, prior.tau, TRUE)
      post.new3 = Timepartsum(mu, sigma_tau.new, senders, timeinc)
      loglike.diff = log(dtruncnorm(sigma_tau, 0, Inf, sigma_tau.new, sqrt(sigma.Q[3])))-
                   log(dtruncnorm(sigma_tau.new, 0, Inf, sigma_tau, sqrt(sigma.Q[3])))+
                   prior.new3+post.new3-prior.old3-post.old3
      if (log(runif(1, 0, 1)) < loglike.diff) {
        sigma_tau = sigma_tau.new
        prior.old3 = prior.new3
        post.old3 = post.new3
     #   accept.rates[3] = accept.rates[3]+1
      }
        sigma_taumat = c(sigma_taumat, sigma_tau)
    }
	
	#very long inner iterations at the end
    if (o == Outer) {
    	Inner = Inner * 1000
    	bmat = list()
 	   etamat = list()
         for (IP in 1:nIP) {
         bmat[[IP]] = matrix(NA, nrow = P, ncol = Inner[1])
         etamat[[IP]] = matrix(NA, nrow = Q, ncol = Inner[2])
         }
  	   deltamat = rep(NA, Inner[1])
 	   sigma_taumat = rep(NA, Inner[3])		
    
    	b.new = matrix(NA, nIP, P)
    	for (inner in 1:Inner[1]) {
      		for (IP in 1:nIP) {
			  b.new[IP, ] = rmvnorm_arma(1, b.old[IP,], sigma.Q[1]*proposal.var1)
	  		}
     		delta.new = rnorm(1, delta, sqrt(sigma.Q[4]))
      		prior.new1 = priorsum(prior.b[[2]], prior.b[[1]], b.new)+
    				 dnorm(delta.new, prior.delta[1], sqrt(prior.delta[2]), TRUE)
      		post.new1 = Edgepartsum(X[[max.edge]], b.new[cd[max.edge],], u[[max.edge]], delta)
      		loglike.diff = prior.new1+post.new1-prior.old1-post.old1
      		if (log(runif(1, 0, 1)) < loglike.diff) {
        		b.old = b.new
        		delta = delta.new
        		prior.old1 = prior.new1
        		post.old1 = post.new1
        	}
        	for (IP in 1:nIP) {
          		bmat[[IP]][,inner] = b.old[IP, ] 
        	}
        	deltamat[inner] = delta
    	}
     
      	eta.new = matrix(NA, nIP, Q)
      	for (inner in 1:Inner[2]) {
          	for (IP in 1:nIP) {
			  eta.new[IP, ] = rmvnorm_arma(1, eta.old[IP,], sigma.Q[2]*proposal.var2)
		  	}
      	mu = mu_mat(timemat, A, eta.new, cd)
      	prior.new2 = priorsum(prior.eta[[2]], prior.eta[[1]], eta.old)
      	post.new2 = Timepartsum(mu, sigma_tau, senders, timeinc)
      	loglike.diff = prior.new2+post.new2-prior.old2-post.old2
      	if (log(runif(1, 0, 1)) < loglike.diff) {
        	eta.old = eta.new
        	prior.old2 = prior.new2
        	post.old2 = post.new2
        }
      	for (IP in 1:nIP) {
           etamat[[IP]][,inner] = eta.old[IP, ]
        	}
    	}
      
    	for (inner in 1:Inner[3]) {
      		sigma_tau.new = rtruncnorm(1, 0, Inf, sigma_tau, sqrt(sigma.Q[3]))
      		prior.new3 = dhalfcauchy(sigma_tau.new, prior.tau, TRUE)
      		post.new3 = Timepartsum(mu, sigma_tau.new, senders, timeinc)
      		loglike.diff = log(dtruncnorm(sigma_tau, 0, Inf, sigma_tau.new, sqrt(sigma.Q[3])))-
                   log(dtruncnorm(sigma_tau.new, 0, Inf, sigma_tau, sqrt(sigma.Q[3])))+
                   prior.new3+post.new3-prior.old3-post.old3
      		if (log(runif(1, 0, 1)) < loglike.diff) {
        		sigma_tau = sigma_tau.new
        		prior.old3 = prior.new3
        		post.old3 = post.new3
      		}
       		sigma_taumat[inner] = sigma_tau
    	}
    }
    edgepart[o] = post.old1 + post.old3
    topicpart[o] = topicsum                 
	convergence[o] = topicpart[o] + edgepart[o]
  }
  chain.final = list(cd = cd, z = z, b = bmat, eta = etamat, delta = deltamat, sigma_tau = sigma_taumat,
                     u = u, sigma.Q =sigma.Q, edge.trim = edge.trim,
                     edgepart = edgepart, topicpart = topicpart, convergence = convergence)
  return(chain.final)
}	


#' @title IPTM.inference.min
#' @description Iterate Markov Chain Monte Carlo (MCMC) algorithm for the interaction-partitioned topic model
#'
#' @param edge list of tie data with 3 elements (1: author, 2: recipient, 3: timestamp in unix.time format)
#' @param node vector of node id's (ID starting from 1)
#' @param textlist list of text containing the words in each document
#' @param vocab all vocabularies used over the corpus
#' @param nIP total number of interaction patterns
#' @param K total number of topics
#' @param sigma.Q proposal distribution variance parameter
#' @param alphas Dirichlet concentration prior for document-topic distribution (alpha0, alpha1, alpha)
#' @param beta Dirichlet concentration prior for topic-word distribution
#' @param zeta Dirichlet concentration prior for document-interaction-pattern distribution
#' @param prior.b prior mean and covariance of b in multivariate Normal distribution
#' @param prior.delta prior mean and variance of delta in Normal distribution
#' @param prior.eta prior mean and covariance of eta in multivariate Normal distribution
#' @param prior.tau prior shape and scale parameter of sigma_tau in inverse-Gamma distribution
#' @param Outer size of outer iterations 
#' @param Inner size of inner iteration for Metropolis-Hastings updates
#' @param netstat which type of network statistics to use ("dyadic", "triadic", "degree")
#' @param timestat additional statistics to be used for timestamps other than netstat ("sender", "receiver","timeofday", "dayofweek")
#' @param initial list of initial values user wants to assign including (delta, b, eta, cd, z, u, sigma_tau, proposal.var1, proposal.var2)
#' @param timeunit hour (= 3600) or day (=3600*24) and so on
#' @param tz timezone such as EST and PST
#'
#' @return MCMC output containing all parameter estimates
#'
#' @export

IPTM.inference.min = function(edge, node, textlist, vocab, nIP, K, sigma.Q, alphas, beta, zeta,
                          prior.b, prior.delta, prior.eta, prior.tau, Outer, Inner,
                          netstat, timestat, initial = NULL, timeunit = 3600, tz = "America/New_York") {
  
  # trim the edge so that we only model edges after 384 hours
  A = length(node)
  D = length(edge)
  timestamps = vapply(edge, function(d) { d[[3]] }, c(1))
  senders = vapply(edge, function(d) { d[[1]] }, c(1))
  edge.trim = which_num(384*timeunit, timestamps-timestamps[1]):D
  max.edge = max(edge.trim)
  timeinc = c(timestamps[1], timestamps[-1]-timestamps[-length(timestamps)])/timeunit
  timeinc[timeinc==0] = runif(sum(timeinc==0), 0, min(timeinc[timeinc!=0]))
  emptytext = which(sapply(textlist, function(d){length(d)})==0)
  # initialization 
  convergence = c()
  topicpart = c()
  edgepart = c()
  netstat = as.numeric(c("degree", "dyadic", "triadic") %in% netstat)
  timestat = as.numeric(c("dayofweek","timeofday") %in% timestat)
  timemat = matrix(0, nrow = D, ncol = sum(timestat))
  if (sum(timestat) > 0) {
    Sys.setenv(TZ = tz)
    time_ymd = as.POSIXct(timestamps, tz = getOption("tz"), origin = "1970-01-01")
    if (timestat[1] > 0) {
      days = vapply(time_ymd, function(d) {wday(d)}, c(1))
      days[days==1] = 8
      timemat[,1] = as.numeric(cut(days, c(1,6,8), c("weekdays","weekends")))-1
      it = 1
    }
    if (timestat[2] > 0) {
      hours = vapply(time_ymd, function(d) {hour(d)}, c(1))
      timemat[,it+1] = as.numeric(cut(hours, c(-1,12,24), c("AM", "PM")))-1
    }     
  }
  L = 3
  P = L*(2*netstat[1]+2*netstat[2]+4*netstat[3])
  Q = length(prior.eta[[1]])
  proposal.var1 = diag(P)
  proposal.var2 = diag(Q)
  V = length(vocab)
  phi = matrix(NA, K, V)
  for (k in 1:K) {
    phi[k,] = rdirichlet_cpp(1, rep(beta/V, V))
  }
  alpha0 = alphas[1]
  alpha1 = alphas[2]
  alpha = alphas[3]
  m = rdirichlet_cpp(1, rep(alpha0/K, K))
  mc = matrix(NA, nIP, K)
  for (IP in 1:nIP) {
    mc[IP,] = rdirichlet_cpp(1, alpha1*m)
  }
  psi = rdirichlet_cpp(1, rep(zeta/nIP, nIP))
  if (length(initial) == 0) {
    cd = multinom_vec(D, psi) 
    theta = tvapply(seq(along = edge), function(d) {rdirichlet_cpp(1, alpha*mc[cd[d],]) }, rep(0, K))
    delta = rnorm(1, prior.delta[1], sqrt(prior.delta[2]))
    sigma_tau = rhalfcauchy(1, prior.tau)
    b.old = rmvnorm_arma(nIP, prior.b[[1]], prior.b[[2]])
    eta.old = rmvnorm_arma(nIP, prior.eta[[1]], prior.eta[[2]])
    z = lapply(seq(along = edge), function(d) multinom_vec(max(1, length(textlist[[d]])), theta[d, ]))
    sigma.Q = sigma.Q
    u = list()
    for (d in edge.trim) {
      u[[d]] = matrix(rbinom(A^2, 1, 1/A), nrow = A, ncol = A)
      diag(u[[d]]) = 0
      u[[d]][senders[d],] = tabulateC(as.numeric(unlist(edge[[d]][2])), A)
    }
  } else {
    delta = initial$delta
    sigma_tau = initial$sigma_tau
    b.old = initial$b
    eta.old = initial$eta
    cd = initial$cd
    z = initial$z
    sigma.Q = initial$sigma.Q
    u = initial$u
  }						 
  bmat = list()
  etamat = list()
  for (IP in 1:nIP) {
    bmat[[IP]] = matrix(b.old[IP,], nrow = P, ncol = 1)
    etamat[[IP]] = matrix(eta.old[IP,], nrow = Q, ncol = 1)
  }
  deltamat = rep(delta, 1)
  sigma_taumat = rep(sigma_tau, 1)		
  mu = matrix(0, nrow = D, ncol = A)
  textlist.raw = unlist(textlist)
  accept.rates = rep(0, 4)
  hist.d = c()
  for (d in 1:D) {
  if (timestamps[d]+384*timeunit > timestamps[max.edge]) {
        hist.d[d] = max.edge
      } else {
        hist.d[d] = which_num(timestamps[d]+384*timeunit, timestamps)
      }
  }
  timeinterval = timefinder(timestamps, edge.trim, timeunit)
  X = list()
  for (d in edge.trim) {
    X[[d]] = Netstats_cpp(edge, timestamps, timeinterval[[d]], senders, cd, A, timeunit, netstat)
  }
  table.W = tvapply(1:K, function(k) {
      tabulateC(textlist.raw[which(unlist(z[-emptytext]) == k)], V)
        }, rep(0, V))
  zuniq = lapply(z, function(x) {sortuniq(x)})            
  table.cd = vapply(1:nIP, function(IP) {
    if (sum(cd==IP) > 0) {
      tabulateC(unlist(zuniq[which(cd == IP)]), K)
    } else {
      rep(0, K)
    }  
  }, rep(0, K))
  table.dk = vapply(1:D, function(d) {
      tabulateC(unlist(z[[d]]), K)
        }, rep(0, K))
  table.k = rowSums(table.cd > 0) 
  totalN = nIP    
  #start outer iteration
  for (o in 1:Outer) {
    print(o)
    # Data augmentation
    for (d in edge.trim) {
        lambda = MultiplyXB(X[[d]], b.old[cd[d],])
        for (i in node[-senders[d]]) {
            for (j in sample(node[-i], A-1)) {
                probij = u_Gibbs(u[[d]][i, ], lambda[i,], delta, j)
                u[[d]][i, j] = lmultinom(probij)-1
            }
        }
    }
    # Z update	
    topicsum = 0
    for (d in 1:D) {
	   	textlist.d = textlist[[d]]
	   	for (w in 1:length(z[[d]])) {
	   		zw.old = z[[d]][w]
	   		table.dk[zw.old, d] = table.dk[zw.old, d]-1
	   		if (!identical(zuniq[[d]], sortuniq(z[[d]][-w]))) {
	   			zuniq[[d]] = sortuniq(z[[d]][-w])
	   			table.cd = vapply(1:nIP, function(IP) {
	   			  if (sum(cd == IP) > 0) {
	   			    tabulateC(unlist(zuniq[which(cd == IP)]), K)
	   			  } else {
	   			    rep(0, K)
	   			  }  
	   			}, rep(0, K))      
	   			table.k = rowSums(table.cd > 0) 
	   		}
	   	 if (length(textlist.d) > 0) {
          table.W[zw.old, textlist.d[w]] = table.W[zw.old, textlist.d[w]]-1
       	  topicword.d = TopicWord(K, table.dk[,d], table.W[,textlist.d[w]], table.cd[,cd[d]], table.k, totalN, 
       	  						 alphas, beta, V)
          zw.new = lmultinom(topicword.d)
          if (zw.new != zw.old) {
              z[[d]][w] = zw.new
          }
          table.W[z[[d]][w], textlist.d[w]] = table.W[z[[d]][w], textlist.d[w]]+1
        } else {
          topicword.d = TopicWord0(K, table.cd[,cd[d]], table.k, totalN, alphas, beta, V)
          zw.new = lmultinom(topicword.d)
          if (zw.new != zw.old) {
              z[[d]][w] = zw.new
          }
        }
		if (!identical(zuniq[[d]], sortuniq(z[[d]]))) {
	   			zuniq[[d]] = sortuniq(z[[d]])
	   			table.cd = vapply(1:nIP, function(IP) {
	   			  if (sum(cd == IP) > 0) {
	   			    tabulateC(unlist(zuniq[which(cd == IP)]), K)
	   			  } else {
	   			    rep(0, K)
	   			  }  
	   			}, rep(0, K))   
	   			table.k = rowSums(table.cd > 0) 
	   		}
       table.dk[z[[d]][w], d] = table.dk[z[[d]][w], d]+1
       topicsum = topicsum + topicword.d[z[[d]][w]] 
      }
    }

    # C update
    table.C = rep(0, nIP)
    const.C = rep(NA, nIP)
    for (d in 1:D) {
    	IPpart = log(table.C + zeta/nIP)
      for (IP in 1:nIP) {
        cd[d] = IP
       	Xnew =  Netstats_cpp(edge, timestamps, timeinterval[[hist.d[d]]], senders, cd, A, timeunit, netstat)
        Edgepart = Edgepartsum(Xnew, b.old[IP,], u[[hist.d[d]]], delta)
        munew = mu_vec(timemat[d,], A, eta.old[IP,])
        Timepart = Timepart(munew, sigma_tau, senders[d], timeinc[d])
        Topicpart = Topicpart(K, z[[d]], table.cd[,IP], table.k, alphas) 	  	
        const.C[IP] = Edgepart+Timepart+Topicpart+ IPpart[IP]
      }
      cd[d] = lmultinom(const.C)
      table.C[cd[d]] = table.C[cd[d]]+1
	}
	  for (d in edge.trim) {
        X[[d]] = Netstats_cpp(edge, timestamps, timeinterval[[d]], senders, cd, A, timeunit, netstat)
 	  }
    table.cd = vapply(1:nIP, function(IP) {
      if (sum(cd == IP) > 0) {
        tabulateC(unlist(zuniq[which(cd == IP)]), K)
      } else {
        rep(0, K)
      }  
    }, rep(0, K)) 
    table.k = rowSums(table.cd > 0) 
    
  # adaptive M-H   
  #  if (o > 1) {
  #  	accept.rates[1] = accept.rates[1]/5
  #  	accept.rates[2] = accept.rates[2]/5
  #    accept.rates[3] = accept.rates[3]/5
  #    accept.rates[4] = accept.rates[1]
  #  	sigma.Q = adaptive.MH(sigma.Q, accept.rates, update.size = 0.2*sigma.Q)
  #  }
  #  accept.rates = rep(0, 4)
    
    prior.old1 = priorsum(prior.b[[2]], prior.b[[1]], b.old)+
    			 dnorm(delta, prior.delta[1], sqrt(prior.delta[2]), TRUE)
    post.old1 = Edgepartsum(X[[max.edge]], b.old[cd[max.edge],], u[[max.edge]], delta)
    b.new = matrix(NA, nIP, P)
    for (inner in 1:Inner[1]) {
      for (IP in 1:nIP) {
			  b.new[IP, ] = rmvnorm_arma(1, b.old[IP,], sigma.Q[1]*proposal.var1)
	  }
      delta.new = rnorm(1, delta, sqrt(sigma.Q[4]))
      prior.new1 = priorsum(prior.b[[2]], prior.b[[1]], b.new)+
    				 dnorm(delta.new, prior.delta[1], sqrt(prior.delta[2]), TRUE)
      post.new1 = Edgepartsum(X[[max.edge]], b.new[cd[max.edge],], u[[max.edge]], delta)
      loglike.diff = prior.new1+post.new1-prior.old1-post.old1
      if (log(runif(1, 0, 1)) < loglike.diff) {
        b.old = b.new
        delta = delta.new
        prior.old1 = prior.new1
        post.old1 = post.new1
      #  accept.rates[1] = accept.rates[1]+1
      }
        for (IP in 1:nIP) {
          bmat[[IP]] = cbind(bmat[[IP]], b.old[IP,])
        }
        deltamat = c(deltamat, delta)
    }
     
      mu = mu_mat(timemat, A, eta.old, cd)
	  prior.old2 = priorsum(prior.eta[[2]], prior.eta[[1]], eta.old)
	  post.old2 = Timepartsum(mu, sigma_tau, senders, timeinc)
      eta.new = matrix(NA, nIP, Q)
      for (inner in 1:Inner[2]) {
          for (IP in 1:nIP) {
			  eta.new[IP, ] = rmvnorm_arma(1, eta.old[IP,], sigma.Q[2]*proposal.var2)
		  }
      mu = mu_mat(timemat, A, eta.new, cd)
      prior.new2 = priorsum(prior.eta[[2]], prior.eta[[1]], eta.old)
      post.new2 = Timepartsum(mu, sigma_tau, senders, timeinc)
      loglike.diff = prior.new2+post.new2-prior.old2-post.old2
      if (log(runif(1, 0, 1)) < loglike.diff) {
        eta.old = eta.new
        prior.old2 = prior.new2
        post.old2 = post.new2
      #  accept.rates[2] = accept.rates[2]+1
      }
      for (IP in 1:nIP) {
          etamat[[IP]] = cbind(etamat[[IP]], eta.old[IP,])
        }
    }
      
    prior.old3 = dhalfcauchy(sigma_tau, prior.tau, TRUE)
    post.old3 = post.old2
    for (inner in 1:Inner[3]) {
      sigma_tau.new = rtruncnorm(1, 0, Inf, sigma_tau, sqrt(sigma.Q[3]))
      prior.new3 = dhalfcauchy(sigma_tau.new, prior.tau, TRUE)
      post.new3 = Timepartsum(mu, sigma_tau.new, senders, timeinc)
      loglike.diff = log(dtruncnorm(sigma_tau, 0, Inf, sigma_tau.new, sqrt(sigma.Q[3])))-
                   log(dtruncnorm(sigma_tau.new, 0, Inf, sigma_tau, sqrt(sigma.Q[3])))+
                   prior.new3+post.new3-prior.old3-post.old3
      if (log(runif(1, 0, 1)) < loglike.diff) {
        sigma_tau = sigma_tau.new
        prior.old3 = prior.new3
        post.old3 = post.new3
      #  accept.rates[3] = accept.rates[3]+1
      }
        sigma_taumat = c(sigma_taumat, sigma_tau)
    }
	
	#very long inner iterations at the end
    if (o == Outer) {
    	Inner = Inner * 1000
    	bmat = list()
 	    etamat = list()
         for (IP in 1:nIP) {
         bmat[[IP]] = matrix(NA, nrow = P, ncol = Inner[1])
         etamat[[IP]] = matrix(NA, nrow = Q, ncol = Inner[2])
         }
  	   deltamat = rep(NA, Inner[1])
 	   sigma_taumat = rep(NA, Inner[3])		
    
    	b.new = matrix(NA, nIP, P)
    	for (inner in 1:Inner[1]) {
      		for (IP in 1:nIP) {
			  b.new[IP, ] = rmvnorm_arma(1, b.old[IP,], sigma.Q[1]*proposal.var1)
	  		}
     		delta.new = rnorm(1, delta, sqrt(sigma.Q[4]))
      		prior.new1 = priorsum(prior.b[[2]], prior.b[[1]], b.new)+
    				 dnorm(delta.new, prior.delta[1], sqrt(prior.delta[2]), TRUE)
      		post.new1 = Edgepartsum(X[[max.edge]], b.new[cd[max.edge],], u[[max.edge]], delta.new)
      		loglike.diff = prior.new1+post.new1-prior.old1-post.old1
      		if (log(runif(1, 0, 1)) < loglike.diff) {
        		b.old = b.new
        		delta = delta.new
        		prior.old1 = prior.new1
        		post.old1 = post.new1
        	}
        	for (IP in 1:nIP) {
          		bmat[[IP]][,inner] = b.old[IP, ] 
        	}
        	deltamat[inner] = delta
    	}
     
      	eta.new = matrix(NA, nIP, Q)
      	for (inner in 1:Inner[2]) {
          	for (IP in 1:nIP) {
			  eta.new[IP, ] = rmvnorm_arma(1, eta.old[IP,], sigma.Q[2]*proposal.var2)
		  	}
      	mu = mu_mat(timemat, A, eta.new, cd)
      	prior.new2 = priorsum(prior.eta[[2]], prior.eta[[1]], eta.old)
      	post.new2 = Timepartsum(mu, sigma_tau, senders, timeinc)
      	loglike.diff = prior.new2+post.new2-prior.old2-post.old2
      	if (log(runif(1, 0, 1)) < loglike.diff) {
        	eta.old = eta.new
        	prior.old2 = prior.new2
        	post.old2 = post.new2
        }
      	for (IP in 1:nIP) {
           etamat[[IP]][,inner] = eta.old[IP, ]
        	}
    	}
      
    	for (inner in 1:Inner[3]) {
      		sigma_tau.new = rtruncnorm(1, 0, Inf, sigma_tau, sqrt(sigma.Q[3]))
      		prior.new3 = dhalfcauchy(sigma_tau.new, prior.tau, TRUE)
      		post.new3 = Timepartsum(mu, sigma_tau.new, senders, timeinc)
      		loglike.diff = log(dtruncnorm(sigma_tau, 0, Inf, sigma_tau.new, sqrt(sigma.Q[3])))-
                   log(dtruncnorm(sigma_tau.new, 0, Inf, sigma_tau, sqrt(sigma.Q[3])))+
                   prior.new3+post.new3-prior.old3-post.old3
      		if (log(runif(1, 0, 1)) < loglike.diff) {
        		sigma_tau = sigma_tau.new
        		prior.old3 = prior.new3
        		post.old3 = post.new3
      		}
       		sigma_taumat[inner] = sigma_tau
    	}
    }
    edgepart[o] = post.old1 + post.old3
    topicpart[o] = topicsum                 
	convergence[o] = topicpart[o] + edgepart[o]
  }
  chain.final = list(cd = cd, z = z, b = bmat, eta = etamat, delta = deltamat, sigma_tau = sigma_taumat,
                     u = u, sigma.Q =sigma.Q, edge.trim = edge.trim,
                     edgepart = edgepart, topicpart = topicpart, convergence = convergence)
  return(chain.final)
}	

#' @title IPTM.inference.GiR2
#' @description Iterate Markov Chain Monte Carlo (MCMC) algorithm for the interaction-partitioned topic model
#'
#' @param edge list of tie data with 3 elements (1: author, 2: recipient, 3: timestamp in unix.time format)
#' @param node vector of node id's (ID starting from 1)
#' @param textlist list of text containing the words in each document
#' @param vocab all vocabularies used over the corpus
#' @param nIP total number of interaction patterns
#' @param K total number of topics
#' @param sigma.Q proposal distribution variance parameter
#' @param alphas Dirichlet concentration prior for document-topic distribution (alpha0, alpha1, alpha)
#' @param beta Dirichlet concentration prior for topic-word distribution
#' @param zeta Dirichlet concentration prior for document-interaction-pattern distribution
#' @param prior.b prior mean and covariance of b in multivariate Normal distribution
#' @param prior.delta prior mean and variance of delta in Normal distribution
#' @param prior.eta prior mean and covariance of eta in multivariate Normal distribution
#' @param prior.tau prior shape and scale parameter of sigma_tau in inverse-Gamma distribution
#' @param Outer size of outer iterations 
#' @param Inner size of inner iteration for Metropolis-Hastings updates
#' @param netstat which type of network statistics to use ("dyadic", "triadic", "degree")
#' @param timestat additional statistics to be used for timestamps other than netstat ("sender", "receiver","timeofday", "dayofweek")
#' @param initial list of initial values user wants to assign including (delta, b, eta, cd, z, u, sigma_tau, proposal.var1, proposal.var2)
#' @param timeunit hour (= 3600) or day (=3600*24) and so on
#' @param tz timezone such as EST and PST
#'
#' @return MCMC output containing all parameter estimates
#'
#' @export

IPTM.inference.GiR2 = function(edge, node, textlist, vocab, nIP, K, sigma.Q, alphas, beta, zeta,
                          prior.b, prior.delta, prior.eta, prior.tau, Outer, Inner,
                          netstat, timestat, initial = NULL, timeunit = 3600, tz = "America/New_York") {
  
  # trim the edge so that we only model edges after 384 hours
  A = length(node)
  D = length(edge)
  timestamps = vapply(edge, function(d) { d[[3]] }, c(1))
  senders = vapply(edge, function(d) { d[[1]] }, c(1))
  edge.trim = which_num(384*timeunit, timestamps-timestamps[1]):D
  max.edge = max(edge.trim)
  timeinc = c(timestamps[1], timestamps[-1]-timestamps[-length(timestamps)])/timeunit
  timeinc[timeinc==0] = runif(sum(timeinc==0), 0, min(timeinc[timeinc!=0]))
  emptytext = which(sapply(textlist, function(d){length(d)})==0)
  # initialization 
  netstat = as.numeric(c("degree", "dyadic", "triadic") %in% netstat)
  timestat = as.numeric(c("dayofweek","timeofday") %in% timestat)
  timemat = matrix(0, nrow = D, ncol = sum(timestat))
  if (sum(timestat) > 0) {
    Sys.setenv(TZ = tz)
    time_ymd = as.POSIXct(timestamps, tz = getOption("tz"), origin = "1970-01-01")
    if (timestat[1] > 0) {
      days = vapply(time_ymd, function(d) {wday(d)}, c(1))
      days[days==1] = 8
      timemat[,1] = as.numeric(cut(days, c(1,6,8), c("weekdays","weekends")))-1
      it = 1
    }
    if (timestat[2] > 0) {
      hours = vapply(time_ymd, function(d) {hour(d)}, c(1))
      timemat[,it+1] = as.numeric(cut(hours, c(-1,12,24), c("AM", "PM")))-1
    }     
  }
  L = 3
  P = L*(2*netstat[1]+2*netstat[2]+4*netstat[3])
  Q = length(prior.eta[[1]])
  proposal.var1 = diag(P)
  proposal.var2 = diag(Q)
  V = length(vocab)
  phi = matrix(NA, K, V)
  for (k in 1:K) {
    phi[k,] = rdirichlet_cpp(1, rep(beta/V, V))
  }
  alpha0 = alphas[1]
  alpha1 = alphas[2]
  alpha = alphas[3]
  if (length(initial) == 0) {
  	m = rdirichlet_cpp(1, rep(alpha0/K, K))
  	mc = matrix(NA, nIP, K)
  	for (IP in 1:nIP) {
    	mc[IP,] = rdirichlet_cpp(1, alpha1*m)
  	}
    cd = multinom_vec(D, psi) 
    theta = tvapply(seq(along = edge), function(d) {rdirichlet_cpp(1, alpha*mc[cd[d],]) }, rep(0, K))
    delta = rnorm(1, prior.delta[1], sqrt(prior.delta[2]))
    sigma_tau = rhalfcauchy(1, prior.tau)
    b.old = rmvnorm_arma(nIP, prior.b[[1]], prior.b[[2]])
    eta.old = rmvnorm_arma(nIP, prior.eta[[1]], prior.eta[[2]])
    z = lapply(seq(along = edge), function(d) multinom_vec(max(1, length(textlist[[d]])), theta[d, ]))
    sigma.Q = sigma.Q
    psi = rdirichlet_cpp(1, rep(zeta/nIP, nIP))
    u = list()
    for (d in edge.trim) {
      u[[d]] = matrix(rbinom(A^2, 1, 1/A), nrow = A, ncol = A)
      diag(u[[d]]) = 0
      u[[d]][senders[d],] = tabulateC(as.numeric(unlist(edge[[d]][2])), A)
    }
  } else {
    delta = initial$delta
    sigma_tau = initial$sigma_tau
    b.old = initial$b
    eta.old = initial$eta
    cd = initial$cd
    z = initial$z
    sigma.Q = initial$sigma.Q
    u = initial$u
  }						 
  bmat = list()
  etamat = list()
  for (IP in 1:nIP) {
    bmat[[IP]] = matrix(b.old[IP,], nrow = P, ncol = 1)
    etamat[[IP]] = matrix(eta.old[IP,], nrow = Q, ncol = 1)
  }
  deltamat = rep(delta, 1)
  sigma_taumat = rep(sigma_tau, 1)		
  mu = matrix(0, nrow = D, ncol = A)
  textlist.raw = unlist(textlist[edge.trim])
  accept.rates = rep(0, 4)
  hist.d = c()
  for (d in 1:D) {
  if (timestamps[d]+384*timeunit > timestamps[max.edge]) {
        hist.d[d] = max.edge
      } else {
        hist.d[d] = which_num(timestamps[d]+384*timeunit, timestamps)
      }
  }
  timeinterval = timefinder(timestamps, edge.trim, timeunit)
  X = list()
  for (d in edge.trim) {
    X[[d]] = Netstats_cpp(edge, timestamps, timeinterval[[d]], senders, cd, A, timeunit, netstat)
  }
  const.C = rep(NA, nIP)
  table.W = tvapply(1:K, function(k) {
      tabulateC(textlist.raw[which(unlist(z[edge.trim]) == k)], V)
        }, rep(0, V))
  totalN = sum(table.W)-1 
  table.C = tabulateC(cd[edge.trim], nIP)
  table.cd = vapply(1:nIP, function(IP) {
    if (sum(cd[edge.trim]==IP) > 0) {
      tabulateC(unlist(z[edge.trim][which(cd[edge.trim] == IP)]), K)
    } else {
      rep(0, K)
    }  
  }, rep(0, K))
  table.dk = vapply(1:D, function(d) {
      tabulateC(unlist(z[[d]]), K)
        }, rep(0, K))
  table.k = rowSums(table.dk[,edge.trim])   
  #start outer iteration
  for (o in 1:Outer) {
    # Data augmentation
    for (d in edge.trim) {
        lambda = MultiplyXB(X[[d]], b.old[cd[d],])
        for (i in node[-senders[d]]) {
            for (j in sample(node[-i], A-1)) {
                probij = u_Gibbs(u[[d]][i, ], lambda[i,], delta, j)
                u[[d]][i, j] = lmultinom(probij)-1
            }
        }
    }
    # cd update	
     for (d in edge.trim) {
     	table.C[cd[d]] = table.C[cd[d]]-1
     	table.cd[, cd[d]] = table.cd[, cd[d]] - tabulateC(z[[d]], K)
        for (IP in 1:nIP) {
         cd[d] = IP
         IPpart = log(table.C[IP] + 1 + zeta/nIP)
       	 Xnew =  Netstats_cpp(edge, timestamps, timeinterval[[d]], senders, cd, A, timeunit, netstat)
       	 Edgepart = Edgepartsum(Xnew, b.old[IP,], u[[d]], delta)
       	 munew = mu_vec(timemat[d,], A, eta.old[IP,])
         Timepart = Timepart(munew, sigma_tau, senders[d], timeinc[d])
         Topicpart = Topicpart(K, z[[d]], table.cd[,IP], table.k-tabulateC(z[[d]], K), alphas)
         const.C[IP] = Edgepart+Timepart+Topicpart+IPpart
       }
       cd[d] = lmultinom(const.C)
       table.C[cd[d]] = table.C[cd[d]] + 1
       table.cd[, cd[d]] = table.cd[, cd[d]] + tabulateC(z[[d]], K)
    }
    for (d in edge.trim) {
         X[[d]] = Netstats_cpp(edge, timestamps, timeinterval[[d]], senders, cd, A, timeunit, netstat)
	 }
    
    #Z update
     for (d in edge.trim) {   
       textlist.d = textlist[[d]]
 	   for (w in 1:length(z[[d]])) {
 	   		zw.old = z[[d]][w]
 	   		table.cd[zw.old, cd[d]] = table.cd[zw.old, cd[d]]-1
 	   		table.dk[zw.old, d] = table.dk[zw.old, d]-1
 	   		table.k[zw.old] = table.k[zw.old]-1
         if (length(textlist.d) > 0) {
             table.W[zw.old, textlist.d[w]] = table.W[zw.old, textlist.d[w]]-1
        	 topicword.d = TopicWord(K, table.dk[,d], table.W[,textlist.d[w]], table.cd[,cd[d]], table.k, totalN,
        	  						alphas, beta, V)
            zw.new = lmultinom(topicword.d)
            if (zw.new != zw.old) {
               z[[d]][w] = zw.new
           }
           table.W[z[[d]][w], textlist.d[w]] = table.W[z[[d]][w], textlist.d[w]]+1
           } else {
           topicword.d = TopicWord0(K, table.cd[,cd[d]], table.k, totalN, alphas, beta, V)
           zw.new = lmultinom(topicword.d)
           if (zw.new != zw.old) {
               z[[d]][w] = zw.new
           }
         }
        table.cd[z[[d]][w], cd[d]] = table.cd[z[[d]][w], cd[d]]+1
        table.dk[z[[d]][w], d] = table.dk[z[[d]][w], d]+1
 	    table.k[z[[d]][w]] = table.k[z[[d]][w]]+1
	 	}
	 }
 	       
  # M-H       
    prior.old1 = priorsum(prior.b[[2]], prior.b[[1]], b.old)+
    			 dnorm(delta, prior.delta[1], sqrt(prior.delta[2]), TRUE)
    post.old1 = Edgepartsum(X[[max.edge]], b.old[cd[max.edge],], u[[max.edge]], delta)
    b.new = matrix(NA, nIP, P)
    for (inner in 1:Inner[1]) {
      for (IP in 1:nIP) {
		 b.new[IP, ] = rmvnorm_arma(1, b.old[IP,], sigma.Q[1]*proposal.var1)
	  }
      delta.new = rnorm(1, delta, sqrt(sigma.Q[4]))
      prior.new1 = priorsum(prior.b[[2]], prior.b[[1]], b.new)+
    				 dnorm(delta.new, prior.delta[1], sqrt(prior.delta[2]), TRUE)
      post.new1 = Edgepartsum(X[[max.edge]], b.new[cd[max.edge],], u[[max.edge]], delta.new)
      loglike.diff = prior.new1+post.new1-prior.old1-post.old1
      if (log(runif(1, 0, 1)) < loglike.diff) {
        b.old = b.new
        delta = delta.new
        prior.old1 = prior.new1
        post.old1 = post.new1
      }
        for (IP in 1:nIP) {
          bmat[[IP]] = cbind(bmat[[IP]], b.old[IP,])
        }
        deltamat = c(deltamat, delta)
    }
     
      mu = mu_mat(timemat, A, eta.old, cd)
	  prior.old2 = priorsum(prior.eta[[2]], prior.eta[[1]], eta.old)
	  post.old2 = Timepartsum(mu[edge.trim,], sigma_tau, senders[edge.trim], timeinc[edge.trim])
      eta.new = matrix(NA, nIP, Q)
      for (inner in 1:Inner[2]) {
          for (IP in 1:nIP) {
			  eta.new[IP, ] = rmvnorm_arma(1, eta.old[IP,], sigma.Q[2]*proposal.var2)
		  }
      mu = mu_mat(timemat, A, eta.new, cd)
      prior.new2 = priorsum(prior.eta[[2]], prior.eta[[1]], eta.old)
      post.new2 = Timepartsum(mu[edge.trim,], sigma_tau, senders[edge.trim], timeinc[edge.trim])
      loglike.diff = prior.new2+post.new2-prior.old2-post.old2
      if (log(runif(1, 0, 1)) < loglike.diff) {
        eta.old = eta.new
        prior.old2 = prior.new2
        post.old2 = post.new2
      }
      for (IP in 1:nIP) {
          etamat[[IP]] = cbind(etamat[[IP]], eta.old[IP,])
        }
    }
      
    prior.old3 = dhalfcauchy(sigma_tau, prior.tau, TRUE)
    post.old3 = post.old2
    for (inner in 1:Inner[3]) {
      sigma_tau.new = rtruncnorm(1, 0, Inf, sigma_tau, sqrt(sigma.Q[3]))
      prior.new3 = dhalfcauchy(sigma_tau.new, prior.tau, TRUE)
      post.new3 =  Timepartsum(mu[edge.trim,], sigma_tau, senders[edge.trim], timeinc[edge.trim])
      loglike.diff = log(dtruncnorm(sigma_tau, 0, Inf, sigma_tau.new, sqrt(sigma.Q[3])))-
                   log(dtruncnorm(sigma_tau.new, 0, Inf, sigma_tau, sqrt(sigma.Q[3])))+
                   prior.new3+post.new3-prior.old3-post.old3
      if (log(runif(1, 0, 1)) < loglike.diff) {
        sigma_tau = sigma_tau.new
        prior.old3 = prior.new3
        post.old3 = post.new3
      }
        sigma_taumat = c(sigma_taumat, sigma_tau)
    }
  }
  chain.final = list(cd = cd, z = z, b = bmat, eta = etamat, delta = deltamat, sigma_tau = sigma_taumat,
                     u = u, sigma.Q =sigma.Q, edge.trim = edge.trim)
  return(chain.final)
}	


#' @title IPTM.inference.GiR
#' @description Iterate Markov Chain Monte Carlo (MCMC) algorithm for the interaction-partitioned topic model
#'
#' @param edge list of tie data with 3 elements (1: author, 2: recipient, 3: timestamp in unix.time format)
#' @param node vector of node id's (ID starting from 1)
#' @param textlist list of text containing the words in each document
#' @param vocab all vocabularies used over the corpus
#' @param nIP total number of interaction patterns
#' @param K total number of topics
#' @param sigma.Q proposal distribution variance parameter
#' @param alphas Dirichlet concentration prior for document-topic distribution (alpha0, alpha1, alpha)
#' @param beta Dirichlet concentration prior for topic-word distribution
#' @param zeta Dirichlet concentration prior for document-interaction-pattern distribution
#' @param prior.b prior mean and covariance of b in multivariate Normal distribution
#' @param prior.delta prior mean and variance of delta in Normal distribution
#' @param prior.eta prior mean and covariance of eta in multivariate Normal distribution
#' @param prior.tau prior shape and scale parameter of sigma_tau in inverse-Gamma distribution
#' @param Outer size of outer iterations 
#' @param Inner size of inner iteration for Metropolis-Hastings updates
#' @param netstat which type of network statistics to use ("dyadic", "triadic", "degree")
#' @param timestat additional statistics to be used for timestamps other than netstat ("sender", "receiver","timeofday", "dayofweek")
#' @param initial list of initial values user wants to assign including (delta, b, eta, cd, z, u, sigma_tau, proposal.var1, proposal.var2)
#' @param timeunit hour (= 3600) or day (=3600*24) and so on
#' @param tz timezone such as EST and PST
#'
#' @return MCMC output containing all parameter estimates
#'
#' @export

IPTM.inference.GiR = function(edge, node, textlist, vocab, nIP, K, sigma.Q, alphas, beta, zeta,
                          prior.b, prior.delta, prior.eta, prior.tau, Outer, Inner,
                          netstat, timestat, initial = NULL, timeunit = 3600, tz = "America/New_York") {
  
  # trim the edge so that we only model edges after 384 hours
  A = length(node)
  D = length(edge)
  timestamps = vapply(edge, function(d) { d[[3]] }, c(1))
  senders = vapply(edge, function(d) { d[[1]] }, c(1))
  edge.trim = which_num(384*timeunit, timestamps-timestamps[1]):D
  max.edge = max(edge.trim)
  timeinc = c(timestamps[1], timestamps[-1]-timestamps[-length(timestamps)])/timeunit
  timeinc[timeinc==0] = runif(sum(timeinc==0), 0, min(timeinc[timeinc!=0]))
  emptytext = which(sapply(textlist, function(d){length(d)})==0)
  # initialization 
  netstat = as.numeric(c("degree", "dyadic", "triadic") %in% netstat)
  timestat = as.numeric(c("dayofweek","timeofday") %in% timestat)
  timemat = matrix(0, nrow = D, ncol = sum(timestat))
  if (sum(timestat) > 0) {
    Sys.setenv(TZ = tz)
    time_ymd = as.POSIXct(timestamps, tz = getOption("tz"), origin = "1970-01-01")
    if (timestat[1] > 0) {
      days = vapply(time_ymd, function(d) {wday(d)}, c(1))
      days[days==1] = 8
      timemat[,1] = as.numeric(cut(days, c(1,6,8), c("weekdays","weekends")))-1
      it = 1
    }
    if (timestat[2] > 0) {
      hours = vapply(time_ymd, function(d) {hour(d)}, c(1))
      timemat[,it+1] = as.numeric(cut(hours, c(-1,12,24), c("AM", "PM")))-1
    }     
  }
  L = 3
  P = L*(2*netstat[1]+2*netstat[2]+4*netstat[3])
  Q = length(prior.eta[[1]])
  proposal.var1 = diag(P)
  proposal.var2 = diag(Q)
  V = length(vocab)
  phi = matrix(NA, K, V)
  for (k in 1:K) {
    phi[k,] = rdirichlet_cpp(1, rep(beta/V, V))
  }
  alpha0 = alphas[1]
  alpha1 = alphas[2]
  alpha = alphas[3]
  m = rdirichlet_cpp(1, rep(alpha0/K, K))
  mc = matrix(NA, nIP, K)
  for (IP in 1:nIP) {
    mc[IP,] = rdirichlet_cpp(1, alpha1*m)
  }
  psi = rdirichlet_cpp(1, rep(zeta/nIP, nIP))
  if (length(initial) == 0) {
    cd = multinom_vec(D, psi) 
    theta = tvapply(seq(along = edge), function(d) {rdirichlet_cpp(1, alpha*mc[cd[d],]) }, rep(0, K))
    delta = rnorm(1, prior.delta[1], sqrt(prior.delta[2]))
    sigma_tau = rhalfcauchy(1, prior.tau)
    b.old = rmvnorm_arma(nIP, prior.b[[1]], prior.b[[2]])
    eta.old = rmvnorm_arma(nIP, prior.eta[[1]], prior.eta[[2]])
    z = lapply(seq(along = edge), function(d) multinom_vec(max(1, length(textlist[[d]])), theta[d, ]))
    sigma.Q = sigma.Q
    u = list()
    for (d in edge.trim) {
      u[[d]] = matrix(rbinom(A^2, 1, 1/A), nrow = A, ncol = A)
      diag(u[[d]]) = 0
      u[[d]][senders[d],] = tabulateC(as.numeric(unlist(edge[[d]][2])), A)
    }
  } else {
    delta = initial$delta
    sigma_tau = initial$sigma_tau
    b.old = initial$b
    eta.old = initial$eta
    cd = initial$cd
    z = initial$z
    sigma.Q = initial$sigma.Q
    u = initial$u
  }						 
  bmat = list()
  etamat = list()
  for (IP in 1:nIP) {
    bmat[[IP]] = matrix(b.old[IP,], nrow = P, ncol = 1)
    etamat[[IP]] = matrix(eta.old[IP,], nrow = Q, ncol = 1)
  }
  deltamat = rep(delta, 1)
  sigma_taumat = rep(sigma_tau, 1)		
  mu = matrix(0, nrow = D, ncol = A)
  textlist.raw = unlist(textlist[edge.trim])
  accept.rates = rep(0, 4)
  hist.d = c()
  for (d in 1:D) {
  if (timestamps[d]+384*timeunit > timestamps[max.edge]) {
        hist.d[d] = max.edge
      } else {
        hist.d[d] = which_num(timestamps[d]+384*timeunit, timestamps)
      }
  }
  timeinterval = timefinder(timestamps, edge.trim, timeunit)
  X = list()
  for (d in edge.trim) {
    X[[d]] = Netstats_cpp(edge, timestamps, timeinterval[[d]], senders, cd, A, timeunit, netstat)
  }
  table.W = tvapply(1:K, function(k) {
      tabulateC(textlist.raw[which(unlist(z[edge.trim]) == k)], V)
  }, rep(0, V))
  zuniq = lapply(z, function(x) {sortuniq(x)}) 
  table.C = tabulateC(cd[edge.trim], nIP)           
  table.cd = vapply(1:nIP, function(IP) {
    if (sum(cd[edge.trim] == IP) > 0) {
      tabulateC(unlist(zuniq[edge.trim][which(cd[edge.trim] == IP)]), K)
    } else {
      rep(0, K)
    }  
  }, rep(0, K))
  total.cd = length(edge.trim)
  table.dk = vapply(1:D, function(d) {
     tabulateC(unlist(z[[d]]), K)
       }, rep(0, K))
  table.k = rowSums(table.cd > 0) 
  totalN = nIP    
  const.C = rep(NA, nIP)
  #start outer iteration
  for (o in 1:Outer) {
    # Data augmentation
    for (d in edge.trim) {
        lambda = MultiplyXB(X[[d]], b.old[cd[d],])
        for (i in node[-senders[d]]) {
            for (j in sample(node[-i], A-1)) {
                probij = u_Gibbs(u[[d]][i, ], lambda[i,], delta, j)
                u[[d]][i, j] = lmultinom(probij)-1
            }
        }
    }

     # cd update	
     for (d in edge.trim) {
     	table.C[cd[d]] = table.C[cd[d]]-1
     	table.cd[, cd[d]] = table.cd[, cd[d]] - tabulateC(zuniq[[d]], K)
     	table.k = rowSums(table.cd > 0) 
        for (IP in 1:nIP) {
         cd[d] = IP
         IPpart = log(table.C[IP] + 1 + zeta/nIP)
       	 Xnew =  Netstats_cpp(edge, timestamps, timeinterval[[d]], senders, cd, A, timeunit, netstat)
       	 Edgepart = Edgepartsum(Xnew, b.old[IP,], u[[d]], delta)
       	 munew = mu_vec(timemat[d,], A, eta.old[IP,])
         Timepart = Timepart(munew, sigma_tau, senders[d], timeinc[d])
         Topicpart = Topicpart_min(K, z[[d]], table.cd[,IP], total.cd-1, table.k, alphas, nIP)
         const.C[IP] = Edgepart+Timepart+Topicpart+IPpart
       }
       cd[d] = lmultinom(const.C)
       table.C[cd[d]] = table.C[cd[d]] + 1
       table.cd[, cd[d]] = table.cd[, cd[d]] + tabulateC(zuniq[[d]], K)
       table.k = rowSums(table.cd > 0) 
    }  
	 for (d in edge.trim) {
         X[[d]] = Netstats_cpp(edge, timestamps, timeinterval[[d]], senders, cd, A, timeunit, netstat)
	 }
    # # Z update	
    # for (d in edge.trim) {
	   	# textlist.d = textlist[[d]]
	   	# for (w in 1:length(z[[d]])) {
	   		# zw.old = z[[d]][w]
	   		# table.dk[zw.old, d] = table.dk[zw.old, d]-1
	   		# if (!identical(zuniq[[d]], sortuniq(z[[d]][-w]))) {
	   			# zuniq[[d]] = sortuniq(z[[d]][-w])
	   			# table.cd = vapply(1:nIP, function(IP) {
	   			# if (sum(cd[edge.trim] == IP) > 0) {
	   			    # tabulateC(unlist(zuniq[edge.trim][which(cd[edge.trim] == IP)]), K)
	   			  # } else {
	   			    # rep(0, K)
	   			  # }  
	   			# }, rep(0, K))      
	   			# table.k = rowSums(table.cd > 0) 
	   		 # }
	   	 # if (length(textlist.d) > 0) {
          # table.W[zw.old, textlist.d[w]] = table.W[zw.old, textlist.d[w]]-1
       	  # topicword.d = TopicWord_min(K, table.dk[,d], table.W[,textlist.d[w]], table.cd[,cd[d]], total.cd, table.k, totalN, 
       	  						 # alphas, beta, V)
          # zw.new = lmultinom(topicword.d)
          # if (zw.new != zw.old) {
              # z[[d]][w] = zw.new
          # }
          # table.W[z[[d]][w], textlist.d[w]] = table.W[z[[d]][w], textlist.d[w]]+1
        # } else {
          # topicword.d = TopicWord0_min(K, table.cd[,cd[d]], total.cd, table.k, totalN, alphas, beta, V)
          # zw.new = lmultinom(topicword.d)
          # if (zw.new != zw.old) {
              # z[[d]][w] = zw.new
          # }
        # }
		# if (!identical(zuniq[[d]], sortuniq(z[[d]]))) {
	   			# zuniq[[d]] = sortuniq(z[[d]])
	   			# table.cd = vapply(1:nIP, function(IP) {
	   			  # if (sum(cd[edge.trim] == IP) > 0) {
	   			    # tabulateC(unlist(zuniq[edge.trim][which(cd[edge.trim] == IP)]), K)
	   			  # } else {
	   			    # rep(0, K)
	   			  # }  
	   			# }, rep(0, K))   
	   			# table.k = rowSums(table.cd > 0) 
	   		# }
       # table.dk[z[[d]][w], d] = table.dk[z[[d]][w], d]+1
      # }
    # }
    
    prior.old1 = priorsum(prior.b[[2]], prior.b[[1]], b.old)+
    			 dnorm(delta, prior.delta[1], sqrt(prior.delta[2]), TRUE)
    post.old1 = Edgepartsum(X[[max.edge]], b.old[cd[max.edge],], u[[max.edge]], delta)
    b.new = matrix(NA, nIP, P)
    for (inner in 1:Inner[1]) {
      for (IP in 1:nIP) {
			  b.new[IP, ] = rmvnorm_arma(1, b.old[IP,], sigma.Q[1]*proposal.var1)
	  }
      delta.new = rnorm(1, delta, sqrt(sigma.Q[4]))
      prior.new1 = priorsum(prior.b[[2]], prior.b[[1]], b.new)+
    				 dnorm(delta.new, prior.delta[1], sqrt(prior.delta[2]), TRUE)
      post.new1 = Edgepartsum(X[[max.edge]], b.new[cd[max.edge],], u[[max.edge]], delta)
      loglike.diff = prior.new1+post.new1-prior.old1-post.old1
      if (log(runif(1, 0, 1)) < loglike.diff) {
        b.old = b.new
        delta = delta.new
        prior.old1 = prior.new1
        post.old1 = post.new1
      }
        for (IP in 1:nIP) {
          bmat[[IP]] = cbind(bmat[[IP]], b.old[IP,])
        }
        deltamat = c(deltamat, delta)
    }
     
      mu = mu_mat(timemat, A, eta.old, cd)
	  prior.old2 = priorsum(prior.eta[[2]], prior.eta[[1]], eta.old)
	  post.old2 = Timepartsum(mu, sigma_tau, senders, timeinc)
      eta.new = matrix(NA, nIP, Q)
      for (inner in 1:Inner[2]) {
          for (IP in 1:nIP) {
			  eta.new[IP, ] = rmvnorm_arma(1, eta.old[IP,], sigma.Q[2]*proposal.var2)
		  }
      mu = mu_mat(timemat, A, eta.new, cd)
      prior.new2 = priorsum(prior.eta[[2]], prior.eta[[1]], eta.old)
      post.new2 = Timepartsum(mu, sigma_tau, senders, timeinc)
      loglike.diff = prior.new2+post.new2-prior.old2-post.old2
      if (log(runif(1, 0, 1)) < loglike.diff) {
        eta.old = eta.new
        prior.old2 = prior.new2
        post.old2 = post.new2
      }
      for (IP in 1:nIP) {
          etamat[[IP]] = cbind(etamat[[IP]], eta.old[IP,])
        }
    }
      
    prior.old3 = dhalfcauchy(sigma_tau, prior.tau, TRUE)
    post.old3 = post.old2
    for (inner in 1:Inner[3]) {
      sigma_tau.new = rtruncnorm(1, 0, Inf, sigma_tau, sqrt(sigma.Q[3]))
      prior.new3 = dhalfcauchy(sigma_tau.new, prior.tau, TRUE)
      post.new3 = Timepartsum(mu, sigma_tau.new, senders, timeinc)
      loglike.diff = log(dtruncnorm(sigma_tau, 0, Inf, sigma_tau.new, sqrt(sigma.Q[3])))-
                   log(dtruncnorm(sigma_tau.new, 0, Inf, sigma_tau, sqrt(sigma.Q[3])))+
                   prior.new3+post.new3-prior.old3-post.old3
      if (log(runif(1, 0, 1)) < loglike.diff) {
        sigma_tau = sigma_tau.new
        prior.old3 = prior.new3
        post.old3 = post.new3
      }
        sigma_taumat = c(sigma_taumat, sigma_tau)
    }
  }
  chain.final = list(cd = cd, z = z, b = bmat, eta = etamat, delta = deltamat, sigma_tau = sigma_taumat,
                     u = u, sigma.Q =sigma.Q, edge.trim = edge.trim)
  return(chain.final)
}	

#' @title GenerateDocs
#' @description Generate a collection of documents according to the generative process of IPTM
#'
#' @param nDocs number of documents to be generated
#' @param node vector of node id's (ID starting from 1)
#' @param vocab all vocabularies used over the corpus
#' @param nIP total number of interaction patterns
#' @param K total number of topics
#' @param n.d number of words in a document (fixed constant for now)
#' @param alphas Dirichlet concentration prior for document-topic distribution
#' @param beta Dirichlet concentration prior for topic-word distribution
#' @param zeta Dirichlet concentration prior for document-interaction pattern distribution
#' @param b coefficients for recipients
#' @param eta coefficients for timestamps
#' @param delta tuning parameter for the number of recipients
#' @param sigma_tau variance parameter for the timestamps
#' @param cd_assignments document-interaction pattern assignments
#' @param support support of latent recipients
#' @param netstat which type of network statistics to use ("dyadic", "triadic", "degree")
#' @param timestat additional statistics to be used for timestamps other than netstat ("sender", "receiver", "timeofday", "dayofweek")
#' @param base.data edges before 384 hours that is used to calculate initial history of interactions
#' @param topic_token_assignments matrix of topic-token assignments
#' @param backward Logigal indicating whether we are generating backward samples (if FALSE -> forward)
#' @param base Logical indicating whether or not we are generating base edges (< 384)
#' @param timeunit hour (= 3600) or day (=3600*24) and so on
#' @param tz timezone such as EST and PST
#'
#' @return generated data including (author, recipients, timestamp, words)
#'
#' @export
GenerateDocs = function(nDocs, node, vocab, nIP, K, n.d, alphas, beta, zeta,
                        b, eta, delta, sigma_tau, cd_assignments = NULL, support, netstat, timestat,
                        base.data = NULL, topic_token_assignments = NULL,
                        backward = FALSE, base = FALSE, timeunit = 3600, tz = "America/New_York") { 
  A = length(node)
  V = length(vocab)

  netstat = as.numeric(c("degree", "dyadic", "triadic" ) %in% netstat)
  timestat = as.numeric(c("dayofweek","timeofday") %in% timestat)
  L = 3
  P = L*(2*netstat[1]+2*netstat[2]+4*netstat[3])
  t.d = ifelse(base, 0, base.data$edge[[1]][[3]]+384*timeunit)
  edge = base.data$edge
  text = base.data$text
  base.length = length(edge)
  timemat = matrix(0, nrow = nDocs+base.length, ncol = sum(timestat))
  timestamps = rep(NA, nDocs+base.length)
  senders = rep(NA, nDocs+base.length)
  cd = rep(NA, nDocs+base.length)
  timeinterval = list()
  if (base.length > 0) {
  	timestamps[1:base.length] = vapply(edge, function(d) { d[[3]] }, c(1))
  	senders[1:base.length] = vapply(edge, function(d) { d[[1]] }, c(1))
  	timeinterval[1:base.length] = lapply(1:base.length, function(d) timefinder_vec(timestamps[1:d], d, timeunit))
  	cd[1:base.length] = base.data$cd
  }

  if (!base) {
    if (sum(timestat) > 0) {
      Sys.setenv(TZ = tz)
      time_ymd = as.POSIXct(vapply(base.data$edge, function(d) {d[[3]]}, c(1)), tz = getOption("tz"), origin = "1970-01-01")
      if (timestat[1] > 0) {
        days = vapply(time_ymd, function(d) {wday(d)}, c(1))
        days[days==1] = 8
        timemat[1:base.length,1] = as.numeric(cut(days, c(1,6,8), c("weekdays","weekends")))-1
        it = 1
      }
      if (timestat[2] > 0) {
        hours = vapply(time_ymd, function(d) {hour(d)}, c(1))
        timemat[1:base.length,it+1] = as.numeric(cut(hours, c(-1,12,24), c("AM", "PM")))-1
      }
    }
  }  
  if (!backward) {
  phi = matrix(NA, K, V)
  for (k in 1:K) {
    phi[k,] = rdirichlet_cpp(1, rep(beta/V, V))
  } 
  psi = rdirichlet_cpp(1, rep(zeta/nIP, nIP))
  alpha0 = alphas[1]
  alpha1 = alphas[2]
  alpha = alphas[3]
  m = rdirichlet_cpp(1, rep(alpha0/K, K))
  mc = matrix(NA, nIP, K)
  for (IP in 1:nIP) {
    mc[IP,] = rdirichlet_cpp(1, alpha1*m)
  }
  }
  VKmat = matrix(0, V, K)
  Cmat = rep(0, nIP)
  u = list()
  for (d in 1:nDocs) {
    u[[base.length+d]] = matrix(0, A, A)
    text[[base.length+d]] = rep(NA, n.d)
    if (!backward) {
      cd[base.length+d] =  multinom_vec(1, psi)
      theta.d = rdirichlet_cpp(1, alpha*mc[cd[base.length+d], ])	
      topic.d = multinom_vec(max(1, n.d), theta.d)
      for (n in 1:n.d){
        text[[base.length+d]][n] = multinom_vec(1, phi[topic.d[n],])
      }
      names(text[[base.length+d]]) = topic.d
    } else {
      cd[base.length+d] = cd_assignments[base.length+d]
      phi.k = rep(NA, K)
      topic.d = topic_token_assignments[[base.length+d]]
      for (n in 1:n.d){
        for (w in 1:V) {
          phi.k[w] = (VKmat[w, topic.d[n]]+beta/V)/(sum(VKmat[,topic.d[n]])+beta)
        } 
        text[[base.length+d]][n] = multinom_vec(1, phi.k)
        VKmat[text[[base.length+d]][n], topic.d[n]] = VKmat[text[[base.length+d]][n], topic.d[n]]+1
      }
      names(text[[base.length+d]]) = topic.d
    }
    if (t.d >= 384*timeunit) {
       timeinterval = timefinder_vec(timestamps[1:(base.length+d-1)], base.length+d-1, timeunit)
       X = Netstats_cpp(edge, timestamps[1:(base.length+d-1)], timeinterval, senders[1:(base.length+d-1)], cd[1:(base.length+d)], A, timeunit, netstat)
       lambda = MultiplyXB(X, b[cd[base.length+d],])   
    } else {
       lambda = matrix(0, A, A)
    }
    timevec = rep(NA, A)
    for (i in node) {
      u[[base.length+d]][i,-i] = r.gibbs.measure(1, lambda[i,-i], delta, support)
      mu = mu_vec(timemat[base.length+d, ], A, eta[cd[base.length+d],])
      time = rlnorm(1, mu, sigma_tau)
      timevec[i] = time*timeunit
    }
    i.d = which(timevec == min(timevec))
    j.d = which(u[[base.length+d]][i.d,] == 1)
    t.d = t.d+timevec[i.d]
    senders[base.length+d] = i.d
    timestamps[base.length+d] = t.d
    edge[[base.length+d]] = list(author = i.d, recipients = j.d, timestamp = t.d)
    if (t.d <= exp(38.7) & sum(timestat) > 0) {
        Sys.setenv(TZ = tz)
     	it = 0
      time_ymd = as.POSIXct(edge[[base.length+d]][[3]], tz = getOption("tz"), origin = "1970-01-01")
      if (timestat[1] > 0) {
      	it = it + 1
        days = vapply(time_ymd, function(d) {wday(d)}, c(1))
        days[days==1] = 8
        timemat[base.length+d,it] = as.numeric(cut(days, c(1,6,8), c("weekdays","weekends")))-1
      }
      if (timestat[2] > 0) {
      	it = it + 1
        hours = vapply(time_ymd, function(d) {hour(d)}, c(1))
        timemat[base.length+d,it] = as.numeric(cut(hours, c(-1,12,24), c("AM", "PM")))-1
      }
    }
  }
  if (base == TRUE & t.d > 384*timeunit) {
    cutoff = which_num(384*timeunit, vapply(1:length(edge), function(d) {edge[[d]][[3]]}, c(1)))-1
    edge = edge[1:cutoff]
    text = text[1:cutoff]
    cd = cd[1:cutoff]
  }
  return(list(edge = edge, text = text, base = base.length, u = u, z = lapply(text, function(d){as.numeric(names(d))}), 
  		b = b, eta = eta, delta = delta, sigma_tau = sigma_tau, cd = cd))							
} 


#' @title GiR_stats
#' @description Calculate several statistics from samples generated from forward or backward sampling
#'
#' @param GiR_sample one sample from generative process
#' @param V number of unique words
#' @param K number of topics
#' @param timeunit hour (= 3600) or day (=3600*24) and so on
#'
#' @return A vector of statistics calculated from one GiR sample
#'
#' @export

GiR_stats = function(GiR_sample, V, K, timeunit = 3600) {
  edge = GiR_sample$edge
  text = GiR_sample$text
  cd = GiR_sample$cd
  if (GiR_sample$base > 0)  {
  	edge = edge[-(1:GiR_sample$base)]
  	text = text[-(1:GiR_sample$base)]
  	cd = cd[-(1:GiR_sample$base)]
  }
  GiR_stats = c()
  D = length(edge)
  P = ncol(GiR_sample$b)
  Q = ncol(GiR_sample$eta)
  nIP = nrow(GiR_sample$b)
  n.d = length(text[[1]])
  
  GiR_stats[1:P] = colMeans(GiR_sample$b)
  GiR_stats[(P+1):(P+Q)] = colMeans(GiR_sample$eta)
  GiR_stats[P+Q+1] = GiR_sample$delta
  GiR_stats[P+Q+2] = GiR_sample$sigma_tau
  GiR_stats[P+Q+3] = mean(vapply(1:D, function(d) {length(edge[[d]][[2]])}, c(1)))
  GiR_stats[P+Q+4] = mean(vapply(2:D, function(d) {edge[[d]][[3]]-edge[[d-1]][[3]]}, c(1))/timeunit) 			
  GiR_stats[P+Q+5] = mean(GiR_sample$cd)
  Tokens_in_Topic = tabulate(vapply(1:D, function(d){as.numeric(names(text[[d]]))}, rep(0, n.d)), K)
  GiR_stats[(P+Q+6):(P+Q+5+nIP)] = vapply(1:nIP, function(IP) {sum(cd==IP)}, c(1))
  GiR_stats[(P+Q+6+nIP):(P+Q+5+nIP+K)] = Tokens_in_Topic
  GiR_stats[(P+Q+6+nIP+K):(P+Q+5+nIP+K+V)] = tabulate(unlist(text), V)
  return(GiR_stats)
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
          
     
#' @title GiR
#' @description Getting it Right test for the IPTM
#'
#' @param Nsamp number of GiR samples to be generated 
#' @param nDocs number of documents to be generated per each sample
#' @param node vector of node id's (ID starting from 1)
#' @param vocab all vocabularies used over the corpus
#' @param nIP total number of interaction patterns
#' @param K total number of topics
#' @param n.d number of words in a document (fixed constant for now)
#' @param alphas Dirichlet concentration prior for document-topic distribution
#' @param beta Dirichlet concentration prior for topic-word distribution
#' @param zeta Dirichlet concentration prior for document-interaction pattern distribution
#' @param prior.b prior mean and covariance of b in multivariate Normal distribution
#' @param prior.delta prior mean and variance of delta in Normal distribution
#' @param prior.eta prior mean and covariance of eta in multivariate Normal distribution
#' @param prior.tau prior shape and scale parameter of sigma_tau in inverse-Gamma distribution
#' @param sigma.Q proposal distribution variance parameter
#' @param Outer size of outer iterations 
#' @param Inner size of inner iteration for Metropolis-Hastings updates
#' @param netstat which type of network statistics to use ("dyadic", "triadic", "degree")
#' @param timestat additional statistics to be used for timestamps other than netstat ("sender", "receiver", "timeofday", "dayofweek")
#' @param base.data artificial collection of documents to be used as initial state of history
#' @param generate_PP_plots Logical indicating whether to draw PP plots
#'
#' @return Forward and Backward samples and corresponding test results
#'
#' @export
GiR = function(Nsamp, nDocs, node, vocab, nIP, K, n.d, alphas, beta, zeta,
              prior.b, prior.delta, prior.eta, prior.tau, sigma.Q, Outer, Inner,
              netstat = c("dyadic"), timestat = c("timeofday", "dayofweek"),
              base.data, generate_PP_plots = TRUE) {
  
  A = length(node)
  netstat2 = as.numeric(c("degree", "dyadic", "triadic") %in% netstat)
  timestat2 = as.numeric(c("dayofweek","timeofday") %in% timestat)
  L = 3
  P = L*(2*netstat2[1]+2*netstat2[2]+4*netstat2[3])
  Q = length(prior.eta[[1]])
  V = length(vocab)
  support = gibbs.measure.support(A-1)
  
  #Forward sampling
  Forward_stats = matrix(NA, nrow = Nsamp, ncol = P+Q+5+nIP+K+V)
  colnames(Forward_stats) = c(paste0("b_",1:P), paste0("eta_",1:Q), "delta", "sigma_tau", 
                              "Mean_recipients", "Mean_timediff", "Mean_TopicIP", paste0("Documents_in_IP_", 1:nIP), 
                              paste0("Tokens_in_Topic", 1:K), paste0("Tokens_in_Word",1:V))
  for (i in 1:Nsamp) { 
    if (i %% 5000 == 0) {cat("Forward sampling", i, "\n")}
    b = rmvnorm_arma(nIP, prior.b[[1]], prior.b[[2]])
    eta = rmvnorm_arma(nIP, prior.eta[[1]], prior.eta[[2]])
    delta = rnorm(1, prior.delta[1], sqrt(prior.delta[2]))
    sigma_tau = rhalfcauchy(1, prior.tau)
    psi = rdirichlet_cpp(1, rep(zeta/nIP, nIP))
    cd_assignments = vapply(1:nDocs+length(base.data$edge), function(i) which(rmultinom(1,1,psi)==1), c(1))

    Forward_sample = GenerateDocs(nDocs, node, vocab, nIP, K, n.d, alphas, beta, zeta, b, eta, delta, sigma_tau, 
                    cd_assignments = cd_assignments, support, netstat, timestat, base.data = base.data,
                    topic_token_assignments = NULL, backward = FALSE, base = FALSE)
    Forward_stats[i, ] = GiR_stats(Forward_sample, V)
  }
  #Backward sampling
  Backward_stats = matrix(NA, nrow = Nsamp, ncol = ncol(Forward_stats))
  Backward_sample = GenerateDocs(nDocs, node, vocab, nIP, K, n.d, alphas, beta, zeta, b, eta, delta, sigma_tau, 
                    cd_assignments = cd_assignments, support, netstat, timestat, base.data = base.data, 
                    topic_token_assignments = NULL, backward = FALSE, base = FALSE)
  for (i in 1:Nsamp) { 
    if (i %% 1 == 0) {cat("Backward Sampling", i, "\n")}
    Inference_samp = IPTM.inference.GiR(Backward_sample$edge, node, Backward_sample$text, vocab, nIP, K, sigma.Q, 
                     alphas, beta, zeta, prior.b, prior.delta,prior.eta, prior.tau, Outer, Inner,
                     netstat, timestat, initial = NULL)
    b = tvapply(1:nIP, function(IP) {Inference_samp$b[[IP]][,ncol(Inference_samp$b[[IP]])]}, rep(0, P))
    eta = tvapply(1:nIP, function(IP) {Inference_samp$eta[[IP]][,ncol(Inference_samp$eta[[IP]])]}, rep(0, Q))
    delta = Inference_samp$delta[length(Inference_samp$delta)]
    sigma_tau = Inference_samp$sigma_tau[length(Inference_samp$sigma_tau)]
    cd_assignments = Inference_samp$cd
    z = Inference_samp$z
    for (d in 1:length(z)) {
      names(z[[d]]) = Forward_sample$text[[d]]
    }
    Backward_sample = GenerateDocs(nDocs, node, vocab, nIP, K, n.d, alphas, beta, zeta, b, eta, delta, sigma_tau, 
                      cd_assignments = cd_assignments, support, netstat, timestat, base.data = base.data,
                      topic_token_assignments = z, backward = TRUE, base = FALSE)
    Backward_stats[i, ] = GiR_stats(Backward_sample, V)
  }
  				
  if (generate_PP_plots) {
    par(mfrow=c(5,6), oma = c(3,3,3,3), mar = c(2,1,1,1))
    GiR_PP_Plots(Forward_stats, Backward_stats)
  }			
  return(list(Forward = Forward_stats, Backward = Backward_stats))
}                         	          
      

     
#' @title Schein
#' @description Schein (preliminary GiR) test for the IPTM
#'
#' @param Nsamp number of GiR samples to be generated 
#' @param nDocs number of documents to be generated per each sample
#' @param node vector of node id's (ID starting from 1)
#' @param vocab all vocabularies used over the corpus
#' @param nIP total number of interaction patterns
#' @param K total number of topics
#' @param n.d number of words in a document (fixed constant for now)
#' @param alphas Dirichlet concentration prior for document-topic distribution
#' @param beta Dirichlet concentration prior for topic-word distribution
#' @param zeta Dirichlet concentration prior for document-interaction patterndistribution
#' @param prior.b prior mean and covariance of b in multivariate Normal distribution
#' @param prior.delta prior mean and variance of delta in Normal distribution
#' @param prior.eta prior mean and covariance of eta in multivariate Normal distribution
#' @param prior.tau prior shape and scale parameter of sigma_tau in inverse-Gamma distribution
#' @param sigma.Q proposal distribution variance parameter
#' @param Outer size of outer iterations 
#' @param Inner size of inner iteration for Metropolis-Hastings updates
#' @param netstat which type of network statistics to use ("dyadic", "triadic", "degree")
#' @param timestat additional statistics to be used for timestamps other than netstat ("sender", "receiver", "timeofday", "dayofweek")
#' @param base.data artificial collection of documents to be used as initial state of history
#' @param generate_PP_plots Logical indicating whether to draw PP plots
#'
#' @return Forward and Backward samples and corresponding test results
#'
#' @export
Schein = function(Nsamp, nDocs, node, vocab, nIP, K, n.d, alphas, beta, zeta,
              prior.b, prior.delta, prior.eta, prior.tau, sigma.Q, Outer, Inner, 
              netstat = c("dyadic"), timestat = c("timeofday", "dayofweek"),
              base.data, generate_PP_plots = TRUE) {
  
  A = length(node)
  netstat2 = as.numeric(c("degree", "dyadic", "triadic") %in% netstat)
  timestat2 = as.numeric(c("dayofweek","timeofday") %in% timestat)
  L = 3
  P = L*(2*netstat2[1]+2*netstat2[2]+4*netstat2[3])
  Q = length(prior.eta[[1]])
  V = length(vocab)
  support = gibbs.measure.support(A-1)
  
  #Forward sampling
  Forward_stats = matrix(NA, nrow = Nsamp, ncol = P+Q+5+nIP+K+V)
  colnames(Forward_stats) = c(paste0("b_",1:P), paste0("eta_",1:Q), "delta", "sigma_tau",
                            "Mean_recipients", "Mean_timediff", "Mean_TopicIP", paste0("Documents_in_IP_", 1:nIP), 
                            paste0("Tokens_in_Topic", 1:K), paste0("Tokens_in_Word",1:V))
  #Backward sampling
  Backward_stats = matrix(NA, nrow = Nsamp, ncol = ncol(Forward_stats))
					 
  for (i in 1:Nsamp) { 
  	if (i %% 100 == 0) {cat("Sampling", i, "\n")}
    b = rmvnorm_arma(nIP, prior.b[[1]], prior.b[[2]])
    eta = rmvnorm_arma(nIP, prior.eta[[1]], prior.eta[[2]])
    delta = rnorm(1, prior.delta[1], sqrt(prior.delta[2]))
    sigma_tau = rhalfcauchy(1, prior.tau)
	while (sigma_tau > 10) {sigma_tau = rhalfcauchy(1, prior.tau)}
	psi = rdirichlet_cpp(1, rep(zeta/nIP, nIP))
    cd_assignments = vapply(1:nDocs+length(base.data$edge), function(i) which(rmultinom(1,1,psi)==1), c(1))
    Forward_sample = GenerateDocs(nDocs, node, vocab, nIP, K, n.d, alphas, beta, zeta, b, eta, delta, sigma_tau, 
                     cd_assignments = cd_assignments, support, netstat, timestat, base.data = base.data, 
                     topic_token_assignments = NULL, backward = FALSE, base = FALSE)
    Forward_stats[i, ] = GiR_stats(Forward_sample, V, K)
  	initial = list(delta = delta, sigma_tau = sigma_tau, b = b, eta = eta,
  	cd = Forward_sample$cd, z = Forward_sample$z, sigma.Q = sigma.Q, u = Forward_sample$u)
    Inference_samp = IPTM.inference.GiR2(Forward_sample$edge, node, Forward_sample$text, vocab, nIP, K, sigma.Q,
                     alphas, beta, zeta, prior.b, prior.delta, prior.eta, prior.tau, Outer, Inner, netstat, timestat, initial = initial)
    b = tvapply(1:nIP, function(IP) {Inference_samp$b[[IP]][,ncol(Inference_samp$b[[IP]])]}, rep(0, P))
    eta = tvapply(1:nIP, function(IP) {Inference_samp$eta[[IP]][,ncol(Inference_samp$eta[[IP]])]}, rep(0, Q))
    delta = Inference_samp$delta[length(Inference_samp$delta)]
    sigma_tau = Inference_samp$sigma_tau[length(Inference_samp$sigma_tau)]
    cd_assignments = Inference_samp$cd
    z = Inference_samp$z
    for (d in 1:length(z)) {
      names(z[[d]]) = Forward_sample$text[[d]]
    }
    Backward_sample = GenerateDocs(nDocs, node, vocab, nIP, K, n.d, alphas, beta, zeta, b, eta, delta, sigma_tau,
                      cd_assignments = cd_assignments, support, netstat, timestat, base.data = base.data, 
                      topic_token_assignments = z, backward = TRUE, base = FALSE)
    Backward_stats[i, ] = GiR_stats(Backward_sample, V, K)
 }
  if (generate_PP_plots) {
    par(mfrow=c(5,6), oma = c(3,3,3,3), mar = c(2,1,1,1))
    GiR_PP_Plots(Forward_stats, Backward_stats)
  } 
  return(list(Forward = Forward_stats, Backward = Backward_stats))
}                         	          
      
   
