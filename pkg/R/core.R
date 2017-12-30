#' @useDynLib IPTM
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
		exp(sum((delta+lambda.i)*support[s,]))
		}, c(1))		
	logitNumerator[logitNumerator==Inf] = exp(700)
	logitNumerator[logitNumerator==0] = exp(-745)
	samp = multinom_vec(1, logitNumerator)	
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
adaptive.MH = function(sigma.Q, accept.rates, target = 0.25, update.size, tol = 0.15) {
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
		nK.word.list[d, ] = tabulateC(z[[d]], K)
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
#' @param alpha Dirichlet concentration prior for document-topic distribution
#' @param mvec Dirichlet base prior for document-topic distribution
#' @param beta Dirichlet concentration prior for topic-word distribution
#' @param prior.b prior mean and covariance of b in multivariate Normal distribution
#' @param prior.delta prior mean and variance of delta in Normal distribution
#' @param prior.eta prior mean and covariance of eta in multivariate Normal distribution
#' @param prior.tau prior shape and scale parameter of sigma_tau in inverse-Gamma distribution
#' @param Outer size of outer iterations 
#' @param Inner size of inner iteration for Metropolis-Hastings updates
#' @param burn iterations to be discarded at the beginning of Metropolis-Hastings chains
#' @param thin the thinning interval of Metropolis-Hastings chains
#' @param netstat which type of network statistics to use ("dyadic", "triadic", "degree")
#' @param timestat additional statistics to be used for timestamps other than netstat ("sender", "receiver","timeofday", "dayofweek")
#' @param optimize logical to optimize alpha (Dirichlet concentration prior for document-topic distribution)
#' @param initial list of initial values user wants to assign including (alpha, mvec, delta, b, eta, l, z, u, sigma_tau, proposal.var1, proposal.var2)
#' @param timeunit hour (= 3600) or day (=3600*24) and so on
#' @param tz timezone such as EST and PST
#'
#' @return MCMC output containing all parameter estimates
#'
#' @export

IPTM.inference = function(edge, node, textlist, vocab, nIP, K, sigma.Q, alpha, mvec, beta, 
                          prior.b, prior.delta, prior.eta, prior.tau, Outer, Inner, burn, thin, 
                          netstat, timestat, optimize = FALSE, initial = NULL, timeunit = 3600, tz = "America/New_York") {
  
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
  V = length(vocab)
  phi = matrix(NA, K, V)
  for (k in 1:K) {
    phi[k,] = rdirichlet_cpp(1, rep(beta/V, V))
  }
  if (length(initial) == 0) {
    theta = rdirichlet_cpp(D, alpha*mvec)
    delta = rnorm(1, prior.delta[1], sqrt(prior.delta[2]))
    sigma_tau = rhalfcauchy(1, prior.tau)
    b.old = rmvnorm_arma(nIP, prior.b[[1]], prior.b[[2]])
    eta.old = rmvnorm_arma(nIP, prior.eta[[1]], prior.eta[[2]])
    l = sample(1:nIP, K, replace = TRUE) 
    z = lapply(seq(along = edge), function(d) multinom_vec(max(1, length(textlist[[d]])), theta[d, ]))
    p.d = pdmat(z, l, nIP) 
    proposal.var1 = diag(P)
    proposal.var2 = diag(Q)
    sigma.Q = sigma.Q
    u = list()
    for (d in edge.trim) {
      u[[d]] = matrix(rbinom(A^2, 1, 1/A), nrow =A, ncol = A)
      diag(u[[d]]) = 0
    } 
  } else {
    theta = rdirichlet_cpp(D, initial$alpha*initial$mvec)
    delta = initial$delta
    sigma_tau = initial$sigma_tau
    b.old = initial$b
    eta.old = initial$eta
	  l = initial$l
 	  z = initial$z
  	p.d = pdmat(z, l, nIP) 
    proposal.var1 = initial$proposal.var1
    proposal.var2 = initial$proposal.var2
    sigma.Q = initial$sigma.Q
    u = initial$u
  }						 
  bmat = list()
  etamat = list()
  for (IP in 1:nIP) {
    bmat[[IP]] = matrix(b.old[IP,], nrow = P, ncol = (Inner[1]-burn[1])/thin[1])
    etamat[[IP]] = matrix(eta.old[IP,], nrow = Q, ncol = (Inner[2]-burn[2])/thin[2])
  }
  deltamat = rep(delta, (Inner[1]-burn[1])/thin[1])
  sigma_taumat = rep(sigma_tau, (Inner[3]-burn[3])/thin[3])		
  xi = xi_all(timemat, eta.old[,node], eta.old[,-node], edge.trim)
  mu = matrix(0, nrow = D, ncol = A)
  textlist.raw = unlist(textlist)
  alphavec = c()
  mvecmat = matrix(NA, nrow = 0, ncol = K)
  accept.rates = rep(0, 4)
  hist.d = c()
  for (d in edge.trim) {
  if (timestamps[d]+384*timeunit > timestamps[max.edge]) {
        hist.d[d] = max.edge
      } else {
        hist.d[d] = which_num(timestamps[d]+384*timeunit, timestamps)-1
      }
  }
  X = list()
  for (d in edge.trim) {
        history.t = History(edge, p.d, node, timestamps[d-1]+exp(-745), timeunit)
        X[[d]] = Netstats_cpp(history.t, node, netstat)
  }    

  #start outer iteration
  for (o in 1:Outer) {
    print(o)
    if (o == Outer) {
      Inner = Inner * 10
      burn = burn * 10
      for (IP in 1:nIP) {
        bmat[[IP]] = matrix(b.old[IP,], nrow = P, ncol = (Inner[1]-burn[1])/thin[1])
        etamat[[IP]] = matrix(eta.old[IP,], nrow = Q, ncol = (Inner[2]-burn[2])/thin[2])
      }
      deltamat = rep(delta, (Inner[1]-burn[1])/thin[1])
      sigma_taumat = rep(sigma_tau, (Inner[3]-burn[3])/thin[3])						 
    }
    if (optimize & o > 1) {
      #update the hyperparameter alpha and mvec
      vec = AlphamvecOpt(K, z, alpha, mvec, 5)
      alpha = sum(vec)
      mvec = vec/alpha
      alphavec = c(alphavec, alpha)
      mvecmat = rbind(mvecmat, mvec)
    }
      
    # Data augmentation
    for (d in edge.trim) {
        vu = MultiplyXB(X[[d]], b.old)
        lambda = lambda_cpp(p.d[d,], vu)
        for (i in node[-senders[d]]) {
            for (j in sample(node[-i], A-1)) {
                probij = u_Gibbs(u[[d]][i, ], lambda[i,], delta, j)
                u[[d]][i, j] = multinom_vec(1, expconst(probij))-1
            }
        }
        u[[d]][senders[d],] = tabulateC(as.numeric(unlist(edge[[d]][2])), A)
    }
    # Z update	
    table.W = lapply(1:K, function(k) tabulateC(textlist.raw[which(unlist(z[-emptytext]) == k)], V))
    for (d in 1:(edge.trim[1]-1)) {
	   	textlist.d = textlist[[d]] 
	   	for (w in 1:length(z[[d]])) {
	   		zw.old = z[[d]][w]
        if (length(textlist.d) > 0) {
          table.W[[zw.old]][textlist.d[w]] = table.W[[zw.old]][textlist.d[w]]-1
       	  topicword.d = TopicWord(K, z[[d]][-w], textlist.d, table.W, alpha, mvec, beta, V)
        } else {
          topicword.d = matrix(0, nrow = length(z[[d]]), ncol = K)
        }
        const.Z = topicword.d[w, ]
        zw.new = multinom_vec(1, expconst(const.Z))
        if (zw.new != zw.old) {
          z[[d]][w] = zw.new
      	  p.d[d, ] = pdmat(list(z[[d]]), l, nIP)
        }
        if (length(textlist.d) > 0) {
           	table.W[[z[[d]][w]]][textlist.d[w]] = table.W[[z[[d]][w]]][textlist.d[w]]+1
        } 
      }
    }
    for (d in edge.trim) {
      textlist.d = textlist[[d]]
      edgetime.d = rep(NA, K)
      for (w in 1:length(z[[d]])) {
       	zw.old = z[[d]][w]
       	if (length(textlist.d) > 0) {
       	  table.W[[zw.old]][textlist.d[w]] = table.W[[zw.old]][textlist.d[w]]-1
          topicword.d = TopicWord(K, z[[d]][-w], textlist.d, table.W, alpha, mvec, beta, V)
        } else {
          topicword.d = matrix(0, nrow = length(z[[d]]), ncol = K)
        }
	      for (IP in unique(l)) {
	  	    lK = which(l == IP)
       	    z[[d]][w] = min(lK)
       	    p.d[d, ] = pdmat(list(z[[d]]), l, nIP)           
            history.t = History(edge, p.d, node, timestamps[hist.d[d]-1]+exp(-745), timeunit)
    	  	 	Xnew = Netstats_cpp(history.t, node, netstat)
		    mu[d, ] = mu_vec(p.d[d,], xi[[d]])
           	edgetime.d[lK] = Edgepartsum(Xnew, p.d[hist.d[d], ], b.old, u[[hist.d[d]]], delta)+
                           Timepart(mu[d,], sigma_tau, senders[d], timeinc[d])
	      }
        const.Z = edgetime.d+topicword.d[w, ]
        zw.new = multinom_vec(1, expconst(const.Z))
        if (zw.new != zw.old) {
          z[[d]][w] = zw.new
        } else {
          z[[d]][w] = zw.old
        }
        if (length(textlist.d) > 0) {
           	table.W[[z[[d]][w]]][textlist.d[w]] = table.W[[z[[d]][w]]][textlist.d[w]]+1
        }
      }
    }

    # C update 
    for (k in sort(unique(unlist(z)))) {
    	const.C = rep(NA, nIP)
      for (IP in 1:nIP) {
        l[k] = IP
        p.dnew = pdmat(z, l, nIP) 
        	history.t = History(edge, p.dnew, node, timestamps[max.edge-1]+exp(-745), timeunit)
       	Xnew = Netstats_cpp(history.t, node, netstat)
        Edgepartsum = Edgepartsum(Xnew, p.dnew[max.edge, ], b.old, u[[max.edge]], delta)
        mu = mu_mat(p.dnew, xi, edge.trim)
        Timepartsum = Timepartsum(mu, sigma_tau, senders, timeinc, edge.trim)
        prob = Edgepartsum+Timepartsum
        const.C[IP] = prob
      }
      l[k] = multinom_vec(1, expconst(const.C))
	}
 	  p.d = pdmat(z, l, nIP)  
 	  mu = mu_mat(p.d, xi, edge.trim)
 	  Timepartsum = Timepartsum(mu, sigma_tau, senders, timeinc, edge.trim)
 	  for (d in edge.trim) {
        history.t = History(edge, p.d, node, timestamps[d-1]+exp(-745), timeunit)
       	X[[d]] = Netstats_cpp(history.t, node, netstat)
      }
      
  # adaptive M-H   
    if (o > 1) {
    		accept.rates[1] = accept.rates[1]/Inner[1]
    		accept.rates[2] = accept.rates[2]/Inner[2]
        accept.rates[3] = accept.rates[3]/Inner[3]
        accept.rates[4] = accept.rates[1]
    		sigma.Q = adaptive.MH(sigma.Q, accept.rates, update.size = 0.2*sigma.Q)
    }
    accept.rates = rep(0, 4)
    
    prior.old1 = priorsum(prior.b[[2]], prior.b[[1]], b.old)+
    				 dnorm(delta, prior.delta[1], sqrt(prior.delta[2]), TRUE)
    post.old1 = Edgepartsum(X[[max.edge]], p.d[max.edge, ], b.old, u[[max.edge]], delta)
    b.new = matrix(NA, nIP, P)
    for (inner in 1:Inner[1]) {
      for (IP in 1:nIP) {
			  b.new[IP, ] = rmvnorm_arma(1, b.old[IP,], sigma.Q[1]*proposal.var1)
	  }
      delta.new = rnorm(1, delta, sqrt(sigma.Q[4]))
      prior.new1 = priorsum(prior.b[[2]], prior.b[[1]], b.new)+
    				 dnorm(delta.new, prior.delta[1], sqrt(prior.delta[2]), TRUE)
      post.new1 = Edgepartsum(X[[max.edge]], p.d[max.edge, ], b.new, u[[max.edge]], delta.new)
      loglike.diff = prior.new1+post.new1-prior.old1-post.old1
      if (log(runif(1, 0, 1)) < loglike.diff) {
        b.old = b.new
        delta = delta.new
        prior.old1 = prior.new1
        post.old1 = post.new1
        accept.rates[1] = accept.rates[1]+1
      }
      if (inner > burn[1] & inner %% (thin[1]) == 0) {
        for (IP in 1:nIP) {
          bmat[[IP]][ ,(inner-burn[1])/thin[1]] = b.old[IP,]
        }
        deltamat[(inner-burn[1])/thin[1]] = delta
      }
    }
	
	  prior.old2 = priorsum(prior.eta[[2]], prior.eta[[1]], eta.old)
	  post.old2 = Timepartsum
      eta.new = matrix(NA, nIP, Q)
      for (inner in 1:Inner[2]) {
          for (IP in 1:nIP) {
			  eta.new[IP, ] = rmvnorm_arma(1, eta.old[IP,], sigma.Q[2]*proposal.var2)
		  }
      xi = xi_all(timemat, eta.new[,node], eta.new[,-node], edge.trim)
      mu = mu_mat(p.d, xi, edge.trim)
      Timepartsum = Timepartsum(mu, sigma_tau, senders, timeinc, edge.trim)
      prior.new2 = priorsum(prior.eta[[2]], prior.eta[[1]], eta.old)
      post.new2 = Timepartsum
      loglike.diff = prior.new2+post.new2-prior.old2-post.old2
      if (log(runif(1, 0, 1)) < loglike.diff) {
        eta.old = eta.new
        prior.old2 = prior.new2
        post.old2 = post.new2
        accept.rates[2] = accept.rates[2]+1
      }
      if (inner > burn[2] & inner %% (thin[2]) == 0) {
        for (IP in 1:nIP) {
          etamat[[IP]][ ,(inner-burn[2])/thin[2]] = eta.old[IP,]
        }
      }
    }
	xi = xi_all(timemat, eta.old[,node], eta.old[,-node], edge.trim)
    mu = mu_mat(p.d, xi, edge.trim)
    Timepartsum = Timepartsum(mu, sigma_tau, senders, timeinc, edge.trim)
      
    prior.old3 = dhalfcauchy(sigma_tau, prior.tau, TRUE)
    post.old3 = Timepartsum
    for (inner in 1:Inner[3]) {
      sigma_tau.new = rtruncnorm(1, 0, Inf, sigma_tau, sqrt(sigma.Q[3]))
      prior.new3 = dhalfcauchy(sigma_tau.new, prior.tau, TRUE)
      post.new3 = Timepartsum(mu, sigma_tau.new, senders, timeinc, edge.trim)
      loglike.diff = log(dtruncnorm(sigma_tau, 0, Inf, sigma_tau.new, sqrt(sigma.Q[3])))-
                   log(dtruncnorm(sigma_tau.new, 0, Inf, sigma_tau, sqrt(sigma.Q[3])))+
                   prior.new3+post.new3-prior.old3-post.old3
      if (log(runif(1, 0, 1)) < loglike.diff) {
        sigma_tau = sigma_tau.new
        prior.old3 = prior.new3
        post.old3 = post.new3
        accept.rates[3] = accept.rates[3]+1
      }
      if (inner > burn[3] & inner %% (thin[3]) == 0) {
        sigma_taumat[(inner-burn[3])/thin[3]] = sigma_tau
      }
    }

    convergence[o] = Edgepartsum(X[[max.edge]], p.d[max.edge, ], b.old, u[[max.edge]], delta)+
                     Timepartsum(mu, sigma_tau, senders, timeinc, edge.trim)
  }
 
  chain.final = list(l = l, z = z, b = bmat, eta = etamat, delta = deltamat, sigma_tau = sigma_taumat,
                     u = u, sigma.Q =sigma.Q, alpha = alphavec, mvec = mvecmat, edge.trim = edge.trim,
                     convergence = convergence)
  return(chain.final)
}	


#' @title IPTM.inference.noIP
#' @description Iterate Markov Chain Monte Carlo (MCMC) algorithm for the interaction-partitioned topic model
#'
#' @param edge list of tie data with 3 elements (1: author, 2: recipient, 3: timestamp in unix.time format)
#' @param node vector of node id's (ID starting from 1)
#' @param textlist list of text containing the words in each document
#' @param vocab all vocabularies used over the corpus
#' @param nIP total number of interaction patterns
#' @param K total number of topics
#' @param sigma.Q proposal distribution variance parameter
#' @param alpha Dirichlet concentration prior for document-topic distribution
#' @param mvec Dirichlet base prior for document-topic distribution
#' @param beta Dirichlet concentration prior for topic-word distribution
#' @param prior.b prior mean and covariance of b in multivariate Normal distribution
#' @param prior.delta prior mean and variance of delta in Normal distribution
#' @param prior.eta prior mean and covariance of eta in multivariate Normal distribution
#' @param prior.tau prior shape and scale parameter of sigma_tau in inverse-Gamma distribution
#' @param Outer size of outer iterations 
#' @param Inner size of inner iteration for Metropolis-Hastings updates
#' @param burn iterations to be discarded at the beginning of Metropolis-Hastings chains
#' @param thin the thinning interval of Metropolis-Hastings chains
#' @param netstat which type of network statistics to use ("dyadic", "triadic", "degree")
#' @param timestat additional statistics to be used for timestamps other than netstat ("sender", "receiver","timeofday", "dayofweek")
#' @param optimize logical to optimize alpha (Dirichlet concentration prior for document-topic distribution)
#' @param initial list of initial values user wants to assign including (alpha, mvec, delta, b, eta, l, z, u, sigma_tau, proposal.var1, proposal.var2)
#' @param timeunit hour (= 3600) or day (=3600*24) and so on
#' @param tz timezone such as EST and PST
#'
#' @return MCMC output containing all parameter estimates
#'
#' @export

IPTM.inference.noIP = function(edge, node, textlist, vocab, nIP, K, sigma.Q, alpha, mvec, beta, 
                          prior.b, prior.delta, prior.eta, prior.tau, Outer, Inner, burn, thin, 
                          netstat, timestat, optimize = FALSE, initial = NULL, timeunit = 3600, tz = "America/New_York") {
  
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
  V = length(vocab)
  phi = matrix(NA, K, V)
  for (k in 1:K) {
    phi[k,] = rdirichlet_cpp(1, rep(beta/V, V))
  }
  if (length(initial) == 0) {
    theta = rdirichlet_cpp(D, alpha*mvec)
    delta = rnorm(1, prior.delta[1], sqrt(prior.delta[2]))
    sigma_tau = rhalfcauchy(1, prior.tau)
    b.old = rmvnorm_arma(nIP, prior.b[[1]], prior.b[[2]])
    eta.old = rmvnorm_arma(nIP, prior.eta[[1]], prior.eta[[2]])
    l = sample(1:nIP, K, replace = TRUE) 
    z = lapply(seq(along = edge), function(d) multinom_vec(max(1, length(textlist[[d]])), theta[d, ]))
    p.d = pdmat(z, l, nIP) 
    proposal.var1 = diag(P)
    proposal.var2 = diag(Q)
    sigma.Q = sigma.Q
    u = list()
    for (d in edge.trim) {
      u[[d]] = matrix(rbinom(A^2, 1, 1/A), nrow =A, ncol = A)
      diag(u[[d]]) = 0
    } 
  } else {
    theta = rdirichlet_cpp(D, initial$alpha*initial$mvec)
    delta = initial$delta
    sigma_tau = initial$sigma_tau
    b.old = initial$b
    eta.old = initial$eta
	l = initial$l
 	z = initial$z
  	p.d = pdmat(z, l, nIP) 
    proposal.var1 = initial$proposal.var1
    proposal.var2 = initial$proposal.var2
    sigma.Q = initial$sigma.Q
    u = initial$u
  }						 
  bmat = list()
  etamat = list()
  for (IP in 1:nIP) {
    bmat[[IP]] = matrix(b.old[IP,], nrow = P, ncol = (Inner[1]-burn[1])/thin[1])
    etamat[[IP]] = matrix(eta.old[IP,], nrow = Q, ncol = (Inner[2]-burn[2])/thin[2])
  }
  deltamat = rep(delta, (Inner[1]-burn[1])/thin[1])
  sigma_taumat = rep(sigma_tau, (Inner[3]-burn[3])/thin[3])		
  xi = xi_all(timemat, matrix(eta.old[,node], nrow = 1), matrix(eta.old[,-node], nrow = 1), edge.trim)
  mu = matrix(0, nrow = D, ncol = A)
  textlist.raw = unlist(textlist)
  alphavec = c()
  mvecmat = matrix(NA, nrow = 0, ncol = K)
  accept.rates = rep(0, 4)
  hist.d = c()
  for (d in edge.trim) {
  if (timestamps[d]+384*timeunit > timestamps[max.edge]) {
        hist.d[d] = max.edge
      } else {
        hist.d[d] = which_num(timestamps[d]+384*timeunit, timestamps)-1
      }
  }
  X = list()
  for (d in edge.trim) {
        history.t = History(edge, p.d, node, timestamps[d-1]+exp(-745), timeunit)
        X[[d]] = Netstats_cpp(history.t, node, netstat)
  }    
  #start outer iteration
  for (o in 1:Outer) {
    print(o)
    if (o == Outer) {
      Inner = Inner * 10
      burn = burn * 10
      for (IP in 1:nIP) {
        bmat[[IP]] = matrix(b.old[IP,], nrow = P, ncol = (Inner[1]-burn[1])/thin[1])
        etamat[[IP]] = matrix(eta.old[IP,], nrow = Q, ncol = (Inner[2]-burn[2])/thin[2])
      }
      deltamat = rep(delta, (Inner[1]-burn[1])/thin[1])
      sigma_taumat = rep(sigma_tau, (Inner[3]-burn[3])/thin[3])						 
    }
    if (optimize & o > 1) {
      #update the hyperparameter alpha and mvec
      vec = AlphamvecOpt(K, z, alpha, mvec, 5)
      alpha = sum(vec)
      mvec = vec/alpha
      alphavec = c(alphavec, alpha)
      mvecmat = rbind(mvecmat, mvec)
    }
      
    # Data augmentation
    for (d in edge.trim) {
        vu = MultiplyXB(X[[d]], b.old)
        lambda = lambda_cpp(p.d[d,], vu)
        for (i in node[-senders[d]]) {
            for (j in sample(node[-i], A-1)) {
                probij = u_Gibbs(u[[d]][i, ], lambda[i,], delta, j)
                u[[d]][i, j] = multinom_vec(1, expconst(probij))-1
            }
        }
        u[[d]][senders[d],] = tabulateC(as.numeric(unlist(edge[[d]][2])), A)
    }
    # Z update	
    table.W = lapply(1:K, function(k) tabulateC(textlist.raw[which(unlist(z[-emptytext]) == k)], V))
    for (d in 1:(edge.trim[1]-1)) {
	   	textlist.d = textlist[[d]] 
	   	for (w in 1:length(z[[d]])) {
	   		zw.old = z[[d]][w]
        if (length(textlist.d) > 0) {
          table.W[[zw.old]][textlist.d[w]] = table.W[[zw.old]][textlist.d[w]]-1
       	  topicword.d = TopicWord(K, z[[d]][-w], textlist.d, table.W, alpha, mvec, beta, V)
        } else {
          topicword.d = matrix(0, nrow = length(z[[d]]), ncol = K)
        }
        const.Z = topicword.d[w, ]
        zw.new = multinom_vec(1, expconst(const.Z))
        if (zw.new != zw.old) {
          z[[d]][w] = zw.new
        }
        if (length(textlist.d) > 0) {
           	table.W[[z[[d]][w]]][textlist.d[w]] = table.W[[z[[d]][w]]][textlist.d[w]]+1
        } 
      }
    }
    for (d in edge.trim) {
      textlist.d = textlist[[d]]
      edgetime.d = rep(NA, K)
      for (w in 1:length(z[[d]])) {
       	zw.old = z[[d]][w]
       	if (length(textlist.d) > 0) {
       	  table.W[[zw.old]][textlist.d[w]] = table.W[[zw.old]][textlist.d[w]]-1
          topicword.d = TopicWord(K, z[[d]][-w], textlist.d, table.W, alpha, mvec, beta, V)
        } else {
          topicword.d = matrix(0, nrow = length(z[[d]]), ncol = K)
        }
	      for (IP in unique(l)) {
	  	    lK = which(l == IP)
       	    z[[d]][w] = min(lK)
            history.t = History(edge, p.d, node, timestamps[hist.d[d]-1]+exp(-745), timeunit)
    	  	 	Xnew = Netstats_cpp(history.t, node, netstat)
		    mu[d, ] = mu_vec(p.d[d,], xi[[d]])
           	edgetime.d[lK] = Edgepartsum(Xnew, p.d[hist.d[d], ], b.old, u[[hist.d[d]]], delta)+
                           Timepart(mu[d,], sigma_tau, senders[d], timeinc[d])
	      }
        const.Z = edgetime.d+topicword.d[w, ]
        zw.new = multinom_vec(1, expconst(const.Z))
        if (zw.new != zw.old) {
          z[[d]][w] = zw.new
        } else {
          z[[d]][w] = zw.old
        }
        if (length(textlist.d) > 0) {
           	table.W[[z[[d]][w]]][textlist.d[w]] = table.W[[z[[d]][w]]][textlist.d[w]]+1
        }
      }
    }
 	  mu = mu_mat(p.d, xi, edge.trim)
 	  Timepartsum = Timepartsum(mu, sigma_tau, senders, timeinc, edge.trim)
 	  for (d in edge.trim) {
        history.t = History(edge, p.d, node, timestamps[d-1]+exp(-745), timeunit)
       	X[[d]] = Netstats_cpp(history.t, node, netstat)
      }
      
  # adaptive M-H   
    if (o > 1) {
    		accept.rates[1] = accept.rates[1]/Inner[1]
    		accept.rates[2] = accept.rates[2]/Inner[2]
        accept.rates[3] = accept.rates[3]/Inner[3]
        accept.rates[4] = accept.rates[1]
    		sigma.Q = adaptive.MH(sigma.Q, accept.rates, update.size = 0.2*sigma.Q)
    }
    accept.rates = rep(0, 4)
    
    prior.old1 = priorsum(prior.b[[2]], prior.b[[1]], b.old)+
    				 dnorm(delta, prior.delta[1], sqrt(prior.delta[2]), TRUE)
    post.old1 = Edgepartsum(X[[max.edge]], p.d[max.edge, ], b.old, u[[max.edge]], delta)
    b.new = matrix(NA, nIP, P)
    for (inner in 1:Inner[1]) {
      for (IP in 1:nIP) {
			  b.new[IP, ] = rmvnorm_arma(1, b.old[IP,], sigma.Q[1]*proposal.var1)
	  }
      delta.new = rnorm(1, delta, sqrt(sigma.Q[4]))
      prior.new1 = priorsum(prior.b[[2]], prior.b[[1]], b.new)+
    				 dnorm(delta.new, prior.delta[1], sqrt(prior.delta[2]), TRUE)
      post.new1 = Edgepartsum(X[[max.edge]], p.d[max.edge, ], b.new, u[[max.edge]], delta.new)
      loglike.diff = prior.new1+post.new1-prior.old1-post.old1
      if (log(runif(1, 0, 1)) < loglike.diff) {
        b.old = b.new
        delta = delta.new
        prior.old1 = prior.new1
        post.old1 = post.new1
        accept.rates[1] = accept.rates[1]+1
      }
      if (inner > burn[1] & inner %% (thin[1]) == 0) {
        for (IP in 1:nIP) {
          bmat[[IP]][ ,(inner-burn[1])/thin[1]] = b.old[IP,]
        }
        deltamat[(inner-burn[1])/thin[1]] = delta
      }
    }
	
	  prior.old2 = priorsum(prior.eta[[2]], prior.eta[[1]], eta.old)
	  post.old2 = Timepartsum
      eta.new = matrix(NA, nIP, Q)
      for (inner in 1:Inner[2]) {
          for (IP in 1:nIP) {
			  eta.new[IP, ] = rmvnorm_arma(1, eta.old[IP,], sigma.Q[2]*proposal.var2)
		  }
      xi = xi_all(timemat, matrix(eta.new[,node], nrow = 1), matrix(eta.new[,-node], nrow = 1), edge.trim)
      mu = mu_mat(p.d, xi, edge.trim)
      Timepartsum = Timepartsum(mu, sigma_tau, senders, timeinc, edge.trim)
      prior.new2 = priorsum(prior.eta[[2]], prior.eta[[1]], eta.old)
      post.new2 = Timepartsum
      loglike.diff = prior.new2+post.new2-prior.old2-post.old2
      if (log(runif(1, 0, 1)) < loglike.diff) {
        eta.old = eta.new
        prior.old2 = prior.new2
        post.old2 = post.new2
        accept.rates[2] = accept.rates[2]+1
      }
      if (inner > burn[2] & inner %% (thin[2]) == 0) {
        for (IP in 1:nIP) {
          etamat[[IP]][ ,(inner-burn[2])/thin[2]] = eta.old[IP,]
        }
      }
    }
	xi = xi_all(timemat, matrix(eta.old[,node], nrow = 1), matrix(eta.old[,-node], nrow = 1), edge.trim)
    mu = mu_mat(p.d, xi, edge.trim)
    Timepartsum = Timepartsum(mu, sigma_tau, senders, timeinc, edge.trim)
      
    prior.old3 = dhalfcauchy(sigma_tau, prior.tau, TRUE)
    post.old3 = Timepartsum
    for (inner in 1:Inner[3]) {
      sigma_tau.new = rtruncnorm(1, 0, Inf, sigma_tau, sqrt(sigma.Q[3]))
      prior.new3 = dhalfcauchy(sigma_tau.new, prior.tau, TRUE)
      post.new3 = Timepartsum(mu, sigma_tau.new, senders, timeinc, edge.trim)
      loglike.diff = log(dtruncnorm(sigma_tau, 0, Inf, sigma_tau.new, sqrt(sigma.Q[3])))-
                   log(dtruncnorm(sigma_tau.new, 0, Inf, sigma_tau, sqrt(sigma.Q[3])))+
                   prior.new3+post.new3-prior.old3-post.old3
      if (log(runif(1, 0, 1)) < loglike.diff) {
        sigma_tau = sigma_tau.new
        prior.old3 = prior.new3
        post.old3 = post.new3
        accept.rates[3] = accept.rates[3]+1
      }
      if (inner > burn[3] & inner %% (thin[3]) == 0) {
        sigma_taumat[(inner-burn[3])/thin[3]] = sigma_tau
      }
    }

    convergence[o] = Edgepartsum(X[[max.edge]], p.d[max.edge, ], b.old, u[[max.edge]], delta)+
                     Timepartsum(mu, sigma_tau, senders, timeinc, edge.trim)
  }
 
  chain.final = list(l = l, z = z, b = bmat, eta = etamat, delta = deltamat, sigma_tau = sigma_taumat,
                     u = u, sigma.Q =sigma.Q, alpha = alphavec, mvec = mvecmat, edge.trim = edge.trim,
                     convergence = convergence)
  return(chain.final)
}	


#' @title IPTM.inference.PPE
#' @description Iterate Markov Chain Monte Carlo (MCMC) algorithm for the interaction-partitioned topic model
#'
#' @param missing D x 3 binary matrix where 1 denotes missing
#' @param edge list of tie data with 3 elements (1: author, 2: recipient, 3: timestamp in unix.time format)
#' @param node vector of node id's (ID starting from 1)
#' @param textlist list of text containing the words in each document
#' @param vocab all vocabularies used over the corpus
#' @param nIP total number of interaction patterns
#' @param K total number of topics
#' @param sigma.Q proposal distribution variance parameter
#' @param alpha Dirichlet concentration prior for document-topic distribution
#' @param mvec Dirichlet base prior for document-topic distribution
#' @param beta Dirichlet concentration prior for topic-word distribution
#' @param prior.b prior mean and covariance of b in multivariate Normal distribution
#' @param prior.delta prior mean and variance of delta in Normal distribution
#' @param prior.eta prior mean and covariance of eta in multivariate Normal distribution
#' @param prior.tau prior shape and scale parameter of sigma_tau in inverse-Gamma distribution
#' @param Outer size of outer iterations 
#' @param Inner size of inner iteration for Metropolis-Hastings updates
#' @param burn iterations to be discarded at the beginning of Metropolis-Hastings chains
#' @param thin the thinning interval of Metropolis-Hastings chains
#' @param netstat which type of network statistics to use ("dyadic", "triadic", "degree")
#' @param timestat additional statistics to be used for timestamps other than netstat ("sender", "receiver","timeofday", "dayofweek")
#' @param optimize logical to optimize alpha (Dirichlet concentration prior for document-topic distribution)
#' @param initial list of initial values user wants to assign including (alpha, mvec, delta, b, eta, l, z, u, sigma_tau, proposal.var1, proposal.var2)
#' @param timeunit hour (= 3600) or day (= 3600*24) and so on
#' @param tz timezone such as EST and PST
#'
#' @return MCMC output containing all parameter estimates
#'
#' @export

IPTM.inference.PPE = function(missing, edge, node, textlist, vocab, nIP, K, sigma.Q, alpha, mvec, beta, 
                          prior.b, prior.delta, prior.eta, prior.tau, Outer, Inner, burn, thin, 
                          netstat, timestat, optimize = FALSE, initial = NULL, timeunit = 3600, tz = "America/New_York") {
  
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
  V = length(vocab)
  phi = matrix(NA, K, V)
  for (k in 1:K) {
    phi[k,] = rdirichlet_cpp(1, rep(beta/V, V))
  }
  if (length(initial) == 0) {
    theta = rdirichlet_cpp(D, alpha*mvec)
    delta = rnorm(1, prior.delta[1], sqrt(prior.delta[2]))
    sigma_tau = rhalfcauchy(1, prior.tau)
    b.old = rmvnorm_arma(nIP, prior.b[[1]], prior.b[[2]])
    eta.old = rmvnorm_arma(nIP, prior.eta[[1]], prior.eta[[2]])
    l = sample(1:nIP, K, replace = TRUE) 
    z = lapply(seq(along = edge), function(d) multinom_vec(max(1, length(textlist[[d]])), theta[d, ]))
    p.d = pdmat(z, l, nIP) 
    proposal.var1 = diag(P)
    proposal.var2 = diag(Q)
    sigma.Q = sigma.Q
    u = list()
    for (d in edge.trim) {
      u[[d]] = matrix(rbinom(A^2, 1, 1/A), nrow =A, ncol = A)
      diag(u[[d]]) = 0
    } 
  } else {
    theta = rdirichlet_cpp(D, initial$alpha*initial$mvec)
    delta = initial$delta
    sigma_tau = initial$sigma_tau
    b.old = initial$b
    eta.old = initial$eta
	l = initial$l
 	z = initial$z
  	p.d = pdmat(z, l, nIP) 
    proposal.var1 = initial$proposal.var1
    proposal.var2 = initial$proposal.var2
    sigma.Q = initial$sigma.Q
    u = initial$u
  }						 
  bmat = list()
  etamat = list()
  for (IP in 1:nIP) {
    bmat[[IP]] = matrix(b.old[IP,], nrow = P, ncol = (Inner[1]-burn[1])/thin[1])
    etamat[[IP]] = matrix(eta.old[IP,], nrow = Q, ncol = (Inner[2]-burn[2])/thin[2])
  }
  deltamat = rep(delta, (Inner[1]-burn[1])/thin[1])
  sigma_taumat = rep(sigma_tau, (Inner[3]-burn[3])/thin[3])						 
  mu = matrix(0, nrow = D, ncol = A)
  textlist.raw = unlist(textlist)
  alphavec = c()
  mvecmat = matrix(NA, nrow = 0, ncol = K)
  accept.rates = rep(0, 4)
  hist.d = c()
  for (d in edge.trim) {
  if (timestamps[d]+384*timeunit > timestamps[max.edge]) {
        hist.d[d] = max.edge
      } else {
        hist.d[d] = which_num(timestamps[d]+384*timeunit, timestamps)-1
      }
  }
  X = list()
  senderpredict = matrix(NA, nrow = nrow(missing), ncol = Outer)
  receiverpredict = lapply(1:nrow(missing), function(d) {c()})
  timepredict = matrix(NA, nrow = nrow(missing), ncol = Outer)
  xi = xi_all(timemat, eta.old[,node], eta.old[,-node], edge.trim)
  mu = mu_mat(p.d, xi, edge.trim)
  #start outer iteration
  for (o in 1:Outer) {
    print(o)
    if (optimize & o > 1) {
      #update the hyperparameter alpha and mvec
      vec = AlphamvecOpt(K, z, alpha, mvec, 5)
      alpha = sum(vec)
      mvec = vec/alpha
      alphavec = c(alphavec, alpha)
      mvecmat = rbind(mvecmat, mvec)
    }
    
    #imputation
    iter1 = 1
    iter2 = 1
    iter3 = 1
    for (d in edge.trim) {
        if (missing[d,1] == 1) {
            probi = Timepartindiv(mu[d,], sigma_tau, timeinc[d])
            senders[d] = multinom_vec(1, expconst(probi))
            senderpredict[iter1, o] = senders[d]
            iter1 = iter1+1
        }
        if (missing[d,2] == 1) {
            edge[[d]][[2]] = which(u[[d]][senders[d], ] > 0)
            receiverpredict[[iter2]] = rbind(receiverpredict[[iter2]], u[[d]][senders[d], ])
            iter2 = iter2+1
         }
        if (missing[d,3] == 1) {
            timeinc[d] = rlnorm(1, mu[d, senders[d]], sigma_tau)
            timepredict[iter3, o] = timeinc[d]
            iter3 = iter3+1
        }
    }
    timeinc[timeinc==0] = runif(sum(timeinc==0), 0, min(timeinc[timeinc!=0]))
    timestamps[-1] = timestamps[1]+cumsum(timeinc[-1])*timeunit

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
    
    xi = xi_all(timemat, eta.old[,node], eta.old[,-node], edge.trim)
    mu = mu_mat(p.d, xi, edge.trim)
    
    #start inference
    for (d in edge.trim) {
        history.t = History(edge, p.d, node, timestamps[d-1]+exp(-745), timeunit)
        X[[d]] = Netstats_cpp(history.t, node, netstat)
        if (timestamps[d]+384*timeunit > timestamps[max.edge]) {
            hist.d[d] = max.edge
        } else {
            hist.d[d] = which_num(timestamps[d]+384*timeunit, timestamps)-1
        }
    }
    
    # Data augmentation
    for (d in edge.trim) {
        vu = MultiplyXB(X[[d]], b.old)
        lambda = lambda_cpp(p.d[d,], vu)
        for (i in node[-senders[d]]) {
            for (j in sample(node[-i], A-1)) {
                probij = u_Gibbs(u[[d]][i, ], lambda[i,], delta, j)
                u[[d]][i, j] = multinom_vec(1, expconst(probij))-1
            }
        }
        u[[d]][senders[d],] = tabulateC(as.numeric(unlist(edge[[d]][2])), A)
    }
    # Z update	
    table.W = lapply(1:K, function(k) tabulateC(textlist.raw[which(unlist(z[-emptytext]) == k)], V))
    for (d in 1:(edge.trim[1]-1)) {
	   	textlist.d = textlist[[d]] 
	   	for (w in 1:length(z[[d]])) {
	   		zw.old = z[[d]][w]
        if (length(textlist.d) > 0) {
          table.W[[zw.old]][textlist.d[w]] = table.W[[zw.old]][textlist.d[w]]-1
       	  topicword.d = TopicWord(K, z[[d]][-w], textlist.d, table.W, alpha, mvec, beta, V)
        } else {
          topicword.d = matrix(0, nrow = length(z[[d]]), ncol = K)
        }
        const.Z = topicword.d[w, ]
        zw.new = multinom_vec(1, expconst(const.Z))
        if (zw.new != zw.old) {
          z[[d]][w] = zw.new
      	  p.d[d, ] = pdmat(list(z[[d]]), l, nIP)
        }
        if (length(textlist.d) > 0) {
           	table.W[[z[[d]][w]]][textlist.d[w]] = table.W[[z[[d]][w]]][textlist.d[w]]+1
        } 
      }
    }
    for (d in edge.trim) {
      textlist.d = textlist[[d]]
      edgetime.d = rep(NA, K)
      for (w in 1:length(z[[d]])) {
       	zw.old = z[[d]][w]
       	if (length(textlist.d) > 0) {
       	  table.W[[zw.old]][textlist.d[w]] = table.W[[zw.old]][textlist.d[w]]-1
          topicword.d = TopicWord(K, z[[d]][-w], textlist.d, table.W, alpha, mvec, beta, V)
        } else {
          topicword.d = matrix(0, nrow = length(z[[d]]), ncol = K)
        }
	      for (IP in unique(l)) {
	  	    lK = which(l == IP)
       	    z[[d]][w] = min(lK)
       	    p.d[d, ] = pdmat(list(z[[d]]), l, nIP)           
            history.t = History(edge, p.d, node, timestamps[hist.d[d]-1]+exp(-745), timeunit)
    	   	Xnew = Netstats_cpp(history.t, node, netstat)
            mu[d, ] = mu_vec(p.d[d,], xi[[d]])
           	edgetime.d[lK] = Edgepartsum(Xnew, p.d[hist.d[d], ], b.old, u[[hist.d[d]]], delta)+
                           Timepart(mu[d,], sigma_tau, senders[d], timeinc[d])
	      }
        const.Z = edgetime.d+topicword.d[w, ]
        zw.new = multinom_vec(1, expconst(const.Z))
        if (zw.new != zw.old) {
          z[[d]][w] = zw.new
        } else {
          z[[d]][w] = zw.old
        }
        if (length(textlist.d) > 0) {
           	table.W[[z[[d]][w]]][textlist.d[w]] = table.W[[z[[d]][w]]][textlist.d[w]]+1
        }
      }
    }

    # C update 
    for (k in sort(unique(unlist(z)))) {
    	const.C = rep(NA, nIP)
      for (IP in 1:nIP) {
        l[k] = IP
        p.dnew = pdmat(z, l, nIP) 
		history.t = History(edge, p.dnew, node, timestamps[max.edge-1]+exp(-745), timeunit)
        Xnew = Netstats_cpp(history.t, node, netstat)
        Edgepartsum = Edgepartsum(Xnew, p.dnew[max.edge, ], b.old, u[[max.edge]], delta)
        mu = mu_mat(p.dnew, xi, edge.trim)
        Timepartsum = Timepartsum(mu, sigma_tau, senders, timeinc, edge.trim)
        prob = Edgepartsum+Timepartsum
        const.C[IP] = prob
      }
      l[k] = multinom_vec(1, expconst(const.C))
	  }
 	  p.d = pdmat(z, l, nIP)  
 	  mu = mu_mat(p.d, xi, edge.trim)
 	  Timepartsum = Timepartsum(mu, sigma_tau, senders, timeinc, edge.trim)
 	  for (d in edge.trim) {
        history.t = History(edge, p.d, node, timestamps[d-1]+exp(-745), timeunit)
       	X[[d]] = Netstats_cpp(history.t, node, netstat)
      }
      
  # adaptive M-H   
    if (o > 1) {
    		accept.rates[1] = accept.rates[1]/Inner[1]
    		accept.rates[2] = accept.rates[2]/Inner[2]
        accept.rates[3] = accept.rates[3]/Inner[3]
        accept.rates[4] = accept.rates[1]
    		sigma.Q = adaptive.MH(sigma.Q, accept.rates, update.size = 0.2*sigma.Q)
    }
    accept.rates = rep(0, 4)
    
    prior.old1 = priorsum(prior.b[[2]], prior.b[[1]], b.old)+
    				 dnorm(delta, prior.delta[1], sqrt(prior.delta[2]), TRUE)
    post.old1 = Edgepartsum(X[[max.edge]], p.d[max.edge, ], b.old, u[[max.edge]], delta)
    b.new = matrix(NA, nIP, P)
    for (inner in 1:Inner[1]) {
      for (IP in 1:nIP) {
			  b.new[IP, ] = rmvnorm_arma(1, b.old[IP,], sigma.Q[1]*proposal.var1)
	  }
      delta.new = rnorm(1, delta, sqrt(sigma.Q[4]))
      prior.new1 = priorsum(prior.b[[2]], prior.b[[1]], b.new)+
    				 dnorm(delta.new, prior.delta[1], sqrt(prior.delta[2]), TRUE)
      post.new1 = Edgepartsum(X[[max.edge]], p.d[max.edge, ], b.new, u[[max.edge]], delta.new)
      loglike.diff = prior.new1+post.new1-prior.old1-post.old1
      if (log(runif(1, 0, 1)) < loglike.diff) {
        b.old = b.new
        delta = delta.new
        prior.old1 = prior.new1
        post.old1 = post.new1
        accept.rates[1] = accept.rates[1]+1
      }
      if (inner > burn[1] & inner %% (thin[1]) == 0) {
        for (IP in 1:nIP) {
          bmat[[IP]][ ,(inner-burn[1])/thin[1]] = b.old[IP,]
        }
        deltamat[(inner-burn[1])/thin[1]] = delta
      }
    }
	
	  prior.old2 = priorsum(prior.eta[[2]], prior.eta[[1]], eta.old)
	  post.old2 = Timepartsum
      eta.new = matrix(NA, nIP, Q)
      for (inner in 1:Inner[2]) {
          for (IP in 1:nIP) {
			  eta.new[IP, ] = rmvnorm_arma(1, eta.old[IP,], sigma.Q[2]*proposal.var2)
		  }
      xi = xi_all(timemat, eta.new[,node], eta.new[,-node], edge.trim)
      mu = mu_mat(p.d, xi, edge.trim)
      Timepartsum = Timepartsum(mu, sigma_tau, senders, timeinc, edge.trim)
      prior.new2 = priorsum(prior.eta[[2]], prior.eta[[1]], eta.old)
      post.new2 = Timepartsum
      loglike.diff = prior.new2+post.new2-prior.old2-post.old2
      if (log(runif(1, 0, 1)) < loglike.diff) {
        eta.old = eta.new
        prior.old2 = prior.new2
        post.old2 = post.new2
        accept.rates[2] = accept.rates[2]+1
      }
      if (inner > burn[2] & inner %% (thin[2]) == 0) {
        for (IP in 1:nIP) {
          etamat[[IP]][ ,(inner-burn[2])/thin[2]] = eta.old[IP,]
        }
      }
    }
	  xi = xi_all(timemat, eta.old[,node], eta.old[,-node], edge.trim)
    mu = mu_mat(p.d, xi, edge.trim)
    Timepartsum = Timepartsum(mu, sigma_tau, senders, timeinc, edge.trim)

    prior.old3 = dhalfcauchy(sigma_tau, prior.tau, TRUE)
    post.old3 = Timepartsum
    for (inner in 1:Inner[3]) {
      sigma_tau.new = rtruncnorm(1, 0, Inf, sigma_tau, sqrt(sigma.Q[3]))
      prior.new3 = dhalfcauchy(sigma_tau.new, prior.tau, TRUE)
      post.new3 = Timepartsum(mu, sigma_tau.new, senders, timeinc, edge.trim)
      loglike.diff = log(dtruncnorm(sigma_tau, 0, Inf, sigma_tau.new, sqrt(sigma.Q[3])))-
                   log(dtruncnorm(sigma_tau.new, 0, Inf, sigma_tau, sqrt(sigma.Q[3])))+
                   prior.new3+post.new3-prior.old3-post.old3
      if (log(runif(1, 0, 1)) < loglike.diff) {
        sigma_tau = sigma_tau.new
        prior.old3 = prior.new3
        post.old3 = post.new3
        accept.rates[3] = accept.rates[3]+1
      }
      if (inner > burn[3] & inner %% (thin[3]) == 0) {
        sigma_taumat[(inner-burn[3])/thin[3]] = sigma_tau
      }
    }

    convergence[o] = Edgepartsum(X[[max.edge]], p.d[max.edge, ], b.old, u[[max.edge]], delta)+
                     Timepartsum(mu, sigma_tau, senders, timeinc, edge.trim)
  }
 
  chain.final = list(l = l, z = z, b = bmat, eta = etamat, delta = deltamat, sigma_tau = sigma_taumat,
                     u = u, sigma.Q = sigma.Q, alpha = alphavec, mvec = mvecmat, edge.trim = edge.trim,
                     convergence = convergence, senderpredict = senderpredict, receiverpredict = receiverpredict, 
                     timepredict = timepredict)
  return(chain.final)
}	


#' @title IPTM.PPE
#' @description Iterate Markov Chain Monte Carlo (MCMC) algorithm for the interaction-partitioned topic model
#'
#' @param missing D x 3 binary matrix where 1 denotes missing
#' @param edge list of tie data with 3 elements (1: author, 2: recipient, 3: timestamp in unix.time format)
#' @param node vector of node id's (ID starting from 1)
#' @param textlist list of text containing the words in each document
#' @param vocab all vocabularies used over the corpus
#' @param nIP total number of interaction patterns
#' @param K total number of topics
#' @param Outer size of outer iterations 
#' @param netstat which type of network statistics to use ("dyadic", "triadic", "degree")
#' @param timestat additional statistics to be used for timestamps other than netstat ("sender", "receiver","timeofday", "dayofweek")
#' @param initial list of initial values user wants to assign including (alpha, mvec, delta, b, eta, l, z, u, sigma_tau, proposal.var1, proposal.var2)
#' @param timeunit hour (= 3600) or day (= 3600*24) and so on
#' @param tz timezone such as EST and PST
#'
#' @return MCMC output containing all parameter estimates
#'
#' @export

IPTM.PPE = function(missing, edge, node, textlist, vocab, nIP, K, 
                    Outer, netstat, timestat, initial = NULL, timeunit = 3600, tz = "America/New_York") {
  
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
  Q = length(node)+sum(timestat)
  V = length(vocab)
    theta = rdirichlet_cpp(D, initial$alpha*initial$mvec)
    delta = initial$delta
    sigma_tau = initial$sigma_tau
    b.old = initial$b
    eta.old = initial$eta
    l = initial$l
    z = initial$z
    p.d = pdmat(z, l, nIP) 
    u = initial$u
  X = list()
  if (nIP == 1) {
  	  xi = xi_all(timemat, matrix(eta.old[,node], nrow = 1), matrix(eta.old[,-node], nrow = 1), edge.trim)
  } else {
  	  xi = xi_all(timemat, eta.old[,node], eta.old[,-node], edge.trim)
  }
  mu = mu_mat(p.d, xi, edge.trim)
  senderpredict = matrix(NA, nrow = nrow(missing), ncol = Outer)
  receiverpredict = lapply(1:nrow(missing), function(d) {c()})
  timepredict = matrix(NA, nrow = nrow(missing), ncol = Outer)

  #start outer iteration
  for (o in 1:Outer) {
    print(o)
    #imputation
    iter1 = 1
    iter2 = 1
    iter3 = 1
    for (d in edge.trim) {
      if (missing[d,3] == 1) {
        timeinc[d] = rlnorm(1, mu[d, senders[d]], sigma_tau)
        timepredict[iter3, o] = timeinc[d]
        iter3 = iter3+1
      }
      if (missing[d,1] == 1) {
        probi = Timepartindiv(mu[d,], sigma_tau, timeinc[d])
        senders[d] = multinom_vec(1, expconst(probi))
        senderpredict[iter1, o] = senders[d]
        iter1 = iter1+1
      }
    }
    timeinc[timeinc==0] = runif(sum(timeinc==0), 0, min(timeinc[timeinc!=0]))
    timestamps[-1] = timestamps[1]+cumsum(timeinc[-1])*timeunit
    for (d in edge.trim) {
      history.t = History(edge, p.d, node, timestamps[d-1]+exp(-745), timeunit)
      X[[d]] = Netstats_cpp(history.t, node, netstat)
    }
    for (d in edge.trim) {
      if (missing[d,2] == 1) {
        vu = MultiplyXB(X[[d]], b.old)
        lambda = lambda_cpp(p.d[d,], vu)
        i = senders[d]
          for (j in sample(node[-i], A-1)) {
            probij = u_Gibbs(u[[d]][i, ], lambda[i,], delta, j)
            u[[d]][i, j] = multinom_vec(1, expconst(probij))-1
          }
        receiverpredict[[iter2]] = rbind(receiverpredict[[iter2]], u[[d]][i, ])
        iter2 = iter2+1
      }
    }
  }
  chain.final = list(senderpredict = senderpredict, receiverpredict = receiverpredict, 
                     timepredict = timepredict)
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
#' @param alpha Dirichlet concentration prior for document-topic distribution
#' @param mvec Dirichlet base prior for document-topic distribution
#' @param beta Dirichlet concentration prior for topic-word distribution
#' @param prior.b prior mean and covariance of b in multivariate Normal distribution
#' @param prior.delta prior mean and variance of delta in Normal distribution
#' @param prior.eta prior mean and covariance of eta in multivariate Normal distribution
#' @param prior.tau prior shape and scale parameter of sigma_tau in inverse-Gamma distribution
#' @param Outer size of outer iterations
#' @param Inner size of inner iteration for Metropolis-Hastings updates
#' @param burn iterations to be discarded at the beginning of Metropolis-Hastings chains
#' @param thin the thinning interval of Metropolis-Hastings chains
#' @param netstat which type of network statistics to use ("dyadic", "triadic", "degree")
#' @param timestat additional statistics to be used for timestamps other than netstat ("sender", "receiver", "timeofday", "dayofweek")
#' @param optimize logical to optimize alpha (Dirichlet concentration prior for document-topic distribution)
#' @param initial list of initial values user wants to assign including (alpha, mvec, delta, b, eta, l, z, u, sigma_tau, proposal.var1, proposal.var2)
#' @param timeunit hour (= 3600) or day (=3600*24) and so on
#' @param tz timezone such as EST and PST
#'
#' @return MCMC output containing all parameter estimates
#'
#' @export

IPTM.inference.GiR = function(edge, node, textlist, vocab, nIP, K, sigma.Q, alpha, mvec, beta, 
                              prior.b, prior.delta, prior.eta, prior.tau, Outer, Inner, burn, thin, 
                              netstat, timestat, optimize = FALSE, initial = NULL, timeunit = 3600, tz = "America/New_York") {
    
  # trim the edge so that we only model edges after 384 hours
  A = length(node)
  D = length(edge)
  timestamps = vapply(edge, function(d) { d[[3]]}, c(1))
  senders = vapply(edge, function(d) { d[[1]] }, c(1))
  edge.trim = which_num(384*timeunit, timestamps-timestamps[1]):D
  max.edge = max(edge.trim)
  timeinc = c(timestamps[1], timestamps[-1]-timestamps[-length(timestamps)])/timeunit
  timeinc[timeinc==0] = runif(sum(timeinc==0), 0, min(timeinc[timeinc!=0]))
  # initialization
  netstat = as.numeric(c("degree", "dyadic", "triadic") %in% netstat)
  timestat = as.numeric(c("dayofweek","timeofday") %in% timestat)
  timemat = matrix(0, nrow = D, ncol = sum(timestat))
  if (sum(timestat) > 0) {
    Sys.setenv(TZ = tz)
    unixtime = timestamps
    time_ymd = as.POSIXct(unixtime[which(unixtime<exp(38.7))], tz = getOption("tz"), origin = "1970-01-01")
    if (timestat[1] > 0) {
      days = vapply(time_ymd, function(d) {wday(d)}, c(1))
      days[days==1] = 8
      timemat[1:length(time_ymd),1] = as.numeric(cut(days, c(1,6,8), c("weekdays","weekends")))-1
      it = 1
    }
    if (timestat[2] > 0) {
      hours = vapply(time_ymd, function(d) {hour(d)}, c(1))
      timemat[1:length(time_ymd),it+1] = as.numeric(cut(hours, c(-1,12,24), c("AM", "PM")))-1
    }
  }
  L = 3
  P = L*(2*netstat[1]+2*netstat[2]+4*netstat[3])
  Q = length(prior.eta[[1]])
  V = length(vocab)
  phi = matrix(NA, K, V)
  for (k in 1:K) {
    phi[k,] = rdirichlet_cpp(1, rep(beta/V, V))
  }
  if (length(initial) == 0) {
    theta = rdirichlet_cpp(D, alpha*mvec)
    delta = rnorm(1, prior.delta[1], sqrt(prior.delta[2]))
    sigma_tau = rhalfcauchy(1, prior.tau)
    b.old = rmvnorm_arma(nIP, prior.b[[1]], prior.b[[2]])
    eta.old = rmvnorm_arma(nIP, prior.eta[[1]], prior.eta[[2]])
    l = sample(1:nIP, K, replace = TRUE)
    z = lapply(seq(along = edge), function(d) multinom_vec(max(1, length(textlist[[d]])), theta[d, ]))
    p.d = pdmat(z, l, nIP)
    proposal.var1 = diag(P)
    proposal.var2 = diag(Q)
    sigma.Q = sigma.Q
    u = list()
    for (d in edge.trim) {
      u[[d]] = matrix(rbinom(A^2, 1, 1/A), nrow =A, ncol = A)
      diag(u[[d]]) = 0
    }
  } else {
    theta = rdirichlet_cpp(D, initial$alpha*initial$mvec)
    delta = initial$delta
    sigma_tau = initial$sigma_tau
    b.old = initial$b
    eta.old = initial$eta
    l = initial$l
    z = initial$z
    p.d = pdmat(z, l, nIP)
    proposal.var1 = initial$proposal.var1
    proposal.var2 = initial$proposal.var2
    sigma.Q = initial$sigma.Q
    u = initial$u
  }
  bmat = list()
  etamat = list()
  for (IP in 1:nIP) {
    bmat[[IP]] = matrix(b.old[IP,], nrow = P, ncol = (Inner[1]-burn[1])/thin[1])
    etamat[[IP]] = matrix(eta.old[IP,], nrow = Q, ncol = (Inner[2]-burn[2])/thin[2])
  }
  deltamat = rep(delta, (Inner[1]-burn[1])/thin[1])
  sigma_taumat = rep(sigma_tau, (Inner[3]-burn[3])/thin[3])
  xi = xi_all(timemat, eta.old[,node], eta.old[,-node], edge.trim)
  mu = matrix(0, nrow = D, ncol = A)
  textlist.raw = unlist(textlist[edge.trim])
  alphavec = c()
  mvecmat = matrix(NA, nrow = 0, ncol = K)
  accept.rates = rep(0, 4)
  hist.d = c()
  for (d in edge.trim) {
  if (timestamps[d]+384*timeunit >= timestamps[max.edge]) {
        hist.d[d] = max.edge
      } else {
        hist.d[d] = which_num(timestamps[d]+384*timeunit, timestamps)-1
        if (hist.d[d] < min(edge.trim)) {hist.d[d] = d}
      }
  }
  X = list()
  for (d in edge.trim) {
        history.t = History(edge, p.d, node, timestamps[d-1]+exp(-745), timeunit)
        X[[d]] = Netstats_cpp(history.t, node, netstat)
  }
  #start outer iteration
  for (o in 1:Outer) {
    if (optimize & o > 1) {
    #update the hyperparameter alpha and mvec
      vec = AlphamvecOpt(K, z[edge.trim], alpha, mvec, 5)
      alpha = sum(vec)
      mvec = vec/alpha
      alphavec = c(alphavec, alpha)
      mvecmat = rbind(mvecmat, mvec)
    }
    
    # Data augmentation
    for (d in edge.trim) {
      vu = MultiplyXB(X[[d]], b.old)
      lambda = lambda_cpp(p.d[d,], vu)
	  for (i in node[-senders[d]]) {
        for (j in sample(node[-i], A-1)) {
          probij = u_Gibbs(u[[d]][i, ], lambda[i,], delta, j)
          u[[d]][i, j] = multinom_vec(1, expconst(probij))-1
        }
      }
      u[[d]][senders[d],] = tabulateC(as.numeric(unlist(edge[[d]][2])), A)
    }
    # Z update
    table.W = lapply(1:K, function(k) {tabulateC(textlist.raw[which(unlist(z[edge.trim]) == k)], V)})
    for (d in edge.trim) {
      textlist.d = textlist[[d]]
      edgetime.d = rep(NA, K)
      for (w in 1:length(z[[d]])) {
        zw.old = z[[d]][w]
        if (length(textlist.d) > 0) {
          table.W[[zw.old]][textlist.d[w]] = table.W[[zw.old]][textlist.d[w]]-1
          topicword.d = TopicWord(K, z[[d]][-w], textlist.d, table.W, alpha, mvec, beta, V)
        } else {
          topicword.d = matrix(0, nrow = length(z[[d]]), ncol = K)
        }
        for (IP in unique(l)) {
          lK = which(l == IP)
          z[[d]][w] = min(lK)
          p.d[d, ] = pdmat(list(z[[d]]), l, nIP)           
          history.t = History(edge, p.d, node, timestamps[hist.d[d]-1]+exp(-745), timeunit)
          Xnew = Netstats_cpp(history.t, node, netstat)
          mu[d, ] = mu_vec(p.d[d,], xi[[d]])
          edgetime.d[lK] = Edgepartsum(Xnew, p.d[hist.d[d], ], b.old, u[[hist.d[d]]], delta)+
          Timepart(mu[d,], sigma_tau, senders[d], timeinc[d])
        }
        const.Z = edgetime.d+topicword.d[w, ]
        zw.new = multinom_vec(1, expconst(const.Z))
        if (zw.new != zw.old) {
          z[[d]][w] = zw.new
        } else {
          z[[d]][w] = zw.old
        }
        if (length(textlist.d) > 0) {
           	table.W[[z[[d]][w]]][textlist.d[w]] = table.W[[z[[d]][w]]][textlist.d[w]]+1
        }
      }
    }
    
    # C update 
    for (k in sort(unique(unlist(z)))) {
    	const.C = rep(NA, nIP)
      	for (IP in 1:nIP) {
        l[k] = IP
        p.dnew = pdmat(z, l, nIP)
        history.t = History(edge, p.dnew, node, timestamps[max.edge-1]+exp(-745), timeunit)
        Xnew = Netstats_cpp(history.t, node, netstat)
        Edgepartsum = Edgepartsum(Xnew, p.dnew[max.edge, ], b.old, u[[max.edge]], delta)
        mu = mu_mat(p.dnew, xi, edge.trim)
        Timepartsum = Timepartsum(mu, sigma_tau, senders, timeinc, edge.trim)
        prob = Edgepartsum+Timepartsum
        const.C[IP] = prob
      }
      l[k] = multinom_vec(1, expconst(const.C))
    }
 	  p.d = pdmat(z, l, nIP)  
 	  mu = mu_mat(p.d, xi, edge.trim)
 	  Timepartsum = Timepartsum(mu, sigma_tau, senders, timeinc, edge.trim)
 	  for (d in edge.trim) {
        history.t = History(edge, p.d, node, timestamps[d-1]+exp(-745), timeunit)
       	X[[d]] = Netstats_cpp(history.t, node, netstat)
      }
	     
    # adaptive M-H
    if (o > 1) {
      accept.rates[1] = accept.rates[1]/Inner[1]
      accept.rates[2] = accept.rates[2]/Inner[2]
      accept.rates[3] = accept.rates[3]/Inner[3]
      accept.rates[4] = accept.rates[1]
      sigma.Q = adaptive.MH(sigma.Q, accept.rates, update.size = 0.2*sigma.Q)
    }
    accept.rates = rep(0, 4)

    prior.old1 = priorsum(prior.b[[2]], prior.b[[1]], b.old)+
        dnorm(delta, prior.delta[1], sqrt(prior.delta[2]), TRUE)
    post.old1 = Edgepartsum(X[[max.edge]], p.d[max.edge, ], b.old, u[[max.edge]], delta)
    b.new = matrix(NA, nIP, P)
    for (inner in 1:Inner[1]) {
      for (IP in 1:nIP) {
 			 b.new[IP, ] = rmvnorm_arma(1, b.old[IP,], sigma.Q[1]*proposal.var1)
 	  }
      delta.new = rnorm(1, delta, sqrt(sigma.Q[4]))
      prior.new1 = priorsum(prior.b[[2]], prior.b[[1]], b.new)+
          dnorm(delta.new, prior.delta[1], sqrt(prior.delta[2]), TRUE)
      post.new1 = Edgepartsum(X[[max.edge]], p.d[max.edge, ], b.new, u[[max.edge]], delta.new)
      loglike.diff = prior.new1+post.new1-prior.old1-post.old1
      if (log(runif(1, 0, 1)) < loglike.diff) {
         b.old = b.new
         delta = delta.new
         prior.old1 = prior.new1
         post.old1 = post.new1
         accept.rates[1] = accept.rates[1]+1
       }
       if (inner > burn[1] & inner %% (thin[1]) == 0) {
         for (IP in 1:nIP) {
           bmat[[IP]][ ,(inner-burn[1])/thin[1]] = b.old[IP,]
         }
         deltamat[(inner-burn[1])/thin[1]] = delta
       }
     }

     prior.old2 = priorsum(prior.eta[[2]], prior.eta[[1]], eta.old)
     post.old2 = Timepartsum
     eta.new = matrix(NA, nIP, Q)
     for (inner in 1:Inner[2]) {
       for (IP in 1:nIP) {
 			  eta.new[IP, ] = rmvnorm_arma(1, eta.old[IP,], sigma.Q[2]*proposal.var2)
 		  }
      xi = xi_all(timemat, eta.new[,node], eta.new[,-node], edge.trim)
      mu = mu_mat(p.d, xi, edge.trim)
      Timepartsum = Timepartsum(mu, sigma_tau, senders, timeinc, edge.trim)
      prior.new2 = priorsum(prior.eta[[2]], prior.eta[[1]], eta.new)
      post.new2 = Timepartsum
      loglike.diff = prior.new2+post.new2-prior.old2-post.old2
      if (log(runif(1, 0, 1)) < loglike.diff) {
         eta.old = eta.new
         prior.old2 = prior.new2
         post.old2 = post.new2
         accept.rates[2] = accept.rates[2]+1
       }
       if (inner > burn[2] & inner %% (thin[2]) == 0) {
         for (IP in 1:nIP) {
           etamat[[IP]][ ,(inner-burn[2])/thin[2]] = eta.old[IP,]
         }
       }
     }
	  xi = xi_all(timemat, eta.old[,node], eta.old[,-node], edge.trim)
	  mu = mu_mat(p.d, xi, edge.trim)
    Timepartsum = Timepartsum(mu, sigma_tau, senders, timeinc, edge.trim)

    prior.old3 = dhalfcauchy(sigma_tau, prior.tau, TRUE)
    post.old3 = Timepartsum
    for (inner in 1:Inner[3]) {
      sigma_tau.new = rtruncnorm(1, 0, Inf, sigma_tau, sqrt(sigma.Q[3]))
      prior.new3 = dhalfcauchy(sigma_tau.new, prior.tau, TRUE)
      post.new3 = Timepartsum(mu, sigma_tau.new, senders, timeinc, edge.trim)
      loglike.diff = log(dtruncnorm(sigma_tau, 0, Inf, sigma_tau.new, sqrt(sigma.Q[3])))-
                     log(dtruncnorm(sigma_tau.new, 0, Inf, sigma_tau, sqrt(sigma.Q[3])))+
                     prior.new3+post.new3-prior.old3-post.old3
      if (log(runif(1, 0, 1)) < loglike.diff) {
        sigma_tau = sigma_tau.new
        prior.old3 = prior.new3
        post.old3 = post.new3
        accept.rates[3] = accept.rates[3]+1
      }
      if (inner > burn[3] & inner %% (thin[3]) == 0) {
        sigma_taumat[(inner-burn[3])/thin[3]] = sigma_tau
      }
    }
  }
    
  chain.final = list(l = l, z = z, b = bmat, eta = etamat, delta = deltamat, sigma_tau = sigma_taumat,
                    u = u, sigma.Q = sigma.Q, edge.trim = edge.trim)
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
#' @param alpha Dirichlet concentration prior for document-topic distribution
#' @param mvec Dirichlet base prior for document-topic distribution
#' @param beta Dirichlet concentration prior for topic-word distribution
#' @param b coefficients for recipients
#' @param eta coefficients for timestamps
#' @param delta tuning parameter for the number of recipients
#' @param sigma_tau variance parameter for the timestamps
#' @param l topic-interaction pattern assignment
#' @param support support of latent recipients
#' @param netstat which type of network statistics to use ("dyadic", "triadic", "degree")
#' @param timestat additional statistics to be used for timestamps other than netstat ("sender", "receiver", "timeofday", "dayofweek")
#' @param base.edge edges before 384 hours that is used to calculate initial history of interactions
#' @param base.text texts corresponding to base.edge
#' @param topic_token_assignments matrix of topic-token assignments
#' @param backward Logigal indicating whether we are generating backward samples (if FALSE -> forward)
#' @param base Logical indicating whether or not we are generating base edges (< 384)
#' @param timeunit hour (= 3600) or day (=3600*24) and so on
#' @param tz timezone such as EST and PST
#'
#' @return generated data including (author, recipients, timestamp, words)
#'
#' @export
GenerateDocs = function(nDocs, node, vocab, nIP, K, n.d, alpha, mvec, beta,
                        b, eta, delta, sigma_tau, l, support, netstat, timestat,
                        base.edge = NULL, base.text = NULL, topic_token_assignments = NULL,
                        backward = FALSE, base = FALSE, timeunit = 3600, tz = "America/New_York") { 
  A = length(node)
  V = length(vocab)
  phi = matrix(NA, K, V)
  for (k in 1:K) {
    phi[k,] = rdirichlet_cpp(1, rep(beta/V, V))
  }  
  netstat = as.numeric(c("degree", "dyadic", "triadic" ) %in% netstat)
  timestat = as.numeric(c("dayofweek","timeofday") %in% timestat)
  L = 3
  P = L*(2*netstat[1]+2*netstat[2]+4*netstat[3])
  t.d = ifelse(base, 0, base.edge[[1]][[3]]+384*timeunit)
  edge = base.edge
  text = base.text
  base.length = length(edge)
  p.d = matrix(NA, nrow = nDocs+base.length, ncol = nIP)
  timemat = matrix(0, nrow = nDocs+base.length, ncol = sum(timestat))

  if (!base) {
    for (d in 1:base.length) {
  	  p.d[d, ] = pdmat(list(as.integer(names(text[[d]]))), l, nIP)
    }
    if (sum(timestat) > 0) {
      Sys.setenv(TZ = tz)
      time_ymd = as.POSIXct(vapply(base.edge, function(d) {d[[3]]}, c(1)), tz = getOption("tz"), origin = "1970-01-01")
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
  history.t = lapply(1:nIP, function(IP) {
        	   lapply(1:3, function(x) matrix(0, A, A))})  
  VKmat = matrix(0, V, K)
  u = list()
  timestamps = matrix(0, nrow = nDocs, ncol = A)
  for (d in 1:nDocs) {
    u[[base.length+d]] = matrix(0, A, A)
    text[[base.length+d]] = rep(NA, n.d)
    if (!backward) {
      theta.d = rdirichlet_cpp(1, alpha*mvec)	
      topic.d = multinom_vec(max(1, n.d), theta.d)
      for (n in 1:n.d){
        text[[base.length+d]][n] = multinom_vec(1, phi[topic.d[n],])
      }
      names(text[[base.length+d]]) = topic.d
    } else {
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
    p.d[base.length+d, ] = pdmat(list(as.integer(names(text[[base.length+d]]))), l, nIP)
    if (t.d >= 384*timeunit) {
      history.t = History(edge, p.d, node, t.d+exp(-745), timeunit)
    }
    X = Netstats_cpp(history.t, node, netstat)
    vu = MultiplyXB(X, b)   
    lambda = lambda_cpp(p.d[base.length+d,], vu)
    for (i in node) {
      sendervec = rep(0, A)
      u[[base.length+d]][i,-i] = r.gibbs.measure(1, lambda[i,-i], delta, support)
      sendervec[i] = 1
      Y = c(sendervec, timemat[base.length+d-1,])
      xi = MultiplyYeta(Y, eta)
      mu = mu_cpp(p.d[base.length+d,], xi)
      time = rlnorm(1, mu, sigma_tau)
      timestamps[d,i] = time*timeunit
    }
    i.d = which(timestamps[d,] == min(timestamps[d,]))
    j.d = which(u[[base.length+d]][i.d,] == 1)
    t.d = t.d+timestamps[d,i.d]
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
  }
  return(list(edge = edge, text = text, base = base.length, u = u, z = lapply(text, function(d){as.numeric(names(d))}), 
  		b = b, eta = eta, delta = delta, sigma_tau = sigma_tau, l = l))							
} 


#' @title GenerateDocs.PPC
#' @description Generate a collection of documents according to the generative process of IPTM using Gibbs measure
#'
#' @param nDocs number of documents to be generated
#' @param node vector of node id's (ID starting from 1)
#' @param vocab all vocabularies used over the corpus
#' @param nIP total number of interaction patterns
#' @param K total number of topics
#' @param alpha Dirichlet concentration prior for document-topic distribution
#' @param mvec Dirichlet base prior for document-topic distribution
#' @param beta Dirichlet concentration prior for topic-word distribution
#' @param b coefficients for recipients
#' @param eta coefficients for timestamps
#' @param delta tuning parameter for the number of recipients
#' @param sigma_tau variance parameter for the timestamps
#' @param l topic-interaction pattern assignment
#' @param u inferred latent receivers
#' @param netstat which type of network statistics to use ("dyadic", "triadic", "degree")
#' @param timestat additional statistics to be used for timestamps other than netstat ("sender", "receiver","timeofday", "dayofweek")
#' @param base.edge edges before 384 hours that is used to calculate initial history of interactions
#' @param base.text texts corresponding to base.edge
#' @param z list topic-token assignments
#' @param word_type_topic_counts word_type_topic_counts from inference
#' @param text.length number of words in each document
#' @param timeunit hour (= 3600) or day (=3600*24) and so on
#' @param tz timezone such as EST and PST
#'
#' @return generated edge and text
#'
#' @export
GenerateDocs.PPC = function(nDocs, node, vocab, nIP, K, alpha, mvec, beta, b, eta, delta, sigma_tau,
 						    l, u, netstat, timestat, base.edge, base.text, z, word_type_topic_counts, text.length, timeunit = 3600, tz = "America/New_York") {
  A = length(node)
  V = length(vocab)
  netstat = as.numeric(c("degree", "dyadic", "triadic" ) %in% netstat)
  timestat = as.numeric(c("dayofweek","timeofday") %in% timestat)
  L = 3
  P = L*(2*netstat[1]+2*netstat[2]+4*netstat[3])
  base.length = length(base.edge)
  t.d = base.edge[[base.length]][[3]]
  edge = base.edge
  text = base.text
  p.d = pdmat(z, l, nIP)  
  timemat = matrix(0, nrow = nDocs+base.length, ncol = sum(timestat))
  if (sum(timestat) > 0) {
      Sys.setenv(TZ = tz)
      time_ymd = as.POSIXct(vapply(edge, function(d) {d[[3]]}, c(1)), tz = getOption("tz"), origin = "1970-01-01")
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
  
  for (d in 1:nDocs) {
    N.d = text.length[base.length+d]
    text[[base.length+d]] = rep(NA, N.d)
    phi.k = rep(NA, V)
    topic.d = z[[base.length+d]]
    if (N.d > 0) {
    		for (n in 1:N.d){
          for (w in 1:V) {
            phi.k[w] = (word_type_topic_counts[w, topic.d[n]]+beta/V) / (sum(word_type_topic_counts[, topic.d[n]])+beta)
          }
          text[[base.length+d]][n] = multinom_vec(1, phi.k)
          word_type_topic_counts[text[[base.length+d]][n], topic.d[n]] = word_type_topic_counts[text[[base.length + d]][n], topic.d[n]]+1
   		}
   		names(text[[base.length+d]]) = topic.d
   	} 	
    history.t = History(edge, p.d, node, t.d+exp(-745), timeunit)
    X = Netstats_cpp(history.t, node, netstat)
    vu = MultiplyXB(X, b)     
    lambda = lambda_cpp(p.d[base.length+d,], vu)
    	for (i in node) {
        for (j in sample(node[-i], A-1)) {
          probij = u_Gibbs(u[[base.length+d]][i, ], lambda[i,], delta, j)
          u[[base.length+d]][i, j] = multinom_vec(1, expconst(probij))-1
        }
    	}
    xi = ximat(timemat[base.length+d-1,], eta[,node], eta[,-node])
    mu = mu_vec(p.d[base.length+d,], xi)
    timestamps = rlnorm(1, mu, sigma_tau)*timeunit
    i.d = which(timestamps == min(timestamps))
    j.d = which(u[[base.length+d]][i.d,] == 1)
    t.d = t.d+timestamps[i.d]
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
  return(list(edge = edge, text = text, u = u, timemat = timemat))							
}


#' @title IPTM.PPC
#' @description posterior predictive experiments
#'
#' @param Out number of outer iterations of inference from which to generate predictions
#' @param edge list of document information with 3 elements (element 1 sender, element 2 receiver, element 3 time in unix.time format)
#' @param node nodelist containing the ID of nodes (ID starting from 1)
#' @param textlist list of text (length=number of documents in total) containing the words in each document
#' @param vocab all vocabularies used over the corpus
#' @param nIP total number of interaction patterns
#' @param K total number of topics
#' @param netstat which type of network statistics to use ("dyadic", "triadic", "degree")
#' @param timestat additional statistics to be used for timestamps other than netstat ("timeofday", "dayofweek")
#' @param inference inference result
#' @param timeunit hour (= 3600) or day (=3600*24) and so on
#' @param tz timezone such as EST and PST
#'
#' @return prediction output
#'
#' @export
IPTM.PPC = function(Out, edge, node, textlist, vocab, nIP, K, netstat, timestat, inference, timeunit = 3600, tz = "America/New_York") {
    New_sample = list()
    word_type_topic_counts = matrix(0, length(vocab), K)
    textlist.raw = unlist(textlist)
    text.length = sapply(textlist, function(d){length(d)})
    emptytext = which(text.length==0)
    z = inference$z
    edge.trim = inference$edge.trim
    alpha = inference$alpha[length(inference$alpha)]
    mvec = inference$mvec[nrow(inference$mvec),]
    b = tvapply(1:nIP, function(IP) {inference$b[[IP]][,ncol(inference$b[[IP]])]}, rep(0, nrow(inference$b[[1]])))
    eta = tvapply(1:nIP, function(IP) {inference$eta[[IP]][,ncol(inference$eta[[IP]])]}, rep(0, nrow(inference$eta[[1]])))
    delta = inference$delta[length(inference$delta)]
    sigma_tau = inference$sigma_tau[length(inference$sigma_tau)]
    l = inference$l
 	  u = inference$u
    for (k in 1:K) {
      word_type_topic_counts[,k] = tabulate(textlist.raw[which(unlist(z[-emptytext])==k)], length(vocab))
    }
    base.edge = edge[-edge.trim]
    base.text = textlist[-edge.trim]
    for (o in 1:Out) {
        print(o)
        New_sample[[o]] = GenerateDocs.PPC(length(edge.trim), node, vocab, nIP, K, alpha, mvec, beta = 2, b, eta, delta, sigma_tau,
 						  l, u, netstat, timestat, base.edge, base.text, z, word_type_topic_counts, text.length, timeunit = 3600, tz = tz)
    }
    return(New_sample)
}


#' @title GiR_stats
#' @description Calculate several statistics from samples generated from forward or backward sampling
#'
#' @param GiR_sample one sample from generative process
#' @param V number of unique words
#' @param timeunit hour (= 3600) or day (=3600*24) and so on
#'
#' @return A vector of statistics calculated from one GiR sample
#'
#' @export

GiR_stats = function(GiR_sample, V, timeunit = 3600) {
  edge = GiR_sample$edge
  text = GiR_sample$text
  edge = edge[-(1:GiR_sample$base)]
  text = text[-(1:GiR_sample$base)]
  GiR_stats = c()
  D = length(edge)
  P = ncol(GiR_sample$b)
  Q = ncol(GiR_sample$eta)
  K = length(GiR_sample$l)
  nIP = nrow(GiR_sample$b)
  n.d = length(text[[1]])
  
  GiR_stats[1:P] = colMeans(GiR_sample$b)
  GiR_stats[(P+1):(P+Q)] = colMeans(GiR_sample$eta) 
  GiR_stats[P+Q+1] = GiR_sample$delta
  GiR_stats[P+Q+2] = GiR_sample$sigma_tau
  GiR_stats[P+Q+3] = mean(vapply(1:D, function(d) {length(edge[[d]][[2]])}, c(1)))
  GiR_stats[P+Q+4] = mean(vapply(2:D, function(d) {edge[[d]][[3]]-edge[[d-1]][[3]]}, c(1))/timeunit) 			
  GiR_stats[P+Q+5] = mean(GiR_sample$l)
  Tokens_in_Topic = tabulate(vapply(1:D, function(d){as.numeric(names(text[[d]]))}, rep(0, n.d)), K)
  GiR_stats[(P+Q+6):(P+Q+5+nIP)] = vapply(1:nIP, function(IP) {Tokens_in_Topic %*% (GiR_sample$l == IP)}, c(1))
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
#' @param alpha Dirichlet concentration prior for document-topic distribution
#' @param mvec Dirichlet base prior for document-topic distribution
#' @param beta Dirichlet concentration prior for topic-word distribution
#' @param prior.b prior mean and covariance of b in multivariate Normal distribution
#' @param prior.delta prior mean and variance of delta in Normal distribution
#' @param prior.eta prior mean and covariance of eta in multivariate Normal distribution
#' @param prior.tau prior shape and scale parameter of sigma_tau in inverse-Gamma distribution
#' @param sigma.Q proposal distribution variance parameter
#' @param Outer size of outer iterations 
#' @param Inner size of inner iteration for Metropolis-Hastings updates
#' @param burn iterations to be discarded at the beginning of Metropolis-Hastings chains
#' @param thin the thinning interval of Metropolis-Hastings chains
#' @param netstat which type of network statistics to use ("dyadic", "triadic", "degree")
#' @param timestat additional statistics to be used for timestamps other than netstat ("sender", "receiver", "timeofday", "dayofweek")
#' @param base.edge artificial collection of documents to be used as initial state of history
#' @param base.text artificial collection of documents to be used as initial state of history
#' @param generate_PP_plots Logical indicating whether to draw PP plots
#'
#' @return Forward and Backward samples and corresponding test results
#'
#' @export
GiR = function(Nsamp, nDocs, node, vocab, nIP, K, n.d, alpha, mvec, beta, 
              prior.b, prior.delta, prior.eta, prior.tau, sigma.Q, Outer, Inner, burn, thin,
              netstat = c("dyadic"), timestat = c("timeofday", "dayofweek"),
              base.edge, base.text, generate_PP_plots = TRUE) {
  
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
                              "Mean_recipients", "Mean_timediff", "Mean_TopicIP", paste0("Tokens_in_IP_", 1:nIP), 
                              paste0("Tokens_in_Topic", 1:K), paste0("Tokens_in_Word",1:V))
  for (i in 1:Nsamp) { 
    if (i %% 5000 == 0) {cat("Forward sampling", i, "\n")}
    b = rmvnorm_arma(nIP, prior.b[[1]], prior.b[[2]])
    eta = rmvnorm_arma(nIP, prior.eta[[1]], prior.eta[[2]])
    delta = rnorm(1, prior.delta[1], sqrt(prior.delta[2]))
    sigma_tau = rhalfcauchy(1, prior.tau)
    l = sample(1:nIP, K, replace = TRUE)
    Forward_sample = GenerateDocs(nDocs, node, vocab, nIP, K, n.d, alpha, mvec, beta, b, eta, delta, sigma_tau, 
                    l, support, netstat, timestat, base.edge = base.edge, base.text = base.text, 
                    topic_token_assignments = NULL, backward = FALSE, base = FALSE)
    Forward_stats[i, ] = GiR_stats(Forward_sample, V)
  }
  #Backward sampling
  Backward_stats = matrix(NA, nrow = Nsamp, ncol = ncol(Forward_stats))
  Backward_sample = GenerateDocs(nDocs, node, vocab, nIP, K, n.d, alpha, mvec, beta, b, eta, delta, sigma_tau, 
                    l, support, netstat, timestat, base.edge = base.edge, base.text = base.text, 
                    topic_token_assignments = NULL, backward = FALSE, base = FALSE)
  for (i in 1:Nsamp) { 
    if (i %% 1 == 0) {cat("Backward Sampling", i, "\n")}
    Inference_samp = IPTM.inference.GiR(Backward_sample$edge, node, Backward_sample$text, vocab, nIP, K, sigma.Q, 
                     alpha, mvec, beta, prior.b, prior.delta,prior.eta, prior.tau, Outer, Inner, burn, thin, 
                     netstat, timestat, optimize = FALSE, initial = NULL)
    b = tvapply(1:nIP, function(IP) {Inference_samp$b[[IP]][,ncol(Inference_samp$b[[IP]])]}, rep(0, P))
    eta = tvapply(1:nIP, function(IP) {Inference_samp$eta[[IP]][,ncol(Inference_samp$eta[[IP]])]}, rep(0, Q))
    delta = Inference_samp$delta[length(Inference_samp$delta)]
    sigma_tau = Inference_samp$sigma_tau[length(Inference_samp$sigma_tau)]
    l = Inference_samp$l
    z = Inference_samp$z
    for (d in 1:length(z)) {
      names(z[[d]]) = Forward_sample$text[[d]]
    }
    Backward_sample = GenerateDocs(nDocs, node, vocab, nIP, K, n.d, alpha, mvec, beta, b, eta, delta, sigma_tau, 
                      l, support, netstat, timestat, base.edge = base.edge, base.text = base.text, 
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
#' @param alpha Dirichlet concentration prior for document-topic distribution
#' @param mvec Dirichlet base prior for document-topic distribution
#' @param beta Dirichlet concentration prior for topic-word distribution
#' @param prior.b prior mean and covariance of b in multivariate Normal distribution
#' @param prior.delta prior mean and variance of delta in Normal distribution
#' @param prior.eta prior mean and covariance of eta in multivariate Normal distribution
#' @param prior.tau prior shape and scale parameter of sigma_tau in inverse-Gamma distribution
#' @param sigma.Q proposal distribution variance parameter
#' @param Outer size of outer iterations 
#' @param Inner size of inner iteration for Metropolis-Hastings updates
#' @param burn iterations to be discarded at the beginning of Metropolis-Hastings chains
#' @param thin the thinning interval of Metropolis-Hastings chains
#' @param netstat which type of network statistics to use ("dyadic", "triadic", "degree")
#' @param timestat additional statistics to be used for timestamps other than netstat ("sender", "receiver", "timeofday", "dayofweek")
#' @param base.edge artificial collection of documents to be used as initial state of history
#' @param base.text artificial collection of documents to be used as initial state of history
#' @param generate_PP_plots Logical indicating whether to draw PP plots
#'
#' @return Forward and Backward samples and corresponding test results
#'
#' @export
Schein = function(Nsamp, nDocs, node, vocab, nIP, K, n.d, alpha, mvec, beta, 
              prior.b, prior.delta, prior.eta, prior.tau, sigma.Q, Outer, Inner, burn, thin,
              netstat = c("dyadic"), timestat = c("timeofday", "dayofweek"),
              base.edge, base.text, generate_PP_plots = TRUE) {
  
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
                            "Mean_recipients", "Mean_timediff", "Mean_TopicIP", paste0("Tokens_in_IP_", 1:nIP), 
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
    l = sample(1:nIP, K, replace = TRUE)
    Forward_sample = GenerateDocs(nDocs, node, vocab, nIP, K, n.d, alpha, mvec, beta, b, eta, delta, sigma_tau, 
                     l, support, netstat, timestat, base.edge = base.edge, base.text = base.text, 
                     topic_token_assignments = NULL, backward = FALSE, base = FALSE)
    Forward_stats[i, ] = GiR_stats(Forward_sample, V)
  	initial = list(alpha = alpha, mvec = mvec, delta = delta, sigma_tau = sigma_tau, b = b, eta = eta, l = l, 
                   z = Forward_sample$z, proposal.var1 = diag(P), proposal.var2 = diag(Q), sigma.Q = sigma.Q, 
                   u = Forward_sample$u)
    Inference_samp = IPTM.inference.GiR(Forward_sample$edge, node, Forward_sample$text, vocab, nIP, K, sigma.Q, 
                     alpha, mvec, beta, prior.b, prior.delta,prior.eta, prior.tau, Outer, Inner, burn, thin, 
                     netstat, timestat, optimize = FALSE, initial = initial)
    b = tvapply(1:nIP, function(IP) {Inference_samp$b[[IP]][,ncol(Inference_samp$b[[IP]])]}, rep(0, P))
    eta = tvapply(1:nIP, function(IP) {Inference_samp$eta[[IP]][,ncol(Inference_samp$eta[[IP]])]}, rep(0, Q))
    delta = Inference_samp$delta[length(Inference_samp$delta)]
    sigma_tau = Inference_samp$sigma_tau[length(Inference_samp$sigma_tau)]
    l = Inference_samp$l
    z = Inference_samp$z
    for (d in 1:length(z)) {
      names(z[[d]]) = Forward_sample$text[[d]]
    }
    Backward_sample = GenerateDocs(nDocs, node, vocab, nIP, K, n.d, alpha, mvec, beta, b, eta, delta, sigma_tau, 
                      l, support, netstat, timestat, base.edge = base.edge, base.text = base.text, 
                      topic_token_assignments = z, backward = TRUE, base = FALSE)
    Backward_stats[i, ] = GiR_stats(Backward_sample, V)
 }
  if (generate_PP_plots) {
    par(mfrow=c(5,6), oma = c(3,3,3,3), mar = c(2,1,1,1))
    GiR_PP_Plots(Forward_stats, Backward_stats)
  }			
  return(list(Forward = Forward_stats, Backward = Backward_stats))
}                         	          
      
