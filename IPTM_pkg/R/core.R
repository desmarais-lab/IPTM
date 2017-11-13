#' @useDynLib IPTM
#' @import stats
#' @import grDevices
#' @import graphics
#' @import RcppArmadillo
#' @importFrom Rcpp sourceCpp
#' @importFrom reshape melt
#' @importFrom coda mcmc geweke.diag
#' @importFrom MCMCpack dinvgamma
#' @importFrom combinat permn
#' @importFrom mgcv uniquecombs
#' @importFrom psych geometric.mean
#' @importFrom lubridate wday hour
#' @importFrom FastGP rcpp_rmvnorm rcpp_log_dmvnorm

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
	gibbsNormalizer = prod(exp(delta+log(lambda.i))+1)-1
	if (gibbsNormalizer == 0) {gibbsNormalizer = exp(-745)}
	logitNumerator = vapply(1:nrow(support), function(s) {
					 	exp(sum((delta+log(lambda.i))*support[s,]))
					  }, c(1))					  
	samp <- sample(1:nrow(support), nsamp, replace = TRUE, prob = logitNumerator/gibbsNormalizer)	
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
#' @param prior.tau prior shape and scale parameter of sigma2_tau in inverse-Gamma distribution
#' @param Outer size of outer iterations 
#' @param Inner size of inner iteration for Metropolis-Hastings updates
#' @param burn iterations to be discarded at the beginning of Metropolis-Hastings chains
#' @param thin the thinning interval of Metropolis-Hastings chains
#' @param netstat which type of network statistics to use ("intercept", dyadic", "triadic", "degree")
#' @param timestat additional statistics to be used for timestamps other than netstat ("timeofday", "dayofweek")
#' @param optimize logical to optimize alpha (Dirichlet concentration prior for document-topic distribution)
#' @param initial list of initial values user wants to assign including (alpha, mvec, delta, b, eta, l, z, u, sigma2_tau, proposal.var1, proposal.var2)
#'
#' @return MCMC output containing all parameter estimates
#'
#' @export

IPTM.inference = function(edge, node, textlist, vocab, nIP, K, sigma.Q, alpha, mvec, beta, prior.b, prior.delta,
				 prior.eta, prior.tau, Outer, Inner, burn, thin, netstat, timestat, optimize = FALSE, initial = NULL) {
  
  # trim the edge so that we only model edges after 384 hours
  timestamps = vapply(edge, function(d) { d[[3]]/3600 }, c(1))
  edge.trim = which_int(384, timestamps):length(edge)
  max.edge = max(edge.trim)
  timeinc = c(timestamps[1], timestamps[-1]-timestamps[-length(timestamps)])
  timeinc[timeinc==0] = exp(-745)
  # initialization
  convergence = c()
  netstat = as.numeric(c("intercept", "degree", "dyadic", "triadic") %in% netstat)
  timestat = as.numeric(c("dayofweek","timeofday") %in% timestat)
  timemat = matrix(0, nrow = length(edge), ncol = sum(timestat))
  if (sum(timestat) > 0) {
      time_ymd = as.POSIXct(vapply(edge, function(d) {d[[3]]}, c(1)), tz = "America/New_York", origin="1970-01-01")
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
  P = netstat[1]+L*(2*netstat[2]+2*netstat[3]+4*netstat[4])
  Q = length(timestat)
  V = length(vocab)
  phi = lapply(1:K, function(k) {rdirichlet_cpp(1, rep(beta/V, V))})						 	
  if (length(initial) == 0) {
  theta = rdirichlet_cpp(length(edge), alpha*mvec)
  delta = rnorm(1, prior.delta[1], sqrt(prior.delta[2]))
  sigma2_tau = 1/rgamma(1, prior.tau[1], prior.tau[2])
  b.old = lapply(1:nIP, function(IP) {c(rcpp_rmvnorm(1, prior.b[[2]], prior.b[[1]]))}) 
  eta.old = lapply(1:nIP, function(IP) {c(rcpp_rmvnorm(1, prior.eta[[2]], prior.eta[[1]]))})
  l = sample(1:nIP, K, replace = TRUE) 
  z = lapply(seq(along = edge), function(d) {
  	multinom_vec(max(1, length(textlist[[d]])), theta[d, ])})
  p.d = pdmat(z, l, nIP) 
  proposal.var1 = lapply(1:nIP, function(IP){diag(P)})
  proposal.var2 = lapply(1:nIP, function(IP){diag(P+Q)})
  sigma.Q = sigma.Q
  u = lapply(seq(along = edge), function(d) {
    matrix(0, nrow = length(node), ncol = length(node))
  })
  for (d in edge.trim) {
    u[[d]] = matrix(rbinom(length(node)^2, 1, 1/length(node)), nrow =length(node), ncol = length(node))
    diag(u[[d]]) = 0
  } 
  } else {
    theta = rdirichlet_cpp(length(edge), initial$alpha*initial$mvec)
    delta = initial$delta
    sigma2_tau = initial$sigma2_tau
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
    bmat[[IP]] = matrix(b.old[[IP]], nrow = P, ncol = (Inner[1]-burn[1])/thin[1])
    etamat[[IP]] = matrix(eta.old[[IP]], nrow = P+Q, ncol = (Inner[2]-burn[2])/thin[2])
  }
  deltamat = rep(delta, (Inner[1]-burn[1])/thin[1])
  sigma2_taumat = rep(sigma2_tau, (Inner[2]-burn[2])/thin[2])
						 
  lambda = list()
  mu = matrix(0, nrow = length(edge), ncol = length(node))

  textlist.raw = unlist(textlist)
  alphavec = c()
  mvecmat = matrix(NA, nrow = 0, ncol = K)
  accept.rates = rep(0, 2)

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
      	history.t = History(edge, p.d, node, timestamps[d-1]+exp(-745))
      	X = Netstats_cpp(history.t, node, netstat)
      	vu = MultiplyXBList(X, b.old)
     	lambda[[d]] = lambda_cpp(p.d[d,], vu)
     	for (i in node[-edge[[d]][[1]]]) {
     		for (j in sample(node[-i], length(node)-1)) {
          		probij = u_Gibbs(u[[d]][i, ], lambda[[d]][i,], delta, j)
          		u[[d]][i, j] = multinom_vec(1, probij)-1
        		}
      	}
      	u[[d]][edge[[d]][[1]],] = tabulateC(as.numeric(unlist(edge[[d]][2])), length(node))
      }	 
       
    # Z update	
      table.W = lapply(1:K, function(k) {tabulateC(textlist.raw[which(unlist(z) == k)], V)})
      for (d in 1:(edge.trim[1]-1)) {
	   	textlist.d = textlist[[d]]
	   	for (w in 1:length(z[[d]])) {
	   		zw.old = z[[d]][w]
          	if (length(textlist.d) > 0) {
          		table.W[[zw.old]][textlist.d[w]] = table.W[[zw.old]][textlist.d[w]]-1
        			topicpart.d = TopicInEqZ(K, z[[d]][-w], alpha, mvec)
       		 	wordpart.d = WordInEqZ(K, textlist.d, table.W, beta, V)
        		 } else {
        			topicpart.d = 0
        			wordpart.d = matrix(0, nrow = length(z[[d]]), ncol = K)
        		 }
          		 const.Z = topicpart.d+wordpart.d[w, ]
          		 zw.new = multinom_vec(1, expconst(const.Z))
          		 if (zw.new != zw.old) {
                     z[[d]][w] = zw.new
           			 table.W[[zw.new]][textlist.d[w]] = table.W[[zw.new]][textlist.d[w]]+1
      				 p.d[d, ] = pdmat(list(z[[d]]), l, nIP) 	
               	 } else {
               	     table.W[[zw.old]][textlist.d[w]] = table.W[[zw.old]][textlist.d[w]]+1 	
               	 }
        	 	}
     }
     for (d in edge.trim) {
       textlist.d = textlist[[d]]
       if (timestamps[d]+384 > timestamps[max.edge]) {
      	 hist.d = max.edge
       } else {
      	 hist.d = which_num(timestamps[d]+384, timestamps)
       }
       edgetime.d = rep(NA, K)
       for (w in 1:length(z[[d]])) {
       	 zw.old = z[[d]][w]
       	 if (length(textlist.d) > 0) {
       	 table.W[[zw.old]][textlist.d[w]] = table.W[[zw.old]][textlist.d[w]]-1
         topicpart.d = TopicInEqZ(K, z[[d]][-w], alpha, mvec)
         wordpart.d = WordInEqZ(K, textlist.d, table.W, beta, V)
       } else {
         topicpart.d = 0
         wordpart.d = matrix(0, nrow = length(z[[d]]), ncol = K)
       }
	  for (IP in unique(l)) {
	  	lK = which(l == IP)
       	z[[d]][w] = min(lK)
       	p.d[d, ] = pdmat(list(z[[d]]), l, nIP)           
        history.t = History(edge, p.d, node, timestamps[hist.d-1]+exp(-745))
    	   	X = Netstats_cpp(history.t, node, netstat)
    	    XB = MultiplyXBList(X, b.old)
    	    lambda[[hist.d]] = lambda_cpp(p.d[hist.d,], XB)
    	    for (i in node) {
            X_u = lapply(X[[i]], function(X_IP) {X_IP[which(u[[hist.d]][i,]==1),]})
            if (length(X_u[[1]]) > P) {
                X_u = lapply(X_u, function(X_u_IP) {geometric.mean(X_u_IP)})
            }
            Y = lapply(X_u, function(X_u_IP) {c(X_u_IP, timemat[hist.d-1,])})
            xi = MultiplyYeta(Y, eta.old)
            mu[hist.d, i] = mu_cpp(p.d[hist.d,], xi)
      	}
        	edgetime.d[lK] = Edgepart(u[[hist.d]], lambda[[hist.d]], delta)+
        					 Timepart(mu[hist.d,], sigma2_tau, edge[[hist.d]][[1]], timeinc[hist.d])
        }
        const.Z = edgetime.d+topicpart.d+wordpart.d[w, ]
        zw.new = multinom_vec(1, expconst(const.Z))
        if (zw.new != zw.old) {
           z[[d]][w] = zw.new
           table.W[[zw.new]][textlist.d[w]] = table.W[[zw.new]][textlist.d[w]]+1
        } else {
         	z[[d]][w] = zw.old
            table.W[[zw.old]][textlist.d[w]] = table.W[[zw.old]][textlist.d[w]]+1
        }
      }
     }
     
    # C update 
    for (k in sort(unique(unlist(z[edge.trim])))) { 
    		const.C = rep(NA, nIP)
        for (IP in 1:nIP) {
          	l[k] = IP
          	p.d = pdmat(z, l, nIP) 
 		 	history.t = History(edge, p.d, node, timestamps[max.edge-1]+exp(-745))
    	   		X = Netstats_cpp(history.t, node, netstat)
    	   		XB = MultiplyXBList(X, b.old)    
           	lambda[[max.edge]] = lambda_cpp(p.d[max.edge,], XB)
		    for (i in node) {
            	X_u = lapply(X[[i]], function(X_IP) {X_IP[which(u[[max.edge]][i,]==1),]})
            	if (length(X_u[[1]]) > P) {
                	X_u = lapply(X_u, function(X_u_IP) {geometric.mean(X_u_IP)})
            	}
            	Y = lapply(X_u, function(X_u_IP) {c(X_u_IP, timemat[max.edge-1,])})
            	xi = MultiplyYeta(Y, eta.old)
            	mu[max.edge, i] = mu_cpp(p.d[max.edge,], xi)
      		}
           	prob = Edgepart(u[[max.edge]], lambda[[max.edge]], delta)+
        		   Timepart(mu[max.edge,], sigma2_tau, edge[[max.edge]][[1]], timeinc[max.edge])
           	const.C[IP] = prob
      	}
        l[k] = multinom_vec(1, expconst(const.C))
	}    
 	p.d = pdmat(z, l, nIP)  
    for (d in max.edge) {
        	history.t = History(edge, p.d, node, timestamps[d-1]+exp(-745))
    	    X = Netstats_cpp(history.t, node, netstat)
    	    XB = MultiplyXBList(X, b.old)   
    	    lambda[[d]] = lambda_cpp(p.d[d,], XB)
    	     for (i in node) {
            	X_u = lapply(X[[i]], function(X_IP) {X_IP[which(u[[d]][i,]==1),]})
            	if (length(X_u[[1]]) > P) {
                	X_u = lapply(X_u, function(X_u_IP) {geometric.mean(X_u_IP)})
            	}
            	Y = lapply(X_u, function(X_u_IP) {c(X_u_IP, timemat[d-1,])})
            	xi = MultiplyYeta(Y, eta.old)
            	mu[d, i] = mu_cpp(p.d[d,], xi)
      	}
	}
	# adaptive M-H   
    if (o > 1) {
    	accept.rates[1] = accept.rates[1]/Inner[1]
    	accept.rates[2] = accept.rates[2]/Inner[2]
    sigma.Q = adaptive.MH(sigma.Q, accept.rates, update.size = 0.2*sigma.Q)
    #  if (accept.rates[1] > 1/Inner[1]) { 
    #    for (IP in 1:nIP) {
    #      proposal.var1[[IP]] = var(uniquecombs(t(bmat[[IP]])))
    #      proposal.var1[[IP]] = proposal.var1[[IP]]/max(abs(proposal.var1[[IP]]))
    #      }
    #  }
    #  if (accept.rates[2] > 1/Inner[2]) { 
    #    for (IP in 1:nIP) {
    #      proposal.var2[[IP]] = var(uniquecombs(t(etamat[[IP]])))
    #      proposal.var2[[IP]] = proposal.var2[[IP]]/max(abs(proposal.var2[[IP]]))
    #    }
    #  }
    }
    accept.rates = rep(0, 2)
    
    prior.old1 = sum(vapply(1:nIP, function(IP) {rcpp_log_dmvnorm(prior.b[[2]], prior.b[[1]], b.old[[IP]], FALSE)}, c(1)))+
    				 dnorm(delta, prior.delta[1], sqrt(prior.delta[2]), TRUE)
    post.old1 = Edgepart(u[[max.edge]], lambda[[max.edge]], delta)   
    for (inner in 1:Inner[1]) {
      b.new = lapply(1:nIP, function(IP) {
        c(rcpp_rmvnorm(1, sigma.Q[1]*proposal.var1[[IP]], b.old[[IP]]))
      }) 
      delta.new = rnorm(1, delta, sqrt(sigma.Q[1]))
      for (d in max.edge) {
        XB = MultiplyXBList(X, b.new)
        lambda[[d]] = lambda_cpp(p.d[d,], XB)    
      }
    prior.new1 = sum(vapply(1:nIP, function(IP) {rcpp_log_dmvnorm(prior.b[[2]], prior.b[[1]], b.new[[IP]], FALSE)}, c(1)))+
    				 dnorm(delta.new, prior.delta[1], sqrt(prior.delta[2]), TRUE)
    post.new1 = Edgepart(u[[max.edge]], lambda[[max.edge]], delta.new)
    loglike.diff = prior.new1+post.new1-prior.old1-post.old1
      if (log(runif(1, 0, 1)) < loglike.diff) {
        for (IP in 1:nIP) {
          b.old[[IP]] = b.new[[IP]]
        }
        delta = delta.new
        prior.old1 = prior.new1
        post.old1 = post.new1
        accept.rates[1] = accept.rates[1]+1
      }
      if (inner > burn[1] & inner %% (thin[1]) == 0) {
        for (IP in 1:nIP) {
          bmat[[IP]][ ,(inner-burn[1])/thin[1]] = b.old[[IP]]
        }
        deltamat[(inner-burn[1])/thin[1]] = delta
      }
    }
	
	prior.old2 = sum(vapply(1:nIP, function(IP) {rcpp_log_dmvnorm(prior.eta[[2]], prior.eta[[1]], eta.old[[IP]], FALSE)}, c(1)))+
    				 log(dinvgamma(sigma2_tau, prior.tau[1], prior.tau[2]))
    post.old2 = Timepart(mu[max.edge,], sigma2_tau, edge[[max.edge]][[1]], timeinc[max.edge]) 
    for (inner in 1:Inner[2]) {
      eta.new = lapply(1:nIP, function(IP) {
        c(rcpp_rmvnorm(1, sigma.Q[2]*proposal.var2[[IP]], eta.old[[IP]]))
      }) 
      sigma2_tau.new = exp(rnorm(1, log(sigma2_tau), sqrt(sigma.Q[2])))
      for (d in max.edge) {
       for (i in node) {
            	X_u = lapply(X[[i]], function(X_IP) {X_IP[which(u[[d]][i,]==1),]})
            	if (length(X_u[[1]]) > P) {
                	X_u = lapply(X_u, function(X_u_IP) {geometric.mean(X_u_IP)})
            	}
            	Y = lapply(X_u, function(X_u_IP) {c(X_u_IP, timemat[d-1,])})
            	xi = MultiplyYeta(Y, eta.old)
            	mu[d, i] = mu_cpp(p.d[d,], xi)
      	}   
      }
    prior.new2 = sum(vapply(1:nIP, function(IP) {rcpp_log_dmvnorm(prior.eta[[2]], prior.eta[[1]], eta.new[[IP]], FALSE)}, c(1)))+
    			 log(dinvgamma(sigma2_tau.new, prior.tau[1], prior.tau[2]))
    post.new2 =Timepart(mu[max.edge,], sigma2_tau.new, edge[[max.edge]][[1]], timeinc[max.edge])
    loglike.diff = prior.new2+post.new2-prior.old2-post.old2
      if (log(runif(1, 0, 1)) < loglike.diff) {
        for (IP in 1:nIP) {
          eta.old[[IP]] = eta.new[[IP]]
        }
       	sigma2_tau = sigma2_tau.new
        prior.old2 = prior.new2
        post.old2 = post.new2
        accept.rates[2] = accept.rates[2]+1
      }
      if (inner > burn[2] & inner %% (thin[2]) == 0) {
        for (IP in 1:nIP) {
          etamat[[IP]][ ,(inner-burn[2])/thin[2]] = eta.old[[IP]]
        }
        sigma2_taumat[(inner-burn[2])/thin[2]] = sigma2_tau
      }
    }
   
    convergence[o] = Edgepart(u[[max.edge]], lambda[[max.edge]], delta)+
    					 Timepart(mu[max.edge,], sigma2_tau, edge[[max.edge]][[1]], timeinc[max.edge]) 
   }
   
   for (d in max.edge) {
        	history.t = History(edge, p.d, node, timestamps[d-1]+exp(-745))
    	    X = Netstats_cpp(history.t, node, netstat)
    	    XB = MultiplyXBList(X, b.old)   
    	    lambda[[d]] = lambda_cpp(p.d[d,], XB)
    	     for (i in node) {
            	X_u = lapply(X[[i]], function(X_IP) {X_IP[which(u[[d]][i,]==1),]})
            	if (length(X_u[[1]]) > P) {
                	X_u = lapply(X_u, function(X_u_IP) {geometric.mean(X_u_IP)})
            	}
            	Y = lapply(X_u, function(X_u_IP) {c(X_u_IP, timemat[d-1,])})
            	xi = MultiplyYeta(Y, eta.old)
            	mu[d, i] = mu_cpp(p.d[d,], xi)
      	}
	}
 
  chain.final = list(l = l, z = z, b = bmat, eta = etamat, delta = deltamat, sigma2_tau = sigma2_taumat,
                     u = u, sigma.Q =sigma.Q, alpha = alphavec, mvec = mvecmat, edge.trim = edge.trim,
                     proposal.var1= proposal.var1, proposal.var2= proposal.var2, convergence = convergence)
  return(chain.final)
}	


#' @title GenerateDocs
#' @description Generate a collection of documents according to the generative process of IPTM
#'
#' @param D number of documents to be generated
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
#' @param sigma2_tau variance parameter for the timestamps
#' @param l topic-interaction pattern assignment
#' @param support support of latent recipients
#' @param netstat which type of network statistics to use ("intercept", "dyadic", "triadic", "degree")
#' @param timestat additional statistics to be used for timestamps other than netstat ("timeofday", "dayofweek")
#' @param base.edge edges before 384 hours that is used to calculate initial history of interactions
#' @param base.text texts corresponding to base.edge
#' @param topic_token_assignments matrix of topic-token assignments
#' @param backward Logigal indicating whether we are generating backward samples (if FALSE -> forward)
#' @param base Logical indicating whether or not we are generating base edges (< 384)
#'
#' @return generated data including (author, recipients, timestamp, words)
#'
#' @export
GenerateDocs = function(D, node, vocab, nIP, K, n.d, alpha, mvec, beta,
                        b, eta, delta, sigma2_tau, l, support, netstat, timestat,
                        base.edge = NULL, base.text = NULL, topic_token_assignments = NULL,
                        backward = FALSE, base = FALSE) { 
  V = length(vocab)
  phi = lapply(1:K, function(k) {rdirichlet_cpp(1, rep(beta/V, V))})
  netstat = as.numeric(c("intercept", "degree", "dyadic", "triadic" ) %in% netstat)
  timestat = as.numeric(c("dayofweek","timeofday") %in% timestat)
  L = 3
  P = netstat[1]+L*(2*netstat[2]+2*netstat[3]+4*netstat[4])
  t.d = ifelse(base, 0, 384*3600)
  edge = base.edge
  text = base.text
  base.length = length(edge)
  p.d = matrix(NA, nrow = D+base.length, ncol = nIP)
  timemat = matrix(0, nrow = D+base.length, ncol = sum(timestat))

  if (!base) {
  for (d in 1:base.length) {
  	 p.d[d, ] = pdmat(list(as.integer(names(text[[d]]))), l, nIP)
  }
    if (sum(timestat) > 0) {
      time_ymd = as.POSIXct(vapply(base.edge, function(d) {d[[3]]}, c(1)), tz = "America/New_York", origin = "1970-01-01")
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
        	   lapply(1:3, function(x){
         	   matrix(0, length(node), length(node))
       	   })})  
  VKmat = matrix(0, V, K)
  u = list()
  timestamps = matrix(0, nrow = D, ncol = length(node))
  for (d in 1:D) {
    u[[d]] = matrix(0, length(node), length(node))
    text[[base.length+d]] = rep(NA, n.d)
    if (!backward) {
      theta.d = rdirichlet_cpp(1, alpha*mvec)	
      topic.d = multinom_vec(max(1, n.d), theta.d)
        for (n in 1:n.d){
          text[[base.length+d]][n] = multinom_vec(1, phi[[topic.d[n]]])
        }
        names(text[[base.length+d]]) = topic.d
    } else {
      phi.k = rep(NA, K)
      topic.d = topic_token_assignments[[d]]
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
    if (t.d >= 384*3600) {
      history.t = History(edge, p.d, node, t.d+exp(-745))
    }
    X = Netstats_cpp(history.t, node, netstat)
    vu = MultiplyXBList(X, b)     
    lambda = lambda_cpp(p.d[base.length+d,], vu)
    	for (i in node) {
    		u[[d]][i,-i] = r.gibbs.measure(1, lambda[i,-i], delta, support)
    		X_u = lapply(X[[i]], function(X_IP) {X_IP[which(u[[d]][i,]==1),]})
           if (length(X_u[[1]]) > P) {
           		X_u = lapply(X_u, function(X_u_IP) {geometric.mean(X_u_IP)})
           } 
           Y = lapply(X_u, function(X_u_IP) {c(X_u_IP, timemat[base.length+d-1,])})
           xi = MultiplyYeta(Y, eta)
           mu = mu_cpp(p.d[base.length+d,], xi)
           timestamps[d,i] = exp(rnorm(1, mu, sqrt(sigma2_tau)))*3600
	}
    i.d = which(timestamps[d,] == min(timestamps[d,]))
    j.d = which(u[[d]][i.d,] == 1)
    t.d = t.d+timestamps[d,i.d]
    edge[[base.length+d]] = list(author = i.d, recipients = j.d, timestamp = t.d)
    if (sum(timestat) > 0) {
    	it = 0
      time_ymd = as.POSIXct(edge[[base.length+d]][[3]], tz = "America/New_York", origin="1970-01-01")
      if (!is.na(time_ymd)) {
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
  }
   if (base == TRUE & t.d > 384*3600) {
    cutoff = which_int(384*3600, vapply(1:length(edge), function(d) {edge[[d]][[3]]}, c(1)))-1
    edge = edge[1:cutoff]
    text = text[1:cutoff]
   }
  return(list(edge = edge, text = text, base = base.length, 
  		b = b, eta = eta, delta = delta, sigma2_tau = sigma2_tau, l = l))							
} 



#' @title GiR_stats
#' @description Calculate several statistics from samples generated from forward or backward sampling
#'
#' @param GiR_sample one sample from generative process
#' @param V number of unique words
#'
#' @return A vector of statistics calculated from one GiR sample
#'
#' @export

GiR_stats = function(GiR_sample, V) {
  edge = GiR_sample$edge
  text = GiR_sample$text
  edge = edge[-(1:GiR_sample$base)]
  text = text[-(1:GiR_sample$base)]
  
  GiR_stats = c()
  D = length(edge)
  P = length(GiR_sample$b[[1]])
  Q = length(GiR_sample$eta[[1]])
  K = length(GiR_sample$l)
  nIP = length(GiR_sample$b)
  n.d = length(text[[1]])
  
  GiR_stats[1:P] = Reduce('+', GiR_sample$b) / nIP 
  GiR_stats[(P+1):(P+Q)] = Reduce('+', GiR_sample$eta) / nIP 
  GiR_stats[P+Q+1] = GiR_sample$delta
  GiR_stats[P+Q+2] = GiR_sample$sigma2_tau
  GiR_stats[P+Q+3] = mean(vapply(1:D, function(d) {length(edge[[d]][[2]])}, c(1)))
  GiR_stats[P+Q+4] = median(vapply(2:D, function(d) {edge[[d]][[3]]-edge[[d-1]][[3]]}, c(1))) 			
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
    if (grepl("sigma2_tau", nms[i]) ) {
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
           x = 0.65,
           y = 0.15,
           cex = 0.4)
  }
}      
     
     
#' @title GiR
#' @description Getting it Right test for the IPTM
#'
#' @param Nsamp number of GiR samples to be generated 
#' @param D number of documents to be generated per each sample
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
#' @param prior.tau prior shape and scale parameter of sigma2_tau in inverse-Gamma distribution
#' @param sigma.Q proposal distribution variance parameter
#' @param Outer size of outer iterations 
#' @param Inner size of inner iteration for Metropolis-Hastings updates
#' @param burn iterations to be discarded at the beginning of Metropolis-Hastings chains
#' @param thin the thinning interval of Metropolis-Hastings chains
#' @param netstat which type of network statistics to use ("intercept", dyadic", "triadic", "degree")
#' @param timestat additional statistics to be used for timestamps other than netstat ("timeofday", "dayofweek")
#' @param base.edge artificial collection of documents to be used as initial state of history
#' @param base.text artificial collection of documents to be used as initial state of history
#' @param generate_PP_plots Logical indicating whether to draw PP plots
#'
#' @return Forward and Backward samples and corresponding test results
#'
#' @export
GiR = function(Nsamp, D, node, vocab, nIP, K, n.d, alpha, mvec, beta, 
              prior.b, prior.delta, prior.eta, prior.tau, sigma.Q, Outer, Inner, burn, thin,
              netstat = c("intercept", "dyadic"), timestat = c("timeofday", "dayofweek"),
              base.edge, base.text, generate_PP_plots = TRUE) {
  
  netstat2 = as.numeric(c("intercept", "degree", "dyadic", "triadic") %in% netstat)
  timestat2 = as.numeric(c("dayofweek","timeofday") %in% timestat)
  L = 3
  P = netstat2[1]+L*(2*netstat2[2]+2*netstat2[3]+4*netstat2[4])
  Q = length(timestat2)
  V = length(vocab)
  support = gibbs.measure.support(length(node)-1)
  
  #Forward sampling
  Forward_stats = matrix(NA, nrow = Nsamp, ncol = P+(P+Q)+5+nIP+K+V)
  colnames(Forward_stats) = c(paste0("b_",1:P), paste0("eta_",1:(P+Q)), "delta", "simga2_tau", "Mean_recipients", "Mean_timediff", 
  					 "Mean_TopicIP", paste0("Tokens_in_IP_", 1:nIP), paste0("Tokens_in_Topic", 1:K), paste0("Tokens_in_Word",1:V))
  for (i in 1:Nsamp) { 
    if (i %% 5000 == 0) {cat("Forward sampling", i, "\n")}
    b = lapply(1:nIP, function(IP) {c(rcpp_rmvnorm(1, prior.b[[2]], prior.b[[1]]))}) 
    eta = lapply(1:nIP, function(IP) {c(rcpp_rmvnorm(1, prior.eta[[2]], prior.eta[[1]]))})
    delta = rnorm(1, prior.delta[1], sqrt(prior.delta[2]))
    sigma2_tau = 1/rgamma(1, prior.tau[1], prior.tau[2])

    l = sample(1:nIP, K, replace = TRUE)
    Forward_sample = GenerateDocs(D, node, vocab, nIP, K, n.d, alpha, mvec, beta, b, eta, delta, sigma2_tau, l, support, netstat, timestat,
                        base.edge = base.edge, base.text = base.text, topic_token_assignments = NULL, 
                        	     backward = FALSE, base = FALSE)
    Forward_stats[i, ] = GiR_stats(Forward_sample, V)
  }
  #Backward sampling
  Backward_stats = matrix(NA, nrow = Nsamp, ncol = ncol(Forward_stats))
  Backward_sample = GenerateDocs(D, node, vocab, nIP, K, n.d, alpha, mvec, beta, b, eta, delta, sigma2_tau, l, support, netstat, timestat,
                        base.edge = base.edge, base.text = base.text, topic_token_assignments = NULL, backward = FALSE, base = FALSE)
  for (i in 1:Nsamp) { 
     if (i %% 100 == 0) {cat("Backward Sampling", i, "\n")}
    Inference_samp = IPTM.inference(Backward_sample$edge, node, Backward_sample$text, vocab, nIP, K, sigma.Q, alpha, mvec, beta, 
    			    prior.b, prior.delta,prior.eta, prior.tau, Outer, Inner, burn, thin, netstat, timestat, optimize = FALSE, initial = NULL)
    b = lapply(1:nIP, function(IP) {Inference_samp$b[[IP]][,ncol(Inference_samp$b[[IP]])]})
    eta = lapply(1:nIP, function(IP) {Inference_samp$eta[[IP]][,ncol(Inference_samp$eta[[IP]])]})
    delta = Inference_samp$delta[length(Inference_samp$delta)]
    sigma2_tau = Inference_samp$sigma2_tau[length(Inference_samp$sigma2_tau)]
    l = Inference_samp$l
    z = Inference_samp$z
    for (d in 1:length(z)) {
      names(z[[d]]) = Backward_sample$text[[d]]
    }
    
    Backward_sample = GenerateDocs(D, node, vocab, nIP, K, n.d, alpha, mvec, beta, b, eta, delta, sigma2_tau, l, support, netstat, timestat,
                        base.edge = base.edge, base.text = base.text, topic_token_assignments = z, 
                        	     backward = TRUE, base = FALSE)
    Backward_stats[i, ] = GiR_stats(Backward_sample, V)
    }
  				
  if (generate_PP_plots) {
    par(mfrow=c(4,8), oma = c(3,3,3,3), mar = c(2,1,1,1))
    GiR_PP_Plots(Forward_stats, Backward_stats)
  }			
  return(list(Forward = Forward_stats, Backward = Backward_stats))
}                         	          
      