#' @useDynLib IPTM
#' @import stats
#' @import grDevices
#' @import graphics
#' @import RcppArmadillo
#' @importFrom Rcpp sourceCpp
#' @importFrom entropy entropy.empirical
#' @importFrom lda top.topic.words
#' @importFrom reshape melt
#' @importFrom coda mcmc geweke.diag
#' @importFrom combinat permn
#' @importFrom FastGP rcpp_rmvnorm rcpp_log_dmvnorm

NULL

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
		gibbs.mat.i = do.call('rbind',permn(c(rep(1, i), rep(0, n - i))))
		gibbs.support = rbind(gibbs.support, gibbs.mat.i)
	}
	as.matrix(unique(gibbs.support))
}

#' @title r.gibbs.measure
#' @description List out the support of Gibbs measure
#'
#' @param nsamp the number of binary vector samples to draw
#' @param lambdai a vector of coefficients according to which each element
#' @param delta a positive real valued parameter that controls the density penalty 
#' @param support support of Gibbs measure
#'
#' @return nsamp number of samples with each row denoting binary vector
#'
#' @export
r.gibbs.measure <- function(nsamp, lambdai, delta, support) {
	gibbsNormalizer = prod(exp(delta + log(lambdai)) + 1) - 1
	if (gibbsNormalizer == 0) {gibbsNormalizer = exp(-745)}
	logitNumerator = vapply(1:nrow(support), function(s) {
					 	exp(sum((delta + log(lambdai)) * support[s,]))
					  }, c(1))
					  
	samp <- sample(1:nrow(support), nsamp, replace = TRUE, prob = logitNumerator / gibbsNormalizer)
	
	support[samp,]	
}

#' @title adaptive_MH 
#' @description adaptive Metropolis Hastings to maintain target acceptance rate
#'
#' @param sigma_Q proposal distribution variance parameter for beta and delta
#' @param accept_rates acceptance rate from previous iteration
#' @param target target acceptance rate
#' @param update_size size of update to be adjusted 
#' @param tol tolerance level to determine if acceptance rate is too high
#'
#' @return nsamp number of samples with each row denoting binary vector
#'
#' @export
adaptive_MH = function(sigma_Q, accept_rates, target = 0.25, update_size, tol = 0.15) {
	for (i in 1:length(sigma_Q)) {
		if (accept_rates[i] < target) {
				sigma_Q[i] = sigma_Q[i] - update_size[i]
		}
		if (accept_rates[i] > (target + tol)) {
			sigma_Q[i] = sigma_Q[i] + update_size[i]
		}		
	}
	return(sigma_Q)
}

#' @title AlphamvecOpt
#' @description Optimize the hyperparmeter vector (= alpha * mvec) given the current assignment of C and Z
#'
#' @param K total number of topics specified by the user
#' @param currentZ current state of the assignment of topics 
#' @param alpha Dirichlet concentration prior for topic distribution
#' @param mvec Dirichlet base prior for topic distribution
#' @param niter number of iterations to perfom
#'
#' @return The optimized value of the vector (= alpha * mvec)
#'
#' @export
AlphamvecOpt =  function(K, currentZ, alpha, mvec, niter) {
	# Optimize the hyperparmeter vector given the current assignment of Z
	#
	# Args 
	#  K total number of topics specified by the user
	#  currentZ current state of the assignment of topics 
	#  alpha Dirichlet concentration prior for topic distribution
	#  mvec Dirichlet base prior for topic distribution
	#
	# Returns
	#  The optimized value of the vector (= alpha * mvec)
	
	final.vec = list()
	
	current.vec = alpha * mvec
	n.word = mapply(length, currentZ)
	n.word.table = tabulate(n.word)
	nK.word.list = matrix(NA, nrow = length(currentZ), ncol = K)
	for (d in seq(along = currentZ)) {
		nK.word.list[d, ] = tabulateC(currentZ[[d]], K)
	}
	nK.word.table = lapply(1:K, function(k){
			  tabulate(nK.word.list[,k])
			  })
	for (i in 1:niter) {
	alpha = sum(current.vec)
	S = UpdateDenom(alpha, n.word.table)		
	s = UpdateNum(current.vec, nK.word.table)
	current.vec = current.vec * s / S
	}
	final.vec = current.vec	
	return(final.vec)	
}	


#' @title IPTM_inference.Gibbs
#' @description Iterate Markov Chain Monte Carlo (MCMC) algorithm using Gibbs measure to sequentially update the assignments of Z, C and B
#'
#' @param edge list of document information with 3 elements (element 1 sender, element 2 receiver, element 3 time in unix.time format)
#' @param node nodelist containing the ID of nodes (ID starting from 1)
#' @param textlist list of text (length=number of documents in total) containing the words in each document
#' @param vocabulary all vocabularies used over the corpus
#' @param nIP total number of interaction patterns specified by the user
#' @param K total number of topics specified by the user
#' @param sigma_Q proposal distribution variance parameter for beta and delta
#' @param alpha Dirichlet concentration prior for document-topic distribution
#' @param mvec Dirichlet base prior for document-topic distribution
#' @param betas Dirichlet concentration prior for topic-word distribution
#' @param nvec Dirichlet base prior for topic-word distribution
#' @param prior.b.mean mean vector of b in multivariate normal distribution
#' @param prior.b.var covairance matrix of b in multivariate normal distribution
#' @param prior.delta parameter of delta in Normal prior
#' @param out size of outer iterations 
#' @param n_B size of third inner iteration for updates of B
#' @param n_d size of third inner iteration for updates of delta
#' @param burn iterations to be discarded at the beginning of beta chain
#' @param thinning the thinningning interval of beta chain
#' @param netstat which type of network statistics to use ("intercept", dyadic", "triadic", "degree")
#' @param optimize to optimize alpha (Dirichlet concentration prior for document-topic distribution) or not (TRUE/FALSE)
#'
#' @return MCMC output containing IP assignment, topic assignment, and (beta, mu, delta) chain 
#'
#' @export
IPTM_inference.Gibbs = function(edge, node, textlist, vocabulary, nIP, K, sigma_Q, alpha, mvec, betas, nvec, prior.b.mean, prior.b.var, prior.delta, 
					out, n_B, n_d, burn, thinning, netstat, optimize = FALSE) {
   
  # trim the edge so that we only model edges after 384 hours
	timestamps = vapply(edge, function(d) {
  			  d[[3]]
 			  }, c(1))

    edge2 = which_int(384, timestamps) : length(edge)
    maxedge2 = max(edge2)
    timeinc = c(timestamps[1], timestamps[-1] - timestamps[-length(timestamps)])
   
  # initialize alpha, mvec, delta, nvec, delta, lvec, and gammas
 	W = length(vocabulary)
  	phi = lapply(1:K, function(k) {
		rdirichlet_cpp(1, betas * nvec)
	})
	delta = rnorm(1, prior.delta[1], sqrt(prior.delta[2]))
	beta.old = lapply(1:nIP, function(IP) {
      			  c(rcpp_rmvnorm(1,  prior.b.var, prior.b.mean))
      			 })
	# initialize C, theta and Z
	 currentC = sample(1:nIP, K, replace = TRUE)  	 
	 theta = rdirichlet_cpp(length(edge), alpha * mvec)
     currentZ = lapply(seq(along = edge), function(d) {
     			 multinom_vec(max(1, length(textlist[[d]])), theta[d, ])
 				 })
 	p.d = pdmat(currentZ, currentC, nIP)
	accept.rates = rep(0, 2)
    # initialize beta
    netstat = as.numeric(c("intercept", "degree", "dyadic", "triadic" ) %in% netstat)
    L = 3
    P = netstat[1] + L * (2 * netstat[2] + 2 * netstat[3] + 4 * netstat[4])
    bmat = list()
	for (IP in 1:nIP) {
		bmat[[IP]] = matrix(beta.old[[IP]], nrow = P, ncol = (n_B - burn) / thinning)
  	}
  	deltamat = rep(delta, n_d)
    proposal.var = lapply(1:nIP, function(IP){diag(P)})

  	#initialize the latent sender-receiver pairs
  	iJi = lapply(seq(along = edge), function(d) {
  		matrix(0, nrow = length(node), ncol = length(node))
  	})
  	lambda = list()
    LambdaiJi = list()
	observediJi = list()
	support = gibbs.measure.support(length(node) - 1)
	for (d in edge2) {
   	 	history.t = History(edge, p.d, node, edge[[d-1]][[3]] + exp(-745))
   	 	X = Netstats_cpp(history.t, node, netstat)
   	 	XB = MultiplyXBList(X, beta.old)     
		lambda[[d]] = lambda_cpp(p.d[d,], XB)
    	for (i in node) {
    		iJi[[d]][i, -i] = r.gibbs.measure(1, lambda[[d]][i, -i], delta, support)
    		}
    	}
     textlist.raw = unlist(textlist[edge2])
     table.W = lapply(1:K, function(k) {
      			 tabulateC(textlist.raw[which(unlist(currentZ[edge2]) == k)], W)
      			 })	
	alphamat = matrix(NA, nrow = 0, ncol = 1)
	mvecmat = matrix(NA, nrow = 0, ncol = K)
    #start outer iteration
    for (o in 1:out) {
      
      if (optimize) {
      #update the hyperparameter alpha and mvec
      vec = AlphamvecOpt(K, currentZ[edge2], alpha, mvec, 5)
      alpha = sum(vec)
      mvec = vec / alpha
      alphamat = rbind(alphamat, alpha)
      mvecmat = rbind(mvecmat, mvec)
      }
     # Data augmentation
      for (d in edge2) {
   	 	history.t = History(edge, p.d, node, edge[[d-1]][[3]] + exp(-745))
   	 	X = Netstats_cpp(history.t, node, netstat)
   	 	XB = MultiplyXBList(X, beta.old)   
		lambda[[d]] = lambda_cpp(p.d[d,], XB)
		#calculate the resampling probability	
		for (i in node[-edge[[d]][[1]]]) {
			XB_IP = lapply(XB, function(IP) {IP[i,]}) 
			for (j in sample(node[-i], length(node) - 1)) {
				probij = DataAug_cpp_Gibbs(iJi[[d]][i, ], lambda[[d]][i,], XB_IP, p.d[d, ], delta, timeinc[d], j)
				iJi[[d]][i, j] = multinom_vec(1, probij) - 1		
				}
		}
		iJi[[d]][edge[[d]][[1]],] = tabulateC(as.numeric(unlist(edge[[d]][2])), length(node))
		LambdaiJi[[d]] = lambdaiJi(p.d[d,], XB, iJi[[d]])
		observediJi[[d]] = LambdaiJi[[d]][edge[[d]][[1]]]
		}	 

      # C update
       for (k in sort(unique(unlist(currentZ[edge2])))) { 
         const.C = rep(NA, nIP)
         for (IP in 1:nIP) {
          	currentC[k] = IP
          	p.d = pdmat(currentZ, currentC, nIP)
 		 	      history.t = History(edge, p.d, node, edge[[maxedge2-1]][[3]] + exp(-745))
    	   	X = Netstats_cpp(history.t, node, netstat)
    	    XB = MultiplyXBList(X, beta.old)    
           	lambda[[maxedge2]] = lambda_cpp(p.d[maxedge2,], XB)
		    LambdaiJi[[maxedge2]] = lambdaiJi(p.d[maxedge2,], XB, iJi[[maxedge2]])
           	observediJi[[maxedge2]] = LambdaiJi[[maxedge2]][edge[[maxedge2]][[1]]]
           	prob = EdgeTime(iJi[[maxedge2]], lambda[[maxedge2]], delta, LambdaiJi[[maxedge2]], timeinc[maxedge2], observediJi[[maxedge2]]) 
           	const.C[IP] = prob
      	 }	
         currentC[k] = multinom_vec(1, expconst(const.C))
	}
     
     # Z update
   	 for (d in edge2) {
       textlist.d = textlist[[d]]
       if (edge[[d]][[3]] + 384 > edge[[maxedge2]][[3]]) {
      	 hist.d = maxedge2
       } else {
      	 hist.d = which_num(edge[[d]][[3]] + 384, timestamps)
       }
       edgetime.d = rep(NA, K)
       for (w in 1:length(currentZ[[d]])) {
       	 zw.old = currentZ[[d]][w]
       	 if (length(textlist.d) > 0) {
       	 table.W[[zw.old]][textlist.d[w]] = table.W[[zw.old]][textlist.d[w]] - 1
         topicpart.d = TopicInEqZ(K, currentZ[[d]][-w], alpha, mvec)
         wordpart.d = WordInEqZ(K, textlist.d, table.W, betas, nvec)
       } else {
         topicpart.d = 0
         wordpart.d = matrix(0, nrow = length(currentZ[[d]]), ncol = K)
       }
       	  for (IP in unique(currentC)) {
       	  	currentCK = which(currentC == IP)
       		currentZ[[d]][w] = min(currentCK)
       	 	p.d[d, ] = pdmat(list(currentZ[[d]]), currentC, nIP)           
            history.t = History(edge, p.d, node, edge[[hist.d-1]][[3]] + exp(-745))
    	    		X = Netstats_cpp(history.t, node, netstat)
    	    		XB = MultiplyXBList(X, beta.old)   
    	    		lambda[[hist.d]] = lambda_cpp(p.d[hist.d,], XB)
	        LambdaiJi[[hist.d]] = lambdaiJi(p.d[hist.d, ], XB, iJi[[hist.d]])
        		observediJi[[hist.d]] = LambdaiJi[[hist.d]][edge[[hist.d]][[1]]]
        		edgetime.d[currentCK] = EdgeTime(iJi[[hist.d]], lambda[[hist.d]], delta, LambdaiJi[[hist.d]], timeinc[hist.d], observediJi[[hist.d]])
        }
         const.Z = edgetime.d + topicpart.d + wordpart.d[w, ]
         zw.new = multinom_vec(1, expconst(const.Z))
      
         if (zw.new != zw.old) {
           currentZ[[d]][w] = zw.new
           table.W[[zw.new]][textlist.d[w]] = table.W[[zw.new]][textlist.d[w]] + 1
         } else {
         	currentZ[[d]][w] = zw.old
         	table.W[[zw.old]][textlist.d[w]] = table.W[[zw.old]][textlist.d[w]] + 1
         }
       }
     }

      # C update given Z and B - withinning each document d
       for (k in sort(unique(unlist(currentZ[edge2])))) { 
         const.C = rep(NA, nIP)
         for (IP in 1:nIP) {
          	currentC[k] = IP
          	p.d = pdmat(currentZ, currentC, nIP)
 		 	      history.t = History(edge, p.d, node, edge[[maxedge2-1]][[3]] + exp(-745))
    	   	X = Netstats_cpp(history.t, node, netstat)
    	    XB = MultiplyXBList(X, beta.old)    
           	lambda[[maxedge2]] = lambda_cpp(p.d[maxedge2,], XB)
		        LambdaiJi[[maxedge2]] = lambdaiJi(p.d[maxedge2,], XB, iJi[[maxedge2]])
           	observediJi[[maxedge2]] = LambdaiJi[[maxedge2]][edge[[maxedge2]][[1]]]
           	prob = EdgeTime(iJi[[maxedge2]], lambda[[maxedge2]], delta, LambdaiJi[[maxedge2]], timeinc[maxedge2], observediJi[[maxedge2]]) 
           	const.C[IP] = prob
      	 }	
         currentC[k] = multinom_vec(1, expconst(const.C))
	}
     
 	p.d = pdmat(currentZ, currentC, nIP)
 
    for (d in maxedge2) {
        	history.t = History(edge, p.d, node, edge[[d-1]][[3]] + exp(-745))
    	    X = Netstats_cpp(history.t, node, netstat)
    	    XB = MultiplyXBList(X, beta.old)   
    	    lambda[[d]] = lambda_cpp(p.d[d,], XB)
	    LambdaiJi[[d]] = lambdaiJi(p.d[d,], XB, iJi[[d]])
        	observediJi[[d]] = LambdaiJi[[d]][edge[[d]][[1]]]
	}
	
    # beta update
    prior.old1 = sum(vapply(1:nIP, function(IP) {rcpp_log_dmvnorm(prior.b.var, prior.b.mean, beta.old[[IP]], FALSE)}, c(1)))
    post.old1 = EdgeTime(iJi[[maxedge2]], lambda[[maxedge2]], delta, LambdaiJi[[maxedge2]], timeinc[maxedge2], observediJi[[maxedge2]])
    if (o != 1) {
    	accept.rates[1] = accept.rates[1] / n_B
    	accept.rates[2] = accept.rates[2] / n_d
      sigma_Q = adaptive_MH(sigma_Q, accept.rates, update_size = 0.2 * sigma_Q)
      if (accept.rates[1] > 0.15) { 
        for (IP in 1:nIP) {
          proposal.var[[IP]] = cor(t(bmat[[IP]]))
        }
      }
    }
    accept.rates = rep(0, 2)
    for (i3 in 1:n_B) {
      beta.new = lapply(1:nIP, function(IP) {
        c(rcpp_rmvnorm(1, sigma_Q[1] * proposal.var[[IP]], beta.old[[IP]]))
      }) 
      for (d in maxedge2) {
        XB = MultiplyXBList(X, beta.new)
        lambda[[d]] = lambda_cpp(p.d[d,], XB)    
        LambdaiJi[[d]] = lambdaiJi(p.d[d,], XB, iJi[[d]])
        observediJi[[d]] = LambdaiJi[[d]][edge[[d]][[1]]]
      }
      prior.new1 = sum(vapply(1:nIP, function(IP) {rcpp_log_dmvnorm(prior.b.var, prior.b.mean, beta.new[[IP]], FALSE)}, c(1)))
      post.new1 = EdgeTime(iJi[[maxedge2]], lambda[[maxedge2]], delta, LambdaiJi[[maxedge2]], timeinc[maxedge2], observediJi[[maxedge2]])
      loglike.diff = prior.new1 + post.new1 - prior.old1 - post.old1
      if (log(runif(1, 0, 1)) < loglike.diff) {
        for (IP in 1:nIP) {
          beta.old[[IP]] = beta.new[[IP]]
        }
        prior.old1 = prior.new1
        post.old1 = post.new1
        accept.rates[1] = accept.rates[1] + 1
      }
      if (i3 > burn[1] && i3 %% (thinning[1]) == 0) {
        for (IP in 1:nIP) {
          bmat[[IP]][ , (i3 - burn[1]) / thinning[1]] = beta.old[[IP]]
        }
      }
    }
    
    #delta update
    for (d in maxedge2) {
      XB = MultiplyXBList(X, beta.old)
      lambda[[d]] = lambda_cpp(p.d[d,], XB)  
    }
    prior.old2 = dnorm(delta, prior.delta[1], sqrt(prior.delta[2]), log = TRUE)
    post.old2 = EdgeInEqZ_Gibbs(iJi[[maxedge2]], lambda[[maxedge2]], delta)
    for (i4 in 1:n_d) {
      delta.new = rnorm(1, delta, sqrt(sigma_Q[2]))
      prior.new2 = dnorm(delta.new, prior.delta[1], sqrt(prior.delta[2]), log = TRUE)
      post.new2 = EdgeInEqZ_Gibbs(iJi[[maxedge2]], lambda[[maxedge2]], delta.new)
      loglike.diff2 = prior.new2 + post.new2 - prior.old2 - post.old2 
      if (log(runif(1, 0, 1)) < loglike.diff2) {
        delta = delta.new
        prior.old2 = prior.new2
        post.old2 = post.new2
        accept.rates[2] = accept.rates[2] + 1
      } 
      if (i4 > burn[2] && i4 %% (thinning[2]) == 0) {
        deltamat[(i4 - burn[2]) / thinning[2]] = delta
      }
    }
 }
     
  chain.final = list(C = currentC, Z = lapply(edge2, function(d) {currentZ[[d]]}), B = bmat, D = deltamat,
                     iJi = iJi, sigma_Q =sigma_Q, alpha = alphamat, mvec = mvecmat, 
                     proposal.var= proposal.var)
  return(chain.final)
}

#' @title IPTM_inference.data
#' @description Iterate Markov Chain Monte Carlo (MCMC) algorithm using Gibbs measure to sequentially update the assignments of Z, C and B
#'
#' @param edge list of document information with 3 elements (element 1 sender, element 2 receiver, element 3 time in unix.time format)
#' @param node nodelist containing the ID of nodes (ID starting from 1)
#' @param textlist list of text (length=number of documents in total) containing the words in each document
#' @param vocabulary all vocabularies used over the corpus
#' @param nIP total number of interaction patterns specified by the user
#' @param K total number of topics specified by the user
#' @param sigma_Q proposal distribution variance parameter for beta and delta
#' @param alpha Dirichlet concentration prior for document-topic distribution
#' @param mvec Dirichlet base prior for document-topic distribution
#' @param betas Dirichlet concentration prior for topic-word distribution
#' @param nvec Dirichlet base prior for topic-word distribution
#' @param prior.b.mean mean vector of b in multivariate normal distribution
#' @param prior.b.var covairance matrix of b in multivariate normal distribution
#' @param prior.delta parameter of delta in Normal prior
#' @param out size of outer iterations 
#' @param n_B size of third inner iteration for updates of B
#' @param n_d size of third inner iteration for updates of delta
#' @param burn iterations to be discarded at the beginning of beta and delta chain 
#' @param thinning the thinningning interval of beta and delta chain
#' @param netstat which type of network statistics to use ("intercept", dyadic", "triadic", "degree")
#' @param optimize to optimize alpha (Dirichlet concentration prior for document-topic distribution) or not (TRUE/FALSE)
#' @param initial list of initial values user wants to assign
#'
#' @return MCMC output containing IP assignment, topic assignment, and (beta, mu, delta) chain 
#'
#' @export
IPTM_inference.data = function(edge, node, textlist, vocabulary, nIP, K, sigma_Q, alpha, mvec, betas, nvec, prior.b.mean, prior.b.var, prior.delta, 
                               out, n_B, n_d, burn, thinning, netstat, optimize = FALSE, initial = NULL) {
  # trim the edge so that we only model edges after 384 hours
  timestamps = vapply(edge, function(d) {
    d[[3]]
  }, c(1))
  
  edge2 = which_int(384, timestamps) : length(edge)
  maxedge2 = max(edge2)
  timeinc = c(timestamps[1], timestamps[-1] - timestamps[-length(timestamps)])
  
  # initialize alpha, mvec, delta, nvec, delta, lvec, and gammas
  W = length(vocabulary)
  phi = lapply(1:K, function(k) {
    rdirichlet_cpp(1, betas * nvec)
  })
  theta = rdirichlet_cpp(length(edge), alpha * mvec)
  if (length(initial) == 0) {
  delta = rnorm(1, prior.delta[1], sqrt(prior.delta[2]))
  beta.old = lapply(1:nIP, function(IP) {
    prior.b.mean
  })
  # initialize C, theta and Z
  currentC = sample(1:nIP, K, replace = TRUE) 
  currentZ = lapply(seq(along = edge), function(d) {
    multinom_vec(max(1, length(textlist[[d]])), theta[d, ])
  })
  } else {
  delta = initial$D
  beta.old = initial$B
  currentC = initial$C
  currentZ = initial$Z  
  }
  p.d = pdmat(currentZ, currentC, nIP) 
  convergence = c()
  netstat = as.numeric(c("intercept", "degree", "dyadic", "triadic" ) %in% netstat)
  L = 3
  P = netstat[1] + L * (2 * netstat[2] + 2 * netstat[3] + 4 * netstat[4])

  bmat = list()
  for (IP in 1:nIP) {
    bmat[[IP]] = matrix(beta.old[[IP]], nrow = P, ncol = (n_B - burn[1]) / thinning[1])
  }
  deltamat = rep(delta,  (n_d - burn[2]) / thinning[2])
  proposal.var = lapply(1:nIP, function(IP){diag(P)})
    
  #initialize the latent sender-receiver pairs
  iJi = lapply(seq(along = edge), function(d) {
    matrix(0, nrow = length(node), ncol = length(node))
  })
  lambda = list()
  LambdaiJi = list()
  observediJi = list()
  for (d in edge2) {
    iJi[[d]] = matrix(rbinom(length(node)^2, 1, 0.5), nrow =length(node), ncol = length(node))
    diag(iJi[[d]]) = 0
  }
  textlist.raw = unlist(textlist)
  alphamat = matrix(NA, nrow = 0, ncol = 1)
  mvecmat = matrix(NA, nrow = 0, ncol = K)
  n_B2 = n_B / 5
  n_d2 = n_d / 5
   accept.rates = rep(0, 2)
    #start outer iteration
    for (o in 1:out) {
    	if (o == out) {
    		n_B2 = n_B
    		n_d2 = n_d
    	}
      print(o)
      if (optimize) {
      #update the hyperparameter alpha and mvec
      vec = AlphamvecOpt(K, currentZ[edge2], alpha, mvec, 5)
      alpha = sum(vec)
      mvec = vec / alpha
      alphamat = rbind(alphamat, alpha)
      mvecmat = rbind(mvecmat, mvec)
      }
      # Data augmentation
    for (d in edge2) {
      history.t = History(edge, p.d, node, edge[[d-1]][[3]] + exp(-745))
      X = Netstats_cpp(history.t, node, netstat)
      XB = MultiplyXBList(X, beta.old)     
      lambda[[d]] = lambda_cpp(p.d[d,], XB)
      #calculate the resampling probability	
      for (i in node[-edge[[d]][[1]]]) {
      	XB_IP = lapply(XB, function(IP) {IP[i,]})
        for (j in sample(node[-i], length(node) - 1)) {
          probij = DataAug_cpp_Gibbs(iJi[[d]][i, ], lambda[[d]][i,], XB_IP, p.d[d, ], delta, timeinc[d], j)
          iJi[[d]][i, j] = multinom_vec(1, probij) - 1
        }
      }
      iJi[[d]][edge[[d]][[1]],] = tabulateC(as.numeric(unlist(edge[[d]][2])), length(node))
      LambdaiJi[[d]] = lambdaiJi(p.d[d,], XB, iJi[[d]])
      observediJi[[d]] = LambdaiJi[[d]][edge[[d]][[1]]]
    }	 
       
    # Z update	
      table.W = lapply(1:K, function(k) {
      			tabulateC(textlist.raw[which(unlist(currentZ) == k)], W)
      			})
	   for (d in 1:(edge2[1] - 1)) {
      		 textlist.d = textlist[[d]]        		 
          	 for (w in 1:length(currentZ[[d]])) {
          	 	 zw.old = currentZ[[d]][w]
          	 	if (length(textlist.d) > 0) {
       	 		 table.W[[zw.old]][textlist.d[w]] = table.W[[zw.old]][textlist.d[w]] - 1
        			 topicpart.d = TopicInEqZ(K, currentZ[[d]][-w], alpha, mvec)
       		 	 wordpart.d = WordInEqZ(K, textlist.d, table.W, betas, nvec)
        		 } else {
        			topicpart.d = 0
        			wordpart.d = matrix(0, nrow = length(currentZ[[d]]), ncol = K)
        		 }
          		 const.Z = topicpart.d + wordpart.d[w, ]
          		 zw.new = multinom_vec(1, expconst(const.Z))
          		 if (zw.new != zw.old) {
            			 currentZ[[d]][w] = zw.new
           			 table.W[[zw.new]][textlist.d[w]] = table.W[[zw.new]][textlist.d[w]] + 1
      				 p.d[d, ] = pdmat(list(currentZ[[d]]), currentC, nIP) 	
               	 } else {
               	     table.W[[zw.old]][textlist.d[w]] = table.W[[zw.old]][textlist.d[w]] + 1 	
               	 }
        	 }
       	 }
 	 for (d in edge2) {
       textlist.d = textlist[[d]]
       if (edge[[d]][[3]] + 384 > edge[[maxedge2]][[3]]) {
      	 hist.d = maxedge2
       } else {
      	 hist.d = which_num(edge[[d]][[3]] + 384, timestamps)
       }
       edgetime.d = rep(NA, K)
       for (w in 1:length(currentZ[[d]])) {
       	 zw.old = currentZ[[d]][w]
       	 if (length(textlist.d) > 0) {
       	 table.W[[zw.old]][textlist.d[w]] = table.W[[zw.old]][textlist.d[w]] - 1
         topicpart.d = TopicInEqZ(K, currentZ[[d]][-w], alpha, mvec)
         wordpart.d = WordInEqZ(K, textlist.d, table.W, betas, nvec)
       } else {
         topicpart.d = 0
         wordpart.d = matrix(0, nrow = length(currentZ[[d]]), ncol = K)
       }
	  for (IP in unique(currentC)) {
       	  	currentCK = which(currentC == IP)
       		currentZ[[d]][w] = min(currentCK)
       	 	p.d[d, ] = pdmat(list(currentZ[[d]]), currentC, nIP)           
            history.t = History(edge, p.d, node, edge[[hist.d-1]][[3]]+ exp(-745))
    	    		X = Netstats_cpp(history.t, node, netstat)
    	    		XB = MultiplyXBList(X, beta.old)   
    	    		lambda[[hist.d]] = lambda_cpp(p.d[hist.d,], XB)
	    		LambdaiJi[[hist.d]] = lambdaiJi(p.d[hist.d, ], XB, iJi[[hist.d]])
       	 	observediJi[[hist.d]] = LambdaiJi[[hist.d]][edge[[hist.d]][[1]]]
        		edgetime.d[currentCK] = EdgeTime(iJi[[hist.d]], lambda[[hist.d]], delta, LambdaiJi[[hist.d]], timeinc[hist.d], observediJi[[hist.d]])
        }
         const.Z = edgetime.d + topicpart.d + wordpart.d[w, ]
         zw.new = multinom_vec(1, expconst(const.Z))
         if (zw.new != zw.old) {
           currentZ[[d]][w] = zw.new
           table.W[[zw.new]][textlist.d[w]] = table.W[[zw.new]][textlist.d[w]] + 1
         } else {
         	currentZ[[d]][w] = zw.old
            table.W[[zw.old]][textlist.d[w]] = table.W[[zw.old]][textlist.d[w]] + 1
         }
       }
     }
     
    # C update 
       for (k in sort(unique(unlist(currentZ[edge2])))) { 
         const.C = rep(NA, nIP)
         for (IP in 1:nIP) {
          	currentC[k] = IP
          	p.d = pdmat(currentZ, currentC, nIP) 
 		 	history.t = History(edge, p.d, node, edge[[maxedge2-1]][[3]] + exp(-745))
    	   		X = Netstats_cpp(history.t, node, netstat)
    	   		XB = MultiplyXBList(X, beta.old)    
           	lambda[[maxedge2]] = lambda_cpp(p.d[maxedge2,], XB)
		    LambdaiJi[[maxedge2]] = lambdaiJi(p.d[maxedge2,], XB, iJi[[maxedge2]])
           	observediJi[[maxedge2]] = LambdaiJi[[maxedge2]][edge[[maxedge2]][[1]]]
           	prob = EdgeTime(iJi[[maxedge2]], lambda[[maxedge2]], delta, LambdaiJi[[maxedge2]], timeinc[maxedge2], observediJi[[maxedge2]])
           	const.C[IP] = prob
      	 }
         currentC[k] = multinom_vec(1, expconst(const.C))
	}    
 	p.d = pdmat(currentZ, currentC, nIP)  
    for (d in maxedge2) {
        	history.t = History(edge, p.d, node, edge[[d-1]][[3]] + exp(-745))
    	    X = Netstats_cpp(history.t, node, netstat)
    	    XB = MultiplyXBList(X, beta.old)   
    	    lambda[[d]] = lambda_cpp(p.d[d,], XB)
	    	LambdaiJi[[d]] = lambdaiJi(p.d[d,], XB, iJi[[d]])
        	observediJi[[d]] = LambdaiJi[[d]][edge[[d]][[1]]]
	}
	        
    # beta update
    prior.old1 = sum(vapply(1:nIP, function(IP) {rcpp_log_dmvnorm(prior.b.var, prior.b.mean, beta.old[[IP]], FALSE)}, c(1)))
    post.old1 = EdgeTime(iJi[[maxedge2]], lambda[[maxedge2]], delta, LambdaiJi[[maxedge2]], timeinc[maxedge2], observediJi[[maxedge2]])
    if (o != 1) {
    	accept.rates[1] = accept.rates[1] / n_B2
    	accept.rates[2] = accept.rates[2] / n_d2
      sigma_Q = adaptive_MH(sigma_Q, accept.rates, update_size = 0.2 * sigma_Q)
      if (accept.rates[1] > 0.15) { 
        for (IP in 1:nIP) {
          proposal.var[[IP]] = cor(t(bmat[[IP]]))
          proposal.var[[IP]][is.na(proposal.var[[IP]])] = 0
        }
      }
    }
    accept.rates = rep(0, 2)
    for (i3 in 1:n_B2) {
      beta.new = lapply(1:nIP, function(IP) {
        c(rcpp_rmvnorm(1, sigma_Q[1] * proposal.var[[IP]], beta.old[[IP]]))
      }) 
      for (d in maxedge2) {
        XB = MultiplyXBList(X, beta.new)
        lambda[[d]] = lambda_cpp(p.d[d,], XB)    
        LambdaiJi[[d]] = lambdaiJi(p.d[d,], XB, iJi[[d]])
        observediJi[[d]] = LambdaiJi[[d]][edge[[d]][[1]]]
      }
      prior.new1 = sum(vapply(1:nIP, function(IP) {rcpp_log_dmvnorm(prior.b.var, prior.b.mean, beta.new[[IP]], FALSE)}, c(1)))
      post.new1 = EdgeTime(iJi[[maxedge2]], lambda[[maxedge2]], delta, LambdaiJi[[maxedge2]], timeinc[maxedge2], observediJi[[maxedge2]])
      loglike.diff = prior.new1 + post.new1 - prior.old1 - post.old1
      if (log(runif(1, 0, 1)) < loglike.diff) {
        for (IP in 1:nIP) {
          beta.old[[IP]] = beta.new[[IP]]
        }
        prior.old1 = prior.new1
        post.old1 = post.new1
        accept.rates[1] = accept.rates[1] + 1
      }
      if (i3 > burn[1] && i3 %% (thinning[1]) == 0) {
        for (IP in 1:nIP) {
          bmat[[IP]][ , (i3 - burn[1]) / thinning[1]] = beta.old[[IP]]
        }
      }
    }
    
    #delta update
    for (d in maxedge2) {
      XB = MultiplyXBList(X, beta.old)
      lambda[[d]] = lambda_cpp(p.d[d,], XB)  
    }
    prior.old2 = dnorm(delta, prior.delta[1], sqrt(prior.delta[2]), log = TRUE)
    post.old2 = EdgeInEqZ_Gibbs(iJi[[maxedge2]], lambda[[maxedge2]], delta)
    for (i4 in 1:n_d2) {
      delta.new = rnorm(1, delta, sqrt(sigma_Q[2]))
      prior.new2 = dnorm(delta.new, prior.delta[1], sqrt(prior.delta[2]), log = TRUE)
      post.new2 = EdgeInEqZ_Gibbs(iJi[[maxedge2]], lambda[[maxedge2]], delta.new)
      loglike.diff2 = prior.new2 + post.new2 - prior.old2 - post.old2 
      if (log(runif(1, 0, 1)) < loglike.diff2) {
        delta = delta.new
        prior.old2 = prior.new2
        post.old2 = post.new2
        accept.rates[2] = accept.rates[2] + 1
      } 
      if (i4 > burn[2] && i4 %% (thinning[2]) == 0) {
        deltamat[(i4 - burn[2]) / thinning[2]] = delta
      }
    }
   	  for (d in maxedge2) {
        XB = MultiplyXBList(X, beta.old)
        lambda[[d]] = lambda_cpp(p.d[d,], XB)    
        LambdaiJi[[d]] = lambdaiJi(p.d[d,], XB, iJi[[d]])
        observediJi[[d]] = LambdaiJi[[d]][edge[[d]][[1]]]
      }
    convergence[o] = EdgeTime(iJi[[maxedge2]], lambda[[maxedge2]], delta, LambdaiJi[[maxedge2]], timeinc[maxedge2], observediJi[[maxedge2]])
      }
 
  chain.final = list(C = currentC, Z = currentZ, B = bmat, D = deltamat, 
                     iJi = iJi, sigma_Q =sigma_Q, alpha = alphamat, mvec = mvecmat, edge2 = edge2,
                     proposal.var= proposal.var, convergence = convergence)
  return(chain.final)
}

#' @title IPTM_inference.Schein
#' @description Iterate Markov Chain Monte Carlo (MCMC) algorithm using Gibbs measure to sequentially update the assignments of Z, C and B
#'
#' @param edge list of document information with 3 elements (element 1 sender, element 2 receiver, element 3 time in unix.time format)
#' @param node nodelist containing the ID of nodes (ID starting from 1)
#' @param textlist list of text (length=number of documents in total) containing the words in each document
#' @param vocabulary all vocabularies used over the corpus
#' @param nIP total number of interaction patterns specified by the user
#' @param K total number of topics specified by the user
#' @param sigma_Q proposal distribution variance parameter for beta and delta
#' @param alpha Dirichlet concentration prior for document-topic distribution
#' @param mvec Dirichlet base prior for document-topic distribution
#' @param betas Dirichlet concentration prior for topic-word distribution
#' @param nvec Dirichlet base prior for topic-word distribution
#' @param prior.b.mean mean vector of b in multivariate normal distribution
#' @param prior.b.var covairance matrix of b in multivariate normal distribution
#' @param prior.delta parameter of delta in Normal prior
#' @param out size of outer iterations 
#' @param n_B size of third inner iteration for updates of B
#' @param n_d size of third inner iteration for updates of delta
#' @param burn iterations to be discarded at the beginning of beta chain
#' @param thinning the thinningning interval of beta chain
#' @param netstat which type of network statistics to use ("intercept", dyadic", "triadic", "degree")
#' @param optimize to optimize alpha (Dirichlet concentration prior for document-topic distribution) or not (TRUE/FALSE)
#' @param initial initial values to be used from Forward sampling
#'
#' @return MCMC output containing IP assignment, topic assignment, and (beta, mu, delta) chain 
#'
#' @export
IPTM_inference.Schein = function(edge, node, textlist, vocabulary, nIP, K, sigma_Q, alpha, mvec, betas, nvec, prior.b.mean, prior.b.var, prior.delta, 
					out, n_B, n_d, burn, thinning, netstat, optimize = FALSE, initial) {
   
  # trim the edge so that we only model edges after 384 hours
	timestamps = vapply(edge, function(d) {
  			  d[[3]]
 			  }, c(1))

    edge2 = which_int(384, timestamps) : length(edge)
    maxedge2 = max(edge2)
    timeinc = c(timestamps[1], timestamps[-1] - timestamps[-length(timestamps)])
   
  # initialize alpha, mvec, delta, nvec, delta, lvec, and gammas
 	W = length(vocabulary)
 	delta = initial$D
	beta.old = initial$B
	# initialize C, theta and Z
	currentC = initial$C
	currentZ = initial$Z
 	p.d = pdmat(currentZ, currentC, nIP) 

    # initialize beta
    netstat = as.numeric(c("intercept", "degree", "dyadic", "triadic" ) %in% netstat)
    L = 3
    P = netstat[1] + L * (2 * netstat[2] + 2 * netstat[3] + 4 * netstat[4])
    bmat = list()
	for (IP in 1:nIP) {
		bmat[[IP]] = matrix(beta.old[[IP]], nrow = P, ncol = (n_B - burn) / thinning)
  	}
  	deltamat = rep(delta, n_d)
    proposal.var = lapply(1:nIP, function(IP){diag(P)})

    #initialize the latent sender-receiver pairs
  iJi = lapply(seq(along = edge), function(d) {
    matrix(0, nrow = length(node), ncol = length(node))
  })
  lambda = list()
  LambdaiJi = list()
  observediJi = list()
  for (d in edge2) {
    iJi[[d]] = initial$iJi[[d]]
  }
  textlist.raw = unlist(textlist[edge2])

    #start outer iteration
    for (o in 1:out) {
      if (optimize) {
      #update the hyperparameter alpha and mvec
      vec = AlphamvecOpt(K, currentZ[edge2], alpha, mvec, 5)
      alpha = sum(vec)
      mvec = vec / alpha
      }
      # Data augmentation
    for (d in edge2) {
      history.t = History(edge, p.d, node, edge[[d-1]][[3]] + exp(-745))
      X = Netstats_cpp(history.t, node, netstat)
      XB = MultiplyXBList(X, beta.old)     
      lambda[[d]] = lambda_cpp(p.d[d,], XB)
      #calculate the resampling probability	
      for (i in node[-edge[[d]][[1]]]) {
      	XB_IP = lapply(XB, function(IP) {IP[i,]})
        for (j in sample(node[-i], length(node) - 1)) {
          probij = DataAug_cpp_Gibbs(iJi[[d]][i, ], lambda[[d]][i,], XB_IP, p.d[d, ], delta, timeinc[d], j)
          iJi[[d]][i, j] = multinom_vec(1, probij) - 1
        }
      }
      iJi[[d]][edge[[d]][[1]],] = tabulateC(as.numeric(unlist(edge[[d]][2])), length(node))
      LambdaiJi[[d]] = lambdaiJi(p.d[d,], XB, iJi[[d]])
      observediJi[[d]] = LambdaiJi[[d]][edge[[d]][[1]]]
    }
    	 	    	
     # Z update
     table.W = lapply(1:K, function(k) {
      			 tabulateC(textlist.raw[which(unlist(currentZ[edge2]) == k)], W)
      			 })
   	 for (d in edge2) {
       textlist.d = textlist[[d]]
       if (edge[[d]][[3]] + 384 > edge[[maxedge2]][[3]]) {
      	 hist.d = maxedge2
       } else {
      	 hist.d = which_num(edge[[d]][[3]] + 384, timestamps)
       }
       edgetime.d = rep(NA, K)
       for (w in 1:length(currentZ[[d]])) {
       	 zw.old = currentZ[[d]][w]
       	 if (length(textlist.d) > 0) {       	 
       	 table.W[[zw.old]][textlist.d[w]] = table.W[[zw.old]][textlist.d[w]] - 1
         topicpart.d = TopicInEqZ(K, currentZ[[d]][-w], alpha, mvec)
         wordpart.d = WordInEqZ(K, textlist.d, table.W, betas, nvec)
       } else {
         topicpart.d = 0
         wordpart.d = matrix(0, nrow = length(currentZ[[d]]), ncol = K)
       }
 	 for (IP in unique(currentC)) {
       	  	currentCK = which(currentC == IP)
       		currentZ[[d]][w] = min(currentCK)
       	 	p.d[d, ] = pdmat(list(currentZ[[d]]), currentC, nIP)           
            history.t = History(edge, p.d, node, edge[[hist.d-1]][[3]] + exp(-745))
    	    		X = Netstats_cpp(history.t, node, netstat)
    	    		XB = MultiplyXBList(X, beta.old)   
    	    		lambda[[hist.d]] = lambda_cpp(p.d[hist.d,], XB)
	   	 	LambdaiJi[[hist.d]] = lambdaiJi(p.d[hist.d, ], XB, iJi[[hist.d]])
        		observediJi[[hist.d]] = LambdaiJi[[hist.d]][edge[[hist.d]][[1]]]
        		edgetime.d[currentCK] = EdgeTime(iJi[[hist.d]], lambda[[hist.d]], delta, LambdaiJi[[hist.d]], timeinc[hist.d], observediJi[[hist.d]])
        }       
         const.Z = edgetime.d + topicpart.d + wordpart.d[w, ] 
         zw.new = multinom_vec(1, expconst(const.Z))
         if (zw.new != zw.old) {
           currentZ[[d]][w] = zw.new
           table.W[[zw.new]][textlist.d[w]] = table.W[[zw.new]][textlist.d[w]] + 1 
         } else {
         	currentZ[[d]][w] = zw.old
         	table.W[[zw.old]][textlist.d[w]] = table.W[[zw.old]][textlist.d[w]] + 1
         }
       }
     }

      # C update
       for (k in sort(unique(unlist(currentZ[edge2])))) { 
         const.C = rep(NA, nIP)
         for (IP in 1:nIP) {
          	currentC[k] = IP
          	p.d = pdmat(currentZ, currentC, nIP) 
 		 	history.t = History(edge, p.d, node, edge[[maxedge2-1]][[3]] + exp(-745))
    	   		X = Netstats_cpp(history.t, node, netstat)
    	   		XB = MultiplyXBList(X, beta.old)    
           	lambda[[maxedge2]] = lambda_cpp(p.d[maxedge2,], XB)
		    LambdaiJi[[maxedge2]] = lambdaiJi(p.d[maxedge2,], XB, iJi[[maxedge2]])
           	observediJi[[maxedge2]] = LambdaiJi[[maxedge2]][edge[[maxedge2]][[1]]]
           	prob = EdgeTime(iJi[[maxedge2]], lambda[[maxedge2]], delta, LambdaiJi[[maxedge2]], timeinc[maxedge2], observediJi[[maxedge2]])
           	const.C[IP] = prob
      	 }
         currentC[k] = multinom_vec(1, expconst(const.C))
	}
     
 	p.d = pdmat(currentZ, currentC, nIP) 
 
    for (d in maxedge2) {
        	history.t = History(edge, p.d, node, edge[[d-1]][[3]] + exp(-745))
    	    X = Netstats_cpp(history.t, node, netstat)
    	    XB = MultiplyXBList(X, beta.old)   
    	    lambda[[d]] = lambda_cpp(p.d[d,], XB)
	    	LambdaiJi[[d]] = lambdaiJi(p.d[d,], XB, iJi[[d]])
        	observediJi[[d]] = LambdaiJi[[d]][edge[[d]][[1]]]
	}
	 
	 # beta update
 	 prior.old1 = sum(vapply(1:nIP, function(IP) {rcpp_log_dmvnorm(prior.b.var, prior.b.mean, beta.old[[IP]], FALSE)}, c(1)))
 	 post.old1 = EdgeTime(iJi[[maxedge2]], lambda[[maxedge2]], delta, LambdaiJi[[maxedge2]], timeinc[maxedge2], observediJi[[maxedge2]])
   
      if (o != 1) {
      accept.rates = c(length(unique(bmat[[1]][1,])) / ncol(bmat[[1]]), length(unique(deltamat)) / n_d)
      sigma_Q = adaptive_MH(sigma_Q, accept.rates, update_size = 0.1 * sigma_Q)
      if (accept.rates[1] > 1 / ncol(bmat[[1]])) {
     		 for (IP in 1:nIP) {
      		 proposal.var[[IP]] = var(t(bmat[[IP]]))
       	 }
        }
     }
      for (i3 in 1:n_B) {
     		 beta.new = lapply(1:nIP, function(IP) {
           		    c(rcpp_rmvnorm(1, sigma_Q[1] * proposal.var[[IP]], beta.old[[IP]]))
          		    }) 
          for (d in maxedge2) {
             XB = MultiplyXBList(X, beta.new)
             lambda[[d]] = lambda_cpp(p.d[d,], XB)    
             LambdaiJi[[d]] = lambdaiJi(p.d[d,], XB, iJi[[d]])
 	         observediJi[[d]] = LambdaiJi[[d]][edge[[d]][[1]]]
          }
          prior.new1 = sum(vapply(1:nIP, function(IP) {rcpp_log_dmvnorm(prior.b.var, prior.b.mean, beta.new[[IP]], FALSE)}, c(1)))
          post.new1 = EdgeTime(iJi[[maxedge2]], lambda[[maxedge2]], delta, LambdaiJi[[maxedge2]], timeinc[maxedge2], observediJi[[maxedge2]])     			   
     	  loglike.diff = prior.new1 + post.new1 - prior.old1 - post.old1
          if (log(runif(1, 0, 1)) < loglike.diff) {
         	for (IP in 1:nIP) {
          		 beta.old[[IP]] = beta.new[[IP]]
          	 }
          	 prior.old1 = prior.new1
          	 post.old1 = post.new1
          }
           if (i3 > burn && i3 %% (thinning) == 0) {
         		 for (IP in 1:nIP) {
            		 bmat[[IP]][ , (i3 - burn) / thinning] = beta.old[[IP]]
            	 }
     		}
      }
     
     #delta update
      for (d in maxedge2) {
     	 XB = MultiplyXBList(X, beta.old)
       lambda[[d]] = lambda_cpp(p.d[d,], XB)  
      }
      prior.old2 = dnorm(delta, prior.delta[1], sqrt(prior.delta[2]), log = TRUE)
      post.old2 = EdgeInEqZ_Gibbs(iJi[[maxedge2]], lambda[[maxedge2]], delta)
 
  	 for (i4 in 1:n_d) {
          delta.new = rnorm(1, delta, sqrt(sigma_Q[2]))
          prior.new2 = dnorm(delta.new, prior.delta[1], sqrt(prior.delta[2]), log = TRUE)
          post.new2 = EdgeInEqZ_Gibbs(iJi[[maxedge2]], lambda[[maxedge2]], delta.new)
          loglike.diff2 = prior.new2 + post.new2 - prior.old2 - post.old2 
          if (log(runif(1, 0, 1)) < loglike.diff2) {
         	 delta = delta.new
          	 prior.old2 = prior.new2
             post.old2 = post.new2
          } 
             deltamat[i4] = delta
      }
 }
         
   chain.final = list(C = currentC, Z = currentZ[edge2], B = bmat, D = deltamat,
                     iJi = iJi, sigma_Q =sigma_Q, alpha = alpha, mvec = mvec, 
                     proposal.var= proposal.var)
  return(chain.final)
}


#' @title IPTM_inference.LDA
#' @description Iterate Markov Chain Monte Carlo (MCMC) algorithm using Gibbs measure to sequentially update the assignments of Z, C and B
#'
#' @param edge list of document information with 3 elements (element 1 sender, element 2 receiver, element 3 time in unix.time format)
#' @param node nodelist containing the ID of nodes (ID starting from 1)
#' @param textlist list of text (length=number of documents in total) containing the words in each document
#' @param vocabulary all vocabularies used over the corpus
#' @param nIP total number of interaction patterns specified by the user
#' @param K total number of topics specified by the user
#' @param sigma_Q proposal distribution variance parameter for beta and delta
#' @param alpha Dirichlet concentration prior for document-topic distribution
#' @param mvec Dirichlet base prior for document-topic distribution
#' @param betas Dirichlet concentration prior for topic-word distribution
#' @param nvec Dirichlet base prior for topic-word distribution
#' @param prior.b.mean mean vector of b in multivariate normal distribution
#' @param prior.b.var covairance matrix of b in multivariate normal distribution
#' @param prior.delta parameter of delta in Normal prior
#' @param out size of outer iterations 
#' @param n_B size of third inner iteration for updates of B
#' @param n_d size of third inner iteration for updates of delta
#' @param burn iterations to be discarded at the beginning of beta and delta chain 
#' @param thinning the thinningning interval of beta and delta chain
#' @param netstat which type of network statistics to use ("intercept", dyadic", "triadic", "degree")
#' @param optimize to optimize alpha (Dirichlet concentration prior for document-topic distribution) or not (TRUE/FALSE)
#'
#' @return MCMC output containing IP assignment, topic assignment, and (beta, mu, delta) chain 
#'
#' @export
IPTM_inference.LDA = function(edge, node, textlist, vocabulary, nIP, K, sigma_Q, alpha, mvec, betas, nvec, prior.b.mean, prior.b.var, prior.delta, 
                               out, n_B, n_d, burn, thinning, netstat, optimize = FALSE) {
  
  # trim the edge so that we only model edges after 384 hours
  timestamps = vapply(edge, function(d) {
    d[[3]]
  }, c(1))
  
  edge2 = which_int(384, timestamps) : length(edge)
  maxedge2 = max(edge2)
  timeinc = c(timestamps[1], timestamps[-1] - timestamps[-length(timestamps)])
  
  # initialize alpha, mvec, delta, nvec, delta, lvec, and gammas
  W = length(vocabulary)
  phi = lapply(1:K, function(k) {
    rdirichlet_cpp(1, betas * nvec)
  })
  delta = rnorm(1, prior.delta[1], sqrt(prior.delta[2]))
  beta.old = lapply(1:nIP, function(IP) {
    prior.b.mean
  })
  # initialize C, theta and Z
  currentC = sample(1:nIP, K, replace = TRUE)
  theta = rdirichlet_cpp(length(edge), alpha * mvec)
  currentZ = lapply(seq(along = edge), function(d) {
    multinom_vec(max(1, length(textlist[[d]])), theta[d, ])
  })
  p.d = pdmat(currentZ, currentC, nIP) 
  # initialize beta
  netstat = as.numeric(c("intercept", "degree", "dyadic", "triadic" ) %in% netstat)
  L = 3
  P = netstat[1] + L * (2 * netstat[2] + 2 * netstat[3] + 4 * netstat[4])
  
  bmat = list()
  for (IP in 1:nIP) {
    bmat[[IP]] = matrix(beta.old[[IP]], nrow = P, ncol = (n_B - burn[1]) / thinning[1])
  }
  deltamat = rep(delta,  (n_d - burn[2]) / thinning[2])
  proposal.var = lapply(1:nIP, function(IP){diag(P)})
  
  #initialize the latent sender-receiver pairs
  iJi = lapply(seq(along = edge), function(d) {
    matrix(0, nrow = length(node), ncol = length(node))
  })
  lambda = list()
  LambdaiJi = list()
  observediJi = list()
  for (d in edge2) {
    iJi[[d]] = matrix(rbinom(length(node)^2, 1, 0.5), nrow =length(node), ncol = length(node))
    diag(iJi[[d]]) = 0
  }
    textlist.raw = unlist(textlist)
  #start outer iteration
  for (o in 1:out) {
    print(o)
    if (optimize) {
      #update the hyperparameter alpha and mvec
      vec = AlphamvecOpt(K, currentZ[edge2], alpha, mvec, 5)
      alpha = sum(vec)
      mvec = vec / alpha
    }
     table.W = lapply(1:K, function(k) {
      tabulateC(textlist.raw[which(unlist(currentZ) == k)], W)
    })

    # Z update  
    for (d in 1:maxedge2) {
      textlist.d = textlist[[d]]        		 
      for (w in 1:length(currentZ[[d]])) {         	 	
        zw.old = currentZ[[d]][w]             	 	      	 	
        if (length(textlist.d) > 0) {
          table.W[[zw.old]][textlist.d[w]] = table.W[[zw.old]][textlist.d[w]] - 1
          topicpart.d = TopicInEqZ(K, currentZ[[d]][-w], alpha, mvec)
          wordpart.d = WordInEqZ(K, textlist.d, table.W, betas, nvec)
        } else {
          topicpart.d = 0
          wordpart.d = matrix(0, nrow = length(currentZ[[d]]), ncol = K)
        }
        const.Z = topicpart.d + wordpart.d[w, ]
        zw.new = multinom_vec(1, expconst(const.Z))
        if (zw.new != zw.old) {
           currentZ[[d]][w] = zw.new
           table.W[[zw.new]][textlist.d[w]] = table.W[[zw.new]][textlist.d[w]] + 1 
        } else {
        	   table.W[[zw.old]][textlist.d[w]] = table.W[[zw.old]][textlist.d[w]] + 1
        }
      }
    }
  }
  chain.final = list(Z = currentZ, alpha = alpha, mvec = mvec, edge2 = edge2)
  return(chain.final)
}

#' @title TablebetaIP
#' @description Generate a table summary of the MCMC chain of network statistics coefficients (beta) for each interaction pattern
#'
#' @param MCMCchain a chain obtained using MCMC function
#'
#' @return List of table containing the posterior summary of beta for each interaction pattern
#'
#' @export
TablebetaIP = function(MCMCchain) {
  # Generate a table summary of the MCMC chain of beta for each IP
  #
  # Args 
  #  MCMCchain a chain obtained using MCMC function
  #
  # Returns
  #  List of table containing the posterior summary of beta for each IP
  
  table.beta = lapply(MCMCchain$B, function(x) {
    summary(mcmc(t(x)))$quantiles[, c(3, 1, 5)]
  })
  return(table.beta)
}

#' @title PlotbetaIP
#' @description Draw a boxplot of the MCMC chain of network statistics (beta) for each interaction pattern
#'
#' @param Bchain a chain of B obtained using MCMC function
#'
#' @return Joint boxplot of the posterior summary of beta (should click for the locator)
#'
#' @export
PlotbetaIP = function(Bchain) {
  # Draw a boxplot of the MCMC chain of beta for each IP
  #
  # Args 
  #  Bchain a chain of B obtained using MCMC function
  #
  # Returns
  #  Joint boxplot of the posterior summary of beta (should click for the locator)
  combined.data = list()
  P = nrow(Bchain[[1]])
  nIP = length(Bchain)
  for(b in 1:P){
    combined.data[[b]] = sapply(1:nIP, function(c){
      cbind(Bchain[[c]][b,])
    })
  }
  forbox = melt(combined.data)
  boxplot = boxplot(forbox$value ~ forbox$X2 + forbox$L1,
                    at = c(sapply(1:P, function(x){
                      ((nIP + 1) * x - nIP):((nIP + 1) * x - 1)
                    })),
                    col = gray.colors(nIP), axes = FALSE, 
                    main = "Comparison of beta coefficients for different IPs")
  abline(h = 0, lty = 1, col = "red")
  axis(2, labels = TRUE)
  box()
  axis(side = 1, line = 0.5, lwd = 0, 
       at = c(sapply(1:P, function(x){
         median(((nIP + 1) * x - nIP):((nIP + 1) * x - 1))
       })), labels = 1:25)
  legend(locator(1), c(paste("IP", 1:nIP)), col = gray.colors(nIP), pch = 15)
}

#' @title PlotTopic
#' @description Draw a barplot of the topic distributions without considering interaction patterns
#'
#' @param Zchain summary of Z obtained using MCMC function
#' @param K total number of topics specified by the user
#'
#' @return Barplot of the topic distribution
#'
#' @export	
PlotTopic = function(Zchain, K) {
  # Draw a barplot of the topic distributions without considering IPs
  #
  # Args 
  #  MCMCchain a chain obtained using MCMC function
  #  K total number of topics specified by the user
  #
  # Returns
  #  Barplot of the topic distribution
  Zsummary = list()
  for (d in seq(along = Zchain)) {
    Zsummary[[d]] = Zchain[[d]]
  }
  topic.dist = t(tabulateC(unlist(Zsummary), K) / length(unlist(Zsummary)))
  colnames(topic.dist) = c(1:K)
  barplot(topic.dist, col = gray.colors(K), beside = TRUE, xlab = "Topic", ylab = "Proportion", 
          main = "Topic Distritubitions")
}

#' @title TableWord
#' @description Generate a table that summarizes token-topic assignments with high probabilities
#'
#' @param Zchain summary of Z obtained using MCMC function
#' @param K total number of topics specified by the user
#' @param textlist list of text (length=number of documents in total) containing the words in each document
#' @param vocabulary all vocabularies used over the corpus
#'
#' @return List of table that summarize token-topic assignments with high probabilities (topic proportion > 0.1) for each IP
#'
#' @export
TableWord = function(Zchain, K, textlist, vocabulary) {
  # Generate a table of token-topic assignments with high probabilities for each IP
  #
  # Args 
  #  Zchain summary of Z obtained using MCMC function
  #  K total number of topics specified by the user
  #  textlist list of text containing the words in each document
  #  vocabulary all vocabularies used over the corpus
  #
  # Returns
  #  List of table that summarize token-topic assignments for each IP
  W = length(vocabulary)
    Zsummary = list()
    topic.word = matrix(0, nrow = K, ncol = W)
    colnames(topic.word) = vocabulary
    iter = 1
    for (d in seq(along = textlist)) {
      if (length(Zchain[[d]]) > 0){
        Zsummary[[iter]] = Zchain[[d]]
        names(Zsummary[[iter]])<- vocabulary[textlist[[d]]]
        iter = iter+1
      }
    }
    topic.dist = t(tabulateC(unlist(Zsummary), K)/length(unlist(Zsummary)))
    colnames(topic.dist) = c(1:K)
    top.topic = topic.dist[, order(topic.dist, decreasing = TRUE)]
    all.word = unlist(Zsummary)
    for (i in seq(along = all.word)){
      matchWZ = which(colnames(topic.word) == names(all.word[i]))
      topic.word[all.word[i], matchWZ] = topic.word[all.word[i], matchWZ] + 1
    }
    table.word = top.topic.words(topic.word, num.words = 10)
    colnames(table.word) = names(top.topic)
  return(table.word)
}

#' @title GenerateDocs.Gibbs
#' @description Generate a collection of documents according to the generative process of IPTM using Gibbs measure
#'
#' @param nDocs number of documents to be generated
#' @param node nodelist containing the ID of nodes (ID starting from 1)
#' @param vocabulary all vocabularies used over the corpus
#' @param nIP total number of interaction patterns specified by the user
#' @param K total number of topics specified by the user
#' @param nwords number of words in a document (fixed constant for now)
#' @param alpha Dirichlet concentration prior for document-topic distribution
#' @param mvec Dirichlet base prior for document-topic distribution
#' @param betas Dirichlet concentration prior for topic-word distribution
#' @param nvec Dirichlet base prior for topic-word distribution
#' @param b List of coefficients estimating the history effect for each IP 
#' @param delta tuning parameter controlling the number of recipients
#' @param currentC topic-IP assignment
#' @param netstat which type of network statistics to use ("intercept", "dyadic", "triadic", "degree")
#' @param base.edge edges before 384 hours that is used to calculate initial history of interactions
#' @param base.text texts corresponding to base.edge
#' @param topic_token_assignments matrix of topic-token assignments
#' @param forward Logigal indicating whether we are generating forward samples
#' @param backward_init Logigal indicating whether we are generating initial backward samples
#' @param backward Logigal indicating whether we are generating backward samples
#' @param base Logical indicating whether or not we are generating base edges (< 384)
#' @param support support of latent edge vector Ji
#'
#' @return generated edge and text, parameter b used to generate those, and base (if base == TRUE)
#'
#' @export
GenerateDocs.Gibbs = function(nDocs, node, vocabulary, nIP, K, nwords, alpha, mvec, betas, nvec,
                        b, delta, currentC, netstat, base.edge, base.text, topic_token_assignments = NULL,
                        forward = FALSE, backward_init = FALSE, backward = FALSE, base = FALSE, support) {
  W = length(vocabulary)
  phi = lapply(1:K, function(k) {
    rdirichlet_cpp(1, betas * nvec)
  })
    netstat = as.numeric(c("intercept", "degree", "dyadic", "triadic" ) %in% netstat)
    L = 3
    P = netstat[1] + L * (2 * netstat[2] + 2 * netstat[3] + 4 * netstat[4])
  if (base) {
    t.d = 0
  } else {
    t.d = 384
  }
  edge = base.edge
  text = base.text
  base.length = length(edge)
  p.d = matrix(NA, nrow = nDocs + base.length, ncol = nIP)
  if (!base) {
  for (d in 1:base.length) {
  	 p.d[d, ] = pdmat(list(as.integer(names(text[[d]]))), currentC, nIP)
  }
  }
 history.t = lapply(1:nIP, function(IP) {
        	 	 lapply(1:3, function(l){
         		matrix(0, length(node), length(node))
       			})
     		 })  
  word_type_topic_counts = matrix(0, W, K)
  iJi = matrix(0, length(node), length(node))
  iJis = list()
  for (d in 1:nDocs) {
    N.d = nwords
    text[[base.length + d]] = rep(NA, N.d)
    
    if (!backward) {
      theta.d = rdirichlet_cpp(1, alpha * mvec)	
      topic.d = multinom_vec(max(1, N.d), theta.d)
        for (n in 1:N.d){
          text[[base.length + d]][n] = multinom_vec(1, phi[[topic.d[n]]])
        }
        names(text[[base.length + d]]) = topic.d
    } else {
      phi.k = rep(NA, K)
      topic.d = topic_token_assignments[[d]]
        for (n in 1:N.d){
          for (w in 1:W) {
            phi.k[w] = (word_type_topic_counts[w, topic.d[n]] + betas * nvec[w]) / (sum(word_type_topic_counts[,topic.d[n]]) + betas)
          } 
          text[[base.length + d]][n] = multinom_vec(1, phi.k)
          word_type_topic_counts[text[[base.length + d]][n], topic.d[n]] = word_type_topic_counts[text[[base.length + d]][n], topic.d[n]] + 1
        }
        names(text[[base.length + d]]) = topic.d
    }
    
    p.d[base.length + d, ] = pdmat(list(as.integer(names(text[[base.length + d]]))), currentC, nIP)
    if (t.d >= 384) {
      history.t = History(edge, p.d, node, t.d + exp(-745))
    }
    X = Netstats_cpp(history.t, node, netstat)
    XB = MultiplyXBList(X, b)     
    lambda = lambda_cpp(p.d[base.length + d,], XB)
    	for (i in node) {
    		iJi[i, -i] = r.gibbs.measure(1, lambda[i, -i], delta, support)
    	}
    	iJis[[d]] = iJi
    LambdaiJi = lambdaiJi(p.d[base.length + d,], XB, iJi)
    i.d = multinom_vec(1, LambdaiJi)
    j.d = which(iJi[i.d,] == 1)
    t.d = t.d + rexp(1, sum(LambdaiJi))
    edge[[base.length + d]] = list(sender = i.d, receiver = j.d, timestamp = t.d)		
  }
  if (forward) {
    edge = edge[-(1:base.length)]
    text = text[-(1:base.length)]
  }
   if (base == TRUE && t.d > 384) {
   cutoff = which_int(384, vapply(1:length(edge), function(d) {edge[[d]][[3]]}, c(1))) - 1
    edge = edge[1:cutoff]
    text = text[1:cutoff]
   }
  return(list(edge = edge, text = text, base = base.length, b = b, d = delta, X = X, iJis = iJis))							
} 

#' @title GenerateDocs.Schein
#' @description Generate a collection of documents according to the generative process of IPTM using Gibbs measure
#'
#' @param nDocs number of documents to be generated
#' @param node nodelist containing the ID of nodes (ID starting from 1)
#' @param vocabulary all vocabularies used over the corpus
#' @param nIP total number of interaction patterns specified by the user
#' @param K total number of topics specified by the user
#' @param nwords number of words in a document (fixed constant for now)
#' @param alpha Dirichlet concentration prior for document-topic distribution
#' @param mvec Dirichlet base prior for document-topic distribution
#' @param betas Dirichlet concentration prior for topic-word distribution
#' @param nvec Dirichlet base prior for topic-word distribution
#' @param b List of coefficients estimating the history effect for each IP 
#' @param delta tuning parameter controlling the number of recipients
#' @param currentC topic-IP assignment
#' @param netstat which type of network statistics to use ("intercept", "dyadic", "triadic", "degree")
#' @param base.edge edges before 384 hours that is used to calculate initial history of interactions
#' @param base.text texts corresponding to base.edge
#' @param topic_token_assignments matrix of topic-token assignments
#' @param forward Logigal indicating whether we are generating forward samples
#' @param backward_init Logigal indicating whether we are generating initial backward samples
#' @param backward Logigal indicating whether we are generating backward samples
#' @param base Logical indicating whether or not we are generating base edges (< 384)
#' @param support support of latent edge vector Ji
#'
#' @return generated edge and text, parameter b used to generate those, and base (if base == TRUE)
#'
#' @export
GenerateDocs.Schein = function(nDocs, node, vocabulary, nIP, K, nwords, alpha, mvec, betas, nvec,
                        b, delta, currentC, netstat, base.edge, base.text, topic_token_assignments = NULL,
                        forward = FALSE, backward_init = FALSE, backward = FALSE, base = FALSE, support) {
  W = length(vocabulary)
  phi = lapply(1:K, function(k) {
    rdirichlet_cpp(1, betas * nvec)
  })
    netstat = as.numeric(c("intercept", "degree", "dyadic", "triadic" ) %in% netstat)
    L = 3
    P = netstat[1] + L * (2 * netstat[2] + 2 * netstat[3] + 4 * netstat[4])
  
  if (base) {
    t.d = 0
  } else {
    t.d = 384
  }
  edge = base.edge
  text = base.text
  base.length = length(edge)
  p.d = matrix(NA, nrow = nDocs + base.length, ncol = nIP)
  if (!base) {
  for (d in 1:base.length) {
  	 p.d[d, ] = pdmat(list(as.integer(names(text[[d]]))), currentC, nIP)
  	   }
  }	   
 history.t = lapply(1:nIP, function(IP) {
        	 	 lapply(1:3, function(l){
         		matrix(0, length(node), length(node))
       			})
     		 })  		  
  word_type_topic_counts = matrix(0, W, K)
  iJi = list()
  for (d in 1:nDocs) {
    N.d = nwords
    text[[base.length + d]] = rep(NA, N.d)
    iJi[[base.length + d]] = matrix(0, length(node), length(node))
    if (!backward) {
      theta.d = rdirichlet_cpp(1, alpha * mvec)	
      topic.d = multinom_vec(max(1, N.d), theta.d)
        for (n in 1:N.d){
          text[[base.length + d]][n] = multinom_vec(1, phi[[topic.d[n]]])
        }
        names(text[[base.length + d]]) = topic.d
    } else {
      phi.k = rep(NA, K)
      topic.d = topic_token_assignments[[d]]
        for (n in 1:N.d){
          for (w in 1:W) {
            phi.k[w] = (word_type_topic_counts[w, topic.d[n]] + betas * nvec[w]) / (sum(word_type_topic_counts[, topic.d[n]]) + betas)
          } 
          text[[base.length + d]][n] = multinom_vec(1, phi.k)
          word_type_topic_counts[text[[base.length + d]][n], topic.d[n]] = word_type_topic_counts[text[[base.length + d]][n], topic.d[n]] + 1
        }
        names(text[[base.length + d]]) = topic.d
    }
    
    p.d[base.length + d, ] = pdmat(list(as.integer(names(text[[base.length + d]]))), currentC, nIP)
     if (t.d >= 384) {
      history.t = History(edge, p.d, node, t.d + exp(-745))
    }
    X = Netstats_cpp(history.t, node, netstat)
    XB = MultiplyXBList(X, b)     
    lambda = lambda_cpp(p.d[base.length + d,], XB)
    
    	for (i in node) {
    		iJi[[base.length + d]][i, -i] = r.gibbs.measure(1, lambda[i, -i], delta, support)
    	}
    LambdaiJi = lambdaiJi(p.d[base.length + d,], XB, iJi[[base.length + d]])
    i.d = multinom_vec(1, LambdaiJi)
    j.d = which(iJi[[base.length + d]][i.d,] == 1)
    t.d = t.d + rexp(1, sum(LambdaiJi))
    edge[[base.length + d]] = list(sender = i.d, receiver = j.d, timestamp = t.d)		
  }
  if (base == TRUE && t.d > 384) {
    cutoff = which_int(384, vapply(1:length(edge), function(d) {edge[[d]][[3]]}, c(1))) - 1
    edge = edge[1:cutoff]
    text = text[1:cutoff]
  }
  return(list(edge = edge, text = text, base = base.length, b = b, d = delta, X = X, iJi = iJi))							
} 


#' @title GiR_stats
#' @description Calculate several statistics from samples generated from forward or backward sampling
#'
#' @param GiR_sample one sample from generative process
#' @param K total number of topics specified by the user
#' @param currentC topic-IP assignment
#' @param vocabulary all vocabularies used over the corpus
#' @param forward Logigal indicating whether the sample is from forward sampling
#' @param backward Logigal indicating whether the sample is from backward sampling
#'
#' @return A vector of statistics calculated from one GiR sample
#'
#' @export

GiR_stats = function(GiR_sample, K, currentC, vocabulary, forward = FALSE, backward = FALSE) {
  edge = GiR_sample$edge
  text = GiR_sample$text
  if (backward) {
    edge = edge[-(1:GiR_sample$base)]
    text = text[-(1:GiR_sample$base)]
  }
  
  GiR_stats = c()
  nDocs = length(edge)
  P = length(GiR_sample$b[[1]])
  nIP = length(GiR_sample$b)
  nwords = length(GiR_sample$text[[1]])
  W = length(vocabulary)
  A = nrow(GiR_sample$X[[1]][[1]])
  
  GiR_stats[1:P] = Reduce('+', GiR_sample$b) / nIP 
  GiR_stats[(P + 1):(P + 3)] = nIP * Reduce('+', lapply(1:A, function(i) {colMeans(GiR_sample$X[[i]][[1]][,2:4])})) / A
  GiR_stats[P + 4] = GiR_sample$d
  GiR_stats[P + 5] = mean(vapply(1:nDocs, function(d) {
    length(edge[[d]][[2]])
  }, c(1)))
  GiR_stats[P + 6] = mean(c(edge[[1]][[3]] - 384, vapply(2:nDocs, function(d) {
    edge[[d]][[3]] - edge[[d-1]][[3]]
  }, c(1)))) 			
  GiR_stats[P + 7] = mean(currentC)
  Tokens_in_Topic = tabulate(vapply(1:nDocs, function(d){
    as.numeric(names(text[[d]]))
  }, rep(0, nwords)), K)
  GiR_stats[(P + 8):(P + 7 + nIP)] = vapply(1:nIP, function(IP) {
    Tokens_in_Topic %*% (currentC == IP)
  }, c(1))
  GiR_stats[(P + 8 + nIP):(P + 7 + nIP + K)] = Tokens_in_Topic
  Tokens_in_Word = tabulate(unlist(text), W)
  GiR_stats[(P + 8 + nIP + K):(P + 7 + nIP + K + W)] = Tokens_in_Word
  
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
    if (grepl("B_", nms[i]) ) {
      quantiles = 1000
    }
    if (grepl("delta", nms[i]) ) {
      quantiles = 1000
    }
    if (grepl("Mean_timediff", nms[i]) ) {
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
             	

#' @title GiR.Gibbs
#' @description Getting it Right test for IPTM
#'
#' @param Nsamp number of GiR samples to be generated 
#' @param nDocs number of documents to be generated per each sample
#' @param node nodelist containing the ID of nodes (ID starting from 1)
#' @param vocabulary all vocabularies used over the corpus
#' @param nIP total number of interaction patterns specified by the user
#' @param K total number of topics specified by the user
#' @param nwords number of words in a document (fixed constant for now)
#' @param alpha Dirichlet concentration prior for document-topic distribution
#' @param mvec Dirichlet base prior for document-topic distribution
#' @param betas Dirichlet concentration prior for topic-word distribution
#' @param nvec Dirichlet base prior for topic-word distribution
#' @param prior.b.mean mean vector of b in multivariate normal distribution
#' @param prior.b.var covairance matrix of b in multivariate normal distribution
#' @param prior.delta parameter of delta in Normal prior
#' @param sigma_Q proposal distribution variance parameter for beta and delta
#' @param niters (out, n_B, n_d, burn, thinning) in the inference
#' @param netstat which type of network statistics to use ("dyadic", "triadic", "degree")
#' @param base.edge artificial collection of documents to be used as initial state of history
#' @param base.text artificial collection of documents to be used as initial state of history
#' @param generate_PP_plots Logical indicating whether to draw PP plots
#'
#' @return Forward and Backward samples and corresponding test results
#'
#' @export
GiR.Gibbs = function(Nsamp, nDocs, node, vocabulary, nIP, K, nwords, alpha, mvec, betas, nvec, 
               prior.b.mean, prior.b.var, prior.delta, sigma_Q, niters, netstat = c("dyadic"), 
               base.edge, base.text, generate_PP_plots = TRUE) {
  
    netstat2 = as.numeric(c("intercept", "degree", "dyadic", "triadic" ) %in% netstat)
    L = 3
    P = netstat2[1] + L * (2 * netstat2[2] + 2 * netstat2[3] + 4 * netstat2[4])
  supportD = gibbs.measure.support(length(node) - 1)

  #Forward sampling
  Forward_stats = matrix(NA, nrow = Nsamp, ncol = P + (P - 1) / 2 + 4 + nIP + K + length(vocabulary))
  colnames(Forward_stats) = c(paste0("B_",1:P), paste0("Send_",1:3), "delta", "Mean_recipients", "Mean_timediff", 
  							  "Mean_TopicIP", paste0("Tokens_in_IP_", 1:nIP), paste0("Tokens_in_Topic", 1:K), paste0("Tokens_in_Word",1:length(vocabulary)))
  for (i in 1:Nsamp) { 
    if (i %% 5000 == 0) {cat("Forward sampling", i, "\n")}
    b = lapply(1:nIP, function(IP) {
      c(rcpp_rmvnorm(1, prior.b.var, prior.b.mean))
    })
    delta = rnorm(1, prior.delta[1], sqrt(prior.delta[2]))
    currentC = sample(1:nIP, K, replace = TRUE)
    Forward_sample = GenerateDocs.Gibbs(nDocs, node, vocabulary, nIP, K, nwords, alpha, mvec, betas, nvec, b,
    		 								delta, currentC, netstat, base.edge = base.edge, base.text = base.text,
    		 								forward = TRUE, base = FALSE, support = supportD) 
    Forward_stats[i, ] = GiR_stats(Forward_sample, K, currentC, vocabulary, forward = TRUE, backward = FALSE)
    }
  
  #Backward sampling
  Backward_stats = matrix(NA, nrow = Nsamp, ncol = ncol(Forward_stats))
  Backward_sample = GenerateDocs.Gibbs(nDocs, node, vocabulary, nIP, K, nwords, alpha, mvec, betas, nvec, b, 
  									   delta, currentC, netstat, base.edge = base.edge, base.text = base.text,
  									   backward_init = TRUE, forward = FALSE, backward = FALSE, base = FALSE, support = supportD) 
  for (i in 1:Nsamp) { 
     if (i %% 100 == 0) {cat("Backward Sampling", i, "\n")}
    Inference_samp = IPTM_inference.Gibbs(Backward_sample$edge, node, Backward_sample$text, vocabulary, nIP, K,
    										  sigma_Q, alpha, mvec, betas, nvec, prior.b.mean, prior.b.var, prior.delta,
                               			  out = niters[1], n_B = niters[2], n_d = niters[3], burn = niters[4], 
                               			  thinning = niters[5], netstat)
    b = lapply(1:nIP, function(IP) {
        Inference_samp$B[[IP]][,ncol(Inference_samp$B[[IP]])]
    })
    delta = Inference_samp$D[length(Inference_samp$D)]
    currentC = Inference_samp$C
    topic_token_assignments = Inference_samp$Z
    for (d in 1:length(topic_token_assignments)) {
      names(topic_token_assignments[[d]]) = Backward_sample$text[[d + length(base.text)]]
    }
    
    Backward_sample = GenerateDocs.Gibbs(nDocs, node, vocabulary, nIP, K, nwords, alpha, mvec, betas, nvec, b, 
    										 delta, currentC, netstat, base.edge = base.edge, base.text = base.text,
    										 topic_token_assignments = topic_token_assignments, 
                             			 forward = FALSE, backward = TRUE, base = FALSE, support = supportD)
    Backward_stats[i, ] = GiR_stats(Backward_sample, K, currentC, vocabulary, forward = FALSE, backward = TRUE)
    }
  				
  if (generate_PP_plots) {
    par(mfrow=c(5,5), oma = c(3,3,3,3), mar = c(2,1,1,1))
    GiR_PP_Plots(Forward_stats, Backward_stats)
  }			
  return(list(Forward = Forward_stats, Backward = Backward_stats))
}                         	          


#' @title Schein.Gibbs
#' @description Getting it Right test for IPTM
#'
#' @param Nsamp number of GiR samples to be generated 
#' @param nDocs number of documents to be generated per each sample
#' @param node nodelist containing the ID of nodes (ID starting from 1)
#' @param vocabulary all vocabularies used over the corpus
#' @param nIP total number of interaction patterns specified by the user
#' @param K total number of topics specified by the user
#' @param nwords number of words in a document (fixed constant for now)
#' @param alpha Dirichlet concentration prior for document-topic distribution
#' @param mvec Dirichlet base prior for document-topic distribution
#' @param betas Dirichlet concentration prior for topic-word distribution
#' @param nvec Dirichlet base prior for topic-word distribution
#' @param prior.b.mean mean vector of b in multivariate normal distribution
#' @param prior.b.var covairance matrix of b in multivariate normal distribution
#' @param prior.delta parameter of delta in Normal prior
#' @param sigma_Q proposal distribution variance parameter for beta and delta
#' @param niters (out, n_B, n_d, burn, thinning) in the inference
#' @param netstat which type of network statistics to use ("dyadic", "triadic", "degree")
#' @param generate_PP_plots Logical indicating whether to draw PP plots
#'
#' @return Forward and Backward samples and corresponding test results
#'
#' @export
Schein.Gibbs = function(Nsamp, nDocs, node, vocabulary, nIP, K, nwords, alpha, mvec, betas, nvec, 
               prior.b.mean, prior.b.var, prior.delta, sigma_Q, niters, netstat = c("dyadic"), 
               generate_PP_plots = TRUE) {
  
  netstat2 = as.numeric(c("intercept", "degree", "dyadic", "triadic" ) %in% netstat)
  L = 3
  P = netstat2[1] + L * (2 * netstat2[2] + 2 * netstat2[3] + 4 * netstat2[4])
  supportD = gibbs.measure.support(length(node) - 1)

  #Forward sampling
  Forward_stats = matrix(NA, nrow = Nsamp, ncol = P + (P - 1) / 2 + 4 + nIP + K + length(vocabulary))
  colnames(Forward_stats) = c(paste0("B_",1:P), paste0("Send_",1:3), "delta", "Mean_recipients", "Mean_timediff", 
  							  "Mean_TopicIP", paste0("Tokens_in_IP_", 1:nIP), paste0("Tokens_in_Topic", 1:K), paste0("Tokens_in_Word",1:length(vocabulary)))
  Backward_stats = matrix(NA, nrow = Nsamp, ncol = ncol(Forward_stats))

  for (i in 1:Nsamp) { 
    b = lapply(1:nIP, function(IP) {
      c(rcpp_rmvnorm(1, prior.b.var, prior.b.mean))
    })
    delta = rnorm(1, prior.delta[1], sqrt(prior.delta[2]))
    currentC = sample(1:nIP, K, replace = TRUE)
    base.data = GenerateDocs.Gibbs(200, node, vocabulary, nIP, K, nwords, alpha, mvec, betas, nvec, b, delta, currentC, netstat, base.edge = list(),  base.text = list(), base = TRUE, support = supportD) 
	base.edge = base.data$edge	   
	base.text = base.data$text
    Forward_sample = GenerateDocs.Schein(nDocs, node, vocabulary, nIP, K, nwords, alpha, mvec, betas, nvec, b,
    		 								delta, currentC, netstat, base.edge = base.edge, base.text = base.text,
    		 								forward = TRUE, base = FALSE, support = supportD) 
    Forward_stats[i, ] = GiR_stats(Forward_sample, K, currentC, vocabulary, forward = FALSE, backward = TRUE)

    if (i %% 100 == 0) {cat("Sampling", i, "\n")}
    initial = list(C = currentC, D = delta, B = b, Z = lapply(Forward_sample$text, function(d){as.numeric(names(d))}),
    				iJi = Forward_sample$iJi)
    Inference_samp = IPTM_inference.Schein(Forward_sample$edge, node, Forward_sample$text, vocabulary, nIP, K,
    										  sigma_Q, alpha, mvec, betas, nvec, prior.b.mean, prior.b.var, prior.delta,
                               			  out = niters[1], n_B = niters[2], n_d = niters[3], burn = niters[4], 
                               			  thinning = niters[5], netstat, initial = initial)
    b = lapply(1:nIP, function(IP) {
        Inference_samp$B[[IP]][,ncol(Inference_samp$B[[IP]])]
    })
    delta = Inference_samp$D[length(Inference_samp$D)]
    currentC = Inference_samp$C
    topic_token_assignments = Inference_samp$Z

    Backward_sample = GenerateDocs.Schein(nDocs, node, vocabulary, nIP, K, nwords, alpha, mvec, betas, nvec, b, 
    										 delta, currentC, netstat, base.edge = base.edge, base.text = base.text,
    										 topic_token_assignments = topic_token_assignments, 
                             			 forward = FALSE, backward = TRUE, base = FALSE, support = supportD)
    Backward_stats[i, ] = GiR_stats(Backward_sample, K, currentC, vocabulary, forward = FALSE, backward = TRUE)
    
  }
  				
  if (generate_PP_plots) {
    par(mfrow=c(5,5), oma = c(2,2,2,2), mar = c(1,1,1,1))
    GiR_PP_Plots(Forward_stats, Backward_stats)
  }			
  return(list(Forward = Forward_stats, Backward = Backward_stats))
}                   


#' @title GenerateDocs.PPC
#' @description Generate a collection of documents according to the generative process of IPTM using Gibbs measure
#'
#' @param nDocs number of documents to be generated
#' @param node nodelist containing the ID of nodes (ID starting from 1)
#' @param vocabulary all vocabularies used over the corpus
#' @param nIP total number of interaction patterns specified by the user
#' @param K total number of topics specified by the user
#' @param alpha Dirichlet concentration prior for document-topic distribution
#' @param mvec Dirichlet base prior for document-topic distribution
#' @param betas Dirichlet concentration prior for topic-word distribution
#' @param nvec Dirichlet base prior for topic-word distribution
#' @param iJi inferred latent receivers
#' @param b List of coefficients estimating the history effect for each IP 
#' @param delta tuning parameter controlling the number of recipients
#' @param currentC topic-IP assignment
#' @param netstat which type of network statistics to use ("intercept", "dyadic", "triadic", "degree")
#' @param base.edge edges before 384 hours that is used to calculate initial history of interactions
#' @param base.text texts corresponding to base.edge
#' @param currentZ topic-token assignments
#' @param R number of iterations to make the generated data considered independent of observed data
#'
#' @return generated edge and text, parameter b used to generate those, and base (if base == TRUE)
#'
#' @export
GenerateDocs.PPC = function(nDocs, node, vocabulary, nIP, K, alpha, mvec, betas, nvec, iJi,
                        b, delta, currentC, netstat, base.edge, base.text, currentZ, R) {
  W = length(vocabulary)
  phi = lapply(1:K, function(k) {
    rdirichlet_cpp(1, betas * nvec)
  })
  netstat = as.numeric(c("intercept", "degree", "dyadic", "triadic") %in% netstat)
  L = 3
  P = netstat[1] + L * (2 * netstat[2] + 2 * netstat[3] + 4 * netstat[4])
  
  edge = base.edge
  text = base.text
  base.length = length(edge)
  lambda = list()
  t.d = base.edge[[base.length]][[3]] 
  p.d = pdmat(currentZ, currentC, nIP)  		  
  word_type_topic_counts = matrix(0, W, K)
  maxedge2 = length(currentZ)
  for (d in 1:nDocs) {
    N.d = length(currentZ[[base.length + d]])
    text[[base.length + d]] = rep(NA, N.d)
      phi.k = rep(NA, K)
      topic.d = currentZ[[base.length + d]]
        for (n in 1:N.d){
          for (w in 1:W) {
            phi.k[w] = (word_type_topic_counts[w, topic.d[n]] + betas * nvec[w]) / (sum(word_type_topic_counts[, topic.d[n]]) + betas)
          } 
          text[[base.length + d]][n] = multinom_vec(1, phi.k)
          word_type_topic_counts[text[[base.length + d]][n], topic.d[n]] = word_type_topic_counts[text[[base.length + d]][n], topic.d[n]] + 1
        }
        names(text[[base.length + d]]) = topic.d

    history.t = History(edge, p.d, node, t.d + exp(-745))
    X = Netstats_cpp(history.t, node, netstat)
    XB = MultiplyXBList(X, b)     
    lambda[[d]] = lambda_cpp(p.d[base.length + d,], XB)
    LambdaiJi = lambdaiJi(p.d[base.length + d,], XB, iJi[[base.length + d]])
    i.d = multinom_vec(1, LambdaiJi)
    j.d = which(iJi[[base.length + d]][i.d,] == 1)
    t.d = t.d + rexp(1, sum(LambdaiJi))
    edge[[base.length + d]] = list(sender = i.d, receiver = j.d, timestamp = t.d)		
  }
  for (r in 1:R) {
  	timestamps = vapply(edge, function(d) {
  			  d[[3]]
 			  }, c(1))
    timeinc = c(timestamps[1], timestamps[-1] - timestamps[-length(timestamps)])
    textlist.raw = unlist(text[-1:-base.length])
	table.W = lapply(1:K, function(k) {
      			 tabulateC(textlist.raw[which(unlist(currentZ[-1:-base.length]) == k)], W)
      			 })
  	for (d in 1:nDocs) {
   	 #iJi update
   		for (i in node) {
   			XB_IP = lapply(XB, function(IP) {IP[i,]})
		for (j in sample(node[-i], length(node) - 1)) {
			 probij = DataAug_cpp_Gibbs(iJi[[base.length + d]][i, ], lambda[[d]][i,], XB_IP, p.d[base.length + d, ], delta, timeinc[base.length + d], j) 
			 print(probij)
			 iJi[[base.length + d]][i, j] = multinom_vec(1, probij) - 1		
			}
      	}
     #Z update 	
        textlist.d = text[[base.length + d]]
       if (edge[[base.length + d]][[3]] + 384 > edge[[maxedge2]][[3]]) {
      	 hist.d = maxedge2
       } else {
      	 hist.d = which_num(edge[[base.length + d]][[3]] + 384, timestamps)
       }
       edgetime.d = rep(NA, K)
       for (w in 1:length(currentZ[[base.length + d]])) {
       	 zw.old = currentZ[[base.length + d]][w]
       	 if (length(textlist.d) > 0) { 
       	 table.W[[zw.old]][textlist.d[w]] = table.W[[zw.old]][textlist.d[w]] - 1
         topicpart.d = TopicInEqZ(K, currentZ[[base.length + d]][-w], alpha, mvec)
         wordpart.d = WordInEqZ(K, textlist.d, table.W, betas, nvec)
       } else {
         topicpart.d = 0
         wordpart.d = matrix(0, nrow = length(currentZ[[base.length + d]]), ncol = K)
       }
 	 for (IP in unique(currentC)) {
       	  	currentCK = which(currentC == IP)
       		currentZ[[base.length + d]][w] = min(currentCK)
       	 	p.d[base.length + d, ] = pdmat(list(currentZ[[base.length + d]]), currentC, nIP)          
            history.t = History(edge, p.d, node, edge[[hist.d-1]][[3]] + exp(-745))
    	    		X = Netstats_cpp(history.t, node, netstat)
    	    		XB = MultiplyXBList(X, b)   
    	    		lambda[[hist.d]] = lambda_cpp(p.d[hist.d,], XB)
	   	 	LambdaiJi = lambdaiJi(p.d[hist.d, ], XB, iJi[[hist.d]])
        		observediJi = LambdaiJi[edge[[hist.d]][[1]]]
        		edgetime.d[currentCK] = EdgeTime(iJi[[hist.d]], lambda[[hist.d]], delta, LambdaiJi, timeinc[hist.d], observediJi)   
            }
          	const.Z = edgetime.d + topicpart.d + wordpart.d[w, ] 
          	zw.new = multinom_vec(1, expconst(const.Z))
         if (zw.new != zw.old) {
           currentZ[[base.length + d]][w] = zw.new
           table.W[[zw.new]][textlist.d[w]] = table.W[[zw.new]][textlist.d[w]] + 1       
         } else {
         	currentZ[[base.length + d]][w] = zw.old
         	table.W[[zw.old]][textlist.d[w]] = table.W[[zw.old]][textlist.d[w]] + 1
         }
       }
      } 
   #regenerate data   
   p.d = pdmat(currentZ, currentC, nIP)  		  
  word_type_topic_counts = matrix(0, W, K)
  maxedge2 = length(currentZ)
  for (d in 1:nDocs) {
    N.d = length(currentZ[[base.length + d]])
    text[[base.length + d]] = rep(NA, N.d)
      phi.k = rep(NA, K)
      topic.d = currentZ[[base.length + d]]
        for (n in 1:N.d){
          for (w in 1:W) {
            phi.k[w] = (word_type_topic_counts[w, topic.d[n]] + betas * nvec[w]) / (sum(word_type_topic_counts[, topic.d[n]]) + betas)
          } 
          text[[base.length + d]][n] = multinom_vec(1, phi.k)
          word_type_topic_counts[text[[base.length + d]][n], topic.d[n]] = word_type_topic_counts[text[[base.length + d]][n], topic.d[n]] + 1
        }
        names(text[[base.length + d]]) = topic.d

    history.t = History(edge, p.d, node, t.d + exp(-745))
    X = Netstats_cpp(history.t, node, netstat)
    XB = MultiplyXBList(X, b)     
    lambda[[d]] = lambda_cpp(p.d[base.length + d,], XB)
    LambdaiJi = lambdaiJi(p.d[base.length + d,], XB, iJi[[base.length + d]])
    i.d = multinom_vec(1, LambdaiJi)
    j.d = which(iJi[[base.length + d]][i.d,] == 1)
    t.d = t.d + rexp(1, sum(LambdaiJi))
    edge[[base.length + d]] = list(sender = i.d, receiver = j.d, timestamp = t.d)		
  }	 		
  	}
  return(list(edge = edge, text = text, iJi = iJi))							
} 

#' @title GenerateDocs.predict
#' @description Generate a collection of documents according to the generative process of IPTM using Gibbs measure
#'
#' @param nDocs number of documents to be generated
#' @param node nodelist containing the ID of nodes (ID starting from 1)
#' @param vocabulary all vocabularies used over the corpus
#' @param nIP total number of interaction patterns specified by the user
#' @param K total number of topics specified by the user
#' @param owords observed words in a document
#' @param alpha Dirichlet concentration prior for document-topic distribution
#' @param mvec Dirichlet base prior for document-topic distribution
#' @param betas Dirichlet concentration prior for topic-word distribution
#' @param nvec Dirichlet base prior for topic-word distribution
#' @param b List of coefficients estimating the history effect for each IP 
#' @param delta tuning parameter controlling the number of recipients
#' @param currentC topic-IP assignment
#' @param netstat which type of network statistics to use ("intercept", "dyadic", "triadic", "degree")
#' @param initial_iJi initial value of iJi from emprical summary of inferred latent receivers
#' @param base.edge edges before 384 hours that is used to calculate initial history of interactions
#' @param base.text texts corresponding to base.edge
#' @param topic_token_assignments matrix of topic-token assignments
#' @param edge2 index of documents after initial history
#' @param niter number of iterations for Gibbs update of topic assignments and latent receivers
#'
#' @return generated edge and text, parameter b used to generate those, and base (if base == TRUE)
#'
#' @export
GenerateDocs.predict = function(nDocs = 1, node, vocabulary, nIP, K, owords, alpha, mvec, betas, nvec,
                        b, delta, currentC, netstat, initial_iJi, base.edge, base.text, topic_token_assignments = NULL, edge2, niter) {

  netstat = as.numeric(c("intercept", "degree", "dyadic", "triadic") %in% netstat)
  L = 3
  P = netstat[1] + L * (2 * netstat[2] + 2 * netstat[3] + 4 * netstat[4])
  t.d = base.edge[[length(base.edge)]][[3]]
  W = length(vocabulary)
  
  edge = base.edge
  text = base.text
  base.length = length(edge)
  history.t = lapply(1:nIP, function(IP) {
        	 	 lapply(1:3, function(l){
         		matrix(0, length(node), length(node))
       			})
     		 }) 
  textlist.raw = unlist(text[edge2])
  table.W = lapply(1:K, function(k) {
   tabulateC(textlist.raw[which(unlist(topic_token_assignments[edge2]) == k)], W)
  })	
  iJi = list()
  for (d in 1:nDocs) {
    N.d = length(owords)
    text[[base.length + d]] = owords
    iJi[[base.length + d]] = initial_iJi
    topic_token_assignments[[base.length + d]] = rep(NA, N.d)
        for (n in 1:N.d){ 
          #initial topic assignment from empirical topic token assignments
            theta.d = rdirichlet_cpp(1, alpha * mvec + vapply(table.W, function(x){x[owords[n]]}, c(1)))
            topic_token_assignments[[base.length + d]][n] = multinom_vec(1, theta.d)
        }
    textlist.raw = unlist(text[c(edge2, base.length + d)])
    table.W = lapply(1:K, function(k) {
      tabulateC(textlist.raw[which(unlist(topic_token_assignments[c(edge2, base.length + d)]) == k)], W)
    })	
    for (it in 1:niter[1]) {
      for (w in 1:N.d) {
      	zw.old = topic_token_assignments[[base.length + d]][w]
      	currentW = text[[base.length + d]][w]
       if (N.d > 0) {
       	table.W[[zw.old]][[currentW]] = table.W[[zw.old]][[currentW]] - 1
        topicpart.d = TopicInEqZ(K, topic_token_assignments[[base.length + d]][-w], alpha, mvec)
        wordpart.d = WordInEqZ(K, text[[base.length + d]], table.W, betas, nvec)
     	 } else {
        topicpart.d = 0
        wordpart.d = matrix(0, nrow = N.d, ncol = K)
      }
        const.Z = topicpart.d + wordpart.d[w, ] 
        zw.new = multinom_vec(1, expconst(const.Z))
        if (zw.new != zw.old) {
          topic_token_assignments[[base.length + d]][w] = zw.new
          table.W = lapply(1:K, function(k) {
            tabulateC(textlist.raw[which(unlist(topic_token_assignments[c(edge2, base.length + d)]) == k)], W)
          }) 
        } else {
        		table.W[[zw.old]][[currentW]] = table.W[[zw.old]][[currentW]] + 1
        }  
      }
    }  
    names(text[[base.length + d]]) = topic_token_assignments[[base.length + d]]
    p.d = pdmat(topic_token_assignments[1:(base.length + nDocs)] , currentC, nIP) 
  	if (t.d >= 384) {
      history.t = History(edge, matrix(p.d[1:base.length, ], ncol = nIP), node, t.d + exp(-745))
    }
    X = Netstats_cpp(history.t, node, netstat)
    XB = MultiplyXBList(X, b)
    lambda = lambda_cpp(p.d[base.length + d,], XB)
    for (iter in 1:niter[2]) {
      for (i in node) {
		for (j in sample(node[-i], length(node) - 1)) {
			 probij = DataAug_cpp_Gibbs_noObs(iJi[[base.length + d]][i, ], lambda[i,], lapply(XB, function(IP) {IP[i,]}), p.d[base.length + d, ], delta, j)
			 iJi[[base.length + d]][i, j] = multinom_vec(1, probij) - 1		
			}
      }
    }  
    LambdaiJi = lambdaiJi(p.d[base.length + d,], XB, iJi[[base.length + d]])
    Time.inc = t.d + vapply(LambdaiJi, function(lambda) {
      rexp(1, lambda)
    }, c(1))
    i.d = which(Time.inc == min(Time.inc[!is.na(Time.inc)]))
    j.d = which(iJi[[base.length + d]][i.d,] == 1)
    t.d = Time.inc[i.d]
    edge[[base.length + d]] = list(sender = i.d, receiver = j.d, timestamp = t.d)		
  }
  return(list(edge = edge[[length(edge)]], text = text[[length(text)]], iJi = iJi[[length(iJi)]], time = Time.inc))							
} 

#' @title IPTM_predict.data
#' @description posterior predictive experiments
#'
#' @param D Dth document the user wants to predict
#' @param O number of outer iterations of inference from which to generate predictions
#' @param R the number of iterations to sample predicted data within each outer iteration
#' @param edge list of document information with 3 elements (element 1 sender, element 2 receiver, element 3 time in unix.time format)
#' @param node nodelist containing the ID of nodes (ID starting from 1)
#' @param textlist list of text (length=number of documents in total) containing the words in each document
#' @param vocabulary all vocabularies used over the corpus
#' @param nIP total number of interaction patterns specified by the user
#' @param K total number of topics specified by the user
#' @param sigma_Q proposal distribution variance parameter for beta and delta
#' @param alpha Dirichlet concentration prior for document-topic distribution
#' @param mvec Dirichlet base prior for document-topic distribution
#' @param betas Dirichlet concentration prior for topic-word distribution
#' @param nvec Dirichlet base prior for topic-word distribution
#' @param prior.b.mean mean vector of b in multivariate normal distribution
#' @param prior.b.var covairance matrix of b in multivariate normal distribution
#' @param prior.delta parameter of delta in Normal prior
#' @param out size of outer iterations 
#' @param n_B size of third inner iteration for updates of B
#' @param n_d size of third inner iteration for updates of delta
#' @param burn iterations to be discarded at the beginning of beta and delta chain 
#' @param thinning the thinningning interval of beta and delta chain
#' @param netstat which type of network statistics to use ("intercept", dyadic", "triadic", "degree")
#' @param optimize to optimize alpha (Dirichlet concentration prior for document-topic distribution) or not (TRUE/FALSE)
#' @param niter number of iterations for Gibbs update of topic assignments and latent receivers
#'
#' @return prediction output 
#'
#' @export
IPTM_predict.data = function(D, O, R, edge, node, textlist, vocabulary, nIP, K, sigma_Q, alpha, mvec, betas, nvec, 
							 prior.b.mean, prior.b.var, prior.delta, out, n_B, n_d, burn, thinning, netstat, optimize = FALSE, niter = c(1,1)) {
   New_sample = list()
   for (o in 1:O) {
   	New_sample[[o]] = list()
    Inference_samp = IPTM_inference.data(edge[1:(D-1)], node, textlist[1:(D-1)], vocabulary, nIP, K,
    									 sigma_Q, alpha, mvec, betas, nvec, prior.b.mean, prior.b.var, prior.delta,
                      out, n_B, n_d, burn, thinning, netstat, optimize = optimize)
	b = lapply(1:nIP, function(IP) {
        rowMeans(Inference_samp$B[[IP]])
    })
    delta = mean(Inference_samp$D)
    currentC = Inference_samp$C
    topic_token_assignments = Inference_samp$Z
    initial_iJi = matrix(rbinom(length(node)^2, 1, c(Reduce('+', Inference_samp$iJi) / (length(Inference_samp$iJi)))), nrow =length(node), ncol = length(node))
	alpha = ifelse(optimize == TRUE, Inference_samp$alpha[out,], alpha)
	if (optimize == TRUE) {
		mvec = Inference_samp$mvec[out,]
		} else {
			mvec = tabulate(unlist(Inference_samp$Z), K) / length(unlist(Inference_samp$Z))
			}
	for (r in 1:R) {
	New_sample[[o]][[r]] = GenerateDocs.predict(1, node, vocabulary, nIP, K, owords = textlist[[D]],
										alpha = alpha, mvec = mvec, 
										betas, nvec, b, 
    										 delta, currentC, netstat, initial_iJi = initial_iJi, base.edge = edge[1:(D-1)], base.text = textlist[1:(D-1)],
    										 topic_token_assignments = topic_token_assignments, edge2 = Inference_samp$edge2, niter = niter)
	}
   }  
 return(New_sample)                                	
}

#' @title IPTM_predict.data2
#' @description posterior predictive experiments
#'
#' @param D Dth document the user wants to predict
#' @param O number of outer iterations of inference from which to generate predictions
#' @param R the number of iterations to sample predicted data within each outer iteration
#' @param edge list of document information with 3 elements (element 1 sender, element 2 receiver, element 3 time in unix.time format)
#' @param node nodelist containing the ID of nodes (ID starting from 1)
#' @param textlist list of text (length=number of documents in total) containing the words in each document
#' @param vocabulary all vocabularies used over the corpus
#' @param nIP total number of interaction patterns specified by the user
#' @param K total number of topics specified by the user
#' @param sigma_Q proposal distribution variance parameter for beta and delta
#' @param alpha Dirichlet concentration prior for document-topic distribution
#' @param mvec Dirichlet base prior for document-topic distribution
#' @param betas Dirichlet concentration prior for topic-word distribution
#' @param nvec Dirichlet base prior for topic-word distribution
#' @param prior.b.mean mean vector of b in multivariate normal distribution
#' @param prior.b.var covairance matrix of b in multivariate normal distribution
#' @param prior.delta parameter of delta in Normal prior
#' @param out size of outer iterations 
#' @param n_B size of third inner iteration for updates of B
#' @param n_d size of third inner iteration for updates of delta
#' @param burn iterations to be discarded at the beginning of beta and delta chain 
#' @param thinning the thinningning interval of beta and delta chain
#' @param netstat which type of network statistics to use ("intercept", dyadic", "triadic", "degree")
#' @param optimize to optimize alpha (Dirichlet concentration prior for document-topic distribution) or not (TRUE/FALSE)
#' @param niter number of iterations for Gibbs update of topic assignments and latent receivers
#' @param initial list of initial values user wants to assign
#'
#' @return prediction output 
#'
#' @export
IPTM_predict.data2 = function(D, O, R, edge, node, textlist, vocabulary, nIP, K, sigma_Q, alpha, mvec, betas, nvec, 
							 prior.b.mean, prior.b.var, prior.delta, out, n_B, n_d, burn, thinning, netstat, optimize = FALSE, niter = c(1,1), 
							 initial) {
   New_sample = list()
   for (o in 1:O) {
   	print(o)
   	New_sample[[o]] = list()
    Inference_samp = IPTM_inference.data(edge[1:(D-1)], node, textlist[1:(D-1)], vocabulary, nIP, K,
    									 sigma_Q, alpha, mvec, betas, nvec, prior.b.mean, prior.b.var, prior.delta,
                      					 out, n_B, n_d, burn, thinning, netstat, optimize = optimize, initial = initial)
	b = lapply(1:nIP, function(IP) {
        rowMeans(Inference_samp$B[[IP]])
    })
    delta = mean(Inference_samp$D)
    currentC = Inference_samp$C
    topic_token_assignments = Inference_samp$Z
    initial_iJi = matrix(rbinom(length(node)^2, 1, c(Reduce('+', Inference_samp$iJi) / (length(Inference_samp$iJi)))), nrow =length(node), ncol = length(node))

	for (r in 1:R) {
	New_sample[[o]][[r]] = GenerateDocs.predict(1, node, vocabulary, nIP, K, owords = textlist[[D]],
												alpha = Inference_samp$alpha, mvec = tabulate(unlist(Inference_samp$Z)) / length(unlist(Inference_samp$Z)), 
												betas, nvec, b, delta, currentC, netstat, initial_iJi = initial_iJi, base.edge = edge[1:(D-1)], 
												base.text = textlist[1:(D-1)], topic_token_assignments = topic_token_assignments, edge2 = Inference_samp$edge2, niter = niter)
	}
   }  
 return(New_sample)                                	
}
