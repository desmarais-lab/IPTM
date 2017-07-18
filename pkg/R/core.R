#' @useDynLib IPTM
#' @import stats
#' @import grDevices
#' @import graphics
#' @import RcppArmadillo
#' @importFrom Rcpp sourceCpp
#' @importFrom mvtnorm rmvnorm dmvnorm
#' @importFrom entropy entropy.empirical
#' @importFrom lda top.topic.words
#' @importFrom reshape melt
#' @importFrom coda mcmc geweke.diag
#' @importFrom combinat permn

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
		gibbs.mat.i = do.call('rbind',permn(c(rep(1, i),rep(0, n - i))))
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
adaptive_MH = function(sigma_Q, accept_rates, target = 0.25, update_size, tol = 0.05) {
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

#' @title Netstats
#' @description Calculate the network statistics (dydadic, triadic, or degree) given the history of interactions
#'
#' @param historyIP list containing the weighted history of interactions for each IP
#' @param node nodelist containing the ID of nodes
#' @param sender specific timepoint that we are calculating the time difference from
#' @param netstat which type of network statistics to use ("intercept", "dyadic", "triadic", "degree")
#'
#' @return list of network statistics for specific sender and all possible receivers for each IP
#'
#' @export
Netstats = function(historyIP, node, sender, netstat) {
  # Calculate the chosen network statistics given the history of interactions
  #
  # Args 
  #  historyIP list containing the weighted history of interactions for each IP
  #  node nodelist containing the ID of nodes
  #  sender specific timepoint that we are calculating the time difference from
  #  netstat which type of network statistics to use ("dyadic", "triadic", "degree")
  #  nIP total number of interaction patterns specified by the user
  #
  # Returns 
  #  list of network statistics for specific sender and all possible receivers for each IP
  nIP = length(historyIP)
  if ("intercept" %in% netstat) {
  intercept = rep(1, length(node))
  # matcolnames = c("intercept")
  } else {
  	intercept = NULL
  }
   if ("degree" %in% netstat) {
    degree = Degree(historyIP, node, sender)
    # matcolnames = c(matcolnames, paste0("outdegree", 1:3), paste0("indegree", 1:3))
  } else {
  	degree = NULL
  }
	if ("dyadic" %in% netstat) {
    dyadic = Dyadic(historyIP, node, sender)
    # matcolnames = c(matcolnames, paste0("send", 1:3), paste0("receive", 1:3))
  } else {
  	dyadic = NULL
  }
  if ("triadic" %in% netstat) {
    triadic0 = Triadic(historyIP, node, sender)
    triadic = Triadic_reduced(triadic0)
    # matcolnames = c(matcolnames, paste0("2-send", 1:3), paste0("2-receive", 1:3), 
    #                paste0("sibling", 1:3),  paste0("cosibling", 1:3) )
  } else {
  	triadic = NULL
  }
 
  netstatmat = lapply(1L:nIP, function(IP) {
  	matrix(cbind(intercept, degree[[IP]], dyadic[[IP]], triadic[[IP]][,1:12]), nrow = length(node))
  	})
  # for (IP in 1L:nIP) {
  #   colnames(netstatmat[[IP]]) = matcolnames
  # }
  return(netstatmat)
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
		nK.word.table = lapply(1L:K, function(k){
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

#' @title logWZ
#' @description Calculate the log of unnormalized constant to check the convergence
#'
#' @param K total number of topics specified by the user
#' @param currentZ current state of the assignment of topics 
#' @param textlist list of text (length=number of documents in total) containing the words in each document
#' @param tableW summary table of topic-word assignments
#' @param alpha Dirichlet concentration prior for document-topic distribution
#' @param mvec Dirichlet base prior for document-topic distribution
#' @param betas Dirichlet concentration prior for topic-word distribution
#' @param nvec Dirichlet base prior for topic-word distribution
#'
#' @return The log value of unnormalized constant
#'
#' @export
logWZ = function(K, currentZ, textlist, tableW, alpha, mvec, betas, nvec) {
  # Calculate the log of unnormalized constant to check the convergence
  #
  # Args 
  #  K total number of topics specified by the user
  #  currentZ current state of the assignment of topics 
  #  textlist list of text (length=number of documents in total) containing the words in each document
  #  tableW summary table of topic-word assignments
  #  alpha Dirichlet concentration prior for document-topic distribution
  #  mvec Dirichlet base prior for document-topic distribution
  #  betas Dirichlet concentration prior for topic-word distribution
  #  nvec Dirichlet base prior for topic-word distribution
  #
  # Returns 
  #  The log value of unnormalized constant
  finalsum = 0
  for (d in seq(along = currentZ)) {
    currentZ_d = currentZ[[d]]
    k_table = tabulateC(currentZ_d, K)
    textlist_d = textlist[[d]]
    part4 = log(sum(k_table) - 1 + alpha)
    for (k in 1:length(currentZ_d)) {
      table_Wk = tableW[[currentZ_d[k]]]
      if (length(textlist_d) > 0 ){
        part1 = log(table_Wk[textlist_d[k]] - (table_Wk[textlist_d[k]] > 0) + betas * nvec[textlist_d[k]])
      } else {part1 = 0}
      part2 = log(sum(table_Wk) - sum(table_Wk > 0) + betas)
      part3 = log(k_table[currentZ_d[k]] - 1 + alpha * mvec[currentZ_d[k]])
      finalsum = finalsum + part1 - part2 + part3 - part4
    }
  }
  return(finalsum)
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
#' @param plot to plot the convergence diagnostics or not (TRUE/FALSE)
#' @param optimize to optimize alpha (Dirichlet concentration prior for document-topic distribution) or not (TRUE/FALSE)
#'
#' @return MCMC output containing IP assignment, topic assignment, and (beta, mu, delta) chain 
#'
#' @export
IPTM_inference.Gibbs = function(edge, node, textlist, vocabulary, nIP, K, sigma_Q, alpha, mvec, betas, nvec, prior.b.mean, prior.b.var, prior.delta, 
					out, n_B, n_d, burn, thinning, netstat, plot = FALSE, optimize = FALSE) {
   
  # trim the edge so that we only model edges after 384 hours
	timestamps = vapply(edge, function(d) {
  			  d[[3]]
 			  }, c(1))

    edge2 = which_int(384, timestamps) : length(edge)
    timeinc = c(timestamps[1], timestamps[-1] - timestamps[-length(timestamps)])
   
  # initialize alpha, mvec, delta, nvec, delta, lvec, and gammas
 	W = length(vocabulary)
  	phi = lapply(1L:K, function(k) {
		rdirichlet_cpp(1, betas * nvec)
	})
	delta = rnorm(1, prior.delta[1], sqrt(prior.delta[2]))
	beta.old = lapply(1:nIP, function(IP) {
      			  c(rmvnorm(1,  prior.b.mean, prior.b.var))
      			 })
	# initialize C, theta and Z
	 currentC = sample(1L:nIP, K, replace = TRUE)  	 
	 theta = rdirichlet_cpp(length(edge), alpha * mvec)
     currentZ = lapply(seq(along = edge), function(d) {
     			 multinom_vec(max(1, length(textlist[[d]])), theta[d, ])
 				 })
 	 p.d = t(vapply(seq(along = edge), function(d) {
    	vapply(1L:nIP, function(IP) {
      	sum(currentZ[[d]] %in% which(currentC == IP))
  	 		 }, c(1)) / length(currentZ[[d]])
 		 }, rep(1, nIP)))

    # initialize beta
    L = 3
    P = 1 * ("intercept" %in% netstat) + L * (2 * ("dyadic" %in% netstat) + 4 * ("triadic" %in% netstat) + 2 *("degree" %in% netstat))
    bmat = list()
	for (IP in 1L:nIP) {
		bmat[[IP]] = matrix(beta.old[[IP]], nrow = P, ncol = (n_B - burn) / thinning)
  	}
  	deltamat = rep(delta, n_d)
    proposal.var = lapply(1:nIP, function(IP){diag(P)})

    # to check the convergence  
    if (plot) {
     	logWZ.mat = c()							  
     	alpha.mat = c()
     	entropy.mat = c()
    }

  	#initialize the latent sender-receiver pairs
  	iJi = lapply(seq(along = edge), function(d) {
  		matrix(0, nrow = length(node), ncol = length(node))
  	})
  	lambda = list()
    LambdaiJi = list()
	observediJi = list()
	support = gibbs.measure.support(length(node) - 1)
	for (d in edge2) {
   	 	history.t = History(edge, p.d, node, as.numeric(edge[[d-1]][3]) + 10^(-10))
   	 	X = lapply(node, function(i) {
  	        Netstats(history.t, node, i, netstat)
            })
   	 	XB = MultiplyXBList(X, beta.old)     
		lambda[[d]] = lambda_cpp(p.d[d,], XB)
    	for (i in node) {
    		iJi[[d]][i, -i] = r.gibbs.measure(1, lambda[[d]][i, -i], delta, support)
    		}
    	}

    #start outer iteration
    for (o in 1L:out) {
      
      if (optimize) {
      #update the hyperparameter alpha and mvec
      vec = AlphamvecOpt(K, currentZ[edge2], alpha, mvec, 5)
      alpha = sum(vec)
      mvec = vec / alpha
      }
     # Data augmentation
      for (d in edge2) {
   	 	history.t = History(edge, p.d, node, as.numeric(edge[[d-1]][3]) + 10^(-10))
   	 	X = lapply(node, function(i) {
  	        Netstats(history.t, node, i, netstat)
            })
   	 	XB = MultiplyXBList(X, beta.old)     
		lambda[[d]] = lambda_cpp(p.d[d,], XB)
		#calculate the resampling probability	
		for (i in node[-as.numeric(edge[[d]][1])]) {
			for (j in sample(node[-i], length(node) - 1)) {
				probij = DataAug_cpp_Gibbs(iJi[[d]][i, ], lambda[[d]][i,], lapply(XB, function(IP) {IP[i,]}), p.d[d, ], delta, timeinc[d], j)
				iJi[[d]][i, j] = multinom_vec(1, probij) - 1		
				}
		}
		iJi[[d]][as.numeric(edge[[d]][1]),] = tabulateC(as.numeric(unlist(edge[[d]][2])), length(node))
		LambdaiJi[[d]] = lambdaiJi(p.d[d,], XB, iJi[[d]])
		observediJi[[d]] = LambdaiJi[[d]][as.numeric(edge[[d]][1])]
		}	 
	
	  # Z update
    	  textlist.raw = unlist(textlist[edge2])
    	  table.W = lapply(1L:K, function(k) {
      			tabulateC(textlist.raw[which(unlist(currentZ[edge2]) == k)], W)
      			})    
      	for (d in edge2) {
      		textlist.d = textlist[[d]]
        		if (length(textlist.d) > 0) {
        			topicpart.d = TopicInEqZ(K, currentZ[[d]], alpha, mvec, d)
       			wordpart.d = WordInEqZ(K, textlist.d, table.W, betas, nvec)
        		} else {
        			topicpart.d = 0
        			wordpart.d = matrix(0, nrow = length(currentZ[[d]]), ncol = K)
        		}
        		edgepart.d = EdgeInEqZ_Gibbs(iJi[[d]], lambda[[d]], delta)
        		timepart.d = TimeInEqZ(LambdaiJi[[d]], timeinc[d])
        		observed.d = ObservedInEqZ(observediJi[[d]]) 
        		fixedpart = topicpart.d + edgepart.d + timepart.d + observed.d 
        		for (w in 1L:length(currentZ[[d]])) {
          		const.Z = fixedpart + wordpart.d[w, ]
          		const.Z = const.Z - max(const.Z)
          		zw.old = currentZ[[d]][w]
          		zw.new = multinom_vec(1, exp(const.Z))
          		if (zw.new != zw.old) {
            			currentZ[[d]][w] = zw.new
            			topicpart.d = TopicInEqZ(K, currentZ[[d]], alpha, mvec, d)
            			if (length(textlist.d) > 0) {	
            				wordpart.d = WordInEqZ(K, textlist.d, table.W, betas, nvec)
            			}
            			table.W = lapply(1L:K, function(k) {
      				  tabulateC(textlist.raw[which(unlist(currentZ[edge2]) == k)], W)
      				})
      				p.d[d, ] = vapply(1L:nIP, function(IP) {
	 					sum(currentZ[[d]] %in% which(currentC == IP))
	 				}, c(1)) / length(currentZ[[d]])
      				LambdaiJi[[d]] = lambdaiJi(p.d[d,], XB, iJi[[d]])
         			observediJi[[d]] = LambdaiJi[[d]][as.numeric(edge[[d]][1])]
          		}
        		}
       	}		
      # C update given Z and B - withinning each document d
       for (k in sort(unique(unlist(currentZ[edge2])))) { 
         const.C = rep(NA, nIP)
         for (IP in 1:nIP) {
          	 currentC[k] = IP
          	 p.d = t(vapply(seq(along = edge), function(d) {
    				 vapply(1L:nIP, function(IP) {
      				 sum(currentZ[[d]] %in% which(currentC == IP))
  	 				  }, c(1)) / length(currentZ[[d]])
 					  }, rep(1, nIP)))
          	 for (d in edge2) {
           		 history.t = History(edge, p.d, node, as.numeric(edge[[d-1]][3]) + 10^(-10))
    	   	   		 X = lapply(node, function(i) {
               		 Netstats(history.t, node, i, netstat)
               		 })
    	       		 XB = MultiplyXBList(X, beta.old)    
           		 lambda[[d]] = lambda_cpp(p.d[d,], XB)
		       	 LambdaiJi[[d]] = lambdaiJi(p.d[d,], XB, iJi[[d]])
           		 observediJi[[d]] = LambdaiJi[[d]][as.numeric(edge[[d]][1])]
          	 }
           const.C[IP] = sum(vapply(edge2, function(d) {
          				 EdgeInEqZ_Gibbs(iJi[[d]], lambda[[d]], delta) + 
          				 TimeInEqZ(LambdaiJi[[d]], timeinc[d]) + 
    						 ObservedInEqZ(observediJi[[d]]) 
          				 }, c(1))) / length(edge2)
      	 }
         const.C = const.C - max(const.C)
         currentC[k] = multinom_vec(1, exp(const.C))
      }
      
    p.d = t(vapply(seq(along = edge), function(d) {
    				vapply(1L:nIP, function(IP) {
      				sum(currentZ[[d]] %in% which(currentC == IP))
  	 				 }, c(1)) / length(currentZ[[d]])
 					 }, rep(1, nIP)))  
    for (d in edge2) {
        	history.t = History(edge, p.d, node, as.numeric(edge[[d-1]][3]) + 10^(-10))
    	    X = lapply(node, function(i) {
            Netstats(history.t, node, i, netstat)
       		})
    	    XB = MultiplyXBList(X, beta.old)   
    	    lambda[[d]] = lambda_cpp(p.d[d,], XB)
	    LambdaiJi[[d]] = lambdaiJi(p.d[d,], XB, iJi[[d]])
        observediJi[[d]] = LambdaiJi[[d]][as.numeric(edge[[d]][1])]
	}
	 
	# beta update
	prior.old1 = sum(vapply(1L:nIP, function(IP) {
		  		 dmvnorm(beta.old[[IP]], prior.b.mean, prior.b.var, log = TRUE)
		  		 }, c(1))) 
	post.old1 = sum(vapply(edge2, function(d) {
	     	    EdgeInEqZ_Gibbs(iJi[[d]], lambda[[d]], delta) + TimeInEqZ(LambdaiJi[[d]], timeinc[d]) +
    			 	ObservedInEqZ(observediJi[[d]])
    			 	}, c(1))) / length(edge2)
  
    if (o != 1) {
    		accept.rates = c(length(unique(bmat[[1]][1,])) / ncol(bmat[[1]]), length(unique(deltamat)) / n_d)
    	sigma_Q = adaptive_MH(sigma_Q, accept.rates, update_size = 0.1 * sigma_Q)
    	if (accept.rates[1] > 1 / ncol(bmat[[1]])) {
    		for (IP in 1:nIP) {
     		proposal.var[[IP]] = var(t(bmat[[IP]]))
      	}
      }
   }
    for (i3 in 1L:n_B) {
    		beta.new = lapply(1L:nIP, function(IP) {
          		   rmvnorm(1, beta.old[[IP]], sigma_Q[1] * proposal.var[[IP]])
         		   }) 
        for (d in edge2) {
           XB = MultiplyXBList(X, beta.new)
           lambda[[d]] = lambda_cpp(p.d[d,], XB)    
           LambdaiJi[[d]] = lambdaiJi(p.d[d,], XB, iJi[[d]])
	       observediJi[[d]] = LambdaiJi[[d]][as.numeric(edge[[d]][1])]
        }
        prior.new1 = sum(vapply(1L:nIP, function(IP) {
        				 dmvnorm(beta.new[[IP]], prior.b.mean, prior.b.var, log = TRUE)
        				 }, c(1))) 
        post.new1 = sum(vapply(edge2, function(d) {
    			   		EdgeInEqZ_Gibbs(iJi[[d]], lambda[[d]], delta) + TimeInEqZ(LambdaiJi[[d]], timeinc[d]) + 
    			    	ObservedInEqZ(observediJi[[d]])
    			    	}, c(1))) / length(edge2)
    		loglike.diff = prior.new1 + post.new1 - prior.old1 - post.old1
        if (log(runif(1, 0, 1)) < loglike.diff) {
        		for (IP in 1L:nIP) {
         		beta.old[[IP]] = beta.new[[IP]]
         	}
         	prior.old1 = prior.new1
         	post.old1 = post.new1
        }
         if (i3 > burn & i3 %% (thinning) == 0) {
        		for (IP in 1L:nIP) {
           		bmat[[IP]][ , (i3 - burn) / thinning] = beta.old[[IP]]
           	}
    		}
    }
    
    #delta update
    for (d in edge2) {
    		history.t = History(edge, p.d, node, as.numeric(edge[[d-1]][3]) + 10^(-10))
    	X = lapply(node, function(i) {
            Netstats(history.t, node, i, netstat)
       		})
    	 	XB = MultiplyXBList(X, beta.old)
     	lambda[[d]] = lambda_cpp(p.d[d,], XB)  
    }
    prior.old2 = dnorm(delta, prior.delta[1], sqrt(prior.delta[2]), log = TRUE)
    post.old2 = sum(vapply(edge2, function(d) {
    				EdgeInEqZ_Gibbs(iJi[[d]], lambda[[d]], delta)
    				}, c(1))) / length(edge2)
 	for (i4 in 1L:n_d) {
        delta.new = rnorm(1, delta, sqrt(sigma_Q[2]))
        prior.new2 = dnorm(delta.new, prior.delta[1], sqrt(prior.delta[2]), log = TRUE)
        post.new2 = sum(vapply(edge2, function(d) {
        				EdgeInEqZ_Gibbs(iJi[[d]], lambda[[d]], delta.new)
        				}, c(1))) / length(edge2)
        loglike.diff2 = prior.new2 + post.new2 - prior.old2 - post.old2 
        if (log(runif(1, 0, 1)) < loglike.diff2) {
        		delta = delta.new
         	prior.old2 = prior.new2
         	post.old2 = post.new2
        } 
            deltamat[i4] = delta
    	}
 }
         
  if (plot) {
    par(mfrow = c(1, 2))
	matplot(bmat[[1]][7,], lty = 1, col = 1L:P, type = "l", 
	          main = "Traceplot of beta", xlab = "(Inner) Iterations", ylab = "")
	abline(h = mean(bmat[[1]][7,]), lty = 1, col = 1L)
	plot(deltamat, type = "l", 
	xlab = "(Outer) Iterations", ylab = "")
	abline(h = mean(deltamat), lty = 1)
	title("Traceplot of delta")
  }
     
  chain.final = list(C = currentC, Z = lapply(edge2, function(d) {currentZ[[d]]}), B = bmat, D = deltamat,
                     iJi = iJi, sigma_Q =sigma_Q, alpha = alpha, mvec = mvec, 
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
#' @param plot to plot the convergence diagnostics or not (TRUE/FALSE)
#' @param optimize to optimize alpha (Dirichlet concentration prior for document-topic distribution) or not (TRUE/FALSE)
#'
#' @return MCMC output containing IP assignment, topic assignment, and (beta, mu, delta) chain 
#'
#' @export
IPTM_inference.data = function(edge, node, textlist, vocabulary, nIP, K, sigma_Q, alpha, mvec, betas, nvec, prior.b.mean, prior.b.var, prior.delta, 
                                out, n_B, n_d, burn, thinning, netstat, plot = FALSE, optimize = FALSE) {
  
  # trim the edge so that we only model edges after 384 hours
  timestamps = vapply(edge, function(d) {
    d[[3]]
  }, c(1))
  
  edge2 = which_int(384, timestamps) : length(edge)
  timeinc = c(timestamps[1], timestamps[-1] - timestamps[-length(timestamps)])
  
  # initialize alpha, mvec, delta, nvec, delta, lvec, and gammas
  W = length(vocabulary)
  phi = lapply(1L:K, function(k) {
    rdirichlet_cpp(1, betas * nvec)
  })
  delta = rnorm(1, prior.delta[1], sqrt(prior.delta[2]))
  beta.old = lapply(1:nIP, function(IP) {
    c(rmvnorm(1,  prior.b.mean, prior.b.var))
  })
  # initialize C, theta and Z
  currentC = rep(1:nIP, K / nIP)
  theta = rdirichlet_cpp(length(edge), alpha * mvec)
  currentZ = lapply(seq(along = edge), function(d) {
    multinom_vec(max(1, length(textlist[[d]])), theta[d, ])
  })
  p.d = t(vapply(seq(along = edge), function(d) {
    vapply(1L:nIP, function(IP) {
      sum(currentZ[[d]] %in% which(currentC == IP))
    }, c(1)) / length(currentZ[[d]])
  }, rep(1, nIP)))
  
  # initialize beta
  L = 3
  P = 1 * ("intercept" %in% netstat) + L * (2 * ("dyadic" %in% netstat) + 4 * ("triadic" %in% netstat) + 2 *("degree" %in% netstat))
  bmat = list()
  for (IP in 1L:nIP) {
    bmat[[IP]] = matrix(beta.old[[IP]], nrow = P, ncol = (n_B - burn[1]) / thinning[1])
  }
  deltamat = rep(delta,  (n_d - burn[2]) / thinning[2])
  proposal.var = lapply(1:nIP, function(IP){diag(P)})
  
  # to check the convergence  
  if (plot) {
    logWZ.mat = c()							  
    alpha.mat = c()
    entropy.mat = c()
  }
  
  #initialize the latent sender-receiver pairs
  iJi = lapply(seq(along = edge), function(d) {
    matrix(0, nrow = length(node), ncol = length(node))
  })
  lambda = list()
  LambdaiJi = list()
  observediJi = list()
  for (d in edge2) {
    history.t = History(edge, p.d, node, as.numeric(edge[[d-1]][3]) + 10^(-10))
    X = lapply(node, function(i) {
      Netstats(history.t, node, i, netstat)
    })
    XB = MultiplyXBList(X, beta.old)     
    lambda[[d]] = lambda_cpp(p.d[d,], XB)
    iJi[[d]] = matrix(rbinom(length(node)^2, 1, 0.5), nrow =length(node), ncol = length(node))
    diag(iJi[[d]]) = 0
  }

  #start outer iteration
  for (o in 1L:out) {
    output = list(C = currentC, Z = currentZ, B = bmat, D = deltamat)
    save(output, file = "output.RData")
    print(o)
    if (optimize) {
      #update the hyperparameter alpha and mvec
      vec = AlphamvecOpt(K, currentZ[edge2], alpha, mvec, 5)
      alpha = sum(vec)
      mvec = vec / alpha
    }
    # Data augmentation
    for (d in edge2) {
      history.t = History(edge, p.d, node, as.numeric(edge[[d-1]][3]) + 10^(-10))
      X = lapply(node, function(i) {
        Netstats(history.t, node, i, netstat)
      })
      XB = MultiplyXBList(X, beta.old)     
      lambda[[d]] = lambda_cpp(p.d[d,], XB)
      #calculate the resampling probability	
      for (i in node[-as.numeric(edge[[d]][1])]) {
        for (j in sample(node[-i], length(node) - 1)) {
          probij = DataAug_cpp_Gibbs(iJi[[d]][i, ], lambda[[d]][i,], lapply(XB, function(IP) {IP[i,]}), p.d[d, ], delta, timeinc[d], j)
          iJi[[d]][i, j] = multinom_vec(1, probij) - 1
        }
      }
      iJi[[d]][as.numeric(edge[[d]][1]),] = tabulateC(as.numeric(unlist(edge[[d]][2])), length(node))
      LambdaiJi[[d]] = lambdaiJi(p.d[d,], XB, iJi[[d]])
      observediJi[[d]] = LambdaiJi[[d]][as.numeric(edge[[d]][1])]
    }	 
    
    # Z update

    	  textlist.raw = unlist(textlist[edge2])
    	  table.W = lapply(1L:K, function(k) {
      			tabulateC(textlist.raw[which(unlist(currentZ[edge2]) == k)], W)
      			})    
      	for (d in edge2) {
      		textlist.d = textlist[[d]]
        		if (length(textlist.d) > 0) {
        			topicpart.d = TopicInEqZ(K, currentZ[[d]], alpha, mvec, d)
       			wordpart.d = WordInEqZ(K, textlist.d, table.W, betas, nvec)
        		} else {
        			topicpart.d = 0
        			wordpart.d = matrix(0, nrow = length(currentZ[[d]]), ncol = K)
        		}
        		edgepart.d = EdgeInEqZ_Gibbs(iJi[[d]], lambda[[d]], delta)
        		timepart.d = TimeInEqZ(LambdaiJi[[d]], timeinc[d])
        		observed.d = ObservedInEqZ(observediJi[[d]]) 
        		fixedpart = topicpart.d + edgepart.d + timepart.d + observed.d 
        		for (w in 1L:length(currentZ[[d]])) {
          		const.Z = fixedpart + wordpart.d[w, ]
          		const.Z = const.Z - max(const.Z)
          		zw.old = currentZ[[d]][w]
          		zw.new = multinom_vec(1, exp(const.Z))
          		if (zw.new != zw.old) {
            			currentZ[[d]][w] = zw.new
            			topicpart.d = TopicInEqZ(K, currentZ[[d]], alpha, mvec, d)
            			if (length(textlist.d) > 0) {	
            				wordpart.d = WordInEqZ(K, textlist.d, table.W, betas, nvec)
            			}
            			table.W = lapply(1L:K, function(k) {
      				  tabulateC(textlist.raw[which(unlist(currentZ[edge2]) == k)], W)
      				})
      				p.d[d, ] = vapply(1L:nIP, function(IP) {
	 					sum(currentZ[[d]] %in% which(currentC == IP))
	 				}, c(1)) / length(currentZ[[d]])
      				LambdaiJi[[d]] = lambdaiJi(p.d[d,], XB, iJi[[d]])
         			observediJi[[d]] = LambdaiJi[[d]][as.numeric(edge[[d]][1])]
          		}
        		}
       	}		
    # C update given Z and B - withinning each document d
    # for (k in sort(unique(unlist(currentZ[edge2])))) { 
      # const.C = rep(NA, nIP)
      # for (IP in 1:nIP) {
        # currentC[k] = IP
        # p.d = t(vapply(seq(along = edge), function(d) {
          # vapply(1L:nIP, function(IP) {
            # sum(currentZ[[d]] %in% which(currentC == IP))
          # }, c(1)) / length(currentZ[[d]])
        # }, rep(1, nIP)))
        # for (d in edge2) {
          # history.t = History(edge, p.d, node, as.numeric(edge[[d-1]][3]) + 10^(-10))
          # X = lapply(node, function(i) {
            # Netstats(history.t, node, i, netstat)
          # })
          # XB = MultiplyXBList(X, beta.old)    
          # lambda[[d]] = lambda_cpp(p.d[d,], XB)
          # LambdaiJi[[d]] = lambdaiJi(p.d[d,], XB, iJi[[d]])
          # observediJi[[d]] = LambdaiJi[[d]][as.numeric(edge[[d]][1])]
        # }
        # const.C[IP] = sum(vapply(edge2, function(d) {
          # EdgeInEqZ_Gibbs(iJi[[d]], lambda[[d]], delta) + 
            # TimeInEqZ(LambdaiJi[[d]], timeinc[d]) + 
            # ObservedInEqZ(observediJi[[d]]) 
        # }, c(1))) / length(edge2)
      # }
      # const.C = const.C - max(const.C)
      # expconst.C = exp(const.C)
      # if (Inf %in% expconst.C) {
        # expconst.C[which(expconst.C == Inf)] = exp(700)
          # }
      # currentC[k] = multinom_vec(1, expconst.C)
    # }
    
    if (plot) {
      entropy.mat = c(entropy.mat, entropy.empirical(currentC))
      alpha.mat = rbind(alpha.mat, alpha)
      logWZ.mat = c(logWZ.mat, logWZ(K, currentZ[edge2], textlist[edge2], table.W, alpha, mvec, betas, nvec))
    }
    
    p.d = t(vapply(seq(along = edge), function(d) {
      vapply(1L:nIP, function(IP) {
        sum(currentZ[[d]] %in% which(currentC == IP))
      }, c(1)) / length(currentZ[[d]])
    }, rep(1, nIP)))  
    for (d in edge2) {
      history.t = History(edge, p.d, node, as.numeric(edge[[d-1]][3]) + 10^(-10))
      X = lapply(node, function(i) {
        Netstats(history.t, node, i, netstat)
      })
      XB = MultiplyXBList(X, beta.old)   
      lambda[[d]] = lambda_cpp(p.d[d,], XB)
      LambdaiJi[[d]] = lambdaiJi(p.d[d,], XB, iJi[[d]])
      observediJi[[d]] = LambdaiJi[[d]][as.numeric(edge[[d]][1])]
    }
    
    # beta update
    prior.old1 = sum(vapply(1L:nIP, function(IP) {
      dmvnorm(beta.old[[IP]], prior.b.mean, prior.b.var, log = TRUE)
    }, c(1))) 
    post.old1 = sum(vapply(edge2, function(d) {
      EdgeInEqZ_Gibbs(iJi[[d]], lambda[[d]], delta) + TimeInEqZ(LambdaiJi[[d]], timeinc[d]) +
        ObservedInEqZ(observediJi[[d]])
    }, c(1))) / length(edge2)
    
    if (o != 1) {
      accept.rates = c(length(unique(bmat[[1]][1,])) / ncol(bmat[[1]]), length(unique(deltamat)) / length(deltamat))
      sigma_Q = adaptive_MH(sigma_Q, accept.rates, update_size = 0.2 * sigma_Q)
      if (accept.rates[1] > 1 / ncol(bmat[[1]])) { 
        for (IP in 1:nIP) {
          proposal.var[[IP]] = var(t(bmat[[IP]]))
        }
      }
    }
    for (i3 in 1L:n_B) {
      beta.new = lapply(1L:nIP, function(IP) {
        rmvnorm(1, beta.old[[IP]], sigma_Q[1] * proposal.var[[IP]])
      }) 
      for (d in edge2) {
        XB = MultiplyXBList(X, beta.new)
        lambda[[d]] = lambda_cpp(p.d[d,], XB)    
        LambdaiJi[[d]] = lambdaiJi(p.d[d,], XB, iJi[[d]])
        observediJi[[d]] = LambdaiJi[[d]][as.numeric(edge[[d]][1])]
      }
      prior.new1 = sum(vapply(1L:nIP, function(IP) {
        dmvnorm(beta.new[[IP]], prior.b.mean, prior.b.var, log = TRUE)
      }, c(1))) 
      post.new1 = sum(vapply(edge2, function(d) {
        EdgeInEqZ_Gibbs(iJi[[d]], lambda[[d]], delta) + TimeInEqZ(LambdaiJi[[d]], timeinc[d]) + 
          ObservedInEqZ(observediJi[[d]])
      }, c(1))) / length(edge2)
      loglike.diff = prior.new1 + post.new1 - prior.old1 - post.old1
      if (log(runif(1, 0, 1)) < loglike.diff) {
        for (IP in 1L:nIP) {
          beta.old[[IP]] = beta.new[[IP]]
        }
        prior.old1 = prior.new1
        post.old1 = post.new1
      }
      if (i3 > burn[1] & i3 %% (thinning[1]) == 0) {
        for (IP in 1L:nIP) {
          bmat[[IP]][ , (i3 - burn[1]) / thinning[1]] = beta.old[[IP]]
        }
      }
    }
    
    #delta update
    for (d in edge2) {
      history.t = History(edge, p.d, node, as.numeric(edge[[d-1]][3]) + 10^(-10))
      X = lapply(node, function(i) {
        Netstats(history.t, node, i, netstat)
      })
      XB = MultiplyXBList(X, beta.old)
      lambda[[d]] = lambda_cpp(p.d[d,], XB)  
    }
    prior.old2 = dnorm(delta, prior.delta[1], sqrt(prior.delta[2]), log = TRUE)
    post.old2 = sum(vapply(edge2, function(d) {
      EdgeInEqZ_Gibbs(iJi[[d]], lambda[[d]], delta)
    }, c(1))) / length(edge2)
    for (i4 in 1L:n_d) {
      delta.new = rnorm(1, delta, sqrt(sigma_Q[2]))
      prior.new2 = dnorm(delta.new, prior.delta[1], sqrt(prior.delta[2]), log = TRUE)
      post.new2 = sum(vapply(edge2, function(d) {
        EdgeInEqZ_Gibbs(iJi[[d]], lambda[[d]], delta.new)
      }, c(1))) / length(edge2)
      loglike.diff2 = prior.new2 + post.new2 - prior.old2 - post.old2 
      if (log(runif(1, 0, 1)) < loglike.diff2) {
        delta = delta.new
        prior.old2 = prior.new2
        post.old2 = post.new2
      } 
      if (i4 > burn[2] & i4 %% (thinning[2]) == 0) {
        deltamat[(i4 - burn[2]) / thinning[2]] = delta
      }
    }
  }
  
  if (plot) {
    par(mfrow = c(2, 2))
    burnin = 1:(0.1 * out)
    matplot(alpha.mat[-burnin], lty = 1, type = "l", col = 1L:nIP, 
            xlab = "(Outer) Iterations", ylab = "alpha")
    abline(h = mean(alpha.mat[-burnin]), lty = 1, col = 1L:nIP)
    title("Convergence of Optimized alpha")
    plot(logWZ.mat[-burnin], type = "l", 
         xlab = "(Outer) Iterations", ylab = "logWZ")
    abline(h = mean(logWZ.mat[-burnin]), lty = 1)
    title("Convergence of logWZ")
    matplot(bmat[[1]][1,], lty = 1, col = 1L:P, type = "l", 
            main = "Traceplot of beta", xlab = "(Inner) Iterations", ylab = "")
    abline(h = mean(bmat[[1]][1,]), lty = 1, col = 1L)
    plot(deltamat, type = "l", 
         xlab = "(Outer) Iterations", ylab = "")
    abline(h = mean(deltamat), lty = 1)
    title("Traceplot of delta")
  }
  
  chain.final = list(C = currentC, Z = currentZ, B = bmat, D = deltamat,
                     iJi = iJi, sigma_Q =sigma_Q, alpha = alpha, mvec = mvec, 
                     proposal.var= proposal.var)
  return(chain.final)
}


#' @title IPTM_inference.data2
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
#' @param plot to plot the convergence diagnostics or not (TRUE/FALSE)
#' @param optimize to optimize alpha (Dirichlet concentration prior for document-topic distribution) or not (TRUE/FALSE)
#' @param initial list of initial values user wants to assign
#'
#' @return MCMC output containing IP assignment, topic assignment, and (beta, mu, delta) chain 
#'
#' @export
IPTM_inference.data2 = function(edge, node, textlist, vocabulary, nIP, K, sigma_Q, alpha, mvec, betas, nvec, prior.b.mean, prior.b.var, prior.delta, 
                               out, n_B, n_d, burn, thinning, netstat, plot = FALSE, optimize = FALSE, initial) {
  
  # trim the edge so that we only model edges after 384 hours
  timestamps = vapply(edge, function(d) {
    d[[3]]
  }, c(1))
  
  edge2 = which_int(384, timestamps) : length(edge)
  timeinc = c(timestamps[1], timestamps[-1] - timestamps[-length(timestamps)])
  
  # initialize alpha, mvec, delta, nvec, delta, lvec, and gammas
  W = length(vocabulary)
  phi = lapply(1L:K, function(k) {
    rdirichlet_cpp(1, betas * nvec)
  })
  delta = initial$D
  beta.old = initial$B
  # initialize C, theta and Z
  currentC = initial$C
  theta = rdirichlet_cpp(length(edge), alpha * mvec)
  currentZ = initial$Z
  p.d = t(vapply(seq(along = edge), function(d) {
    vapply(1L:nIP, function(IP) {
      sum(currentZ[[d]] %in% which(currentC == IP))
    }, c(1)) / length(currentZ[[d]])
  }, rep(1, nIP)))
  
  # initialize beta
  L = 3
  P = 1 * ("intercept" %in% netstat) + L * (2 * ("dyadic" %in% netstat) + 4 * ("triadic" %in% netstat) + 2 *("degree" %in% netstat))
  bmat = list()
  for (IP in 1L:nIP) {
    bmat[[IP]] = matrix(beta.old[[IP]], nrow = P, ncol = (n_B - burn[1]) / thinning[1])
  }
  deltamat = rep(delta,  (n_d - burn[2]) / thinning[2])
  proposal.var = lapply(1:nIP, function(IP){diag(P)})
  
  # to check the convergence  
  if (plot) {
    logWZ.mat = c()							  
    alpha.mat = c()
    entropy.mat = c()
  }
  
  #initialize the latent sender-receiver pairs
  iJi = lapply(seq(along = edge), function(d) {
    matrix(0, nrow = length(node), ncol = length(node))
  })
  lambda = list()
  LambdaiJi = list()
  observediJi = list()
  for (d in edge2) {
    history.t = History(edge, p.d, node, as.numeric(edge[[d-1]][3]) + 10^(-10))
    X = lapply(node, function(i) {
      Netstats(history.t, node, i, netstat)
    })
    XB = MultiplyXBList(X, beta.old)     
    lambda[[d]] = lambda_cpp(p.d[d,], XB)
    iJi[[d]] = matrix(rbinom(length(node)^2, 1, 0.5), nrow =length(node), ncol = length(node))
    diag(iJi[[d]]) = 0
  }
  
  #start outer iteration
  for (o in 1L:out) {
    output = list(C = currentC, Z = currentZ, B = bmat, D = deltamat)
    save(output, file = "output1.RData")
    print(o)
    if (optimize) {
      #update the hyperparameter alpha and mvec
      vec = AlphamvecOpt(K, currentZ[edge2], alpha, mvec, 5)
      alpha = sum(vec)
      mvec = vec / alpha
    }
    # Data augmentation
    for (d in edge2) {
      history.t = History(edge, p.d, node, as.numeric(edge[[d-1]][3]) + 10^(-10))
      X = lapply(node, function(i) {
        Netstats(history.t, node, i, netstat)
      })
      XB = MultiplyXBList(X, beta.old)     
      lambda[[d]] = lambda_cpp(p.d[d,], XB)
      #calculate the resampling probability	
      for (i in node[-as.numeric(edge[[d]][1])]) {
        for (j in sample(node[-i], length(node) - 1)) {
          probij = DataAug_cpp_Gibbs(iJi[[d]][i, ], lambda[[d]][i,], lapply(XB, function(IP) {IP[i,]}), p.d[d, ], delta, timeinc[d], j)
          iJi[[d]][i, j] = multinom_vec(1, probij) - 1
        }
      }
      iJi[[d]][as.numeric(edge[[d]][1]),] = tabulateC(as.numeric(unlist(edge[[d]][2])), length(node))
      LambdaiJi[[d]] = lambdaiJi(p.d[d,], XB, iJi[[d]])
      observediJi[[d]] = LambdaiJi[[d]][as.numeric(edge[[d]][1])]
    }	 
    
    # Z update
    textlist.raw = unlist(textlist[edge2])
    table.W = lapply(1L:K, function(k) {
      tabulateC(textlist.raw[which(unlist(currentZ[edge2]) == k)], W)
    })    
    for (d in edge2) {
      textlist.d = textlist[[d]]
      if (length(textlist.d) > 0) {
        topicpart.d = TopicInEqZ(K, currentZ[[d]], alpha, mvec, d)
        wordpart.d = WordInEqZ(K, textlist.d, table.W, betas, nvec)
      } else {
        topicpart.d = 0
        wordpart.d = matrix(0, nrow = length(currentZ[[d]]), ncol = K)
      }
      edgepart.d = EdgeInEqZ_Gibbs(iJi[[d]], lambda[[d]], delta)
      timepart.d = TimeInEqZ(LambdaiJi[[d]], timeinc[d])
      observed.d = ObservedInEqZ(observediJi[[d]]) 
      fixedpart = topicpart.d + edgepart.d + timepart.d + observed.d 
      for (w in 1L:length(currentZ[[d]])) {
        print(c(wordpart.d))
        const.Z = fixedpart + wordpart.d[w, ]
        const.Z = const.Z - max(const.Z)
        zw.old = currentZ[[d]][w]
        zw.new = multinom_vec(1, exp(const.Z))
        if (zw.new != zw.old) {
          currentZ[[d]][w] = zw.new
          topicpart.d = TopicInEqZ(K, currentZ[[d]], alpha, mvec, d)
          if (length(textlist.d) > 0) {	
            wordpart.d = WordInEqZ(K, textlist.d, table.W, betas, nvec)
          }
          table.W = lapply(1L:K, function(k) {
            tabulateC(textlist.raw[which(unlist(currentZ[edge2]) == k)], W)
          })
          p.d[d, ] = vapply(1L:nIP, function(IP) {
            sum(currentZ[[d]] %in% which(currentC == IP))
          }, c(1)) / length(currentZ[[d]])
          LambdaiJi[[d]] = lambdaiJi(p.d[d,], XB, iJi[[d]])
          observediJi[[d]] = LambdaiJi[[d]][as.numeric(edge[[d]][1])]
        }
      }
    }		
    # C update given Z and B - withinning each document d
    # for (k in sort(unique(unlist(currentZ[edge2])))) { 
    # const.C = rep(NA, nIP)
    # for (IP in 1:nIP) {
    # currentC[k] = IP
    # p.d = t(vapply(seq(along = edge), function(d) {
    # vapply(1L:nIP, function(IP) {
    # sum(currentZ[[d]] %in% which(currentC == IP))
    # }, c(1)) / length(currentZ[[d]])
    # }, rep(1, nIP)))
    # for (d in edge2) {
    # history.t = History(edge, p.d, node, as.numeric(edge[[d-1]][3]) + 10^(-10))
    # X = lapply(node, function(i) {
    # Netstats(history.t, node, i, netstat)
    # })
    # XB = MultiplyXBList(X, beta.old)    
    # lambda[[d]] = lambda_cpp(p.d[d,], XB)
    # LambdaiJi[[d]] = lambdaiJi(p.d[d,], XB, iJi[[d]])
    # observediJi[[d]] = LambdaiJi[[d]][as.numeric(edge[[d]][1])]
    # }
    # const.C[IP] = sum(vapply(edge2, function(d) {
    # EdgeInEqZ_Gibbs(iJi[[d]], lambda[[d]], delta) + 
    # TimeInEqZ(LambdaiJi[[d]], timeinc[d]) + 
    # ObservedInEqZ(observediJi[[d]]) 
    # }, c(1))) / length(edge2)
    # }
    # const.C = const.C - max(const.C)
    # expconst.C = exp(const.C)
    # if (Inf %in% expconst.C) {
    # expconst.C[which(expconst.C == Inf)] = exp(700)
    # }
    # currentC[k] = multinom_vec(1, expconst.C)
    # }
    
    if (plot) {
      entropy.mat = c(entropy.mat, entropy.empirical(currentC))
      alpha.mat = rbind(alpha.mat, alpha)
      logWZ.mat = c(logWZ.mat, logWZ(K, currentZ[edge2], textlist[edge2], table.W, alpha, mvec, betas, nvec))
    }
    
    p.d = t(vapply(seq(along = edge), function(d) {
      vapply(1L:nIP, function(IP) {
        sum(currentZ[[d]] %in% which(currentC == IP))
      }, c(1)) / length(currentZ[[d]])
    }, rep(1, nIP)))  
    for (d in edge2) {
      history.t = History(edge, p.d, node, as.numeric(edge[[d-1]][3]) + 10^(-10))
      X = lapply(node, function(i) {
        Netstats(history.t, node, i, netstat)
      })
      XB = MultiplyXBList(X, beta.old)   
      lambda[[d]] = lambda_cpp(p.d[d,], XB)
      LambdaiJi[[d]] = lambdaiJi(p.d[d,], XB, iJi[[d]])
      observediJi[[d]] = LambdaiJi[[d]][as.numeric(edge[[d]][1])]
    }
    
    # beta update
    prior.old1 = sum(vapply(1L:nIP, function(IP) {
      dmvnorm(beta.old[[IP]], prior.b.mean, prior.b.var, log = TRUE)
    }, c(1))) 
    post.old1 = sum(vapply(edge2, function(d) {
      EdgeInEqZ_Gibbs(iJi[[d]], lambda[[d]], delta) + TimeInEqZ(LambdaiJi[[d]], timeinc[d]) +
        ObservedInEqZ(observediJi[[d]])
    }, c(1))) / length(edge2)
    
    if (o != 1) {
      accept.rates = c(length(unique(bmat[[1]][1,])) / ncol(bmat[[1]]), length(unique(deltamat)) / length(deltamat))
      sigma_Q = adaptive_MH(sigma_Q, accept.rates, update_size = 0.2 * sigma_Q)
      if (accept.rates[1] > 1 / ncol(bmat[[1]])) { 
        for (IP in 1:nIP) {
          proposal.var[[IP]] = var(t(bmat[[IP]]))
        }
      }
    }
    for (i3 in 1L:n_B) {
      beta.new = lapply(1L:nIP, function(IP) {
        rmvnorm(1, beta.old[[IP]], sigma_Q[1] * proposal.var[[IP]])
      }) 
      for (d in edge2) {
        XB = MultiplyXBList(X, beta.new)
        lambda[[d]] = lambda_cpp(p.d[d,], XB)    
        LambdaiJi[[d]] = lambdaiJi(p.d[d,], XB, iJi[[d]])
        observediJi[[d]] = LambdaiJi[[d]][as.numeric(edge[[d]][1])]
      }
      prior.new1 = sum(vapply(1L:nIP, function(IP) {
        dmvnorm(beta.new[[IP]], prior.b.mean, prior.b.var, log = TRUE)
      }, c(1))) 
      post.new1 = sum(vapply(edge2, function(d) {
        EdgeInEqZ_Gibbs(iJi[[d]], lambda[[d]], delta) + TimeInEqZ(LambdaiJi[[d]], timeinc[d]) + 
          ObservedInEqZ(observediJi[[d]])
      }, c(1))) / length(edge2)
      loglike.diff = prior.new1 + post.new1 - prior.old1 - post.old1
      if (log(runif(1, 0, 1)) < loglike.diff) {
        for (IP in 1L:nIP) {
          beta.old[[IP]] = beta.new[[IP]]
        }
        prior.old1 = prior.new1
        post.old1 = post.new1
      }
      if (i3 > burn[1] & i3 %% (thinning[1]) == 0) {
        for (IP in 1L:nIP) {
          bmat[[IP]][ , (i3 - burn[1]) / thinning[1]] = beta.old[[IP]]
        }
      }
    }
    
    #delta update
    for (d in edge2) {
      history.t = History(edge, p.d, node, as.numeric(edge[[d-1]][3]) + 10^(-10))
      X = lapply(node, function(i) {
        Netstats(history.t, node, i, netstat)
      })
      XB = MultiplyXBList(X, beta.old)
      lambda[[d]] = lambda_cpp(p.d[d,], XB)  
    }
    prior.old2 = dnorm(delta, prior.delta[1], sqrt(prior.delta[2]), log = TRUE)
    post.old2 = sum(vapply(edge2, function(d) {
      EdgeInEqZ_Gibbs(iJi[[d]], lambda[[d]], delta)
    }, c(1))) / length(edge2)
    for (i4 in 1L:n_d) {
      delta.new = rnorm(1, delta, sqrt(sigma_Q[2]))
      prior.new2 = dnorm(delta.new, prior.delta[1], sqrt(prior.delta[2]), log = TRUE)
      post.new2 = sum(vapply(edge2, function(d) {
        EdgeInEqZ_Gibbs(iJi[[d]], lambda[[d]], delta.new)
      }, c(1))) / length(edge2)
      loglike.diff2 = prior.new2 + post.new2 - prior.old2 - post.old2 
      if (log(runif(1, 0, 1)) < loglike.diff2) {
        delta = delta.new
        prior.old2 = prior.new2
        post.old2 = post.new2
      } 
      if (i4 > burn[2] & i4 %% (thinning[2]) == 0) {
        deltamat[(i4 - burn[2]) / thinning[2]] = delta
      }
    }
  }
  
  if (plot) {
    par(mfrow = c(2, 2))
    burnin = 1:(0.1 * out)
    matplot(alpha.mat[-burnin], lty = 1, type = "l", col = 1L:nIP, 
            xlab = "(Outer) Iterations", ylab = "alpha")
    abline(h = mean(alpha.mat[-burnin]), lty = 1, col = 1L:nIP)
    title("Convergence of Optimized alpha")
    plot(logWZ.mat[-burnin], type = "l", 
         xlab = "(Outer) Iterations", ylab = "logWZ")
    abline(h = mean(logWZ.mat[-burnin]), lty = 1)
    title("Convergence of logWZ")
    matplot(bmat[[1]][1,], lty = 1, col = 1L:P, type = "l", 
            main = "Traceplot of beta", xlab = "(Inner) Iterations", ylab = "")
    abline(h = mean(bmat[[1]][1,]), lty = 1, col = 1L)
    plot(deltamat, type = "l", 
         xlab = "(Outer) Iterations", ylab = "")
    abline(h = mean(deltamat), lty = 1)
    title("Traceplot of delta")
  }
  
  chain.final = list(C = currentC, Z = currentZ, B = bmat, D = deltamat, 
                     iJi = iJi, sigma_Q =sigma_Q, alpha = alpha, mvec = mvec, 
                     proposal.var= proposal.var)
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
#' @param plot to plot the convergence diagnostics or not (TRUE/FALSE)
#' @param optimize to optimize alpha (Dirichlet concentration prior for document-topic distribution) or not (TRUE/FALSE)
#' @param initial initial values to be used from Forward sampling
#'
#' @return MCMC output containing IP assignment, topic assignment, and (beta, mu, delta) chain 
#'
#' @export
IPTM_inference.Schein = function(edge, node, textlist, vocabulary, nIP, K, sigma_Q, alpha, mvec, betas, nvec, prior.b.mean, prior.b.var, prior.delta, 
					out, n_B, n_d, burn, thinning, netstat, plot = FALSE, optimize = FALSE, initial) {
   
  # trim the edge so that we only model edges after 384 hours
	timestamps = vapply(edge, function(d) {
  			  d[[3]]
 			  }, c(1))

    edge2 = which_int(384, timestamps) : length(edge)
    timeinc = c(timestamps[1], timestamps[-1] - timestamps[-length(timestamps)])
   
  # initialize alpha, mvec, delta, nvec, delta, lvec, and gammas
 	W = length(vocabulary)
 	delta = initial$D
	beta.old = initial$B
	# initialize C, theta and Z
	currentC = initial$C
     currentZ = initial$Z
 	 p.d = t(vapply(seq(along = edge), function(d) {
    	vapply(1L:nIP, function(IP) {
      	sum(currentZ[[d]] %in% which(currentC == IP))
  	 		 }, c(1)) / length(currentZ[[d]])
 		 }, rep(1, nIP)))

    # initialize beta
    L = 3
    P = 1 * ("intercept" %in% netstat) + L * (2 * ("dyadic" %in% netstat) + 4 * ("triadic" %in% netstat) + 2 *("degree" %in% netstat))
    bmat = list()
	for (IP in 1L:nIP) {
		bmat[[IP]] = matrix(beta.old[[IP]], nrow = P, ncol = (n_B - burn) / thinning)
  	}
  	deltamat = rep(delta, n_d)
    proposal.var = lapply(1:nIP, function(IP){diag(P)})

    # to check the convergence  
    if (plot) {
     	logWZ.mat = c()							  
     	alpha.mat = c()
     	entropy.mat = c()
    }

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

    #start outer iteration
    for (o in 1L:out) {
      
      if (optimize) {
      #update the hyperparameter alpha and mvec
      vec = AlphamvecOpt(K, currentZ[edge2], alpha, mvec, 5)
      alpha = sum(vec)
      mvec = vec / alpha
      }
     # Data augmentation
      for (d in edge2) {
   	 	history.t = History(edge, p.d, node, as.numeric(edge[[d-1]][3]) + 10^(-10))
   	 	X = lapply(node, function(i) {
  	        Netstats(history.t, node, i, netstat)
            })
   	 	XB = MultiplyXBList(X, beta.old)     
		lambda[[d]] = lambda_cpp(p.d[d,], XB)
		#calculate the resampling probability	
		for (i in node[-as.numeric(edge[[d]][1])]) {
			for (j in sample(node[-i], length(node) - 1)) {
				probij = DataAug_cpp_Gibbs(iJi[[d]][i, ], lambda[[d]][i,], lapply(XB, function(IP) {IP[i,]}), p.d[d, ], delta, timeinc[d], j)
				iJi[[d]][i, j] = multinom_vec(1, probij) - 1		
				}
		}
		iJi[[d]][as.numeric(edge[[d]][1]),] = tabulateC(as.numeric(unlist(edge[[d]][2])), length(node))
		LambdaiJi[[d]] = lambdaiJi(p.d[d,], XB, iJi[[d]])
		observediJi[[d]] = LambdaiJi[[d]][as.numeric(edge[[d]][1])]
		}	 
	
    	  # Z update

   	  textlist.raw = unlist(textlist[edge2])
      table.W = lapply(1L:K, function(k) {
      			 tabulateC(textlist.raw[which(unlist(currentZ[edge2]) == k)], W)
      			 })    
      	 for (d in edge2) {
      		 textlist.d = textlist[[d]]
        		if (length(textlist.d) > 0) {
        			topicpart.d = TopicInEqZ(K, currentZ[[d]], alpha, mvec, d)
       			wordpart.d = WordInEqZ(K, textlist.d, table.W, betas, nvec)
        	 } else {
        			topicpart.d = 0
        			wordpart.d = matrix(0, nrow = length(currentZ[[d]]), ncol = K)
        		}
        		edgepart.d = EdgeInEqZ_Gibbs(iJi[[d]], lambda[[d]], delta)
        		timepart.d = TimeInEqZ(LambdaiJi[[d]], timeinc[d])
        		observed.d = ObservedInEqZ(observediJi[[d]]) 
        		fixedpart = topicpart.d + edgepart.d + timepart.d + observed.d 
        		for (w in 1L:length(currentZ[[d]])) {
          		const.Z = fixedpart + wordpart.d[w, ]
          		const.Z = const.Z - max(const.Z)
          		zw.old = currentZ[[d]][w]
          		zw.new = multinom_vec(1, exp(const.Z))
          		if (zw.new != zw.old) {
            			currentZ[[d]][w] = zw.new
            			topicpart.d = TopicInEqZ(K, currentZ[[d]], alpha, mvec, d)
            			if (length(textlist.d) > 0) {	
            				wordpart.d = WordInEqZ(K, textlist.d, table.W, betas, nvec)
            			}
            			table.W = lapply(1L:K, function(k) {
      				  tabulateC(textlist.raw[which(unlist(currentZ[edge2]) == k)], W)
      				})
      				p.d[d, ] = vapply(1L:nIP, function(IP) {
	 					sum(currentZ[[d]] %in% which(currentC == IP))
	 				}, c(1)) / length(currentZ[[d]])
      				LambdaiJi[[d]] = lambdaiJi(p.d[d,], XB, iJi[[d]])
         			observediJi[[d]] = LambdaiJi[[d]][as.numeric(edge[[d]][1])]
          		}
        	}
       	}		
      # C update given Z and B - withinning each document d
      for (k in sort(unique(unlist(currentZ[edge2])))) { 
        const.C = rep(NA, nIP)
        for (IP in 1:nIP) {
          	currentC[k] = IP
          	p.d = t(vapply(seq(along = edge), function(d) {
    				vapply(1L:nIP, function(IP) {
      			sum(currentZ[[d]] %in% which(currentC == IP))
  	 				 }, c(1)) / length(currentZ[[d]])
 					}, rep(1, nIP)))
          	for (d in edge2) {
           		history.t = History(edge, p.d, node, as.numeric(edge[[d-1]][3]) + 10^(-10))
    	   	   		X = lapply(node, function(i) {
               		Netstats(history.t, node, i, netstat)
               		})
    	       		XB = MultiplyXBList(X, beta.old)    
           		lambda[[d]] = lambda_cpp(p.d[d,], XB)
		       	LambdaiJi[[d]] = lambdaiJi(p.d[d,], XB, iJi[[d]])
           		observediJi[[d]] = LambdaiJi[[d]][as.numeric(edge[[d]][1])]
          	}
          const.C[IP] = sum(vapply(edge2, function(d) {
          				EdgeInEqZ_Gibbs(iJi[[d]], lambda[[d]], delta) + 
          				TimeInEqZ(LambdaiJi[[d]], timeinc[d]) + 
    						  ObservedInEqZ(observediJi[[d]]) 
          				}, c(1))) / length(edge2)
      	}
        const.C = const.C - max(const.C)
        currentC[k] = multinom_vec(1, exp(const.C))
     }
     
    p.d = t(vapply(seq(along = edge), function(d) {
    				vapply(1L:nIP, function(IP) {
      				sum(currentZ[[d]] %in% which(currentC == IP))
  	 				 }, c(1)) / length(currentZ[[d]])
 					}, rep(1, nIP)))  
    for (d in edge2) {
        	history.t = History(edge, p.d, node, as.numeric(edge[[d-1]][3]) + 10^(-10))
    	    X = lapply(node, function(i) {
            Netstats(history.t, node, i, netstat)
       		})
    	    XB = MultiplyXBList(X, beta.old)   
    	    lambda[[d]] = lambda_cpp(p.d[d,], XB)
	    LambdaiJi[[d]] = lambdaiJi(p.d[d,], XB, iJi[[d]])
        observediJi[[d]] = LambdaiJi[[d]][as.numeric(edge[[d]][1])]
	}
	 
	 # beta update
	prior.old1 = sum(vapply(1L:nIP, function(IP) {
		  		 dmvnorm(beta.old[[IP]], prior.b.mean, prior.b.var, log = TRUE)
		  		 }, c(1))) 
	post.old1 = sum(vapply(edge2, function(d) {
	     	    EdgeInEqZ_Gibbs(iJi[[d]], lambda[[d]], delta) + TimeInEqZ(LambdaiJi[[d]], timeinc[d]) +
    			 	ObservedInEqZ(observediJi[[d]])
    			 	}, c(1))) / length(edge2)
  
    if (o != 1) {
    		accept.rates = c(length(unique(bmat[[1]][1,])) / ncol(bmat[[1]]), length(unique(deltamat)) / n_d)
    	sigma_Q = adaptive_MH(sigma_Q, accept.rates, update_size = 0.1 * sigma_Q)
    	if (accept.rates[1] > 1 / ncol(bmat[[1]])) {
    		for (IP in 1:nIP) {
     		proposal.var[[IP]] = var(t(bmat[[IP]]))
      	}
      }
   }
    for (i3 in 1L:n_B) {
    	if (i3 %% 500 == 0 ){print(i3)}
    		beta.new = lapply(1L:nIP, function(IP) {
          		   rmvnorm(1, beta.old[[IP]], sigma_Q[1] * proposal.var[[IP]])
         		   }) 
        for (d in edge2) {
           XB = MultiplyXBList(X, beta.new)
           lambda[[d]] = lambda_cpp(p.d[d,], XB)    
           LambdaiJi[[d]] = lambdaiJi(p.d[d,], XB, iJi[[d]])
	       observediJi[[d]] = LambdaiJi[[d]][as.numeric(edge[[d]][1])]
        }
        prior.new1 = sum(vapply(1L:nIP, function(IP) {
        				 dmvnorm(beta.new[[IP]], prior.b.mean, prior.b.var, log = TRUE)
        				 }, c(1))) 
        post.new1 = sum(vapply(edge2, function(d) {
    			   	EdgeInEqZ_Gibbs(iJi[[d]], lambda[[d]], delta) + TimeInEqZ(LambdaiJi[[d]], timeinc[d]) + 
    			    	ObservedInEqZ(observediJi[[d]])
    			    	}, c(1))) / length(edge2)
    	loglike.diff = prior.new1 + post.new1 - prior.old1 - post.old1
        if (log(runif(1, 0, 1)) < loglike.diff) {
        		for (IP in 1L:nIP) {
         		beta.old[[IP]] = beta.new[[IP]]
         	}
         	prior.old1 = prior.new1
         	post.old1 = post.new1
        }
         if (i3 > burn & i3 %% (thinning) == 0) {
        		for (IP in 1L:nIP) {
           		bmat[[IP]][ , (i3 - burn) / thinning] = beta.old[[IP]]
           	}
    	}
    }
    
    #delta update
    for (d in edge2) {
    	history.t = History(edge, p.d, node, as.numeric(edge[[d-1]][3]) + 10^(-10))
    	X = lapply(node, function(i) {
            Netstats(history.t, node, i, netstat)
       		})
    	XB = MultiplyXBList(X, beta.old)
     	lambda[[d]] = lambda_cpp(p.d[d,], XB)  
    }
    prior.old2 = dnorm(delta, prior.delta[1], sqrt(prior.delta[2]), log = TRUE)
    post.old2 = sum(vapply(edge2, function(d) {
    				EdgeInEqZ_Gibbs(iJi[[d]], lambda[[d]], delta)
    				}, c(1))) / length(edge2)
 	for (i4 in 1L:n_d) {
        delta.new = rnorm(1, delta, sqrt(sigma_Q[2]))
        prior.new2 = dnorm(delta.new, prior.delta[1], sqrt(prior.delta[2]), log = TRUE)
        post.new2 = sum(vapply(edge2, function(d) {
        				EdgeInEqZ_Gibbs(iJi[[d]], lambda[[d]], delta.new)
        				}, c(1))) / length(edge2)
        loglike.diff2 = prior.new2 + post.new2 - prior.old2 - post.old2 
        if (log(runif(1, 0, 1)) < loglike.diff2) {
        	delta = delta.new
         	prior.old2 = prior.new2
            post.old2 = post.new2
        } 
            deltamat[i4] = delta
    }
 }
         
  if (plot) {
    par(mfrow = c(1, 2))
	matplot(bmat[[1]][7,], lty = 1, col = 1L:P, type = "l", 
	          main = "Traceplot of beta", xlab = "(Inner) Iterations", ylab = "")
	abline(h = mean(bmat[[1]][7,]), lty = 1, col = 1L)
	plot(deltamat, type = "l", 
	xlab = "(Outer) Iterations", ylab = "")
	abline(h = mean(deltamat), lty = 1)
	title("Traceplot of delta")
  }
     
  chain.final = list(C = currentC, Z = lapply(edge2, function(d) {currentZ[[d]]}), B = bmat, D = deltamat,
                     iJi = iJi, sigma_Q =sigma_Q, alpha = alpha, mvec = mvec, 
                     proposal.var= proposal.var)
  return(chain.final)
}


#' @title IPTM_inference.Schein2
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
#' @param plot to plot the convergence diagnostics or not (TRUE/FALSE)
#' @param optimize to optimize alpha (Dirichlet concentration prior for document-topic distribution) or not (TRUE/FALSE)
#' @param initial initial values to be used from Forward sampling
#'
#' @return MCMC output containing IP assignment, topic assignment, and (beta, mu, delta) chain 
#'
#' @export
IPTM_inference.Schein2 = function(edge, node, textlist, vocabulary, nIP, K, sigma_Q, alpha, mvec, betas, nvec, prior.b.mean, prior.b.var, prior.delta, 
					out, n_B, n_d, burn, thinning, netstat, plot = FALSE, optimize = FALSE, initial) {
   
  # trim the edge so that we only model edges after 384 hours
	timestamps = vapply(edge, function(d) {
  			  d[[3]]
 			  }, c(1))

    edge2 = which_int(384, timestamps) : length(edge)
    timeinc = c(timestamps[1], timestamps[-1] - timestamps[-length(timestamps)])
   
  # initialize alpha, mvec, delta, nvec, delta, lvec, and gammas
 	W = length(vocabulary)
 	delta = initial$D
	beta.old = initial$B
	# initialize C, theta and Z
	currentC = initial$C
     currentZ = initial$Z
 	 p.d = t(vapply(seq(along = edge), function(d) {
    	vapply(1L:nIP, function(IP) {
      	sum(currentZ[[d]] %in% which(currentC == IP))
  	 		 }, c(1)) / length(currentZ[[d]])
 		 }, rep(1, nIP)))

    # initialize beta
    L = 3
    P = 1 * ("intercept" %in% netstat) + L * (2 * ("dyadic" %in% netstat) + 4 * ("triadic" %in% netstat) + 2 *("degree" %in% netstat))
    bmat = list()
	for (IP in 1L:nIP) {
		bmat[[IP]] = matrix(beta.old[[IP]], nrow = P, ncol = (n_B - burn) / thinning)
  	}
  	deltamat = rep(delta, n_d)
    proposal.var = lapply(1:nIP, function(IP){diag(P)})

    # to check the convergence  
    if (plot) {
     	logWZ.mat = c()							  
     	alpha.mat = c()
     	entropy.mat = c()
    }

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

    #start outer iteration
    for (o in 1L:out) {
      
      if (optimize) {
      #update the hyperparameter alpha and mvec
      vec = AlphamvecOpt(K, currentZ[edge2], alpha, mvec, 5)
      alpha = sum(vec)
      mvec = vec / alpha
      }
     # Data augmentation
      for (d in edge2) {
   	 	history.t = History(edge, p.d, node, as.numeric(edge[[d-1]][3]) + 10^(-10))
   	 	X = lapply(node, function(i) {
  	        Netstats(history.t, node, i, netstat)
            })
   	 	XB = MultiplyXBList(X, beta.old)     
		lambda[[d]] = lambda_cpp(p.d[d,], XB)
		#calculate the resampling probability	
		for (i in node[-as.numeric(edge[[d]][1])]) {
			for (j in sample(node[-i], length(node) - 1)) {
				probij = DataAug_cpp_Gibbs(iJi[[d]][i, ], lambda[[d]][i,], lapply(XB, function(IP) {IP[i,]}), p.d[d, ], delta, timeinc[d], j)
				iJi[[d]][i, j] = multinom_vec(1, probij) - 1		
				}
		}
		iJi[[d]][as.numeric(edge[[d]][1]),] = tabulateC(as.numeric(unlist(edge[[d]][2])), length(node))
		LambdaiJi[[d]] = lambdaiJi(p.d[d,], XB, iJi[[d]])
		observediJi[[d]] = LambdaiJi[[d]][as.numeric(edge[[d]][1])]
		}	 
	
    	  # Z update

# #     	  textlist.raw = unlist(textlist[edge2])
    	  # table.W = lapply(1L:K, function(k) {
      			# tabulateC(textlist.raw[which(unlist(currentZ[edge2]) == k)], W)
      			# })    
      	# for (d in edge2) {
      		# textlist.d = textlist[[d]]
        		# if (length(textlist.d) > 0) {
        			# topicpart.d = TopicInEqZ(K, currentZ[[d]], alpha, mvec, d)
       			# wordpart.d = WordInEqZ(K, textlist.d, table.W, betas, nvec)
        		# } else {
        			# topicpart.d = 0
        			# wordpart.d = matrix(0, nrow = length(currentZ[[d]]), ncol = K)
        		# }
        		# edgepart.d = EdgeInEqZ_Gibbs(iJi[[d]], lambda[[d]], delta)
        		# timepart.d = TimeInEqZ(LambdaiJi[[d]], timeinc[d])
        		# observed.d = ObservedInEqZ(observediJi[[d]]) 
        		# fixedpart = topicpart.d + edgepart.d + timepart.d + observed.d 
        		# for (w in 1L:length(currentZ[[d]])) {
          		# const.Z = fixedpart + wordpart.d[w, ]
          		# const.Z = const.Z - max(const.Z)
          		# zw.old = currentZ[[d]][w]
          		# zw.new = multinom_vec(1, exp(const.Z))
          		# if (zw.new != zw.old) {
            			# currentZ[[d]][w] = zw.new
            			# topicpart.d = TopicInEqZ(K, currentZ[[d]], alpha, mvec, d)
            			# if (length(textlist.d) > 0) {	
            				# wordpart.d = WordInEqZ(K, textlist.d, table.W, betas, nvec)
            			# }
            			# table.W = lapply(1L:K, function(k) {
      				  # tabulateC(textlist.raw[which(unlist(currentZ[edge2]) == k)], W)
      				# })
      				# p.d[d, ] = vapply(1L:nIP, function(IP) {
	 					# sum(currentZ[[d]] %in% which(currentC == IP))
	 				# }, c(1)) / length(currentZ[[d]])
      				# LambdaiJi[[d]] = lambdaiJi(p.d[d,], XB, iJi[[d]])
         			# observediJi[[d]] = LambdaiJi[[d]][as.numeric(edge[[d]][1])]
          		# }
        		# }
       	# }		
      # # C update given Z and B - withinning each document d
      # for (k in sort(unique(unlist(currentZ[edge2])))) { 
        # const.C = rep(NA, nIP)
        # for (IP in 1:nIP) {
          	# currentC[k] = IP
          	# p.d = t(vapply(seq(along = edge), function(d) {
    				# vapply(1L:nIP, function(IP) {
      			# sum(currentZ[[d]] %in% which(currentC == IP))
  	 				 # }, c(1)) / length(currentZ[[d]])
 					 # }, rep(1, nIP)))
          	# for (d in edge2) {
           		# history.t = History(edge, p.d, node, as.numeric(edge[[d-1]][3]) + 10^(-10))
    	   	   		# X = lapply(node, function(i) {
               		# Netstats(history.t, node, i, netstat)
               		# })
    	       		# XB = MultiplyXBList(X, beta.old)    
           		# lambda[[d]] = lambda_cpp(p.d[d,], XB)
		       	# LambdaiJi[[d]] = lambdaiJi(p.d[d,], XB, iJi[[d]])
           		# observediJi[[d]] = LambdaiJi[[d]][as.numeric(edge[[d]][1])]
          	# }
          # const.C[IP] = sum(vapply(edge2, function(d) {
          				# EdgeInEqZ_Gibbs(iJi[[d]], lambda[[d]], delta) + 
          				# TimeInEqZ(LambdaiJi[[d]], timeinc[d]) + 
    						  # ObservedInEqZ(observediJi[[d]]) 
          				# }, c(1))) / length(edge2)
      	# }
        # const.C = const.C - max(const.C)
        # currentC[k] = multinom_vec(1, exp(const.C))
     # }
     
    # p.d = t(vapply(seq(along = edge), function(d) {
    				# vapply(1L:nIP, function(IP) {
      				# sum(currentZ[[d]] %in% which(currentC == IP))
  	 				 # }, c(1)) / length(currentZ[[d]])
 					 # }, rep(1, nIP)))  
    # for (d in edge2) {
        	# history.t = History(edge, p.d, node, as.numeric(edge[[d-1]][3]) + 10^(-10))
    	    # X = lapply(node, function(i) {
            # Netstats(history.t, node, i, netstat)
       		# })
    	    # XB = MultiplyXBList(X, beta.old)   
    	    # lambda[[d]] = lambda_cpp(p.d[d,], XB)
	    # LambdaiJi[[d]] = lambdaiJi(p.d[d,], XB, iJi[[d]])
        # observediJi[[d]] = LambdaiJi[[d]][as.numeric(edge[[d]][1])]
	# }
	 
	# # beta update
	# prior.old1 = sum(vapply(1L:nIP, function(IP) {
		  		 # dmvnorm(beta.old[[IP]], prior.b.mean, prior.b.var, log = TRUE)
		  		 # }, c(1))) 
	# post.old1 = sum(vapply(edge2, function(d) {
	     	    # EdgeInEqZ_Gibbs(iJi[[d]], lambda[[d]], delta) + TimeInEqZ(LambdaiJi[[d]], timeinc[d]) +
    			 	# ObservedInEqZ(observediJi[[d]])
    			 	# }, c(1))) / length(edge2)
  
    # if (o != 1) {
    		# accept.rates = c(length(unique(bmat[[1]][1,])) / ncol(bmat[[1]]), length(unique(deltamat)) / n_d)
    	# sigma_Q = adaptive_MH(sigma_Q, accept.rates, update_size = 0.1 * sigma_Q)
    	# if (accept.rates[1] > 1 / ncol(bmat[[1]])) {
    		# for (IP in 1:nIP) {
     		# proposal.var[[IP]] = var(t(bmat[[IP]]))
      	# }
      # }
   # }
    # for (i3 in 1L:n_B) {
    	# if (i3 %% 500 == 0 ){print(i3)}
    		# beta.new = lapply(1L:nIP, function(IP) {
          		   # rmvnorm(1, beta.old[[IP]], sigma_Q[1] * proposal.var[[IP]])
         		   # }) 
        # for (d in edge2) {
           # XB = MultiplyXBList(X, beta.new)
           # lambda[[d]] = lambda_cpp(p.d[d,], XB)    
           # LambdaiJi[[d]] = lambdaiJi(p.d[d,], XB, iJi[[d]])
	       # observediJi[[d]] = LambdaiJi[[d]][as.numeric(edge[[d]][1])]
        # }
        # prior.new1 = sum(vapply(1L:nIP, function(IP) {
        				 # dmvnorm(beta.new[[IP]], prior.b.mean, prior.b.var, log = TRUE)
        				 # }, c(1))) 
        # post.new1 = sum(vapply(edge2, function(d) {
    			   		# EdgeInEqZ_Gibbs(iJi[[d]], lambda[[d]], delta) + TimeInEqZ(LambdaiJi[[d]], timeinc[d]) + 
    			    	# ObservedInEqZ(observediJi[[d]])
    			    	# }, c(1))) / length(edge2)
    		# loglike.diff = prior.new1 + post.new1 - prior.old1 - post.old1
        # if (log(runif(1, 0, 1)) < loglike.diff) {
        		# for (IP in 1L:nIP) {
         		# beta.old[[IP]] = beta.new[[IP]]
         	# }
         	# prior.old1 = prior.new1
         	# post.old1 = post.new1
        # }
         # if (i3 > burn & i3 %% (thinning) == 0) {
        		# for (IP in 1L:nIP) {
           		# bmat[[IP]][ , (i3 - burn) / thinning] = beta.old[[IP]]
           	# }
    		# }
    # }
    
    # #delta update
    # for (d in edge2) {
    		# history.t = History(edge, p.d, node, as.numeric(edge[[d-1]][3]) + 10^(-10))
    	# X = lapply(node, function(i) {
            # Netstats(history.t, node, i, netstat)
       		# })
    	 	# XB = MultiplyXBList(X, beta.old)
     	# lambda[[d]] = lambda_cpp(p.d[d,], XB)  
    # }
    # prior.old2 = dnorm(delta, prior.delta[1], sqrt(prior.delta[2]), log = TRUE)
    # post.old2 = sum(vapply(edge2, function(d) {
    				# EdgeInEqZ_Gibbs(iJi[[d]], lambda[[d]], delta)
    				# }, c(1))) / length(edge2)
 	# for (i4 in 1L:n_d) {
        # delta.new = rnorm(1, delta, sqrt(sigma_Q[2]))
        # prior.new2 = dnorm(delta.new, prior.delta[1], sqrt(prior.delta[2]), log = TRUE)
        # post.new2 = sum(vapply(edge2, function(d) {
        				# EdgeInEqZ_Gibbs(iJi[[d]], lambda[[d]], delta.new)
        				# }, c(1))) / length(edge2)
        # loglike.diff2 = prior.new2 + post.new2 - prior.old2 - post.old2 
        # if (log(runif(1, 0, 1)) < loglike.diff2) {
        		# delta = delta.new
         	# prior.old2 = prior.new2
         	# post.old2 = post.new2
        # } 
            # deltamat[i4] = delta
    	# }
 }
         
  if (plot) {
    par(mfrow = c(1, 2))
	matplot(bmat[[1]][7,], lty = 1, col = 1L:P, type = "l", 
	          main = "Traceplot of beta", xlab = "(Inner) Iterations", ylab = "")
	abline(h = mean(bmat[[1]][7,]), lty = 1, col = 1L)
	plot(deltamat, type = "l", 
	xlab = "(Outer) Iterations", ylab = "")
	abline(h = mean(deltamat), lty = 1)
	title("Traceplot of delta")
  }
     
  chain.final = list(C = currentC, Z = lapply(edge2, function(d) {currentZ[[d]]}), B = bmat, D = deltamat,
                     iJi = iJi, sigma_Q =sigma_Q, alpha = alpha, mvec = mvec, 
                     proposal.var= proposal.var)
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
  for(b in 1L:P){
    combined.data[[b]] = sapply(1L:nIP, function(c){
      cbind(Bchain[[c]][b,])
    })
  }
  forbox = melt(combined.data)
  boxplot = boxplot(forbox$value ~ forbox$X2 + forbox$L1,
                    at = c(sapply(1L:P, function(x){
                      ((nIP + 1) * x - nIP):((nIP + 1) * x - 1)
                    })),
                    col = gray.colors(nIP), axes = FALSE, 
                    main = "Comparison of beta coefficients for different IPs")
  abline(h = 0, lty = 1, col = "red")
  axis(2, labels = TRUE)
  box()
  axis(side = 1, line = 0.5, lwd = 0, 
       at = c(sapply(1L:P, function(x){
         median(((nIP + 1) * x - nIP):((nIP + 1) * x - 1))
       })), labels = 1:25)
  legend(locator(1), c(paste("IP", 1L:nIP)), col = gray.colors(nIP), pch = 15)
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
  colnames(topic.dist) = c(1L:K)
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
    colnames(topic.dist) = c(1L:K)
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

#' @title GenerateDocs
#' @description Generate a collection of documents according to the generative process of IPTM
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
#' @param latentiJi inferred latent sender-receiver matrix (only used for backward sampling)
#' @param forward Logigal indicating whether we are generating forward samples
#' @param backward_init Logigal indicating whether we are generating initial backward samples
#' @param backward Logigal indicating whether we are generating backward samples
#' @param base Logical indicating whether or not we are generating base edges (< 384)
#' @param backward.edge edge from previous backward sample
#'
#' @return generated edge and text, parameter b used to generate those, and base (if base == TRUE)
#'
#' @export
GenerateDocs = function(nDocs, node, vocabulary, nIP, K, nwords, alpha, mvec, betas, nvec,
                        b, delta, currentC, netstat, base.edge, base.text, topic_token_assignments = NULL, latentiJi = list(),
                        forward = FALSE, backward_init = FALSE, backward = FALSE, base = FALSE, backward.edge = NULL) {
  W = length(vocabulary)
  phi = lapply(1L:K, function(k) {
    rdirichlet_cpp(1, betas * nvec)
  })
  L = 3
  P = 1 * ("intercept" %in% netstat) + L * (2 * ("dyadic" %in% netstat) + 4 * ("triadic" %in% netstat) + 2 *("degree" %in% netstat))
  
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
  	 p.d[d, ] = vapply(1L:nIP, function(IP) {
      sum(names(text[[d]]) %in% which(currentC == IP))
    }, c(1)) / max(1, length(text[[d]]))
  }
  }
  word_type_topic_counts = matrix(0, W, K)
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
    
    p.d[base.length + d, ] = vapply(1L:nIP, function(IP) {
      sum(topic.d %in% which(currentC == IP))
    }, c(1)) / max(1, N.d)
    if (base & t.d < 384) {
      history.t = lapply(1:nIP, function(IP) {
        			  lapply(1:3, function(l){
         		  matrix(0, length(node), length(node))
       			  })
     			  })
    } else {	
      history.t = History(edge, p.d, node, t.d + 10^(-10))
    }
    X = lapply(node, function(i) {
      	Netstats(history.t, node, i, netstat)
   		})
    XB = MultiplyXBList(X, b)     
    lambda = lambda_cpp(p.d[base.length + d,], XB)
    
    if (!backward) {
    iJi = rbinom_mat((delta * lambda) / (delta * lambda + 1))
    while (sum(iJi) == 0) {
      iJi = rbinom_mat((delta * lambda) / (delta * lambda + 1))
    }
    } else {
    		iJi = latentiJi[[d]]
    		observedi = backward.edge[[base.length + d]][[1]]
    		iJi[observedi, ] = rbinom_mat((delta * lambda) / (delta * lambda + 1))[observedi,]
    		while (sum(iJi[observedi, ]) == 0) {
      		iJi[observedi, ] = rbinom_mat((delta * lambda) / (delta * lambda + 1))[observedi,]
    		}
    	}
    LambdaiJi = lambdaiJi(p.d[base.length + d,], XB, iJi)
    Time.inc = vapply(LambdaiJi, function(lambda) {
      rexp(1, lambda)
    }, c(1))
    i.d = which(Time.inc == min(Time.inc[!is.na(Time.inc)]))
    j.d = which(iJi[i.d,] == 1)
    t.d = t.d + Time.inc[i.d]
    edge[[base.length + d]] = list(sender = i.d, receiver = j.d, timestamp = t.d)		
  }
  if (forward) {
    edge = edge[-(1:base.length)]
    text = text[-(1:base.length)]
  }
  if (base == TRUE & t.d > 384) {
    cutoff = which_int(384, vapply(1:length(edge), function(d) {edge[[d]][[3]]}, c(1))) - 1
    edge = edge[1:cutoff]
    text = text[1:cutoff]
  }
  return(list(edge = edge, text = text, b = b, d = delta, base = base.length))							
} 


#' @title CollapsedGenerateDocs
#' @description Generate a collection of documents according to the generative process of IPTM
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
#' @param latentiJi inferred latent sender-receiver matrix (only used for backward sampling)
#' @param forward Logigal indicating whether we are generating forward samples
#' @param backward_init Logigal indicating whether we are generating initial backward samples
#' @param backward Logigal indicating whether we are generating backward samples
#' @param base Logical indicating whether or not we are generating base edges (< 384)
#' @param backward.edge edge from previous backward sample
#'
#' @return generated edge and text, parameter b used to generate those, and base (if base == TRUE)
#'
#' @export
CollapsedGenerateDocs = function(nDocs, node, vocabulary, nIP, K, nwords, alpha, mvec, betas, nvec,
                        b, delta, currentC, netstat, base.edge, base.text, topic_token_assignments = NULL, latentiJi = list(),
                        forward = FALSE, backward_init = FALSE, backward = FALSE, base = FALSE, backward.edge = NULL) {
  W = length(vocabulary)
  phi = lapply(1L:K, function(k) {
    rdirichlet_cpp(1, betas * nvec)
  })
  L = 3
  P = 1 * ("intercept" %in% netstat) + L * (2 * ("dyadic" %in% netstat) + 4 * ("triadic" %in% netstat) + 2 *("degree" %in% netstat))
  
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
  	 p.d[d, ] = vapply(1L:nIP, function(IP) {
      sum(names(text[[d]]) %in% which(currentC == IP))
    }, c(1)) / max(1, length(text[[d]]))
  }
  }
  word_type_topic_counts = matrix(0, W, K)
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
    
    p.d[base.length + d, ] = vapply(1L:nIP, function(IP) {
      sum(topic.d %in% which(currentC == IP))
    }, c(1)) / max(1, N.d)
    if (base & t.d < 384) {
      history.t = lapply(1:nIP, function(IP) {
        			  lapply(1:3, function(l){
         		  matrix(0, length(node), length(node))
       			  })
     			  })
    } else {	
      history.t = History(edge, p.d, node, t.d + 10^(-10))
    }
    X = lapply(node, function(i) {
      	Netstats(history.t, node, i, netstat)
   		})
    XB = MultiplyXBList(X, b)     
    lambda = lambda_cpp(p.d[base.length + d,], XB)
    
    if (!backward) {
    iJi = rbinom_mat((delta * lambda) / (delta * lambda + 1))
    while (sum(iJi) == 0) {
      iJi = rbinom_mat((delta * lambda) / (delta * lambda + 1))
    }
    } else {
    		iJi = latentiJi[[d]]
    		observedi = backward.edge[[base.length + d]][[1]]
    		iJi[observedi, ] = rbinom_mat((delta * lambda) / (delta * lambda + 1))[observedi,]
    		while (sum(iJi[observedi, ]) == 0) {
      		iJi[observedi, ] = rbinom_mat((delta * lambda) / (delta * lambda + 1))[observedi,]
    		}
    	}
    LambdaiJi = lambdaiJi(p.d[base.length + d,], XB, iJi)
   	LambdaiJi[is.na(LambdaiJi)] = 0
    i.d = multinom_vec(1, LambdaiJi)
    j.d = which(iJi[i.d,] == 1)
    t.d = t.d + rexp(1, sum(LambdaiJi))
    edge[[base.length + d]] = list(sender = i.d, receiver = j.d, timestamp = t.d)		
  }
  if (forward) {
    edge = edge[-(1:base.length)]
    text = text[-(1:base.length)]
  }
  if (base == TRUE & t.d > 384) {
    cutoff = which_int(384, vapply(1:length(edge), function(d) {edge[[d]][[3]]}, c(1))) - 1
    edge = edge[1:cutoff]
    text = text[1:cutoff]
  }
  return(list(edge = edge, text = text, b = b, d = delta, base = base.length))							
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
  phi = lapply(1L:K, function(k) {
    rdirichlet_cpp(1, betas * nvec)
  })
  L = 3
  P = 1 * ("intercept" %in% netstat) + L * (2 * ("dyadic" %in% netstat) + 4 * ("triadic" %in% netstat) + 2 *("degree" %in% netstat))
  
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
  	 p.d[d, ] = vapply(1L:nIP, function(IP) {
      sum(names(text[[d]]) %in% which(currentC == IP))
    }, c(1)) / max(1, length(text[[d]]))
  }
  }
 history.t = lapply(1:nIP, function(IP) {
        	 	 lapply(1:3, function(l){
         		matrix(0, length(node), length(node))
       			})
     		 })  
  word_type_topic_counts = matrix(0, W, K)
  iJi = matrix(0, length(node), length(node))
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
    
    p.d[base.length + d, ] = vapply(1L:nIP, function(IP) {
      sum(names(text[[base.length + d]]) %in% which(currentC == IP))
    }, c(1)) / max(1, N.d)
    if (t.d >= 384) {
      history.t = History(edge, p.d, node, t.d + 10^(-10))
    }
    X = lapply(node, function(i) {
      	Netstats(history.t, node, i, netstat)
   		})
    XB = MultiplyXBList(X, b)     
    lambda = lambda_cpp(p.d[base.length + d,], XB)
    
    	for (i in node) {
    		iJi[i, -i] = r.gibbs.measure(1, lambda[i, -i], delta, support)
    	}
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
  if (base == TRUE & t.d > 384) {
    cutoff = which_int(384, vapply(1:length(edge), function(d) {edge[[d]][[3]]}, c(1))) - 1
    edge = edge[1:cutoff]
    text = text[1:cutoff]
  }
  return(list(edge = edge, text = text, base = base.length, b = b, d = delta, X = X))							
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
  phi = lapply(1L:K, function(k) {
    rdirichlet_cpp(1, betas * nvec)
  })
  L = 3
  P = 1 * ("intercept" %in% netstat) + L * (2 * ("dyadic" %in% netstat) + 4 * ("triadic" %in% netstat) + 2 *("degree" %in% netstat))
  
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
  	 p.d[d, ] = vapply(1L:nIP, function(IP) {
      sum(names(text[[d]]) %in% which(currentC == IP))
    }, c(1)) / max(1, length(text[[d]]))
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
            phi.k[w] = (word_type_topic_counts[w, topic.d[n]] + betas * nvec[w]) / (sum(word_type_topic_counts[,topic.d[n]]) + betas)
          } 
          text[[base.length + d]][n] = multinom_vec(1, phi.k)
          word_type_topic_counts[text[[base.length + d]][n], topic.d[n]] = word_type_topic_counts[text[[base.length + d]][n], topic.d[n]] + 1
        }
        names(text[[base.length + d]]) = topic.d
    }
    
    p.d[base.length + d, ] = vapply(1L:nIP, function(IP) {
      sum(names(text[[base.length + d]]) %in% which(currentC == IP))
    }, c(1)) / max(1, N.d)
    if (t.d >= 384) {
      history.t = History(edge, p.d, node, t.d + 10^(-10))
    }
    X = lapply(node, function(i) {
      	Netstats(history.t, node, i, netstat)
   		})
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
  if (base == TRUE & t.d > 384) {
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

#' @title GiR_QQ_Plots
#' @description Generate QQ-plots (quantile-quantile plots) for Gettig it Right test
#'
#' @param Forward_stats statistics obtained from forward sampling
#' @param Backward_stats statistics obtained from backward sampling
#'
#' @return QQ-plots for different GiR statistics of interest
#'
#' @export
GiR_QQ_Plots = function(Forward_stats, Backward_stats) {
  nms = colnames(Forward_stats)
  
  for (i in 1L:ncol(Forward_stats)) {
    if (nrow(Forward_stats) > 20000) {
      thinning = seq(from = floor(nrow(Forward_stats)/10), to = nrow(Forward_stats), length.out = 10000)
      Forward_test = Forward_stats[thinning, i]
      Backward_test = Backward_stats[thinning, i]
    } else {
      Forward_test = Forward_stats[, i]
      Backward_test = Backward_stats[, i]
    }
    
    all = c(Backward_stats[, i], Forward_stats[, i])
    ylims = c(min(all) - 0.1 * max(abs(all)), max(all) + 0.1 * max(abs(all)))
    xlims = ylims
    
    quantiles = 50
    if (grepl("B_", nms[i]) ) {
      quantiles = 1000
    }
    
    qqplot(x = quantile(Forward_stats[, i], seq(0, 1, length = quantiles)),
           y = quantile(Backward_stats[, i], seq(0, 1, length = quantiles)),
           ylim = ylims,
           xlim = xlims,
           ylab = "Backward",
           xlab = "Forward",
           col = "blue",
           pch = 19,
           main = nms[i],
           cex.lab = 0.2,
           cex.axis = 0.2,
           cex.main = 1)
    lines(x = xlims, y = ylims, col = "red", lwd = 3)
    text(paste("Backward Mean:", round(mean(Backward_stats[,i]), 4),
               "\nForward Mean:", round(mean(Forward_stats[,i]), 4),
               "\nt-test p-value:", round(t.test(Backward_test, Forward_test)$p.value, 4),
               "\nMann-Whitney p-value:", round(wilcox.test(Backward_test, Forward_test)$p.value,4)),
         x = xlims[2] - 0.35 * abs(xlims[2] - xlims[1]),
         y = ylims[1] + 0.15 * abs(ylims[2] - ylims[1]),
         cex = 0.3)
  }
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
  
  for (i in 1L:ncol(Forward_stats)) {
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
    
    if (nrow(Forward_stats) > 1000) {
    thinning2 = seq(from = floor(nrow(Forward_stats) / 10), to = nrow(Forward_stats), length.out = 1000)
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
#' @param generate_trace_plots Logical indicating whether to draw trace plots for each inference
#' @param seed an integer value which controls random number generation
#'
#' @return Forward and Backward samples and corresponding test results
#'
#' @export
GiR.Gibbs = function(Nsamp, nDocs, node, vocabulary, nIP, K, nwords, alpha, mvec, betas, nvec, 
               prior.b.mean, prior.b.var, prior.delta, sigma_Q, niters, netstat = c("dyadic"), 
               base.edge, base.text, generate_PP_plots = TRUE, generate_trace_plots = FALSE, seed = 1) {
  
  set.seed(seed)
  P = 1 * ("intercept" %in% netstat) + 3 * (2 * ("dyadic" %in% netstat) + 
  	  4 * ("triadic" %in% netstat) + 2 *("degree" %in% netstat))
  supportD = gibbs.measure.support(length(node) - 1)

  #Forward sampling
  Forward_stats = matrix(NA, nrow = Nsamp, ncol = P + (P - 1) / 2 + 4 + nIP + K + length(vocabulary))
  colnames(Forward_stats) = c(paste0("B_",1:P), paste0("Send_",1:3), "delta", "Mean_recipients", "Mean_timediff", 
  							  "Mean_TopicIP", paste0("Tokens_in_IP_", 1:nIP), paste0("Tokens_in_Topic", 1:K), paste0("Tokens_in_Word",1:length(vocabulary)))
  bmat1 = matrix(NA, nrow = Nsamp, ncol = nIP * P)
  cmat1 = matrix(NA, nrow = Nsamp, ncol = K)
  for (i in 1:Nsamp) { 
    if (i %% 5000 == 0) {cat("Forward sampling", i, "\n")}
    b = lapply(1:nIP, function(IP) {
      c(rmvnorm(1, prior.b.mean, prior.b.var))
    })
    delta = rnorm(1, prior.delta[1], sqrt(prior.delta[2]))
    currentC = sample(1L:nIP, K, replace = TRUE)
    Forward_sample = GenerateDocs.Gibbs(nDocs, node, vocabulary, nIP, K, nwords, alpha, mvec, betas, nvec, b,
    		 								delta, currentC, netstat, base.edge = base.edge, base.text = base.text,
    		 								forward = TRUE, base = FALSE, support = supportD) 
    Forward_stats[i, ] = GiR_stats(Forward_sample, K, currentC, vocabulary, forward = TRUE, backward = FALSE)
    bmat1[i, ] = unlist(b)
    cmat1[i, ] = currentC
    }
  
  #Backward sampling
  Backward_stats = matrix(NA, nrow = Nsamp, ncol = ncol(Forward_stats))
  Backward_sample = GenerateDocs.Gibbs(nDocs, node, vocabulary, nIP, K, nwords, alpha, mvec, betas, nvec, b, 
  									   delta, currentC, netstat, base.edge = base.edge, base.text = base.text,
  									   backward_init = TRUE, forward = FALSE, backward = FALSE, base = FALSE, support = supportD) 
  accept.rate = matrix(NA, nrow = Nsamp, ncol = 2)
  geweke.diag1 = matrix(NA, nrow = Nsamp, ncol = P)
  geweke.diag2 = matrix(NA, nrow = Nsamp, ncol = P)
  bmat2 = matrix(NA, nrow = Nsamp, ncol = nIP * P)		
  cmat2 = matrix(NA, nrow = Nsamp, ncol = K)
  
  for (i in 1:Nsamp) { 
    if (i %% 500 == 0) {cat("Backward sampling", i, "\n")}
    Inference_samp = IPTM_inference.Gibbs(Backward_sample$edge, node, Backward_sample$text, vocabulary, nIP, K,
    										  sigma_Q, alpha, mvec, betas, nvec, prior.b.mean, prior.b.var, prior.delta,
                               			  out = niters[1], n_B = niters[2], n_d = niters[3], burn = niters[4], 
                               			  thinning = niters[5], netstat, plot = generate_trace_plots)
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
    
    accept.rate[i, 1] = length(unique(Inference_samp$B[[1]][1,])) / ((niters[2] - niters[4]) / niters[5])
    accept.rate[i, 2] = length(unique(Inference_samp$D)) / niters[3]    
    geweke.diag1[i, ] = geweke.diag(t(Inference_samp$B[[1]]))[[1]]
    geweke.diag2[i, ] = geweke.diag(t(Inference_samp$B[[2]]))[[1]]
    bmat2[i, ] = unlist(b)
    cmat2[i, ] = currentC
  }
  				
  if (generate_PP_plots) {
    par(mfrow=c(5,5), oma = c(3,3,3,3), mar = c(2,1,1,1))
    GiR_PP_Plots(Forward_stats, Backward_stats)
  }			
  return(list(Forward = Forward_stats, Backward = Backward_stats, b1 = bmat1, b2= bmat2, c1 = cmat1, c2 = cmat2,
  			  accept.rate = accept.rate, geweke1 = geweke.diag1, geweke2 = geweke.diag2))
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
#' @param base.edge artificial collection of documents to be used as initial state of history
#' @param base.text artificial collection of documents to be used as initial state of history
#' @param generate_PP_plots Logical indicating whether to draw PP plots
#' @param generate_trace_plots Logical indicating whether to draw trace plots for each inference
#' @param seed an integer value which controls random number generation
#'
#' @return Forward and Backward samples and corresponding test results
#'
#' @export
Schein.Gibbs = function(Nsamp, nDocs, node, vocabulary, nIP, K, nwords, alpha, mvec, betas, nvec, 
               prior.b.mean, prior.b.var, prior.delta, sigma_Q, niters, netstat = c("dyadic"), 
               base.edge, base.text, generate_PP_plots = TRUE, generate_trace_plots = FALSE, seed = 1) {
  
  set.seed(seed)
  P = 1 * ("intercept" %in% netstat) + 3 * (2 * ("dyadic" %in% netstat) + 
  	  4 * ("triadic" %in% netstat) + 2 *("degree" %in% netstat))
  supportD = gibbs.measure.support(length(node) - 1)

  #Forward sampling
  Forward_stats = matrix(NA, nrow = Nsamp, ncol = P + (P - 1) / 2 + 4 + nIP + K + length(vocabulary))
  colnames(Forward_stats) = c(paste0("B_",1:P), paste0("Send_",1:3), "delta", "Mean_recipients", "Mean_timediff", 
  							  "Mean_TopicIP", paste0("Tokens_in_IP_", 1:nIP), paste0("Tokens_in_Topic", 1:K), paste0("Tokens_in_Word",1:length(vocabulary)))
  bmat1 = matrix(NA, nrow = Nsamp, ncol = nIP * P)
  cmat1 = matrix(NA, nrow = Nsamp, ncol = K)
  
  accept.rate = matrix(NA, nrow = Nsamp, ncol = 2)
  geweke.diag1 = matrix(NA, nrow = Nsamp, ncol = P)
  geweke.diag2 = matrix(NA, nrow = Nsamp, ncol = P)
  bmat2 = matrix(NA, nrow = Nsamp, ncol = nIP * P)		
  cmat2 = matrix(NA, nrow = Nsamp, ncol = K)
  Backward_stats = matrix(NA, nrow = Nsamp, ncol = ncol(Forward_stats))

  for (i in 1:Nsamp) { 
    b = lapply(1:nIP, function(IP) {
      c(rmvnorm(1, prior.b.mean, prior.b.var))
    })
    delta = rnorm(1, prior.delta[1], sqrt(prior.delta[2]))
    currentC = sample(1L:nIP, K, replace = TRUE)
    Forward_sample = GenerateDocs.Schein(nDocs, node, vocabulary, nIP, K, nwords, alpha, mvec, betas, nvec, b,
    		 								delta, currentC, netstat, base.edge = base.edge, base.text = base.text,
    		 								forward = TRUE, base = FALSE, support = supportD) 
    Forward_stats[i, ] = GiR_stats(Forward_sample, K, currentC, vocabulary, forward = FALSE, backward = TRUE)
    bmat1[i, ] = unlist(b)
    cmat1[i, ] = currentC

    if (i %% 100 == 0) {cat("Sampling", i, "\n")}
    initial = list(C = currentC, D = delta, B = b, Z = lapply(Forward_sample$text, function(d){as.numeric(names(d))}),
    				iJi = Forward_sample$iJi)
    Inference_samp = IPTM_inference.Schein(Forward_sample$edge, node, Forward_sample$text, vocabulary, nIP, K,
    										  sigma_Q, alpha, mvec, betas, nvec, prior.b.mean, prior.b.var, prior.delta,
                               			  out = niters[1], n_B = niters[2], n_d = niters[3], burn = niters[4], 
                               			  thinning = niters[5], netstat, plot = generate_trace_plots, initial = initial)
    b = lapply(1:nIP, function(IP) {
        Inference_samp$B[[IP]][,ncol(Inference_samp$B[[IP]])]
    })
    delta = Inference_samp$D[length(Inference_samp$D)]
    currentC = Inference_samp$C
    topic_token_assignments = Inference_samp$Z
    for (d in 1:length(topic_token_assignments)) {
      names(topic_token_assignments[[d]]) = Forward_sample$text[[d + length(base.text)]]
    }
    
    Backward_sample = GenerateDocs.Gibbs(nDocs, node, vocabulary, nIP, K, nwords, alpha, mvec, betas, nvec, b, 
    										 delta, currentC, netstat, base.edge = base.edge, base.text = base.text,
    										 topic_token_assignments = topic_token_assignments, 
                             			 forward = FALSE, backward = TRUE, base = FALSE, support = supportD)
    Backward_stats[i, ] = GiR_stats(Backward_sample, K, currentC, vocabulary, forward = FALSE, backward = TRUE)
    
    accept.rate[i, 1] = length(unique(Inference_samp$B[[1]][1,])) / ((niters[2] - niters[4]) / niters[5])
    accept.rate[i, 2] = length(unique(Inference_samp$D)) / niters[3]    
    geweke.diag1[i, ] = geweke.diag(t(Inference_samp$B[[1]]))[[1]]
    geweke.diag2[i, ] = geweke.diag(t(Inference_samp$B[[2]]))[[1]]
    bmat2[i, ] = unlist(b)
    cmat2[i, ] = currentC
  }
  				
  if (generate_PP_plots) {
    par(mfrow=c(5,5), oma = c(2,2,2,2), mar = c(1,1,1,1))
    GiR_PP_Plots(Forward_stats, Backward_stats)
  }			
  return(list(Forward = Forward_stats, Backward = Backward_stats, b1 = bmat1, b2= bmat2, c1 = cmat1, c2 = cmat2,
  			  accept.rate = accept.rate, geweke1 = geweke.diag1, geweke2 = geweke.diag2))
}                   




#' @title Comparison.Gibbs
#' @description Comparison test for IPTM
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
#' @param generate_trace_plots Logical indicating whether to draw trace plots for each inference
#' @param seed an integer value which controls random number generation
#'
#' @return Forward and Backward samples and corresponding test results
#'
#' @export
Comparison.Gibbs = function(Nsamp, nDocs, node, vocabulary, nIP, K, nwords, alpha, mvec, betas, nvec, 
               prior.b.mean, prior.b.var, prior.delta, sigma_Q, niters, netstat = c("dyadic"), 
               base.edge, base.text, generate_PP_plots = TRUE, generate_trace_plots = FALSE, seed = 1) {
               	
  set.seed(seed)
  P = 1 * ("intercept" %in% netstat) + 3 * (2 * ("dyadic" %in% netstat) + 
  	  4 * ("triadic" %in% netstat) + 2 *("degree" %in% netstat))
  supportD = gibbs.measure.support(length(node) - 1)

  #Forward sampling
  Forward = list()
  Backward = list()
  for (i in 1:Nsamp) { 
    b = lapply(1:nIP, function(IP) {
      c(rmvnorm(1, prior.b.mean, prior.b.var))
    })
    delta = rnorm(1, prior.delta[1], sqrt(prior.delta[2]))
    currentC = sample(1L:nIP, K, replace = TRUE)
    Forward_sample = GenerateDocs.Schein(nDocs, node, vocabulary, nIP, K, nwords, alpha, mvec, betas, nvec, b,
    		 								delta, currentC, netstat, base.edge = base.edge, base.text = base.text,
    		 								forward = TRUE, base = FALSE, support = supportD) 
    Forward[[i]] = Forward_sample$iJi[(length(base.edge)+1):length(Forward_sample$iJi)]
    initial = list(C = currentC, D = delta, B = b, Z = lapply(Forward_sample$text, function(d){as.numeric(names(d))}), 
    					iJi = Forward_sample$iJi)
    Backward_sample = IPTM_inference.Schein2(Forward_sample$edge, node, Forward_sample$text, vocabulary, nIP, K,
    										  sigma_Q, alpha, mvec, betas, nvec, prior.b.mean, prior.b.var, prior.delta,
                               			  out = niters[1], n_B = niters[2], n_d = niters[3], burn = niters[4], 
                               			  thinning = niters[5], netstat, plot = generate_trace_plots, initial = initial)
    Backward[[i]] = Backward_sample$iJi[(length(base.edge)+1):length(Backward_sample$iJi)]
     }
return(list(Forward = Forward, Backward = Backward))
}                   


  			          	
               	