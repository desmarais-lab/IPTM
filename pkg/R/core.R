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
#' @importFrom coda mcmc
NULL

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
#'
#' @return The optimized value of the vector (= alpha * mvec)
#'
#' @export
AlphamvecOpt =  function(K, currentZ, alpha, mvec) {
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
	
	iter = 1
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
	while ((abs(alpha - sum(current.vec)) > 0.001) | (iter == 1)) {
	alpha = sum(current.vec)
	S = UpdateDenom(alpha, n.word.table)		
	s = UpdateNum(current.vec, nK.word.table)	
	current.vec = current.vec * s / S
	iter = iter + 1
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

#' @title IPTM_inference
#' @description Iterate Markov Chain Monte Carlo (MCMC) algorithm to sequentially update the assignments of Z, C and B 
#'
#' @param edge list of document information with 3 elements (element 1 sender, element 2 receiver, element 3 time in unix.time format)
#' @param node nodelist containing the ID of nodes (ID starting from 1)
#' @param textlist list of text (length=number of documents in total) containing the words in each document
#' @param vocabulary all vocabularies used over the corpus
#' @param nIP total number of interaction patterns specified by the user
#' @param K total number of topics specified by the user
#' @param sigma_Q proposal distribution variance parameter for beta and eta
#' @param alpha Dirichlet concentration prior for document-topic distribution
#' @param mvec Dirichlet base prior for document-topic distribution
#' @param betas Dirichlet concentration prior for topic-word distribution
#' @param nvec Dirichlet base prior for topic-word distribution
#' @param prior.b.mean mean vector of b in multivariate normal distribution
#' @param prior.b.var covairance matrix of b in multivariate normal distribution
#' @param prior.delta parameter of delta in gamma prior
#' @param out size of outer iterations 
#' @param n1 size of first inner iteration for updates of interaction patterns Z
#' @param n2 size of second inner iteration for updates of topics C
#' @param n3 size of third inner iteration for updates of beta B
#' @param burn iterations to be discarded at the beginning of beta chain
#' @param thin the thinning interval of beta chain
#' @param netstat which type of network statistics to use ("intercept", dyadic", "triadic", "degree")
#' @param plot to plot the convergence diagnostics or not (TRUE/FALSE)
#' @param optimize to optimize alpha (Dirichlet concentration prior for document-topic distribution) or not (TRUE/FALSE)
#'
#' @return MCMC output containing IP assignment, topic assignment, and (beta, mu, delta) chain 
#'
#' @export
IPTM_inference = function(edge, node, textlist, vocabulary, nIP, K, sigma_Q, alpha, mvec, betas, nvec, prior.b.mean, prior.b.var, prior.delta, 
					out, n1, n2, n3, burn, thin, netstat, plot = FALSE, optimize = FALSE) {
   
  # trim the edge so that we only model edges after 384 hours
	timestamps = vapply(edge, function(d) {
  			  d[[3]]
 			  }, c(1))

    edge2 = which_int(384, timestamps) : length(edge)
    timeinc = c(timestamps[1], timestamps[-1] - timestamps[-length(timestamps)])
   
  # initialize alpha, mvec, delta, nvec, eta, lvec, and gammas
 	W = length(vocabulary)
  	phi = lapply(1L:K, function(k) {
		rdirichlet_cpp(1, betas * nvec)
	})
	delta = rgamma(1, prior.delta[1], prior.delta[2])
	
	# initialize C, theta and Z
	currentC = sample(1:nIP, K, replace = TRUE)
	theta = rdirichlet_cpp(length(edge), alpha * mvec)
	currentZ = lapply(seq(along = edge), function(d) {
    if (length(textlist[[d]]) > 0) {
    		multinom_vec(length(textlist[[d]]), theta[d, ])
    	} else {
    		multinom_vec(1, theta[d, ])
    	}
    })
	p.d = t(vapply(seq(along = edge), function(d) {
	  vapply(1L:nIP, function(IP) {
	    sum(currentZ[[d]] %in% which(currentC == IP))
	  }, c(1)) / length(currentZ[[d]])
	}, rep(1, nIP)))

    # initialize beta
    sigma = 1
    L = 3
    P = 1 * ("intercept" %in% netstat) + L * (2 * ("dyadic" %in% netstat) + 4 * ("triadic" %in% netstat) + 2 *("degree" %in% netstat))
    bmat = list()
	for (IP in 1L:nIP) {
		bmat[[IP]] = matrix(0, nrow = P, ncol = (n3 - burn) / thin)
		bmat[[IP]][, 1:(n3 - burn) / thin] = c(rmvnorm(1, prior.b.mean, prior.b.var))
  	}
  	deltamat = rep(delta, (n3 - burn) / thin)

    # to check the convergence  
    if (plot) {
     	logWZ.mat = c()							  
     	alpha.mat = c()
     	entropy.mat = c()
    }
  	
    iJi = lapply(seq(along = edge), function(d) {
  	      matrix(0, nrow = length(node), ncol = length(node))
  	      })  
    
    #start outer iteration
    for (o in 1L:out) {
      
      if (optimize) {
      #update the hyperparameter alpha and mvec
      vec = AlphamvecOpt(K, currentZ[edge2], alpha, mvec)
      alpha = sum(vec)
      mvec = vec / alpha
      }
      Beta.old = lapply(bmat, function(b) {
      			rowMeans(b)
         		})
      delta = mean(deltamat)
     # Data augmentation
      lambda = list()
      LambdaiJi = list()
      nonemptyiJi = list()
	  observediJi = list()
      for (d in edge2) {
   	 	history.t = History(edge, p.d, node, as.numeric(edge[[d]][3]))
   	 	X = lapply(node, function(i) {
  	        Netstats(history.t, node, i, netstat)
            })
   	 	XB = MultiplyXBList(X, Beta.old)     
		lambda[[d]] = lambda_cpp(p.d[d,], XB)
		#calculate the resampling probability
		probij = DataAug_cpp(iJi[[d]], lambda[[d]], delta, timeinc[d])
		iJi[[d]] = rbinom_mat(probij)
		# diag(iJi[[d]])[which(rowSums(iJi[[d]]) == 0)] = 1
		# while (sum(iJi[[d]]) == 0) {
		#	iJi[[d]] = rbinom_mat(probij)
		# }
		iJi[[d]][as.numeric(edge[[d]][1]),] = tabulateC(as.numeric(unlist(edge[[d]][2])), length(node))
		LambdaiJi[[d]] = lambdaiJi(p.d[d,], XB, iJi[[d]])
		nonemptyiJi[[d]] = LambdaiJi[[d]][!is.na(LambdaiJi[[d]])]
		observediJi[[d]] = LambdaiJi[[d]][as.numeric(edge[[d]][1])]
		}	 
    
    	  textlist.raw = unlist(textlist[edge2])
    	  table.W = lapply(1L:K, function(k) {
      			tabulateC(textlist.raw[which(unlist(currentZ[edge2]) == k)], W)
      			})
    
    	  for (i1 in 1L:n1) {
      	for (d in edge2) {
      		textlist.d = textlist[[d]]
        		if (length(textlist.d) > 0) {
        			topicpart.d = TopicInEqZ(K, currentZ[[d]], alpha, mvec, d)
       			wordpart.d = WordInEqZ(K, textlist.d, table.W, betas, nvec)
        		} else {
        			topicpart.d = 0
        			wordpart.d = matrix(0, nrow = length(currentZ[[d]]), ncol = K)
        		}
        		edgepart.d = EdgeInEqZ(iJi[[d]], lambda[[d]], delta)
        		timepart.d = TimeInEqZ(nonemptyiJi[[d]], timeinc[d])
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
		    			nonemptyiJi[[d]] = LambdaiJi[[d]][!is.na(LambdaiJi[[d]])]
         			observediJi[[d]] = LambdaiJi[[d]][as.numeric(edge[[d]][1])]
          		}
        		}
       	}		
	}

    for (i2 in 1L:n2) { 
      # C update given Z and B - within each document d
      for (k in unique(unlist(currentZ[edge2]))) { 
        document.k = which(vapply(currentZ, function(d){k %in% d}, c(1)) == 1)
        document.k = document.k[document.k %in% edge2]
        const.C = rep(NA, nIP)
        for (IP in 1L:nIP) {
        		if (!currentC[k] == IP){
          		currentC[k] = IP
          	for (d in document.k) {
            		p.d[d, ] =  vapply(1L:nIP, function(IP) {
              				sum(currentZ[[d]] %in% which(currentC == IP))
            					}, c(1)) / length(currentZ[[d]])
          	}
          	for (d in edge2) {
           		history.t = History(edge, p.d, node, as.numeric(edge[[d]][3]))
    	   	   		X = lapply(node, function(i) {
               		Netstats(history.t, node, i, netstat)
               		})
    	       		XB = MultiplyXBList(X, Beta.old)    
           		lambda[[d]] = lambda_cpp(p.d[d,], XB)
		       	LambdaiJi[[d]] = lambdaiJi(p.d[d,], XB, iJi[[d]])
		       	nonemptyiJi[[d]] = LambdaiJi[[d]][!is.na(LambdaiJi[[d]])]
           		observediJi[[d]] = LambdaiJi[[d]][as.numeric(edge[[d]][1])]
          	}
          }
          const.C[IP] = sum(vapply(document.k, function(d) {
          				EdgeInEqZ(iJi[[d]], lambda[[d]], delta) + 
          				TimeInEqZ(nonemptyiJi[[d]], timeinc[d]) + 
    						ObservedInEqZ(observediJi[[d]]) 
          				}, c(1))) / length(document.k)
      	}
        const.C = const.C - max(const.C)
        currentC[k] = multinom_vec(1, exp(const.C))
     }
	}
    
    if (plot) {
      entropy.mat = c(entropy.mat, entropy.empirical(currentC))
      alpha.mat = rbind(alpha.mat, alpha)
      logWZ.mat = c(logWZ.mat, logWZ(K, currentZ[edge2], textlist[edge2], table.W, alpha, mvec, betas, nvec))
    }
      
    for (d in edge2) {
    		p.d[d, ] = vapply(1L:nIP, function(IP) {
              	   sum(currentZ[[d]] %in% which(currentC == IP))
       			   }, c(1)) / length(currentZ[[d]])
     	history.t = History(edge, p.d, node, as.numeric(edge[[d]][3]))
    	    X = lapply(node, function(i) {
            Netstats(history.t, node, i, netstat)
       		})
    	    XB = MultiplyXBList(X, Beta.old)   
    	    lambda[[d]] = lambda_cpp(p.d[d,], XB)
		LambdaiJi[[d]] = lambdaiJi(p.d[d,], XB, iJi[[d]])
		nonemptyiJi[[d]] = LambdaiJi[[d]][!is.na(LambdaiJi[[d]])]
        observediJi[[d]] = LambdaiJi[[d]][as.numeric(edge[[d]][1])]
	}
	 
	# Beta and delta update
	prior.old1 = sum(vapply(1L:nIP, function(IP) {
		  		 dmvnorm(Beta.old[[IP]], prior.b.mean, prior.b.var, log = TRUE)
		  		 }, c(1))) 
	post.old1 = sum(vapply(edge2, function(d) {
	     	    EdgeInEqZ(iJi[[d]], lambda[[d]], delta) + TimeInEqZ(nonemptyiJi[[d]], timeinc[d]) +
    			 	ObservedInEqZ(observediJi[[d]])
    			 	}, c(1))) / length(edge2)
  
	options(warn = -1)
	proposal.var = list()
	for (IP in 1:nIP) {
		proposal.var[[IP]] = diag(P)
	}
    #for (IP in 1L:nIP) {
    # 	proposal.var[[IP]] = cor(t(bmat[[IP]]))
    #    proposal.var[[IP]][is.na(proposal.var[[IP]])] = 0
    #    if (sum(eigen(proposal.var[[IP]])$values < 0 ) > 0) {
    #     	proposal.var[[IP]] = diag(P)
    #    }
    #}
    proposal.var = lapply(1:nIP, function(IP){diag(P)})
    options(warn = 0)
    for (i3 in 1L:n3) {
    		Beta.new = lapply(1L:nIP, function(IP) {
          		   rmvnorm(1, Beta.old[[IP]], sigma_Q[1] * proposal.var[[IP]])
         		   }) 
        for (d in edge2) {
           XB = MultiplyXBList(X, Beta.new)
           lambda[[d]] = lambda_cpp(p.d[d,], XB)    
           LambdaiJi[[d]] = lambdaiJi(p.d[d,], XB, iJi[[d]])
		   nonemptyiJi[[d]] = LambdaiJi[[d]][!is.na(LambdaiJi[[d]])]
		   observediJi[[d]] = LambdaiJi[[d]][as.numeric(edge[[d]][1])]
        }
        prior.new1 = sum(vapply(1L:nIP, function(IP) {
        				 dmvnorm(Beta.new[[IP]], prior.b.mean, prior.b.var, log = TRUE)
        				 }, c(1))) 
        post.new1 = sum(vapply(edge2, function(d) {
    			   		EdgeInEqZ(iJi[[d]], lambda[[d]], delta) + TimeInEqZ(nonemptyiJi[[d]], timeinc[d]) + 
    			    	ObservedInEqZ(observediJi[[d]])
    			    	}, c(1))) / length(edge2)
    		loglike.diff = prior.new1 + post.new1 - prior.old1 - post.old1
      
        if (log(runif(1, 0, 1)) < loglike.diff) {
        		for (IP in 1L:nIP) {
         		Beta.old[[IP]]  = Beta.new[[IP]]
         	}
         	prior.old1 = prior.new1
         	post.old1 = post.new1
        }
         if (i3 > burn & i3 %% (thin) == 0) {
        		for (IP in 1L:nIP) {
           		bmat[[IP]][ , (i3 - burn) / thin] = Beta.old[[IP]]
           	}
    		}
    }
    Beta.old = lapply(bmat, function(b) {
      			rowMeans(b)
         		})
    for (d in edge2) {
    	 XB = MultiplyXBList(X, Beta.old)
         lambda[[d]] = lambda_cpp(p.d[d,], XB)  
    }
    prior.old2 = dgamma(delta, prior.delta[1], prior.delta[2], log = TRUE)
    post.old2 = sum(vapply(edge2, function(d) {
    				EdgeInEqZ(iJi[[d]], lambda[[d]], delta)
    				}, c(1))) / length(edge2)
 	for (i4 in 1L:n3) {
        delta.new = exp(rnorm(1, log(delta), sigma_Q[2]))
        prior.new2 = dgamma(delta.new, prior.delta[1], prior.delta[2], log = TRUE)
        post.new2 = sum(vapply(edge2, function(d) {
        				EdgeInEqZ(iJi[[d]], lambda[[d]], delta.new)
        				}, c(1))) / length(edge2)
        loglike.diff2 = prior.new2 + post.new2 - prior.old2 - post.old2 + log(delta.new) - log(delta) 
        if (log(runif(1, 0, 1)) < loglike.diff2) {
        		delta = delta.new
         	prior.old2 = prior.new2
         	post.old2 = post.new2
        } 
        if (i4 > burn & i4 %% (thin) == 0) {
           	deltamat[(i4 - burn) / thin] = delta
        }
    	}
 }
         
  if (plot) {
    par(mfrow = c(2, 2))
    matplot(alpha.mat, lty = 1, type = "l", col = 1L:nIP, 
  	        xlab = "(Outer) Iterations", ylab = "alpha")
  	abline(h = mean(alpha.mat), lty = 1, col = 1L:nIP)
	title("Convergence of Optimized alpha")
	plot(logWZ.mat, type = "l", 
	       xlab = "(Outer) Iterations", ylab = "logWZ")
	abline(h = mean(logWZ.mat), lty = 1)
	title("Convergence of logWZ")
	matplot(bmat[[1]][1,], lty = 1, col = 1L:P, type = "l", 
	          main = "Traceplot of beta", xlab = "(Inner) Iterations", ylab = "")
	abline(h = mean(bmat[[1]][1,]), lty = 1, col = 1L)
	plot(deltamat, type = "l", 
	xlab = "(Outer) Iterations", ylab = "Entropy of IP")
	abline(h = mean(deltamat), lty = 1)
	title("Convergence of delta")
  }
     
  chain.final = list(C = currentC, Z = lapply(edge2, function(d) {currentZ[[d]]}), B = bmat, D = deltamat)
  return(chain.final)
}

#' @title TableBetaIP
#' @description Generate a table summary of the MCMC chain of network statistics coefficients (beta) for each interaction pattern
#'
#' @param MCMCchain a chain obtained using MCMC function
#'
#' @return List of table containing the posterior summary of beta for each interaction pattern
#'
#' @export
TableBetaIP = function(MCMCchain) {
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

#' @title PlotBetaIP
#' @description Draw a boxplot of the MCMC chain of network statistics (beta) for each interaction pattern
#'
#' @param Bchain a chain of B obtained using MCMC function
#'
#' @return Joint boxplot of the posterior summary of beta (should click for the locator)
#'
#' @export
PlotBetaIP = function(Bchain) {
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
#' @param forward Logigal indicating whether we are generating forward samples
#' @param backward_init Logigal indicating whether we are generating initial backward samples
#' @param backward Logigal indicating whether we are generating backward samples
#' @param base Logical indicating whether or not we are generating base edges (< 384)
#'
#' @return generated edge and text, parameter b used to generate those, and base (if base == TRUE)
#'
#' @export
GenerateDocs = function(nDocs, node, vocabulary, nIP, K, nwords, alpha, mvec, betas, nvec,
                        b, delta, currentC, netstat, base.edge, base.text, topic_token_assignments = NULL, forward = FALSE, 
                        backward_init = FALSE, backward = FALSE, base = FALSE) {
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
  options(warn = -1)
  for (d in 1:nDocs) {
    N.d = nwords
    text[[base.length + d]] = rep(NA, N.d)
    
    if (!backward) {
      theta.d = rdirichlet_cpp(1, alpha * mvec)	
      topic.d = multinom_vec(max(1, N.d), theta.d)
      
      if (N.d > 0) {
        for (n in 1:N.d){
          text[[base.length + d]][n] = multinom_vec(1, phi[[topic.d[n]]])
        }
        names(text[[base.length + d]]) = topic.d
      }
    } else {
      phi.k = rep(NA, K)
      topic.d = topic_token_assignments[[d]]
      if (N.d > 0) {
        for (n in 1:N.d){
          for (w in 1:W) {
            phi.k[w] = (word_type_topic_counts[w, topic.d[n]] + betas * nvec[w]) / (sum(word_type_topic_counts[,topic.d[n]]) + betas)
          } 
          text[[base.length + d]][n] = multinom_vec(1, phi.k)
          word_type_topic_counts[text[[base.length + d]][n], topic.d[n]] = word_type_topic_counts[text[[base.length + d]][n], topic.d[n]] + 1
        }
        names(text[[base.length + d]]) = topic.d
      }
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
    
    iJi = rbinom_mat((delta * lambda) / (delta * lambda + 1))
    #diag(iJi)[which(rowSums(iJi) == 0)] = 1
    while (sum(iJi) == 0) {
      iJi = rbinom_mat((delta * lambda) / (delta * lambda + 1))
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
  options(warn = 0)
  if (forward) {
    edge = edge[-(1:base.length)]
    text = text[-(1:base.length)]
  }
  if (base == TRUE & t.d > 384) {
    cutoff = which_int(384, vapply(1:length(edge), function(d) {edge[[d]][[3]]}, c(1))) - 1
    edge = edge[1:cutoff]
    text = text[1:cutoff]
  }
  return(list(edge = edge, text = text, b = b, base = base.length))							
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
  
  GiR_stats[1:P] = Reduce('+', GiR_sample$b) / nIP
  GiR_stats[P + 1] = mean(vapply(1:nDocs, function(d) {
    length(edge[[d]][[2]])
  }, c(1)))
  GiR_stats[P + 2] = mean(c(edge[[1]][[3]] - 384, vapply(2:nDocs, function(d) {
    edge[[d]][[3]] - edge[[d-1]][[3]]
  }, c(1)))) 			
  GiR_stats[P + 3] = mean(currentC)
  Tokens_in_Topic = tabulate(vapply(1:nDocs, function(d){
    as.numeric(names(text[[d]]))
  }, rep(0, nwords)), K)
  GiR_stats[(P + 4):(P + 3 + nIP)] = vapply(1:nIP, function(IP) {
    Tokens_in_Topic %*% (currentC == IP)
  }, c(1))
  GiR_stats[(P + 4 + nIP):(P + 3 + nIP + K)] = Tokens_in_Topic
  Tokens_in_Word = tabulate(unlist(text), W)
  GiR_stats[(P + 4 + nIP + K):(P + 3 + nIP + K + W)] = Tokens_in_Word
  
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
      thin = seq(from = floor(nrow(Forward_stats)/10), to = nrow(Forward_stats), length.out = 10000)
      Forward_test = Forward_stats[thin, i]
      Backward_test = Backward_stats[thin, i]
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
           cex.lab = 1,
           cex.axis = 1,
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
    if (nrow(Forward_stats) > 20000) {
      thin = seq(from = floor(nrow(Forward_stats)/10), to = nrow(Forward_stats), length.out = 10000)
      Forward_test = Forward_stats[thin, i]
      Backward_test = Backward_stats[thin, i]
    } else {
      Forward_test = Forward_stats[, i]
      Backward_test = Backward_stats[, i]
    }
    
    all = c(Backward_test, Forward_test)
    
    quantiles = 50
    if (grepl("B_", nms[i]) ) {
      quantiles = 1000
    }
    uniqueValues = sort(unique(all))
    qx1 = numeric(length(uniqueValues))
  	qx2 = numeric(length(uniqueValues))
  		
  	for (j in 1:length(uniqueValues)) {
  		qx1[j] = mean(Forward_test <= uniqueValues[j])
  		qx2[j] = mean(Backward_test <= uniqueValues[j])
  	}
    
    qqplot(x = qx1,
           y = qx2,
           ylim = c(0, 1),
           xlim = c(0, 1),
           ylab = "Backward",
           xlab = "Forward",
           col = "blue",
           pch = 19,
           main = nms[i],
           cex.lab = 1,
           cex.axis = 1,
           cex.main = 1)
    abline(0, 1, lty = 1, col = "red", lwd = 2)
    text(paste("Backward Mean:", round(mean(Backward_stats[,i]), 4),
               "\nForward Mean:", round(mean(Forward_stats[,i]), 4),
               "\nt-test p-value:", round(t.test(Backward_test, Forward_test)$p.value, 4),
               "\nMann-Whitney p-value:", round(wilcox.test(Backward_test, Forward_test)$p.value,4)),
         x = 0.65,
         y = 0.15,
         cex = 0.4)
  }
}      

#' @title GiR
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
#' @param prior.delta parameter of delta in gamma prior
#' @param sigma_Q proposal distribution variance parameter for beta and eta
#' @param niters (out, n1, n2, n3, burn, thin) in the inference
#' @param netstat which type of network statistics to use ("dyadic", "triadic", "degree")
#' @param generate_PP_plots Logical indicating whether to draw PP plots
#' @param seed an integer value which controls random number generation
#'
#' @return Forward and Backward samples and corresponding test results
#'
#' @export
GiR = function(Nsamp, nDocs = 10, node = 1:4, vocabulary =  c("hi", "hello","bye", "mine", "what"), 
               nIP = 2, K = 4, nwords = 4, alpha = 2, mvec = rep(1/4, 4), betas = 2, nvec = rep(1/5, 5), 
               prior.b.mean = c(-2, rep(0, 6)), prior.b.var = diag(7), prior.delta = c(4, 8), sigma_Q = c(0.25, 2), 
               niters = c(1, 1, 1, 50, 0, 1), netstat = c("intercept", "dyadic"), generate_PP_plots = TRUE, seed = 1) {
  
  set.seed(seed)
  P = 1 * ("intercept" %in% netstat) + 3 * (2 * ("dyadic" %in% netstat) + 4 * ("triadic" %in% netstat) + 2 *("degree" %in% netstat))
  b = lapply(1L:nIP, function(IP) {
    c(rmvnorm(1, prior.b.mean, prior.b.var))
  })
  delta = rgamma(1, prior.delta[1], prior.delta[2])
  currentC = sample(1L:nIP, K, replace = TRUE)	 
  base.data = GenerateDocs(nDocs, node, vocabulary, nIP, K, nwords, alpha, mvec, betas, nvec, b, delta, currentC, netstat, base.edge = list(),  base.text = list(), base = TRUE)
  base.edge = base.data$edge	   
  base.text = base.data$text
  
  #Forward sampling
  Forward_stats = matrix(NA, nrow = Nsamp, ncol = P + 3 + nIP + K + length(vocabulary))
  colnames(Forward_stats) = c(paste0("B_",1:length(b[[1]])), "Mean_recipients", "Mean_timediff", "Mean_TopicIP", 
                              paste0("Tokens_in_IP_", 1:nIP), paste0("Tokens_in_Topic", 1:K), 
                              paste0("Tokens_in_Word", 1:length(vocabulary)))
  deltamat1 = rep(NA, Nsamp)
  bmat1 = matrix(NA, nrow = Nsamp, ncol = nIP * P)
  entropyC1 = matrix(NA, Nsamp, 2)
  entropyZ1 = rep(NA, Nsamp)
  Zstat1 = rep(NA, Nsamp)
  for (i in 1:Nsamp) { 
    if (i %% 5000 == 0) {cat("Forward sampling", i, "\n")}
    b = lapply(1:nIP, function(IP) {
      c(rmvnorm(1, prior.b.mean, prior.b.var))
    })
    delta = rgamma(1, prior.delta[1], prior.delta[2])
    currentC = sample(1L:nIP, K, replace = TRUE)
    Forward_sample = GenerateDocs(nDocs, node, vocabulary, nIP, K, nwords, alpha, mvec, betas, nvec, b, delta, currentC, netstat, 
                                  base.edge = base.edge, base.text = base.text, forward = TRUE) 
    Forward_stats[i, ] = GiR_stats(Forward_sample, K, currentC, vocabulary, forward = TRUE, backward = FALSE)
    
    topic_token_counts = tabulate(unlist(lapply(Forward_sample$text, function(d){as.numeric(names(d))})), K)
    bmat1[i, ] = unlist(b)
    deltamat1[i] = delta
    entropyC1[i,] = c(entropy.empirical(currentC), entropy.empirical(topic_token_counts))
    entropyZ1[i] = entropy.empirical(tabulate(unlist(Forward_sample$text), length(vocabulary)))
    Zstat1[i] = mean(sapply(Forward_sample$text, function(d){entropy.empirical(as.numeric(names(d)))}))
  }
  
  #Backward sampling
  Backward_stats = matrix(NA, nrow = Nsamp, ncol = ncol(Forward_stats))
  Backward_sample = GenerateDocs(nDocs, node, vocabulary, nIP, K, nwords, alpha, mvec, betas, nvec, b, delta, currentC, netstat, 
                                 base.edge = base.edge, base.text = base.text, backward_init = TRUE, forward = FALSE, backward = FALSE) 
  deltamat2 = rep(NA, Nsamp)
  bmat2 = matrix(NA, nrow = Nsamp, ncol = nIP * P)		
  entropyC2 = matrix(NA, Nsamp, 2)
  entropyZ2 = rep(NA, Nsamp)	
  Zstat2 = rep(NA, Nsamp)	   
  accept.rate = matrix(NA, nrow = Nsamp, ncol = 2)
  for (i in 1:Nsamp) { 
    if (i %% 500 == 0) {cat("Backward sampling", i, "\n")}
    Inference_samp = IPTM_inference(Backward_sample$edge, node, Backward_sample$text, vocabulary, nIP, K, sigma_Q, 
                               alpha, mvec, betas, nvec, prior.b.mean, prior.b.var, prior.delta,
                               out = niters[1], n1 = niters[2], n2 = niters[3], n3 = niters[4], burn = niters[5], thin = niters[6], 
                               netstat)
    b = lapply(1:nIP, function(IP) {
      rowMeans(Inference_samp$B[[IP]])
    })
    delta = mean(Inference_samp$D)
    currentC = Inference_samp$C
    topic_token_assignments = Inference_samp$Z
    for (d in 1:length(topic_token_assignments)) {
      names(topic_token_assignments[[d]]) = Backward_sample$text[[d + length(base.text)]]
    }
    topic_token_counts = tabulate(unlist(Inference_samp$Z), K)
    
    Backward_sample = GenerateDocs(nDocs, node, vocabulary, nIP, K, nwords, alpha, mvec, betas, nvec, b, delta, currentC, netstat, 
                                   base.edge = base.edge, base.text = base.text, topic_token_assignments = topic_token_assignments, 
                                   forward = FALSE, backward = TRUE)
    Backward_stats[i, ] = GiR_stats(Backward_sample, K, currentC, vocabulary, forward = FALSE, backward = TRUE)
    
    bmat2[i, ] = unlist(b)
    deltamat2[i] = delta
    entropyC2[i,] = c(entropy.empirical(currentC), entropy.empirical(topic_token_counts))
    entropyZ2[i] = entropy.empirical(tabulate(unlist(Backward_sample$text[-1:-length(base.text)]), length(vocabulary)))
    Zstat2[i] = mean(sapply(topic_token_assignments, function(d){entropy.empirical(d)}))
    accept.rate[i, 1] = length(unique(Inference_samp$B[[1]][1,])) / niters[4]
    accept.rate[i, 2] = length(unique(Inference_samp$D)) / niters[4]
  }
  
  tstats = rep(0, ncol(Forward_stats))
  wstats = rep(0, ncol(Forward_stats))
  for (j in 1:ncol(Forward_stats)) {
    thin = seq(from = floor(Nsamp / 5), to = Nsamp, length.out = 400)
    Forward_test = Forward_stats[thin, j]
    Backward_test = Backward_stats[thin, j]
    tstats[j] = t.test(Backward_test, Forward_test)$p.value
    wstats[j] = wilcox.test(Backward_test, Forward_test)$p.value
  }
  names(tstats) = names(wstats) = colnames(Forward_stats)						
  if (generate_PP_plots) {
    par(mfrow=c(5,5), oma = c(3,3,3,3), mar = c(2,1,1,1))
    GiR_PP_Plots(Forward_stats, Backward_stats)
  }			
  return(list(Forward = Forward_stats, Backward = Backward_stats, tstats = tstats, wstats = wstats, delta = cbind(deltamat1, deltamat2),
              b1 = bmat1, b2= bmat2, entropyC1 = entropyC1, entropyC2 = entropyC2, entropyZ1 = entropyZ1, entropyZ2 = entropyZ2, Zstat1 = Zstat1, Zstat2 = Zstat2, accept.rate = accept.rate))
}                         	
