#' @useDynLib IPTM
#' @import stats
#' @import grDevices
#' @import graphics
#' @importFrom Rcpp sourceCpp
#' @importFrom MCMCpack rdirichlet
#' @importFrom mvtnorm rmvnorm dmvnorm
#' @importFrom entropy entropy.empirical
#' @importFrom lda top.topic.words
#' @importFrom reshape melt
#' @importFrom coda mcmc
NULL

#' @title Selected
#' @description Convert matrix output of rmultinom into vector of chosen items
#'
#' @param nSample number of samples to draw
#' @param proportion numeric non-negative vector of length K, specifying the probability for the K classes
#'
#' @return The vector of chosen items from the rmultinom sampling
#'
#' @export
Selected = function(nSample, proportion) {
	# Convert matrix output of rmultinom into vector of chosen items
	#
	# Args 
	#  nSample number of samples to draw
	#  proportion non-negative vector, specifying the probability for the K classes
	#
	# Returns
	#  The vector of chosen items from the rmultinom sampling
	
	multinom.chosen  = rmultinom(nSample, 1, proportion)
	selected.item = vapply(1L:nSample, function(s) {
	               which(multinom.chosen[, s] == TRUE)
	               }, c(1))
	return(selected.item)
}

#' @title Netstats
#' @description Calculate the network statistics (dydadic, triadic, or degree) given the history of interactions
#'
#' @param historyIP list containing the weighted history of interactions for each IP
#' @param node nodelist containing the ID of nodes
#' @param sender specific timepoint that we are calculating the time difference from
#' @param netstat which type of network statistics to use ("dyadic", "triadic", "degree")
#' @param nIP total number of interaction patterns specified by the user
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
  intercept = rep(1, length(node))
  # matcolnames = c("intercept")
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
    triadic = Triadic(historyIP, node, sender)
    for (IP in 1L:nIP) {
      triadic[[IP]][, 1] = rowSums(triadic[[IP]][, c(1:4, 7)])
      triadic[[IP]][, 2] = rowSums(triadic[[IP]][, c(5:6, 8)])
      triadic[[IP]][, 3] = (triadic[[IP]][, 9])
      triadic[[IP]][, 4] = rowSums(triadic[[IP]][, c(1:4, 7) + 9])
      triadic[[IP]][, 5] = rowSums(triadic[[IP]][, c(5:6, 8) + 9])
      triadic[[IP]][, 6] = (triadic[[IP]][, 9 + 9])
  	  triadic[[IP]][, 7] = rowSums(triadic[[IP]][, c(1:4, 7) + 18])
      triadic[[IP]][, 8] = rowSums(triadic[[IP]][, c(5:6, 8) + 18])
      triadic[[IP]][, 9] = (triadic[[IP]][, 9 + 18])
      triadic[[IP]][, 10] = rowSums(triadic[[IP]][, c(1:4, 7) + 27])
      triadic[[IP]][, 11] = rowSums(triadic[[IP]][, c(5:6, 8) + 27])
      triadic[[IP]][, 12] = (triadic[[IP]][, 9 + 27])
     }
     # matcolnames = c(matcolnames, paste0("2-send", 1:3), paste0("2-receive", 1:3), 
                    paste0("sibling", 1:3),  paste0("cosibling", 1:3) )
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

#' @title TopicInEqZ
#' @description Calculate topic part of the equation used in multinomial sampling of Z
#'
#' @param K total number of topics specified by the user
#' @param currentZ current state of the assignment of topics 
#' @param alpha Dirichlet concentration prior for document-topic distribution
#' @param mvec Dirichlet base prior for document-topic distribution
#' @param document the id of specific edge of interest
#'
#' @return The vector of constants representing topic part of each K
#'
#' @export
TopicInEqZ = function(K, currentZ, alpha, mvec, document) {
	# Calculate topic-IP part of the equation used in multinomial sampling of Z
	#
	# Args 
 	#  K total number of topics specified by the user
	#  currentZ current state of the assignment of topics 
	#  alpha Dirichlet concentration prior for document-topic distribution
	#  mvec Dirichlet base prior for document-topic distribution
	#  document the id of specific edge of interest
	#
	# Returns 
	#  The vector of constants representing topic part of each K
	topics = currentZ[[document]]
	table.topics = tabulateC(topics, K) 
	const = log(table.topics - as.numeric(table.topics > 0) + alpha * mvec)
	return(const)
}

#' @title DataAug
#' @description 
#'
#' @param iJi_d
#' @param lambda_d 
#' @param delta 
#' @param timeinc
#'
#' @return 
#'
#' @export
DataAug = function(iJi_d, lambda_d, delta, timeinc_d, observedi) {
  prenum = iJi_d * log(delta * lambda_d) - log(delta * lambda_d + 1)
  lambdadiff = matrix(0, nrow = nrow(lambda_d), ncol = ncol(lambda_d))
  num = matrix(0, nrow = nrow(lambda_d), ncol = ncol(lambda_d))
  for (i in (1:nrow(lambda_d))){
    Ji = iJi_d[i,]
    for (j in (1:ncol(lambda_d))){
      num[i, j] = sum(prenum[i, -c(i,j)])
      lambdaplus = Ji
      lambdaplus[j] = 1
      lambdaminus = Ji
      lambdaminus[j] = 0
      lambdadiff[i, j] = mean(lambda_d[i, lambdaplus > 0]) -
      					ifelse(sum(lambdaminus > 0) > 0, mean(lambda_d[i, lambdaminus > 0]), 0)        							
    }
  }
  y = timeinc_d * lambdadiff - log(lambda_d)
  denom = log(1 + exp(y))
  denom[y > 35] = y[y > 35]
  denom[y < -10] = exp(y[y < -10])
  prob = num - denom
  return(prob)
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
    part4 = log(sum(k_table) -1 + alpha)
    for (k in 1:length(currentZ_d)) {
      table_Wk = tableW[[currentZ_d[k]]]
      if (length(textlist_d) > 0 ){
        part1 = log(table_Wk[textlist_d[k]] -(table_Wk[textlist_d[k]] > 0) + betas * nvec[textlist_d[k]])
      } else {part1 = 0}
      part2 = log(sum(table_Wk) - sum(table_Wk > 0) + betas)
      part3 = log(k_table[currentZ_d[k]] -1 + alpha * mvec[currentZ_d[k]])
      finalsum = finalsum + part1 - part2 + part3 - part4
    }
  }
  return(finalsum)
}

#' @title MCMC
#' @description Iterate Markov Chain Monte Carlo (MCMC) algorithm to sequentially update the assignments of C, Z, and B 
#'
#' @param edge list of document information with 3 elements (element 1 sender, element 2 receiver, element 3 time in unix.time format)
#' @param node nodelist containing the ID of nodes (ID starting from 1)
#' @param textlist list of text (length=number of documents in total) containing the words in each document
#' @param vocabulary all vocabularies used over the corpus
#' @param nIP total number of interaction patterns specified by the user
#' @param K total number of topics specified by the user
#' @param delta.B proposal distribution variance parameter for beta 
#' @param outer size of outer iteration 
#' @param n1 size of first inner iteration for updates of interaction patterns C
#' @param n2 size of second inner iteration for updates of topics K
#' @param n3 size of third inner iteration for updates of beta B
#' @param burn iterations to be discarded at the beginning of beta chain
#' @param thin the thinning interval of beta chain
#' @param netstat which type of network statistics to use ("dyadic", "triadic", "degree")
#' @param seed an integer value which controls random number generation
#' @param plot to plot the convergence diagnostics or not (TRUE/FALSE)
#'
#' @return MMCMC output containing IP assignment, topic assignment, and (beta, mu, delta) chain 
#'
#' @export
MCMC = function(edge, node, textlist, vocabulary, nIP, K, delta.B, 
                outer = 1000, n1 = 3, n2 = 3, n3 = 3300, burn = 300, thin = 10, 
                netstat = c("dyadic", "triadic", "degree"), seed = 1, plot = TRUE) {
  # Iterate MCMC algorithm to sequentially update the assignments of C, Z, and B
  #
  # Args 
  #  edge list of document information with 3 elements (sender, receiver, time)
  #  node nodelist containing the ID of nodes (ID starting from 1)
  #  textlist list of text containing the words in each document
  #  vocabulary all vocabularies used over the corpus
  #  nIP total number of interaction patterns specified by the user
  #  K total number of topics specified by the user
  #  delta.B proposal distribution variance parameter for beta 
  #  outer size of outer iteration
  #  n1 size of first inner iteration for updates of interaction patterns C
  #  n2 size of second inner iteration for updates of topics K
  #  n3 size of third inner iteration for updates of beta B
  #  burn iterations to be discarded at the beginning of beta chain
  #  thin the thinning interval of beta chain
  #  netstat which type of network statistics to use ("dyadic", "triadic", "degree")
  #  seed an integer value which controls random number generation
  #  plot to plot the convergence diagnostics or not (TRUE/FALSE)
  #
  # Returns 
  #  MCMC output containing IP assignment, topic assignment, and (beta, mu, delta) chain
  
  set.seed(seed)
  # trim the edge so that we only model edges after 16 days
	timeinc = c(0, vapply(seq(along = edge)[-1], function(d) {
  	as.numeric(edge[[d]][3]) - as.numeric(edge[[d-1]][3])
 	}, c(1)))

   edge2 = which((cumsum(timeinc) > 384) == TRUE)
   
  # initialize alpha, mvec, delta, nvec, eta, lvec, and gammas
	alpha = 50 / K
  	mvec = rep(1 / K, K)
  	betas = 10
 	W = length(vocabulary)
  	nvec = rep(1, W) / W
  	phi = lapply(1L:K, function(k) {
		rdirichlet(1, betas * nvec)
	})
		
	delta = rbeta(1, 1, 10)
	deltamat = delta
	
	# initialize C, theta and Z
	currentC = sample(1:nIP, K, replace = TRUE)
	theta = rdirichlet(length(edge), alpha * mvec)
	currentZ = lapply(seq(along = edge), function(d) {
    if (length(textlist[[d]]) > 0) {
    	Selected(length(textlist[[d]]), theta[d, ])
    	} else {
    		Selected(1, theta[d, ])
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
    P = 1 + L * (2 * ("dyadic" %in% netstat) + 4 * ("triadic" %in% netstat) + 2 *("degree" %in% netstat))
    bmat = list()
	for (IP in 1L:nIP) {
		bmat[[IP]] = matrix(0, nrow = P, ncol = (n3 - burn) / thin)
  	}

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
    for (o in 1L:outer) {
      cat("outer iteration = ", o, "\n")
      Beta.old = lapply(bmat, function(b) {
      			rowMeans(b)
         		})  
     # Update the hyperparameter alpha and mvec
      vec = AlphamvecOpt(K, currentZ, alpha, mvec)
      alpha = sum(vec)
      mvec = vec / alpha
    
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
		lambda[[d]] = Reduce('+', lapply(1L:nIP, function(IP) {
		   			 p.d[d, IP] * exp(XB[[IP]])
		  				}))
		diag(lambda[[d]]) = 0
		
		#calculate the resampling probability
		probij = DataAug(iJi[[d]], lambda[[d]], delta, timeinc[d], as.numeric(edge[[d]][1]))
		iJi[[d]] = matrix(vapply(exp(probij), function(x){
				  rbinom(1, 1, x)
		   		  }, c(1)), nrow = length(node))
		iJi[[d]][as.numeric(edge[[d]][1]),] = tabulateC(as.numeric(unlist(edge[[d]][2])), length(node))
		LambdaiJi[[d]] = rowSums(vapply(1L:nIP, function(IP) {
						p.d[d, IP] * exp(rowSums(XB[[IP]] * iJi[[d]]) / rowSums(iJi[[d]] > 0))
						}, rep(0, length(node))))
		nonemptyiJi[[d]] = LambdaiJi[[d]][!is.na(LambdaiJi[[d]])]
		observediJi[[d]] = LambdaiJi[[d]][as.numeric(edge[[d]][1])]
		}	 
    
    cat("inner iteration 1", "\n")
    textlist.raw = unlist(textlist)
    table.W = lapply(1L:K, function(k) {
      tabulateC(textlist.raw[which(unlist(currentZ) == k)], W)
      })
    
    for (i1 in 1L:n1) {
      for (d in edge2) { 
        textlist.d = textlist[[d]]
        if (length(textlist.d) > 0) {
        	topicpart.d = TopicInEqZ(K, currentZ, alpha, mvec, d)
       	wordpart.d = WordInEqZ(K, textlist.d, table.W, betas, nvec)
        } else {
        topicpart.d = 0
        wordpart.d = matrix(0, nrow = length(currentZ[[d]]), ncol = K)
        }
        edgepart.d = EdgeInEqZ(iJi[[d]], lambda[[d]], delta)
        timepart.d = TimeInEqZ(nonemptyiJi[[d]], timeinc[d]) 
        observed.d = ObservedInEqZ(nonemptyiJi[[d]], observediJi[[d]]) 
        for (w in 1L:length(currentZ[[d]])) {
          const.Z = topicpart.d + edgepart.d + timepart.d + observed.d + wordpart.d[w, ]
          while (sum(exp(const.Z)) == 0) {
          	const.Z = const.Z + 1000
          }
           while (sum(exp(const.Z)) == Inf) {
          	const.Z = const.Z - 1000
          }
          zw.old = currentZ[[d]][w]
          zw.new = Selected(1, exp(const.Z))
          if (zw.new != zw.old) {
            currentZ[[d]][w] = zw.new
            topicpart.d = TopicInEqZ(K, currentZ, alpha, mvec, d)
            if (length(textlist.d) > 0) {	
            wordpart.d = WordInEqZ(K, textlist.d, table.W, betas, nvec)
            }
            table.W = lapply(1L:K, function(k) {
      				  tabulateC(textlist.raw[which(unlist(currentZ) == k)], W)
      				  })
      		p.d[d, ] = vapply(1L:nIP, function(IP) {
	 			sum(currentZ[[d]] %in% which(currentC == IP))
	 			}, c(1)) / length(currentZ[[d]])
      		LambdaiJi[[d]] = rowSums(vapply(1L:nIP, function(IP) {
							p.d[d, IP] * exp(rowSums(XB[[IP]] * iJi[[d]]) / rowSums(iJi[[d]] > 0))
							}, rep(0, length(node))))
			nonemptyiJi[[d]] = LambdaiJi[[d]][!is.na(LambdaiJi[[d]])]
         	observediJi[[d]] = LambdaiJi[[d]][as.numeric(edge[[d]][1])]
          }
        }
       }
      }

    cat("inner iteration 2", "\n")
    for (i2 in 1L:n2) {
      # C update given Z and B - within each document d
      for (k in 1L:K) { 
        document.k = which(vapply(currentZ, function(d){k %in% d}, c(1)) == 1)
        document.k = document.k[document.k %in% edge2]
        const.C = rep(NA, nIP)
        for (IP in 1L:nIP) {
          if (!currentC[k] == IP) {
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
           lambda[[d]] = Reduce('+', lapply(1L:nIP, function(IP) {
		   			    p.d[d, IP] * exp(XB[[IP]])
		  				}))
		  diag(lambda[[d]]) = 0
		  LambdaiJi[[d]] = rowSums(vapply(1L:nIP, function(IP) {
						p.d[d, IP] * exp(rowSums(XB[[IP]] * iJi[[d]]) / rowSums(iJi[[d]] > 0))
						}, rep(0, length(node))))
		  nonemptyiJi[[d]] = LambdaiJi[[d]][!is.na(LambdaiJi[[d]])]
          observediJi[[d]] = LambdaiJi[[d]][as.numeric(edge[[d]][1])]
          }
          }
          const.C[IP] = sum(vapply(document.k, function(d) {
            EdgeInEqZ(iJi[[d]], lambda[[d]], delta) + TimeInEqZ(nonemptyiJi[[d]], timeinc[d]) + 
    			ObservedInEqZ(nonemptyiJi[[d]], observediJi[[d]]) 
          }, c(1))) / length(document.k)
          }
          while (sum(exp(const.C)) == 0) {
            const.C = const.C + 1000
          }
          while (sum(exp(const.Z)) == Inf) {
          	const.Z = const.Z - 1000
          }
          currentC[k] = Selected(1, exp(const.C))
      }
    }
      
    if (plot) {
      entropy.mat = c(entropy.mat, entropy.empirical(currentC))
      alpha.mat = rbind(alpha.mat, alpha)
      logWZ.mat = c(logWZ.mat, 
                    logWZ(K, currentZ, textlist, table.W, alpha, mvec, betas, nvec))
      }
      
     cat("inner iteration 3", "\n")
     for (d in edge2) {
     	 p.d[d, ] =  vapply(1L:nIP, function(IP) {
              		 sum(currentZ[[d]] %in% which(currentC == IP))
                     }, c(1)) / length(currentZ[[d]])
     	 history.t = History(edge, p.d, node, as.numeric(edge[[d]][3]))
    	     X = lapply(node, function(i) {
              Netstats(history.t, node, i, netstat)
             })
    	     XB = MultiplyXBList(X, Beta.old)   
    	     lambda[[d]] = Reduce('+', lapply(1L:nIP, function(IP) {
		   			 p.d[d, IP] * exp(XB[[IP]])
		  			 }))
		 diag(lambda[[d]]) = 0
		 LambdaiJi[[d]] = rowSums(vapply(1L:nIP, function(IP) {
						p.d[d, IP] * exp(rowSums(XB[[IP]] * iJi[[d]]) / rowSums(iJi[[d]] > 0))
						}, rep(0, length(node))))
		 nonemptyiJi[[d]] = LambdaiJi[[d]][!is.na(LambdaiJi[[d]])]
         observediJi[[d]] = LambdaiJi[[d]][as.numeric(edge[[d]][1])]
	 }
	 
	 # Beta and delta update
	 prior.old1 = sum(vapply(1L:nIP, function(IP) {
		  		 dmvnorm(Beta.old[[IP]], rep(0, P), sigma^2 * diag(P), log = TRUE)
		  		 }, c(1))) 
	 post.old1 = sum(vapply(edge2, function(d) {
	     	     EdgeInEqZ(iJi[[d]], lambda[[d]], delta) + TimeInEqZ(nonemptyiJi[[d]], timeinc[d]) + 
    			 ObservedInEqZ(nonemptyiJi[[d]], observediJi[[d]])
    			 }, c(1))) / length(edge2)
    	    			    
     for (i3 in 1L:n3) {
     	Beta.new = lapply(1L:nIP, function(IP) {
           rmvnorm(1, Beta.old[[IP]], (delta.B)^2 * diag(P))
         }) 
        for (d in edge2) {
           XB = MultiplyXBList(X, Beta.new)    
           LambdaiJi[[d]] = rowSums(vapply(1L:nIP, function(IP) {
						p.d[d, IP] * exp(rowSums(XB[[IP]] * iJi[[d]]) / rowSums(iJi[[d]] > 0))
						}, rep(0, length(node))))
		   nonemptyiJi[[d]] = LambdaiJi[[d]][!is.na(LambdaiJi[[d]])]
		   observediJi[[d]] = LambdaiJi[[d]][as.numeric(edge[[d]][1])]
         }
         prior.new1 = sum(vapply(1L:nIP, function(IP) {
        		dmvnorm(Beta.new[[IP]], rep(0, P), sigma^2 * diag(P), log = TRUE)
         }, c(1))) 
         post.new1 = sum(vapply(edge2, function(d) {
    			    EdgeInEqZ(iJi[[d]], lambda[[d]], delta) + TimeInEqZ(nonemptyiJi[[d]], timeinc[d]) + 
    			ObservedInEqZ(nonemptyiJi[[d]], observediJi[[d]])
    			    }, c(1))) / length(edge2)
    		 loglike.diff = prior.new1 + post.new1 - prior.old1 - post.old1
      
        if (log(runif(1, 0, 1)) < loglike.diff) {
         for (IP in 1L:nIP) {
         Beta.old[[IP]]  = Beta.new[[IP]]
         }
         prior.old1 = prior.new1
         post.old1 = post.new1
         }
      
         if (i3 > burn && i3 %% (thin) == 0) {
          for (IP in 1L:nIP) {
           bmat[[IP]][ , (i3 - burn) / thin] = Beta.old[[IP]]
           }
          }
         }
         
         prior.old2 = dbeta(delta, 1, 10) 
         post.old2 = sum(vapply(edge2, function(d) {
    			    EdgeInEqZ(iJi[[d]], lambda[[d]], delta)
    			    }, c(1))) / length(edge2)
    			    
    	 delta.new = rnorm(1, delta, (delta.B)^2)
       	 while (delta.new < 0) {
       	 delta.new = rnorm(1, delta, (delta.B)^2)
       	 }
		 prior.new2 = dbeta(delta.new, 1, 10)
		 post.new2 = sum(vapply(edge2, function(d) {
    			    EdgeInEqZ(iJi[[d]], lambda[[d]], delta.new)
    			    }, c(1))) / length(edge2)
    	 loglike.diff1 = prior.new2 + post.new2 - prior.old2 - post.old2
		if (log(runif(1, 0, 1)) < loglike.diff1) {
			delta = delta.new
			prior.old1 = prior.new1
			post.old1 = post.new1
			}
			deltamat = c(deltamat, delta)
	 }

    
    if (plot) {
    burnin = round(outer / 10)
    par(mfrow = c(2, 2))
  	plot(entropy.mat[-1L:-burnin], type = "l", 
  	     xlab = "(Outer) Iterations", ylab = "Entropy of IP")
  	abline(h = mean(entropy.mat[-1L:-burnin]), lty = 1)
  	title("Convergence of Entropy")
  	matplot(alpha.mat[-1L:-burnin,], lty = 1, type = "l", col = 1L:nIP, 
  	        xlab = "(Outer) Iterations", ylab = "alpha")
  	abline(h = mean(alpha.mat[-1L:-burnin,]), lty = 1, col = 1L:nIP)
	  title("Convergence of Optimized alpha")
	  plot(logWZ.mat[-1L:-burnin], type = "l", 
	       xlab = "(Outer) Iterations", ylab = "logWZ")
	  abline(h = mean(logWZ.mat[-1L:-burnin]), lty = 1)
	  title("Convergence of logWZ")
	  matplot(bmat[[1]][1,], lty = 1, col = 1L:P, type = "l", 
	          main = "Traceplot of beta", xlab = "(Inner) Iterations", ylab = "")
	  abline(h = mean(bmat[[1]][1,]), lty = 1, col = 1L)
	  }
     
  chain.final = list(C = currentC, Z = currentZ, B = bmat, D = deltamat)

  return(chain.final)
}


# replace M-H with slice sampling
MCMC1 = function(edge, node, textlist, vocabulary, nIP, K, delta.B, 
                outer = 1000, n1 = 3, n2 = 3, n3 = 3300, burn = 300, thin = 10, 
                netstat = c("dyadic", "triadic", "degree"), seed = 1, plot = TRUE) {
  # Iterate MCMC algorithm to sequentially update the assignments of C, Z, and B
  #
  # Args 
  #  edge list of document information with 3 elements (sender, receiver, time)
  #  node nodelist containing the ID of nodes (ID starting from 1)
  #  textlist list of text containing the words in each document
  #  vocabulary all vocabularies used over the corpus
  #  nIP total number of interaction patterns specified by the user
  #  K total number of topics specified by the user
  #  delta.B proposal distribution variance parameter for beta 
  #  outer size of outer iteration
  #  n1 size of first inner iteration for updates of interaction patterns C
  #  n2 size of second inner iteration for updates of topics K
  #  n3 size of third inner iteration for updates of beta B
  #  burn iterations to be discarded at the beginning of beta chain
  #  thin the thinning interval of beta chain
  #  netstat which type of network statistics to use ("dyadic", "triadic", "degree")
  #  seed an integer value which controls random number generation
  #  plot to plot the convergence diagnostics or not (TRUE/FALSE)
  #
  # Returns 
  #  MCMC output containing IP assignment, topic assignment, and (beta, mu, delta) chain
  
  set.seed(seed)

  # initialize alpha, mvec, delta, nvec, eta, lvec, and gammas
	alpha = 50 / K
  	mvec = rep(1 / K, K)
  	betas = 10
 	W = length(vocabulary)
  	nvec = rep(1, W) / W
  	phi = lapply(1L:K, function(k) {
		rdirichlet(1, betas * nvec)
	})
		
	delta = rbeta(1, 1, 10)
	deltamat = delta
	
	# initialize C, theta and Z
	currentC = sample(1:nIP, K, replace = TRUE)
	theta = rdirichlet(length(edge), alpha * mvec)
	currentZ = lapply(seq(along = edge), function(d) {
    if (length(textlist[[d]]) > 0) {
    	Selected(length(textlist[[d]]), theta[d, ])
    	} else {
    		Selected(1, theta[d, ])
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
    P = 1 + (2 * L + 2) * ("dyadic" %in% netstat) + (4 * L^2 + 4) * ("triadic" %in% netstat) + (2 * L + 2) *("degree" %in% netstat)
    bmat = list()
	for (IP in 1L:nIP) {
		bmat[[IP]] = matrix(0, nrow = P, ncol = (n3 - burn) / thin)
  	}

  # to check the convergence  
  if (plot) {
  	logWZ.mat = c()							  
  	alpha.mat = c()
  	entropy.mat = c()
  }
  	
  iJi = lapply(seq(along = edge), function(d) {
  	  matrix(0, nrow = length(node), ncol = length(node))
  	})
  timeinc = c(0, vapply(seq(along = edge)[-1], function(d) {
  	as.numeric(edge[[d]][3]) - as.numeric(edge[[d-1]][3])
  }, c(1)))
    
   #start outer iteration
  for (o in 1L:outer) {
    cat("outer iteration = ", o, "\n")
    
    Beta.old = lapply(bmat, function(b) {
      			rowMeans(b)
         		})  
    # Update the hyperparameter alpha and mvec
    vec = AlphamvecOpt(K, currentZ, alpha, mvec)
    alpha = sum(vec)
    mvec = vec / alpha
    
    # Data augmentation
    lambda = list()
    LambdaiJi = list()
    nonemptyiJi = list()
    nonemptyiJi2 = list()

    for (d in seq(along = edge)) {
   	 	history.t = History(edge, p.d, node, as.numeric(edge[[d]][3]))
   	 	X = lapply(node, function(i) {
  	        Netstats(history.t, node, i, netstat)
            })
   	 	XB = MultiplyXBList(X, Beta.old)     
		lambda[[d]] = Reduce('+', lapply(1L:nIP, function(IP) {
		   			 p.d[d, IP] * exp(XB[[IP]])
		  				}))
		diag(lambda[[d]]) = 0
		
		#calculate the resampling probability
		probij = DataAug(iJi[[d]], lambda[[d]], delta, timeinc[d], as.numeric(edge[[d]][1]))
		iJi[[d]] = matrix(vapply(exp(probij), function(x){
				  rbinom(1, 1, x)
		   		  }, c(1)), nrow = length(node))
		observedi = as.numeric(edge[[d]][1])
		observedj = as.numeric(unlist(edge[[d]][2]))
		iJi[[d]][observedi,] = tabulateC(observedj, length(node))
		LambdaiJi[[d]] = rowSums(vapply(1L:nIP, function(IP) {
						p.d[d, IP] * exp(rowSums(XB[[IP]] * iJi[[d]]) / rowSums(iJi[[d]] > 0))
						}, rep(0, length(node))))
		nonemptyiJi[[d]] = LambdaiJi[[d]][!is.na(LambdaiJi[[d]])]
         nonemptyiJi2[[d]] = (LambdaiJi[[d]][-as.numeric(edge[[d]][1])])[!is.na(LambdaiJi[[d]][-as.numeric(edge[[d]][1])])]
		}	 
    
    cat("inner iteration 1", "\n")
    textlist.raw = unlist(textlist)
    table.W = lapply(1L:K, function(k) {
      tabulateC(textlist.raw[which(unlist(currentZ) == k)], W)
      })
    
    for (i1 in 1L:n1) {
      for (d in seq(along = edge)) {
        textlist.d = textlist[[d]]
        if (length(textlist.d) > 0) {
        	topicpart.d = TopicInEqZ(K, currentZ, alpha, mvec, d)
       	wordpart.d = WordInEqZ(K, textlist.d, table.W, betas, nvec)
        } else {
        topicpart.d = 0
        wordpart.d = matrix(0, nrow = length(currentZ[[d]]), ncol = K)
        }
        edgepart.d = EdgeInEqZ2(iJi[[d]], lambda[[d]], delta)
        timepart.d = TimeInEqZ(nonemptyiJi[[d]]) + ObservedInEqZ(sum(nonemptyiJi[[d]]) + sum(nonemptyiJi2[[d]]), timeinc[d], delta) 
        for (w in 1L:length(currentZ[[d]])) {
          const.Z = topicpart.d + edgepart.d + timepart.d + wordpart.d[w, ]
          while (sum(exp(const.Z)) == 0) {
          	const.Z = const.Z + 1000
          }
           while (sum(exp(const.Z)) == Inf) {
          	const.Z = const.Z - 1000
          }
          zw.old = currentZ[[d]][w]
          zw.new = Selected(1, exp(const.Z))
          if (zw.new != zw.old) {
            currentZ[[d]][w] = zw.new
            topicpart.d = TopicInEqZ(K, currentZ, alpha, mvec, d)
            if (length(textlist.d) > 0) {	
            wordpart.d = WordInEqZ(K, textlist.d, table.W, betas, nvec)
            }
            table.W = lapply(1L:K, function(k) {
      				tabulateC(textlist.raw[which(unlist(currentZ) == k)], W)
      				})
      		p.d[d, ] = vapply(1L:nIP, function(IP) {
	 			sum(currentZ[[d]] %in% which(currentC == IP))
	 			}, c(1)) / length(currentZ[[d]])
      		LambdaiJi[[d]] = rowSums(vapply(1L:nIP, function(IP) {
							p.d[d, IP] * exp(rowSums(XB[[IP]] * iJi[[d]]) / rowSums(iJi[[d]] > 0))
							}, rep(0, length(node))))
			nonemptyiJi[[d]] = LambdaiJi[[d]][!is.na(LambdaiJi[[d]])]
         	nonemptyiJi2[[d]] = (LambdaiJi[[d]][-as.numeric(edge[[d]][1])])[!is.na(LambdaiJi[[d]][-as.numeric(edge[[d]][1])])]
				
          }
        }
       }
      }

    cat("inner iteration 2", "\n")
    for (i2 in 1L:n2) {
      # C update given Z and B - within each document d
      for (k in 1L:K) {
        document.k = which(vapply(currentZ, function(d){k %in% d}, c(1)) ==1)
        const.C = rep(NA, nIP)
        for (IP in 1L:nIP) {
          if (!currentC[k] == IP) {
          currentC[k] = IP
          for (d in document.k) {
            p.d[d, ] =  vapply(1L:nIP, function(IP) {
              sum(currentZ[[d]] %in% which(currentC == IP))
            }, c(1)) / length(currentZ[[d]])
          }
          for (d in seq(along = edge)) {
            history.t = History(edge, p.d, node, as.numeric(edge[[d]][3]))
    	   	   X = lapply(node, function(i) {
               Netstats(history.t, node, i, netstat)
               })
    	       XB = MultiplyXBList(X, Beta.old)    
           lambda[[d]] = Reduce('+', lapply(1L:nIP, function(IP) {
		   			    p.d[d, IP] * exp(XB[[IP]])
		  				}))
		  diag(lambda[[d]]) = 0
		  LambdaiJi[[d]] = rowSums(vapply(1L:nIP, function(IP) {
						p.d[d, IP] * exp(rowSums(XB[[IP]] * iJi[[d]]) / rowSums(iJi[[d]] > 0))
						}, rep(0, length(node))))
		  nonemptyiJi[[d]] = LambdaiJi[[d]][!is.na(LambdaiJi[[d]])]
           nonemptyiJi2[[d]] = (LambdaiJi[[d]][-as.numeric(edge[[d]][1])])[!is.na(LambdaiJi[[d]][-as.numeric(edge[[d]][1])])]
          }
          }
          const.C[IP] = sum(vapply(document.k, function(d) {
              EdgeInEqZ2(iJi[[d]], lambda[[d]], delta) + TimeInEqZ(nonemptyiJi[[d]]) + 
    			 ObservedInEqZ(sum(nonemptyiJi[[d]]) + sum(nonemptyiJi2[[d]]), timeinc[d], delta)
          }, c(1))) / length(document.k)
          }
          while (sum(exp(const.C)) == 0) {
            const.C = const.C + 1000
          }
          while (sum(exp(const.Z)) == Inf) {
          	const.Z = const.Z - 1000
          }
          currentC[k] = Selected(1, exp(const.C))
      }
    }
      
    if (plot) {
      entropy.mat = c(entropy.mat, entropy.empirical(currentC))
      alpha.mat = rbind(alpha.mat, alpha)
      logWZ.mat = c(logWZ.mat, 
                    logWZ(K, currentZ, textlist, table.W, alpha, mvec, betas, nvec))
      }
      
    cat("inner iteration 3", "\n")
    Beta.new = list()
    for (d in seq(along = edge)) {
     	 p.d[d, ] =  vapply(1L:nIP, function(IP) {
              sum(currentZ[[d]] %in% which(currentC == IP))
            }, c(1)) / length(currentZ[[d]])
     	history.t = History(edge, p.d, node, as.numeric(edge[[d]][3]))
    	     X = lapply(node, function(i) {
              Netstats(history.t, node, i, netstat)
             })
    	     XB = MultiplyXBList(X, Beta.old)   
    	    lambda[[d]] = Reduce('+', lapply(1L:nIP, function(IP) {
		   			 p.d[d, IP] * exp(XB[[IP]])
		  			 }))
		diag(lambda[[d]]) = 0
		LambdaiJi[[d]] = rowSums(vapply(1L:nIP, function(IP) {
						p.d[d, IP] * exp(rowSums(XB[[IP]] * iJi[[d]]) / rowSums(iJi[[d]] > 0))
						}, rep(0, length(node))))
		nonemptyiJi[[d]] = LambdaiJi[[d]][!is.na(LambdaiJi[[d]])]
         nonemptyiJi2[[d]] = (LambdaiJi[[d]][-as.numeric(edge[[d]][1])])[!is.na(LambdaiJi[[d]][-as.numeric(edge[[d]][1])])]
	 }
	delta.B2 = 0.01
	# below is the f(x0)
    	 post.old1 = sum(vapply(seq(along = edge), function(d) {
    			    EdgeInEqZ2(iJi[[d]], lambda[[d]], delta) + ObservedInEqZ(sum(nonemptyiJi[[d]]) + sum(nonemptyiJi2[[d]]), timeinc[d], delta)
    			    }, c(1))) / length(edge)	  		 
     post.old2 = post.old1 + sum(vapply(seq(along = edge), function(d) {
    			    TimeInEqZ(nonemptyiJi[[d]])
    			    }, c(1))) / length(edge)
   	    
    y1 = runif(1, 0, exp(post.old1))
    y2 = runif(1, 0, exp(post.old2))
    	Beta.L = lapply(1L:nIP, function(IP) {
    	Beta.old[[IP]] - delta.B * runif(P, 0, 1)
    		}) 
    	Beta.R = lapply(1L:nIP, function(IP) {
    	Beta.L[[IP]] + delta.B
    		})
	delta.L = delta - delta.B2 * runif(1, 0, 1)
    	delta.R = delta.L + delta.B2	   

    for (i3 in 1L:n3) {
    	 delta.new = delta.L + runif(1, 0, 1) * (delta.R - delta.L)
        while (delta.new < 0) {
        		delta.new = delta.L + runif(1, 0, 1) * (delta.R - delta.L)
        }
     post.new1 = sum(vapply(seq(along = edge), function(d) {
    			    EdgeInEqZ2(iJi[[d]], lambda[[d]], delta.new) + ObservedInEqZ(sum(nonemptyiJi[[d]]) + sum(nonemptyiJi2[[d]]), timeinc[d], delta.new)
    			    }, c(1))) / length(edge)	  	   
     if (y1 >= exp(post.new1)) {
        delta = delta.new
        post.old1 = post.new1
        y1 = runif(1, 0, exp(post.old1))
        delta.L = delta - delta.B2 * runif(1, 0, 1)
    	    delta.R = delta.L + delta.B2	
		} else {
			if (delta.new < delta) {
				delta.L = delta.new
				} else {
					delta.R = delta.new
				}	
		}
             
    	Beta.new = lapply(1L:nIP, function(IP) {
        	Beta.L[[IP]] + runif(P, 0, 1) * (Beta.R[[IP]] - Beta.L[[IP]])
        }) 
         for (d in seq(along = edge)) {
           XB = MultiplyXBList(X, Beta.new)    
           LambdaiJi[[d]] = rowSums(vapply(1L:nIP, function(IP) {
						p.d[d, IP] * exp(rowSums(XB[[IP]] * iJi[[d]]) / rowSums(iJi[[d]] > 0))
						}, rep(0, length(node))))
		  nonemptyiJi[[d]] = LambdaiJi[[d]][!is.na(LambdaiJi[[d]])]
		  nonemptyiJi2[[d]] = (LambdaiJi[[d]][-as.numeric(edge[[d]][1])])[!is.na(LambdaiJi[[d]][-as.numeric(edge[[d]][1])])]
         }
         post.new2 = post.old1 + sum(vapply(seq(along = edge), function(d) {
    			    TimeInEqZ(nonemptyiJi[[d]])
    			    }, c(1))) / length(edge)
		if (y2 >= exp(post.new2)) {
      	for (IP in 1L:nIP) {
        	Beta.old[[IP]]  = Beta.new[[IP]]
        	}
        post.old2 = post.new2
         y2 = runif(1, 0, exp(post.old2))
    		Beta.L = lapply(1L:nIP, function(IP) {
    		Beta.old[[IP]] - delta.B * runif(P, 0, 1)
    		}) 
    		Beta.R = lapply(1L:nIP, function(IP) {
    		Beta.L[[IP]] + delta.B
    		})
		} else {
        	for (IP in 1L:nIP) {
        		for (p in 1L:P) {
        			if (Beta.new[[IP]][p] < Beta.old[[IP]][p]) {
        				Beta.L[[IP]][p] = Beta.new[[IP]][p]
        			} else {
        				Beta.R[[IP]][p] = Beta.new[[IP]][p]
        			}
        		}
        	}
        }
        
         if (i3 > burn && i3 %% (thin) == 0) {
       	deltamat = c(deltamat, delta)
         for (IP in 1L:nIP) {
           bmat[[IP]][ , (i3 - burn) / thin] = Beta.old[[IP]]
           }
         }
              
        }	
     }
    
    if (plot) {
    burnin = round(outer / 10)
    par(mfrow = c(2, 2))
  	plot(entropy.mat[-1L:-burnin], type = "l", 
  	     xlab = "(Outer) Iterations", ylab = "Entropy of IP")
  	abline(h = mean(entropy.mat[-1L:-burnin]), lty = 1)
  	title("Convergence of Entropy")
  	matplot(alpha.mat[-1L:-burnin,], lty = 1, type = "l", col = 1L:nIP, 
  	        xlab = "(Outer) Iterations", ylab = "alpha")
  	abline(h = mean(alpha.mat[-1L:-burnin,]), lty = 1, col = 1L:nIP)
	  title("Convergence of Optimized alpha")
	  plot(logWZ.mat[-1L:-burnin], type = "l", 
	       xlab = "(Outer) Iterations", ylab = "logWZ")
	  abline(h = mean(logWZ.mat[-1L:-burnin]), lty = 1)
	  title("Convergence of logWZ")
	  matplot(bmat[[1]][1,], lty = 1, col = 1L:P, type = "l", 
	          main = "Traceplot of beta", xlab = "(Inner) Iterations", ylab = "")
	  abline(h = mean(bmat[[1]][1,]), lty = 1, col = 1L)
	  }
     
  chain.final = list(C = currentC, Z = currentZ, B = bmat, D = deltamat)

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
#' @param MCMCchain a chain obtained using MCMC function
#'
#' @return Joint boxplot of the posterior summary of beta (should click for the locator)
#'
#' @export
PlotBetaIP = function(MCMCchain) {
  # Draw a boxplot of the MCMC chain of beta for each IP
  #
  # Args 
  #  MCMCchain a chain obtained using MCMC function
  #
  # Returns
  #  Joint boxplot of the posterior summary of beta (should click for the locator)
  combined.data = list()
  P = nrow(MCMCchain$B[[1]])
  nIP = length(MCMCchain$B)
  for(b in 1L:P){
    combined.data[[b]] = sapply(1L:nIP, function(c){
      cbind(MCMCchain$B[[c]][b,])
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
       })))
  legend(locator(1), c(paste("IP", 1L:nIP)), col = gray.colors(nIP), pch = 15)
}

#' @title PlotTopic
#' @description Draw a barplot of the topic distributions without considering interaction patterns
#'
#' @param MCMCchain MCMCchain a chain obtained using MCMC function
#' @param K total number of topics specified by the user
#'
#' @return Barplot of the topic distribution
#'
#' @export	
PlotTopic = function(MCMCchain, K) {
  # Draw a barplot of the topic distributions without considering IPs
  #
  # Args 
  #  MCMCchain a chain obtained using MCMC function
  #  K total number of topics specified by the user
  #
  # Returns
  #  Barplot of the topic distribution
  Zsummary = list()
  for (d in seq(along = MCMCchain$Z)) {
    Zsummary[[d]] = MCMCchain$Z[[d]]
  }
  topic.dist = t(tabulateC(unlist(Zsummary), K) / length(unlist(Zsummary)))
  colnames(topic.dist) = c(1L:K)
  barplot(topic.dist, beside = TRUE, xlab = "Topic", ylab = "Proportion", 
          main = "Topic Distritubitions without IPs")
}

#' @title TableWordIP
#' @description Generate a table that summarizes token-topic assignments with high probabilities for each interaaction pattern
#'
#' @param MCMCchain MCMCchain a chain obtained using MCMC function
#' @param K total number of topics specified by the user
#' @param textlist list of text (length=number of documents in total) containing the words in each document
#' @param vocabulary all vocabularies used over the corpus
#'
#' @return List of table that summarize token-topic assignments with high probabilities (topic proportion > 0.1) for each IP
#'
#' @export
TableWord = function(MCMCchain, K, textlist, vocabulary) {
  # Generate a table of token-topic assignments with high probabilities for each IP
  #
  # Args 
  #  MCMCchain a chain obtained using MCMC function
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
      if (length(MCMCchain$Z[[d]])>0){
        Zsummary[[iter]] = MCMCchain$Z[[d]]
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
#' @description Generate the documents according to the generative process
#'
#' @param nDocs number of documents to be generated
#' @param node nodelist containing the ID of nodes (ID starting from 1)
#' @param vocabulary all vocabularies used over the corpus
#' @param nIP total number of interaction patterns specified by the user
#' @param K total number of topics specified by the user
#' @param xi Poisson parameter for the number of words
#' @param netstat which type of network statistics to use ("dyadic", "triadic", "degree")
#' @param seed an integer value which controls random number generation
#'
#' @return List of edge and text generated according to the generative process
#'
#' @export
GenerateDocs = function(nDocs, node, vocabulary, nIP, K, xi = 10,
						netstat = c("degree","dyadic", "triadic"), seed = 1) {
  # Generate the documents according to the generative process
  #
  # Args 
  #  nDocs number of documents to be generated
  #  node nodelist containing the ID of nodes (ID starting from 1)
  #  vocabulary all vocabularies used over the corpus
  #  nIP total number of interaction patterns specified by the user
  #  K total number of topics specified by the user
  #  xi Poisson parameter for the number of words
  #  netstat which type of network statistics to use ("dyadic", "triadic", "degree")
  #  seed an integer value which controls random number generation
  #
  # Returns 
  #  List of edge and text generated according to the generative process

	set.seed(seed)
	alpha = 50 / K
  	mvec = rep(1 / K, K)
  	betas = 10
 	W = length(vocabulary)
  	nvec = rep(1, W) / W
  	phi = lapply(1L:K, function(k) {
		rdirichlet(1, betas * nvec)
	})
	L = 3
	P = 1 + (2 * L + 2) * ("dyadic" %in% netstat) + (4 * L^2 + 4) * ("triadic" %in% netstat) + (2 * L + 2) *("degree" %in% netstat)
	b = lapply(1L:nIP, function(IP) {
		rmvnorm(1, rep(0, P), diag(P))
	})
	
	delta = rbeta(1, 1, 10)
	
	currentC = sample(1L:nIP, K, replace = TRUE)

	t.d = 0
	edge = list()
	text = list()
	p.d = matrix(NA, nrow = nDocs, ncol = nIP)
	
	options(warn = -1)
	for (d in 1:nDocs) {
		N.d = rpois(1, xi)
		theta.d = rdirichlet(1, alpha * mvec)	
		topic.d = Selected(max(1, N.d), theta.d)
		text[[d]] = rep(NA, N.d)
		if (N.d > 0) {
			for (n in 1:N.d){
				text[[d]][n] = vocabulary[Selected(1, phi[[topic.d[n]]])]
				}
			names(text[[d]]) = topic.d
		}

		p.d[d, ] = vapply(1L:nIP, function(IP) {
			sum(topic.d %in% which(currentC == IP))
			}, c(1)) / max(1, N.d)
		history.t = History(edge, p.d, node, t.d + 10^(-10))
		X= lapply(node, function(i) {
				Netstats(history.t, node, i, netstat)
				})
		XB = lapply(1L:nIP, function(IP) {
    	  	 t(vapply(node, function(i) {
    	     MultiplyXB(X[[i]][[IP]], Beta.old[[IP]])
    	     }, rep(0, length(node))))
    	     }) 
		lambda = Reduce('+', lapply(1L:nIP, function(IP) {
		   			 p.d[d, IP] * exp(XB[[IP]])
		  		 }))
		diag(lambda) = 0
		
		iJi = matrix(vapply(lambda, function(x){
				   rbinom(1, 1, 1 - exp(-delta * x))
		   			}, c(1)), nrow = length(node))
		while (sum(iJi) == 0) {
			iJi = matrix(vapply(lambda, function(x){
				   rbinom(1, 1, 1 - exp(-delta * x))
		   			}, c(1)), nrow = length(node))
		}
		LambdaiJi = rowSums(vapply(1L:nIP, function(IP) {
			p.d[d, IP] * exp(rowSums(XB[[IP]] * iJi[[d]]) / rowSums(iJi[[d]] > 0))
			}, rep(0, length(node))))	
		Time.inc = vapply(LambdaiJi, function(lambda) {
			rexp(1, lambda)
			}, c(1))
		i.d = which(Time.inc == min(Time.inc[!is.na(Time.inc)]))
		j.d = which(iJi[[i.d]] == 1)
		j.d[j.d >= i.d] = j.d[j.d >= i.d] + 1
		t.d = t.d + Time.inc[i.d]
		edge[[d]] = list(sender = i.d, receiver = j.d, timestamp = t.d)	
	}
	options(warn = 0)
	return(list(edge = edge, text = text))
}


 
      


