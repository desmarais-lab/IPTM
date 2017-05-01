library(mvtnorm)
library(MCMCpack)
library(entropy)
library(Rcpp)
library(RcppArmadillo)
sourceCpp('~/Desktop/IPTM-master/pkg/src/sampler.cpp')
source('~/Desktop/IPTM-master/pkg/R/core.R')

GenerateDocs_base = function(nDocs, node, vocabulary, nIP, K, nwords, alpha, mvec, betas, nvec, prior.b.mean, prior.b.var, prior.delta, netstat, seed) {
	set.seed(seed)
 	W = length(vocabulary)
  	phi = lapply(1L:K, function(k) {
		rdirichlet_cpp(1, betas * nvec)
	})
	L = 3
    	P = 1 + L * (2 * ("dyadic" %in% netstat) + 4 * ("triadic" %in% netstat) + 2 *("degree" %in% netstat))
	b = lapply(1L:nIP, function(IP) {
		rmvnorm(1, prior.b.mean, prior.b.var)
	})
	
	delta = rbeta(1, prior.delta[1], prior.delta[2])
	currentC = sample(1L:nIP, K, replace = TRUE)

	t.d = 0
	edge = list()
	text = list()
	p.d = matrix(NA, nrow = nDocs, ncol = nIP)
	
	options(warn = -1)
	for (d in 1:nDocs) {
		N.d = nwords
		theta.d = rdirichlet_cpp(1, alpha * mvec)	
		topic.d = multinom_vec(max(1, N.d), theta.d)
		text[[d]] = rep(NA, N.d)
		if (N.d > 0) {
			for (n in 1:N.d){
				text[[d]][n] = multinom_vec(1, phi[[topic.d[n]]])
				}
			names(text[[d]]) = topic.d
		}

		p.d[d, ] = vapply(1L:nIP, function(IP) {
			sum(topic.d %in% which(currentC == IP))
			}, c(1)) / max(1, N.d)
		if (t.d < 384) { 
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
		lambda = lambda_cpp(p.d[d,], XB)
		
		iJi = rbinom_mat((delta * lambda) / (delta * lambda + 1))
		while (sum(iJi) == 0) {
			iJi = rbinom_mat((delta * lambda) / (delta * lambda + 1))
		}
		LambdaiJi = lambdaiJi(p.d[d,], XB, iJi)
		Time.inc = vapply(LambdaiJi, function(lambda) {
			rexp(1, lambda)
			}, c(1))
		i.d = which(Time.inc == min(Time.inc[!is.na(Time.inc)]))
		j.d = which(iJi[i.d,] == 1)
		t.d = t.d + Time.inc[i.d]
		edge[[d]] = list(sender = i.d, receiver = j.d, timestamp = t.d)	
	}
	options(warn = 0)
	return(list(edge = edge, text = text, C = currentC, beta = b, delta = delta))
}

GenerateDocs_forward = function(nDocs, node, vocabulary, nIP, K, nwords, alpha, mvec, betas, nvec,
							    b, delta, currentC, netstat, base.edge, seed) {

	set.seed(seed)
 	W = length(vocabulary)
  	phi = lapply(1L:K, function(k) {
		rdirichlet_cpp(1, betas * nvec)
	})
	L = 3
    P = 1 + L * (2 * ("dyadic" %in% netstat) + 4 * ("triadic" %in% netstat) + 2 *("degree" %in% netstat))

	t.d = 384
	edge = base.edge
	text = list()
	p.d = matrix(NA, nrow = nDocs, ncol = nIP)
	
	options(warn = -1)
	for (d in 1:nDocs) {
		N.d = nwords
		theta.d = rdirichlet_cpp(1, alpha * mvec)	
		topic.d = multinom_vec(max(1, N.d), theta.d)
		text[[d]] = rep(NA, N.d)
		if (N.d > 0) {
			for (n in 1:N.d){
				text[[d]][n] = multinom_vec(1, phi[[topic.d[n]]])
				}
			names(text[[d]]) = topic.d
		}

		p.d[d, ] = vapply(1L:nIP, function(IP) {
			sum(topic.d %in% which(currentC == IP))
			}, c(1)) / max(1, N.d)
		history.t = History(edge, p.d, node, t.d + 10^(-10))
		X = lapply(node, function(i) {
  	        Netstats(history.t, node, i, netstat)
            })
   	 	XB = MultiplyXBList(X, b)     
		lambda = lambda_cpp(p.d[d,], XB)
		
		iJi = rbinom_mat((delta * lambda) / (delta * lambda + 1))
		while (sum(iJi) == 0) {
			iJi = rbinom_mat((delta * lambda) / (delta * lambda + 1))
		}
		LambdaiJi = lambdaiJi(p.d[d,], XB, iJi)
		Time.inc = vapply(LambdaiJi, function(lambda) {
			rexp(1, lambda)
			}, c(1))
		i.d = which(Time.inc == min(Time.inc[!is.na(Time.inc)]))
		j.d = which(iJi[i.d,] == 1)
		t.d = t.d + Time.inc[i.d]
		edge[[length(base.edge)+ d]] = list(sender = i.d, receiver = j.d, timestamp = t.d)	
	}
	options(warn = 0)
	return(list(edge = edge[-(1:length(base.edge))], text = text, C = currentC, beta = b, delta = delta))
}

GenerateDocs_backward_init = function(nDocs, node, vocabulary, nIP, K, nwords, alpha, mvec, betas, nvec, b, delta, currentC, 
									  netstat, base.edge, base.text, seed) {

	set.seed(seed)
 	W = length(vocabulary)
  	phi = lapply(1L:K, function(k) {
		rdirichlet_cpp(1, betas * nvec)
	})
	L = 3
    P = 1 + L * (2 * ("dyadic" %in% netstat) + 4 * ("triadic" %in% netstat) + 2 *("degree" %in% netstat))

	t.d = 384
	edge = base.edge
	text = base.text
	p.d = matrix(NA, nrow = nDocs, ncol = nIP)
	
	options(warn = -1)
	for (d in 1:nDocs) {
		N.d = nwords
		theta.d = rdirichlet_cpp(1, alpha * mvec)	
		topic.d = multinom_vec(max(1, N.d), theta.d)
		text[[length(base.text) + d]] = rep(NA, N.d)
		if (N.d > 0) {
			for (n in 1:N.d){
				text[[length(base.text) + d]][n] = multinom_vec(1, phi[[topic.d[n]]])
				}
			names(text[[length(base.text) + d]]) = topic.d
		}

		p.d[d, ] = vapply(1L:nIP, function(IP) {
			sum(topic.d %in% which(currentC == IP))
			}, c(1)) / max(1, N.d)
		history.t = History(edge, p.d, node, t.d + 10^(-10))
		X = lapply(node, function(i) {
  	        Netstats(history.t, node, i, netstat)
            })
   	 	XB = MultiplyXBList(X, b)     
		lambda = lambda_cpp(p.d[d,], XB)
		
		iJi = rbinom_mat((delta * lambda) / (delta * lambda + 1))
		while (sum(iJi) == 0) {
			iJi = rbinom_mat((delta * lambda) / (delta * lambda + 1))
		}
		LambdaiJi = lambdaiJi(p.d[d,], XB, iJi)
		Time.inc = vapply(LambdaiJi, function(lambda) {
			rexp(1, lambda)
			}, c(1))
		i.d = which(Time.inc == min(Time.inc[!is.na(Time.inc)]))
		j.d = which(iJi[i.d,] == 1)
		t.d = t.d + Time.inc[i.d]
		edge[[length(base.edge)+ d]] = list(sender = i.d, receiver = j.d, timestamp = t.d)	
	}
	options(warn = 0)
	return(list(edge = edge, text = text, C = currentC, beta = b, delta = delta))
}

Inference = function(edge, node, textlist, vocabulary, nIP, K, sigma_Q, alpha, mvec, betas, nvec, prior.b.mean, prior.b.var, prior.delta, 
					out, n1, n2, n3, burn, thin, netstat, seed, plot = FALSE) {
  
  set.seed(seed)
  # trim the edge so that we only model edges after 16 days
	timeinc = c(as.numeric(edge[[1]][3]), vapply(seq(along = edge)[-1], function(d) {
  	as.numeric(edge[[d]][3]) - as.numeric(edge[[d-1]][3])
 	}, c(1)))

    edge2 = which_int(384, cumsum(timeinc)) : length(edge)
   
  # initialize alpha, mvec, delta, nvec, eta, lvec, and gammas
 	W = length(vocabulary)
  	phi = lapply(1L:K, function(k) {
		rdirichlet_cpp(1, betas * nvec)
	})
	delta = rbeta(1, prior.delta[1], prior.delta[2])
	deltamat = delta
	eta = qnorm(delta)
	
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
    P = 1 + L * (2 * ("dyadic" %in% netstat) + 4 * ("triadic" %in% netstat) + 2 *("degree" %in% netstat))
    bmat = list()
	for (IP in 1L:nIP) {
		bmat[[IP]] = matrix(0, nrow = P, ncol = (n3 - burn) / thin)
		bmat[[IP]][, 1:(n3 - burn) / thin] = prior.b.mean
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
    for (o in 1L:out) {
      Beta.old = lapply(bmat, function(b) {
      			rowMeans(b)
         		})  
    
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
		iJi[[d]][as.numeric(edge[[d]][1]),] = tabulateC(as.numeric(unlist(edge[[d]][2])), length(node))
		LambdaiJi[[d]] = lambdaiJi(p.d[d,], XB, iJi[[d]])
		nonemptyiJi[[d]] = LambdaiJi[[d]][!is.na(LambdaiJi[[d]])]
		observediJi[[d]] = LambdaiJi[[d]][as.numeric(edge[[d]][1])]
		}	 
    
    textlist.raw = unlist(textlist)
    table.W = lapply(1L:K, function(k) {
      tabulateC(textlist.raw[which(unlist(currentZ) == k)], W)
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
      				  tabulateC(textlist.raw[which(unlist(currentZ) == k)], W)
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
            EdgeInEqZ(iJi[[d]], lambda[[d]], delta) + TimeInEqZ(nonemptyiJi[[d]], timeinc[d]) + 
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
      logWZ.mat = c(logWZ.mat, logWZ(K, currentZ, textlist, table.W, alpha, mvec, betas, nvec))
      }
      
     for (d in edge2) {
     	 p.d[d, ] =  vapply(1L:nIP, function(IP) {
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
		  		 dmvnorm(Beta.old[[IP]], rep(0, P), sigma^2 * diag(P), log = TRUE)
		  		 }, c(1))) 
	 post.old1 = sum(vapply(edge2, function(d) {
	     	     EdgeInEqZ(iJi[[d]], lambda[[d]], delta) + TimeInEqZ(nonemptyiJi[[d]], timeinc[d]) +
    			 	 ObservedInEqZ(observediJi[[d]])
    			 }, c(1))) / length(edge2)
    	 
	 options(warn = -1)
	 proposal.var = list()
     for (IP in 1L:nIP) {
     	 proposal.var[[IP]] = cor(t(bmat[[IP]]))
         proposal.var[[IP]][is.na(proposal.var[[IP]])] = 0
         if (sum(eigen(proposal.var[[IP]])$values < 0 ) > 0) {
        	 proposal.var[[IP]] = diag(P)
         }
         }
     proposal.var = lapply(1:nIP, function(IP){diag(P)})
     options(warn = 0)
     for (i3 in 1L:n3) {
     	Beta.new = lapply(1L:nIP, function(IP) {
           rmvnorm(1, Beta.old[[IP]], sigma_Q * proposal.var[[IP]])
         }) 
        for (d in edge2) {
           XB = MultiplyXBList(X, Beta.new)    
           LambdaiJi[[d]] = lambdaiJi(p.d[d,], XB, iJi[[d]])
		   nonemptyiJi[[d]] = LambdaiJi[[d]][!is.na(LambdaiJi[[d]])]
		   observediJi[[d]] = LambdaiJi[[d]][as.numeric(edge[[d]][1])]
         }
         prior.new1 = sum(vapply(1L:nIP, function(IP) {
        		dmvnorm(Beta.new[[IP]], rep(0, P), sigma^2 * diag(P), log = TRUE)
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
         
        #update delta for every outer iteration
        prior.old2 = dbeta(delta, 1, 10) 
        post.old2 = sum(vapply(edge2, function(d) {
    			    EdgeInEqZ(iJi[[d]], lambda[[d]], delta)
    			    }, c(1))) / length(edge2)
    			    
    	eta.new = rnorm(1, eta, sigma_Q)
    	delta.new = pnorm(eta.new)
		prior.new2 = dbeta(delta.new, 1, 10)
		post.new2 = sum(vapply(edge2, function(d) {
    			    EdgeInEqZ(iJi[[d]], lambda[[d]], delta.new)
    			    }, c(1))) / length(edge2)
    	 	loglike.diff2 = prior.new2 + post.new2 - prior.old2 - post.old2
    	     	if (log(runif(1, 0, 1)) < loglike.diff2) {
			delta = delta.new
			eta = eta.new
			}
			deltamat = c(deltamat, delta)
	 }
    
    if (plot) {
    burnin = round(out / 10)
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
     
  chain.final = list(C = currentC, Z = lapply(edge2, function(d) {currentZ[[d]]}), B = bmat, D = deltamat)

  return(chain.final)
}

GenerateDocs_backward = function(nDocs, node, vocabulary, nIP, K, nwords, alpha, mvec, betas, nvec,
							    b, delta, currentC, netstat, topic_token_assignments, topic_token_counts, word_type_topic_counts, 
							    base.edge, base.text, seed) {

	set.seed(seed)
 	W = length(vocabulary)
	L = 3
    P = 1 + L * (2 * ("dyadic" %in% netstat) + 4 * ("triadic" %in% netstat) + 2 *("degree" %in% netstat))

	t.d = 384
	edge = base.edge
	text = base.text
	p.d = matrix(NA, nrow = nDocs, ncol = nIP)
	phi_k = rep(NA, K)
	options(warn = -1)
	for (d in 1:nDocs) {
		N.d = nwords
		text[[length(base.text) + d]] = rep(NA, N.d)
		topic.d = topic_token_assignments[[d]]
		word.d = as.numeric(names(topic_token_assignments[[d]]))
		if (N.d > 0) {
			for (n in 1:N.d){
				word_type_topic_counts[word.d [n], topic.d[n]] = word_type_topic_counts[word.d [n], topic.d[n]] - 1
				for (w in 1:W) {
					phi_k[w] = (word_type_topic_counts[w, topic.d[n]] + betas * nvec[w])/(topic_token_counts[topic.d[n]] + betas)
				} 
				text[[length(base.text) + d]][n] = multinom_vec(1, phi_k)
				word_type_topic_counts[text[[d]][n], topic.d[n]] = word_type_topic_counts[text[[d]][n], topic.d[n]] + 1
				}
			names(text[[length(base.text) + d]]) = topic.d
		}

		p.d[d, ] = vapply(1L:nIP, function(IP) {
			sum(topic.d %in% which(currentC == IP))
			}, c(1)) / max(1, N.d)
		history.t = History(edge, p.d, node, t.d + 10^(-10))
		X = lapply(node, function(i) {
  	        Netstats(history.t, node, i, netstat)
            })
   	 	XB = MultiplyXBList(X, b)     
		lambda = lambda_cpp(p.d[d,], XB)
		
		iJi = rbinom_mat((delta * lambda) / (delta * lambda + 1))
		while (sum(iJi) == 0) {
			iJi = rbinom_mat((delta * lambda) / (delta * lambda + 1))
		}
		LambdaiJi = lambdaiJi(p.d[d,], XB, iJi)
		Time.inc = vapply(LambdaiJi, function(lambda) {
			rexp(1, lambda)
			}, c(1))
		i.d = which(Time.inc == min(Time.inc[!is.na(Time.inc)]))
		j.d = which(iJi[i.d,] == 1)
		t.d = t.d + Time.inc[i.d]
		edge[[length(base.edge)+ d]] = list(sender = i.d, receiver = j.d, timestamp = t.d)	
	}
	options(warn = 0)
	return(list(edge = edge, text = text, C = currentC, beta = b, delta = delta))
}

 
#PP_Plots
GiR_PP_Plots = function(Forward_stats, Backward_stats) {
	nms = colnames(Forward_stats)
	
	for (i in 1L:ncol(Forward_stats)) {
		if (nrow(Forward_stats) > 20000) {
			thin = seq(from = floor(nrow(Forward_stats)/10), to = nrow(Forward_stats), length.out = 10000)
			Forward_test = Forward_stats[thin, i]
			Backward_test = Backward_stats[thin, i]
		} else {
			Forward_test = Forward_stats
			Backward_test = Backward_stats
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
			   cex.lab = 2,
			   cex.axis = 1.4,
			   cex.main = 2)
		lines(x = xlims, y = ylims, col = "red", lwd = 3)
		# text(paste("Backward Mean:", round(mean(Backward_stats[,i]), 4),
				   # "\nForward Mean:", round(mean(Forward_stats[,i]), 4),
				   # "\nt-test p-value:", round(t.test(Backward_test, Forward_test)$p.value, 4),
				   # "\nMann-Whitney p-value:", round(wilcox.test(Backward_test, Forward_test)$p.value,4)),
				   # x = xlims[2] - 0.3 * abs(xlims[2] - xlims[1]),
				   # y = ylims[1] + 0.2 * abs(ylims[2] - ylims[1]),
				   # cex = 1.5)
	}
}      


# Getting_It_Right
GiR = function(Nsamp = 100, nDocs = 10, node = 1:4, vocabulary =  c("hi", "hello","bye", "mine", "what"), 
			   nIP = 2, K = 4, nwords = 4, alpha = 2, mvec = rep(1/4, 4), betas = 2, nvec = rep(1/5, 5), 
			   prior.b.mean = c(-3, rep(0, 6)), prior.b.var = diag(7), prior.delta = c(1, 1), sigma_Q = 0.5, 
			   niters = c(10, 1, 1, 1100, 100, 10), netstat = "dyadic", generate_PP_plots = TRUE, seed = 123) {
	
	set.seed(seed)
	base.data = GenerateDocs_base(nDocs, node, vocabulary, nIP, K, nwords, alpha, mvec, betas, nvec, prior.b.mean, prior.b.var, prior.delta,
				netstat, seed)
	base.edge = base.data$edge	   
	base.text = base.data$text


	#Forward sampling
	Forward_stats = matrix(NA, nrow = Nsamp, ncol = 21)
	for (i in 1L:Nsamp) {
		set.seed(i)
		W = length(vocabulary)
		L = 3
    		P = 1 + L * (2 * ("dyadic" %in% netstat) + 4 * ("triadic" %in% netstat) + 2 *("degree" %in% netstat))
		b = lapply(1L:nIP, function(IP) {
			c(rmvnorm(1, prior.b.mean, prior.b.var))
			})
		delta = rbeta(1, prior.delta[1], prior.delta[2])
		currentC = sample(1L:nIP, K, replace = TRUE)
		Forward_sample = GenerateDocs_forward(nDocs, node, vocabulary, nIP, K, nwords, alpha, mvec, betas, nvec, 
						 b, delta, currentC, netstat, base.edge, seed) 
		Forward_stats[i, 1:P] = Reduce('+', b) / nIP
		Forward_stats[i, P+1] = mean(vapply(1L:nDocs, function(d) {length(Forward_sample$edge[[d]][[2]])}, c(1)))
		Forward_stats[i, P+2] = mean(vapply(2:nDocs, function(d) {Forward_sample$edge[[d]][[3]]-Forward_sample$edge[[d-1]][[3]]}, c(1))) 			
		Forward_stats[i, P+3] = mean(currentC)
		Tokens_in_Topic = tabulate(vapply(1:nDocs, function(d){as.numeric(names(Forward_sample$text[[d]]))}, rep(0, nwords)), K)
		Forward_stats[i, (P+4):(P+3+nIP)] = vapply(1:nIP, function(IP) {Tokens_in_Topic %*% (currentC == IP)}, c(1))
		Forward_stats[i, (P+4+nIP):(P+3+nIP+K)] = Tokens_in_Topic
		Tokens_in_Word = tabulate(vapply(1:nDocs, function(d){Forward_sample$text[[d]]}, rep(0, nwords)), W)
		Forward_stats[i, (P+4+nIP+K):21] = Tokens_in_Word
		colnames(Forward_stats) = c(paste0("B_",1:P), "Mean_receipients", "Mean_timediff", "Mean_TopicIP", 
									paste0("Tokens_in_IP_", 1:nIP), paste0("Tokens_in_Topic", 1:K), paste0("Tokens_in_Word", 1:W))
	}
	
	#Backward sampling
	Backward_stats = matrix(NA, nrow = Nsamp, ncol = 21)
	Backward_sample = GenerateDocs_backward_init(nDocs, node, vocabulary, nIP, K, nwords, alpha, mvec, betas, nvec, 
					  b, delta, currentC, netstat, base.edge, base.text, seed) 
	for (i in 1L:Nsamp) { print(i)
		seed = seed + 100
		Inference_samp = Inference(Backward_sample$edge, node, Backward_sample$text, vocabulary, nIP, K, sigma_Q, 
						  alpha, mvec, betas, nvec, prior.b.mean, prior.b.var, prior.delta,
						  out = niters[1], n1 = niters[2], n2 = niters[3], n3 = niters[4], burn = niters[5], thin = niters[6], 
						  netstat, seed)
		b = lapply(1L:nIP, function(IP) {
			rowMeans(Inference_samp$B[[IP]])
		})
		delta = mean(Inference_samp$D)
		currentC = Inference_samp$C
		topic_token_assignments = Inference_samp$Z
		for (d in 1L:length(topic_token_assignments)) {
			names(topic_token_assignments[[d]]) = Backward_sample$text[[d + length(base.text)]]
		}
		topic_token_counts = tabulate(unlist(Inference_samp$Z), K)
		word_type_topic_counts = matrix(NA , W, K)
		for (w in 1L:W) {
			word_type_w = which(unlist(Backward_sample$text[-(1:length(base.text))]) == w)
			word_type_topic_counts[w, ] = tabulate(unlist(Inference_samp$Z)[word_type_w], K)
		}
		
		Backward_sample = GenerateDocs_backward(nDocs, node, vocabulary, nIP, K, nwords, alpha, mvec, betas, nvec, 
						 b, delta, currentC, netstat, topic_token_assignments, topic_token_counts, word_type_topic_counts,
						 base.edge, base.text, seed)
		Backward_stats[i, 1:P] = Reduce('+', b) / nIP
		Backward_stats[i, P+1] = mean(vapply(1L:nDocs + length(base.text), function(d) {length(Backward_sample$edge[[d]][[2]])}, c(1)))
		Backward_stats[i, P+2] = mean(vapply(2:nDocs + length(base.text), function(d) {Backward_sample$edge[[d]][[3]]-Backward_sample$edge[[d-1]][[3]]}, c(1))) 			
		Backward_stats[i, P+3] = mean(currentC)
		Tokens_in_Topic = tabulate(vapply(1:nDocs + length(base.text), function(d){as.numeric(names(Backward_sample$text[[d]]))}, rep(0, nwords)), K)
		Backward_stats[i, (P+4):(P+3+nIP)] = vapply(1:nIP, function(IP) {Tokens_in_Topic %*% (currentC == IP)}, c(1))
		Backward_stats[i, (P+4+nIP):(P+3+nIP+K)] = Tokens_in_Topic
		Tokens_in_Word = tabulate(vapply(1:nDocs + length(base.text), function(d){Backward_sample$text[[d]]}, rep(0, nwords)), W)
		Backward_stats[i, (P+4+nIP+K):21] = Tokens_in_Word				 
	}
	
	tstats = rep(0, 21)
	wstats = rep(0, 21)
	for (j in 1L:21) {
		thin = seq(from = floor(Nsamp / 5), to = Nsamp, length.out = 400)
		Forward_test = Forward_stats[thin, j]
		Backward_test = Backward_stats[thin, j]
		tstats[j] = t.test(Backward_test, Forward_test)$p.value
		wstats[j] = wilcox.test(Backward_test, Forward_test)$p.value
	}
	names(tstats) = names(wstats) = colnames(Forward_stats)						
	if (generate_PP_plots) {
		par(mfrow=c(3,7), oma = c(3,3,3,3), mar = c(5,5,4,1))
		GiR_PP_Plots(Forward_stats, Forward_stats)
	}			
	return(list(Forward = Forward_stats, Backward = Backward_stats, tstats = tstats, wstats = wstats))	
}