#' @useDynLib IPTM
#' @importFrom Rcpp sourceCpp
#' @importFrom mvtnorm dmvnorm rmvnorm
#' @importFrom MCMCpack rdirichlet
NULL

#' @title history
#' @description calculate the time difference from previous interactions to certain time 'beforetime'
#'
#' @param edge
#' @param node
#' @param when
#'
#' @return
#'
#' @export
history = function(edge, node, when) {
  histlist = lapply(node, function(v) {
    lapply(node, function(u) {
      numeric(0)
    })
  })
  edge2 = matrix(edge[edge[, 3] < when, ], byrow = FALSE, ncol = 3)
  if (nrow(edge2) > 0) {
    for (d in 1:nrow(edge2)) {
      histlist[[edge2[d,1]]][[edge2[d,2]]] = c(histlist[[edge2[d, 1]]][[edge2[d, 2]]],
        when - (edge2[d, 3]))
    }
  }
  return(histlist)
}

#' @title parupdate
#' @description parameter optimization of alpha and mvec at the same time -> need to be faster
#'
#' @param nIP
#' @param K
#' @param currentC
#' @param currentZ
#' @param alpha
#' @param mvec
#'
#' @return
#'
#' @export
parupdate =  function(nIP, K, currentC, currentZ, alpha, mvec) {
  vec = alpha*mvec
  zcnew = list()
  clist = c()
  corpusCnew = sortedZ(nIP, currentC, currentZ)
  for (IP in 1:nIP) {
    zcnew[[IP]] = tabulate(unlist(corpusCnew[[IP]]), nbins = K)
    clist[IP] = sum(zcnew[[IP]])
  }
  iter = 1
  ctable = rep(0, max(clist))
  ctable[clist] = tabulate(currentC)
  cklist = list()
  cktable = list()
  for (k in 1:K) {
    cklist[[k]] = vapply(1:nIP, function(IP) {zcnew[[IP]][k]}, c(1))
    cktable[[k]] = rep(0, max(cklist[[k]]))
    cktable[[k]][cklist[[k]]] = tabulate(currentC)
  }

  while ((abs(alpha - sum(vec)) > 0.001) | (iter == 1)) {
    alpha = sum(vec)
    S = Supdate(alpha, ctable)
    s = Skupdate(vec, cktable)
    vec = vec * s/S
    iter=iter + 1
  }
  return(vec)
}

#' @title selected
#' @description find out the chosen category from the multinomial distribution
#'
#' @param samples
#' @param proportions
#'
#' @return
#'
#' @export
selected = function(samples, proportions) {
  chosen  = rmultinom(samples, 1, proportions)
  out = vapply(1:samples, function(s) {which(chosen[, s] == TRUE)}, c(1))
  return(out)
}

#' @title betapartC
#' @description calculate the log of beta part used to obtain the constants for multinomial sampling of IP
#'
#' @param nIP
#' @param lambdai
#' @param specificedge
#'
#' @return
#'
#' @export
betapartC = function(nIP, lambdai, specificedge) {
  const = rep(NA, nIP)
  for (IP in 1:nIP) {
    edge = matrix(specificedge, ncol = 3)
    receiver = ifelse(edge[, 1] > edge[, 2], edge[, 2],  edge[, 2] - 1)
    const[IP] = lambdai[[IP]][receiver] - log(sum(exp(lambdai[[IP]])))
  }
  return(const)
}

#' @title topicpartC
#' @description calculate the log of topic-IP part used to obtain the constants for multinomial sampling of IP
#'
#' @param nIP
#' @param currentZ
#' @param zc
#' @param alpha
#' @param mvec
#' @param document
#'
#' @return
#'
#' @export
topicpartC = function(nIP, currentZ, zc, alpha, mvec, document) {
  const = rep(NA, nIP)
  for (IP in 1:nIP) {
    topics = currentZ[[document]][IP, ]
    num = sum(log(zc[[IP]][topics] -
                    ifelse(zc[[IP]][topics] > 0, 1, 0) +
                    alpha * mvec[topics]))
    denom = length(topics) * log(sum(zc[[IP]]) - sum(ifelse(zc[[IP]] > 0, 1, 0)) + alpha)
    const[IP] = num - denom
  }
  return(const)
}

#' @title topicpartZ
#' @description calculate the log of topic-IP part used to obtain the constants for multinomial sampling of K
#'
#' @param currentIP
#' @param zcnew
#' @param alpha
#' @param mvec
#'
#' @return
#'
#' @export
topicpartZ = function(currentIP, zcnew, alpha, mvec) {
  const = log(zcnew[[currentIP]] -
                ifelse(zcnew[[currentIP]] > 0, 1, 0) + alpha * mvec)
  return(const)
}

#' @title wordpartZ
#' @description: calculate the log of word-topic part used to obtain the constants for multinomial sampling of K
#'
#' @param K
#' @param textlistd
#' @param tableW
#' @param delta
#' @param nvec
#'
#' @return
#'
#' @export
wordpartZ = function(K, textlistd, tableW, delta, nvec) {
  const = matrix(NA, nrow=length(textlistd), ncol = K)
  for (k in 1:K) {
    num  = log(tableW[[k]][textlistd] -
                 ifelse(tableW[[k]][textlistd] > 0, 1, 0) + delta * nvec[textlistd])
    denom  = log(sum(tableW[[k]]) - sum(ifelse(tableW[[k]] > 0, 1, 0)) + delta)
    const[, k] = num - denom
  }
  return(const)
}

#' @title betapartB
#' @description calculate the log of beta part used to obtain the constants for multinomial sampling of Beta
#'
#' @param nIP
#' @param lambdai
#' @param edgeC
#'
#' @return
#'
#' @export
betapartB = function(nIP, lambdai, edgeC) {
  const = rep(NA, nIP)
  for (IP in 1:nIP) {
    edge = matrix(edgeC[[IP]], ncol = 3)
    receiver = ifelse(edge[, 1] > edge[, 2], edge[, 2],  edge[, 2] - 1)
    const[IP] = sum(vapply(1:length(receiver), function(r) {
      lambdai[[IP]][r, receiver[r]] - log(sum(exp(lambdai[[IP]][r, ])))
    }, c(1)))
  }
  return(const)
}

#' @title MCMC
#' @description
#'
#' @param edge
#' @param node
#' @param textlist
#' @param vocabulary
#' @param nIP
#' @param K
#' @param delta_B
#' @param outer
#' @param n1
#' @param n2
#' @param n3
#' @param burn
#' @param thin
#' @param opt
#' @param seed
#'
#' @return
#'
#' @export
MCMC = function(edge, node, textlist, vocabulary, nIP, K, delta_B, outer = 200,
  n1 = 3, n2 = 3, n3 = 3300, burn = 300, thin = 3, opt = TRUE, seed = 1) {

  set.seed(seed)

  # initialize alpha, mvec, delta, nvec, eta, lvec, and gammas
  alpha =  50 / K
  mvec = rep(1, K) / K
  delta = 0.01
  W = length(vocabulary)
  nvec = rep(1, W) / W
  eta = 50 / nIP
  lvec = rep(1, nIP) / nIP
  gammas = rdirichlet(1, eta * lvec)

  # initialize theta
  theta = list()
  for (IP in 1:nIP) {
    theta[[IP]] = rdirichlet(1, alpha * mvec)
  }

  # initialize C and Z
  currentC = selected(nrow(edge), gammas)
  currentZ = lapply(1:nrow(edge), function(d) {
    matrix(vapply(1:nIP, function(IP) {
      selected(length(textlist[[d]]), theta[[IP]])
    }, rep(0, length(textlist[[d]]))), byrow = TRUE,nrow = nIP)
  })
  # initialize beta
  P = 4
  sigma = 1
  bmat = list()
  for (IP in 1:nIP) {
    bmat[[IP]] = matrix(NA, nrow = P, ncol = (n3 - burn) / thin)
    bmat[[IP]][,] = rmvnorm(1, rep(0, P), sigma^2 * diag(P))
  }

  logWZmat = c()

  if (opt == TRUE) {
    alphamat = alpha
    mvecmat = mvec
  }

  #start outer iteration
  for(o in 1:outer) {
    cat("outer iteration = ", o, "\n")

    if (opt == TRUE) {
      op = options(warn = (-1))
      vec = parupdate(nIP, K, currentC, currentZ, alpha, mvec)
      alpha = sum(vec)
      alphamat = c(alphamat, alpha)
      mvec = vec / alpha
      mvecmat = rbind(mvecmat, mvec)
      options(op)
    }

    cat("inner iteration 1", "\n")
    lambdai = list()
    corpusC  = sortedZ(nIP, currentC, currentZ)
    zc =  lapply(1:nIP, function(IP) {
      tabulate(unlist(corpusC[[IP]]), nbins = K)
    })
    for (i1 in 1:n1) {
      # C update given Z and B - within each document d
      for (d in 1:nrow(edge)) {
        currentC2 = currentC[-d]
        edgeC = lapply(1:nIP, function(IP) {
          edge[which(currentC2 == IP), ]
        })
        zc[[currentC[d]]] = zc[[currentC[d]]] -
          tabulate(currentZ[[d]][currentC[d], ], nbins = K)
        for (IP in 1:nIP) {
          histlist = history(edgeC[[IP]], node, edge[d, 3])
          allxmatlist = allxmat(edgeC[[IP]], node, histlist, edge[d, 1], 0.05)
          lambdai[[IP]] = multiplyXB(allxmatlist, bmat[[IP]][, (n3 - burn) / thin])
        }
        const = log(gammas)+betapartC(nIP, lambdai, edge[d, ]) +
          topicpartC(nIP, currentZ, zc, alpha, mvec, d)
        currentC[d] = selected(1, exp(const))
        zc[[currentC[d]]] = zc[[currentC[d]]] +
          tabulate(currentZ[[d]][currentC[d], ], nbins = K)
      }
    }

    cat("inner iteration 2", "\n")
    textlist2 = unlist(textlist)
    corpusCnew = sortedZ(nIP, currentC, currentZ)
    zcnew = lapply(1:nIP, function(IP) {
      tabulate(unlist(corpusCnew[[IP]]), nbins = K)
    })
    finalZlist = finalZ(currentC, currentZ)
    finalZlist2 = unlist(finalZlist)
    tableW = lapply(1:K, function(k) {
      tabulate(textlist2[which(finalZlist2 == k)], nbins = W)
    })

    for (i2 in 1:n2) {
      for (d in 1:nrow(edge)) {
        currentCd = currentC[d]
        textlistd = textlist[[d]]
        topicpartd = topicpartZ(currentCd, zcnew, alpha, mvec)
        wordpartd = wordpartZ(K, textlistd, tableW, delta, nvec)
        for (w in 1:nrow(wordpartd)) {
          const2 = topicpartd + wordpartd[w, ]
          zwold = currentZ[[d]][currentCd, w]
          zwnew = selected(1, exp(const2))
          if (zwnew != zwold) {
            zcnew[[currentCd]][zwold] = zcnew[[currentCd]][zwold] - 1
            zcnew[[currentCd]][zwnew] = zcnew[[currentCd]][zwnew] + 1
            tableW[[zwold]][textlistd[w]] = tableW[[zwold]][textlistd[w]] - 1
            tableW[[zwnew]][textlistd[w]] = tableW[[zwnew]][textlistd[w]] + 1
            topicpartd = topicpartZ(currentCd, zcnew, alpha, mvec)
            wordpartd = wordpartZ(K, textlistd, tableW, delta, nvec)
            currentZ[[d]][currentCd, w] = zwnew
          }
        }
      }
    }

    logWZmat = c(logWZmat, logWZ(nIP, K, textlist2, W, zcnew, tableW,
      alpha, mvec, delta, nvec))

    cat("inner iteration 3", "\n")
    bold = list()
    allxmatlist2 = list()
    for (IP in 1:nIP) {
      bold[[IP]] = bmat[[IP]][, (n3 - burn) / thin]
      edgeC[[IP]] = edge[which(currentC == IP), ]
      allxmatlist2[[IP]] = list()

      for (d in 1:nrow(edgeC[[IP]])) {
        histlist2 = history(edgeC[[IP]], node, edgeC[[IP]][d, 3])
        allxmatlist2[[IP]][[d]] = allxmat(edgeC[[IP]], node, histlist2,
          edgeC[[IP]][d, 1], 0.05)
      }
    }

    lambdaiold = list()
    lambdainew = list()
    bnew = list()

    for (i3 in 1:n3) {
      for (IP in 1:nIP) {
        bnew[[IP]] = rmvnorm(1, bold[[IP]], (delta_B)^2 * diag(P))
        lambdaiold[[IP]] = t(vapply(1:nrow(edgeC[[IP]]), function(d) {
          multiplyXB(allxmatlist2[[IP]][[d]], bold[[IP]])
        }, rep(1, length(node) - 1)))
        lambdainew[[IP]] = t(vapply(1:nrow(edgeC[[IP]]), function(d) {
          multiplyXB(allxmatlist2[[IP]][[d]], bnew[[IP]])
        }, rep(1, length(node) - 1)))
      }

      prior = vapply(1:nIP, function(IP) {
        dmvnorm(bnew[[IP]], rep(0, P), sigma^2 * diag(P), log = TRUE) -
          dmvnorm(bold[[IP]], rep(0, P), sigma^2 * diag(P), log = TRUE)
      }, c(1))
      post = betapartB(nIP, lambdainew, edgeC) - betapartB(nIP, lambdaiold, edgeC)
      loglikediff = prior + post

      u = log(runif(nIP, 0, 1))
      for (IP in which((u < loglikediff) == TRUE)) {
        bold[[IP]]  = bnew[[IP]]
      }

      if (i3 > burn && i3 %% (thin) == 0) {
        for (IP in 1:nIP) {
          bmat[[IP]][ ,(i3 - burn) / thin] = bold[[IP]]
        }
      }
    }
  }

  if (opt == TRUE) {
    out = list(C = currentC, Z = finalZ(currentC, currentZ), B = bmat, L = logWZmat,
      A = alphamat, M = mvecmat)
  }
  else {
    out = list(C = currentC, Z = finalZ(currentC, currentZ), B = bmat, L = logWZmat)
  }
  return(out)
}

#' @title logWZ
#' @description calculate the log of unnormalized constant (corresponding to product of Eq.13 in Section 3.2) to check the convergence
#'
#' @param nIP
#' @param K
#' @param textlist2
#' @param W
#' @param zcnew
#' @param tableW
#' @param alpha
#' @param mvec
#' @param delta
#' @param nvec
#'
#' @return
#'
#' @export
logWZ = function(nIP, K, textlist2, W, zcnew, tableW, alpha, mvec, delta, nvec) {
  out1 = 0
  for (IP in 1:nIP) {
    num1 = zcnew[[IP]] + alpha * mvec
    denom1 = sum(zcnew[[IP]]) + alpha
    out1 = out1 + sum(lgamma(num1)) - lgamma(denom1)
  }
  out2 = 0
  for (k in 1:K) {
    num2 = tableW[[k]] + delta * nvec
    denom2 = sum(tableW[[k]]) + delta
    out2 = out2 + sum(lgamma(num2)) - lgamma(denom2)
  }
  return(out1 + out2)
}
