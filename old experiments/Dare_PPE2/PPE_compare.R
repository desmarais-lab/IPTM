library(MLmetrics)
library(pROC)
load('/Users/bomin8319/Desktop/IPTM/enron/Dare.RData')
set.seed(1)
selectD = sample(1105:2210, 200, replace = FALSE)
setwd('/Users/bomin8319/Desktop/IPTM/Dare_PPE2')
tables = list()
commondocs = selectD
for (nIP in 1:3) {
    for (K in c(2, 5, 10, 25, 50, 75, 100)){
        filename = paste0("_",nIP,"_",K,".RData")
        docs = list.files(pattern = filename)
        docnames = as.numeric(sapply(docs, function(x){strsplit(x, "\\_")[[1]][2]}))
        commondocs = intersect(commondocs, docnames)
    }
}

d = 1
for (D in commondocs) {
  observeds = tabulate(Dare$edge[[D]]$sender, 27)
  observedr = tabulate(Dare$edge[[D]]$receiver, 27)
  observedt = Dare$edge[[D]]$unixtime
  tables[[d]] = matrix(0, nrow = 21, ncol = 3)
  rownames(tables[[d]]) = 1:21
  iter = 1
  for (nIP in 1:3) {
    for (K in c(2, 5, 10, 25, 50, 75, 100)){
    filename = paste0("PPE_",D,"_",nIP,"_",K,".RData")
    load(filename)
    rownames(tables[[d]])[iter] = filename
    MSE = c()
    F1score = 0
    AUC = 0
    it = 1
    for (o in 1:10) {
      for (r in 1:50) {
          if (r > 10 && (r-10) %% 4 == 0) {
        AUC = AUC + as.numeric(multiclass.roc(observeds, tabulate(PPE[[o]][[r]]$edge$sender, 27))$auc)
        F1 =  F1_Score(observedr, tabulate(PPE[[o]][[r]]$edge$receiver, 27))
		if (is.na(F1)) {F1 = 0}
        F1score = F1score + F1       
        MSE[it] = abs(PPE[[o]][[r]]$edge$timestamp - observedt)
        it = it + 1
          }
    }
    }
    tables[[d]][iter, 1] = AUC / (it -1)
    tables[[d]][iter, 2] = F1score / (it - 1)
    tables[[d]][iter, 3] = median(MSE)
  iter = iter + 1
  }
  }
  d = d + 1
}

all =  Reduce('+', tables) / length(tables)

22 - rank(all[,1]) + 22 - rank(all[,2]) + rank(all[,3])

