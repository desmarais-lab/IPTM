library(MLmetrics)
library(pROC)
load('/Users/bomin8319/Desktop/IPTM/enron/Enron.RData')
set.seed(1)
selectD = sample(1963:3925, 200, replace = FALSE)
setwd('/Users/bomin8319/Desktop/IPTM/Enron_PPE')
commondocs = selectD
for (nIP in 1:3) {
    for (K in c(2, 5, 10, 25, 50, 75, 100)){
        filename = paste0("_",nIP,"_",K,".RData")
        docs = list.files(pattern = filename)
        docnames = as.numeric(sapply(docs, function(x){strsplit(x, "\\_")[[1]][2]}))
            commondocs = intersect(commondocs, docnames)
        }
    }
    
tables = list()
d = 1
for (D in commondocs) {
  observeds = tabulate(Enron$edge[[D]]$sender, 33)
  observedr = tabulate(Enron$edge[[D]]$receiver, 33)
  observedt = Enron$edge[[D]]$timestamp
  tables[[d]] = matrix(0, nrow = 21, ncol = 3)
  rownames(tables[[d]]) = 1:21
  iter = 1
  for (nIP in 1:3) {
    for (K in c(2, 5, 10, 25, 50, 75, 100)){
    filename = paste0("PPE_",D,"_",nIP,"_",K,".RData")
    load(filename)
    rownames(tables[[d]])[iter] = filename
    MSE = 0
    F1score = 0
    AUC = 0
    it = 1
    for (o in 1:10) {
      for (r in 1:50) {
          if (r > 10 && (r-10) %% 4 == 0) {
        AUC = AUC + as.numeric(multiclass.roc(observeds, tabulate(PPE[[o]][[r]]$edge$sender, 33))$auc)
        F1 =  F1_Score(observedr, tabulate(PPE[[o]][[r]]$edge$receiver, 33))
        if (is.na(F1)) {F1 = 0}
        F1score = F1score + F1
        MSE = MSE + (PPE[[o]][[r]]$edge$timestamp - observedt)^2
          it = it + 1
          }
      }
    }
    tables[[d]][iter, 1] = AUC / (it - 1)
    tables[[d]][iter, 2] = F1score / (it - 1)
    tables[[d]][iter, 3] = MSE / (it - 1)
  iter = iter + 1
  }
  }
  d = d + 1
}

all =  Reduce('+', tables) / length(tables)

22 - rank(all[,1]) + 22 - rank(all[,2]) + rank(all[,3])

