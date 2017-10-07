library(MLmetrics)
library(pROC)
set.seed(1)
selectD = sample(1963:3925, 200, replace = FALSE)
setwd('/Users/bomin8319/Desktop/IPTM/Dare_PPE2')
##################################################

load('/Users/bomin8319/Desktop/IPTM/enron/Dare.RData')
set.seed(1)
selectD = sample(1105:2210, 200, replace = FALSE)
commondocs = selectD
for (nIP in 1:3) {
    for (K in c(2, 5, 10, 25, 50, 75, 100)){
        filename = paste0("_",nIP,"_",K,".RData")
        docs = list.files(pattern = filename)
        docnames = as.numeric(sapply(docs, function(x){strsplit(x, "\\_")[[1]][2]}))
        commondocs = intersect(commondocs, docnames)
        }
    }

 tables = matrix(0, nrow = length(commondocs), ncol = 3)
  colnames(tables) = c("sender", "receiver", "time")
  rownames(tables) = commondocs

d = 1
for (D in commondocs) {
  observeds = tabulate(Dare$edge[[D]]$sender, 27)
  observedr = tabulate(Dare$edge[[D]]$receiver, 27)
  observedt = Dare$edge[[D]]$timeinc
   filename = paste0("DareSimulatedDocumentsD",D,".RData")
  load(filename)
    MSE = c()
    F1score = 0
    AUC = 0
    for (o in 1:100) {
        AUC = AUC + as.numeric(multiclass.roc(observeds, tabulate(simulated_document_matrix[o,1], 27))$auc)
        F1score = F1score + F1_Score(observedr, tabulate(which(simulated_document_matrix[o,2:28] == 1), 27))
        MSE[o] = abs(simulated_document_matrix[o,29]- observedt)
    }
    tables[d, 1] = AUC / 100
    tables[d, 2] = F1score / 100
    tables[d, 3] = median(MSE)
    d = d+1
  }


table_sum = matrix(rep(colMeans(tables), 7), 3)
rownames(table_sum) =  c("sender", "receiver", "time")
colnames(table_sum) = c("2", "5", "10", "25", "50", "75", "100")


tables2 = list()
d = 1
for (D in commondocs) {
  observeds = tabulate(Dare$edge[[D]]$sender, 27)
  observedr = tabulate(Dare$edge[[D]]$receiver, 27)
  observedt = Dare$edge[[D]]$unixtime
  tables2[[d]] = matrix(0, nrow = 21, ncol = 3)
  rownames(tables2[[d]]) = 1:21
  iter = 1
  for (nIP in 1:3) {
    for (K in c(2, 5, 10, 25, 50, 75, 100)){
    filename = paste0("PPE_",D,"_",nIP,"_",K,".RData")
    load(filename)
    rownames(tables2[[d]])[iter] = filename
    MSE =c()
    F1score = 0
    AUC = 0
    it = 1
    for (o in 1:10) {
      for (r in 1:50) {
          if (r > 10 && (r-10) %% 4 == 0) {
        AUC = AUC + as.numeric(multiclass.roc(observeds, tabulate(PPE[[o]][[r]]$edge$sender, 27))$auc)
        F1 = F1_Score(observedr, tabulate(PPE[[o]][[r]]$edge$receiver, 27))
        if (is.na(F1)) {F1 = 0}
        F1score = F1score + F1
        MSE[it] = abs(PPE[[o]][[r]]$edge$timestamp - observedt)
          it = it + 1
          }
      }
    }
    tables2[[d]][iter, 1] = AUC / (it - 1)
    tables2[[d]][iter, 2] = F1score / (it - 1)
    tables2[[d]][iter, 3] = median(MSE)
  iter = iter + 1
  }
  }
  d = d + 1
}

all = Reduce('+', tables2) / length(tables2)

par(mfrow = c(1,3))
plot(colnames(table_sum), table_sum[1,], type = 'b', pch = 19, col = 'red', xlab = "Number of Topics", lty = 2, ylab = "Sender AUC", ylim = c(0.49, 0.55))
lines(colnames(table_sum), all[1:7,1], lty = 3,col = "blue", type = 'b', pch = 19)
lines(colnames(table_sum), all[8:14,1], lty = 4,col = "green",type = 'b', pch = 19)
lines(colnames(table_sum), all[15:21,1], lty = 5,col = "purple",type = 'b', pch = 19)
legend(60, 0.55, pch = 19, bty = "n", lty = 2:5, legend = c("IPTM with C = 3","IPTM with C = 2","IPTM with C = 1", "Regressions"), col = c("purple", "green", "blue", "red"))

plot(colnames(table_sum), table_sum[2,], type = 'b', pch = 19,  col = 'red', xlab = "Number of Topics", lty = 2, ylab = "Receiver F1Score", ylim=c(0.6, 1.1))
lines(colnames(table_sum), all[1:7,2], lty = 3,col = "blue",type = 'b', pch = 19)
lines(colnames(table_sum), all[8:14,2], lty = 4,col = "green",type = 'b', pch = 19)
lines(colnames(table_sum), all[15:21,2], lty = 5,col = "purple",type = 'b', pch = 19)
legend(60, 1.1,  pch = 19,bty = "n", lty = 2:5, legend = c("IPTM with C = 3","IPTM with C = 2","IPTM with C = 1", "Regressions"), col = c("purple", "green", "blue", "red"))

plot(colnames(table_sum), table_sum[3,], type = 'b', pch = 19, col = 'red', xlab = "Number of Topics", lty = 2, ylab = "Time MAE", ylim = c(1.0, 1.25))
lines(colnames(table_sum), all[1:7,3], lty = 3,col = "blue",type = 'b', pch = 19)
lines(colnames(table_sum), all[8:14,3], lty = 4,col = "green",type = 'b', pch = 19)
lines(colnames(table_sum), all[15:21,3], lty = 5,col = "purple",type = 'b', pch = 19)
legend(60, 1.25,  pch = 19, bty = "n",lty = 2:5, legend = c("IPTM with C = 3","IPTM with C = 2","IPTM with C = 1", "Regressions"), col = c("purple", "green", "blue", "red"))


