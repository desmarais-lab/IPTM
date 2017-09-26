library(MLmetrics)
library(pROC)
load('/Users/bomin8319/Desktop/IPTM/enron/Enron.RData')
set.seed(1)
selectD = sample(1963:3925, 200, replace = FALSE)
setwd('/Users/bomin8319/Desktop/IPTM/paper/code/predictiveComparisonModels/results')
 tables = matrix(0, nrow = 200, ncol = 3)
  colnames(tables) = c("sender", "receiver", "time")
  rownames(tables) = selectD
for (d in 1:22) {
  D = selectD[d]
  observeds = tabulate(Enron$edge[[D]]$sender, 33)
  observedr = tabulate(Enron$edge[[D]]$receiver, 33)
  observedt = Enron$edge[[D]]$timestamp - Enron$edge[[D-1]]$timestamp
   filename = paste0("EnronSimulatedDocumentsD",D,".RData")
  load(filename)
    MSE = 0
    F1score = 0
    AUC = 0
    for (o in 1:100) {
        AUC = AUC + as.numeric(multiclass.roc(observeds, tabulate(simulated_document_matrix[o,1], 33))$auc)
        F1score = F1score + F1_Score(observedr, tabulate(which(simulated_document_matrix[o,2:34] == 1), 33))
        MSE = MSE + (simulated_document_matrix[o,35]- observedt)^2
    }
    tables[d, 1] = AUC / 100
    tables[d, 2] = F1score / 100
    tables[d, 3] = MSE / 100
  }

tables = tables[order(selectD[1:22]),]
tables = tables[- which(tables[,3] > 50),]
par(mfrow = c(1,3))
plot(rownames(tables), tables[,1], type = 'l', col = 'black', xlab = "Document", lty = 2, ylab = "Sender AUC", ylim = c(min(tables[,1]), max(tables[,1])), main = "Enron")
#lines(colnames(Dare_topic), Dare_topic[2,], lty = 3,col = "blue")
#lines(colnames(Dare_topic), Dare_topic[3,], lty = 4, col = "green")
#legend(65, -405, pch = 21, legend = c("IPTM with C = 3","IPTM with C = 2","LDA"), col = c("green","blue","red"))

plot(rownames(tables), tables[,2], type = 'l', col = 'red', xlab = "Document", lty = 2, ylab = "Receiver F1Score", ylim = c(min(tables[,2]), max(tables[,2])), main = "Enron")
#lines(colnames(Dare_topic), Dare_topic[2,], lty = 3,col = "blue")
#lines(colnames(Dare_topic), Dare_topic[3,], lty = 4, col = "green")
#legend(65, -405, pch = 21, legend = c("IPTM with C = 3","IPTM with C = 2","LDA"), col = c("green","blue","red"))


plot(rownames(tables), tables[,3], type = 'l', col = 'red', xlab = "Document", lty = 2, ylab = "Time MSE", ylim = c(min(tables[,3]), max(tables[,3])), main = "Enron")
#lines(colnames(Dare_topic), Dare_topic[2,], lty = 3,col = "blue")
#lines(colnames(Dare_topic), Dare_topic[3,], lty = 4, col = "green")
#legend(65, -405, pch = 21, legend = c("IPTM with C = 3","IPTM with C = 2","LDA"), col = c("green","blue","red"))



table_sum = matrix(rep(colMeans(tables), 7), 3)
rownames(table_sum) =  c("sender", "receiver", "time")
colnames(table_sum) = c("2", "5", "10", "25", "50", "75", "100")



setwd('/Users/bomin8319/Desktop/IPTM/Enron_PPE')
tables = list()
for (d in 1:22) {
  D = selectD[d]
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
        F1score = F1score + F1_Score(observedr, tabulate(PPE[[o]][[r]]$edge$receiver, 33))
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
}
tables = tables[order(selectD[1:22])]
all = Reduce('+', tables[-c(1,12)]) / length(tables[-c(1,12)])

par(mfrow = c(1,3))
plot(colnames(table_sum), table_sum[1,], type = 'b', pch = 19, col = 'red', xlab = "Number of Topics", lty = 2, ylab = "Sender AUC", ylim = c(table_sum[1,1] - 0.02, table_sum[1,1] + 0.02))
lines(colnames(table_sum), all[1:7,1], lty = 3,col = "blue", type = 'b', pch = 19)
lines(colnames(table_sum), all[8:14,1], lty = 4,col = "green",type = 'b', pch = 19)
lines(colnames(table_sum), all[15:21,1], lty = 5,col = "purple",type = 'b', pch = 19)
legend(60, 0.545, pch = 19, bty = "n", lty = 2:5, legend = c("IPTM with C = 3","IPTM with C = 2","IPTM with C = 1", "Regressions"), col = c("purple", "green", "blue", "red"))

plot(colnames(table_sum), table_sum[2,], type = 'b', pch = 19,  col = 'red', xlab = "Number of Topics", lty = 2, ylab = "Receiver F1Score", ylim = c(table_sum[2,1] - 0.02, table_sum[2,1] + 0.01))
lines(colnames(table_sum), all[1:7,2], lty = 3,col = "blue",type = 'b', pch = 19)
lines(colnames(table_sum), all[8:14,2], lty = 4,col = "green",type = 'b', pch = 19)
lines(colnames(table_sum), all[15:21,2], lty = 5,col = "purple",type = 'b', pch = 19)
legend(60, 0.96,  pch = 19,bty = "n", lty = 2:5, legend = c("IPTM with C = 3","IPTM with C = 2","IPTM with C = 1", "Regressions"), col = c("purple", "green", "blue", "red"))

plot(colnames(table_sum), table_sum[3,], type = 'b', pch = 19, col = 'red', xlab = "Number of Topics", lty = 2, ylab = "Time MSE", ylim = c(0, 1.6))
lines(colnames(table_sum), all[1:7,3], lty = 3,col = "blue",type = 'b', pch = 19)
lines(colnames(table_sum), all[8:14,3], lty = 4,col = "green",type = 'b', pch = 19)
lines(colnames(table_sum), all[15:21,3], lty = 5,col = "purple",type = 'b', pch = 19)
legend(60, 1.55,  pch = 19, bty = "n",lty = 2:5, legend = c("IPTM with C = 3","IPTM with C = 2","IPTM with C = 1", "Regressions"), col = c("purple", "green", "blue", "red"))


##################################################

load('/Users/bomin8319/Desktop/IPTM/enron/Dare.RData')
set.seed(1)
selectD = sample(1105:2210, 200, replace = FALSE)
setwd('/Users/bomin8319/Desktop/IPTM/paper/code/predictiveComparisonModels/results')
 tables = matrix(0, nrow = 200, ncol = 3)
  colnames(tables) = c("sender", "receiver", "time")
  rownames(tables) = selectD
for (d in 1:3) {
  D = selectD[d]
  observeds = tabulate(Dare$edge[[D]]$sender, 27)
  observedr = tabulate(Dare$edge[[D]]$receiver, 27)
  observedt = Dare$edge[[D]]$timeinc
   filename = paste0("DareSimulatedDocumentsD",D,".RData")
  load(filename)
    MSE = 0
    F1score = 0
    AUC = 0
    for (o in 1:100) {
        AUC = AUC + as.numeric(multiclass.roc(observeds, tabulate(simulated_document_matrix[o,1], 27))$auc)
        F1score = F1score + F1_Score(observedr, tabulate(which(simulated_document_matrix[o,2:28] == 1), 27))
        MSE = MSE + (simulated_document_matrix[o,29]- observedt)^2
    }
    tables[d, 1] = AUC / 100
    tables[d, 2] = F1score / 100
    tables[d, 3] = MSE / 100
  }

tables = tables[order(selectD[1:3]),]
#tables = tables[- which(tables[,3] > 50),]

par(mfrow = c(1,3))
plot(rownames(tables), tables[,1], type = 'l', col = 'black', xlab = "Document", lty = 2, ylab = "Sender AUC", ylim = c(min(tables[,1]), max(tables[,1])), main = "Dare")
#lines(colnames(Dare_topic), Dare_topic[2,], lty = 3,col = "blue")
#lines(colnames(Dare_topic), Dare_topic[3,], lty = 4, col = "green")
#legend(65, -405, pch = 21, legend = c("IPTM with C = 3","IPTM with C = 2","LDA"), col = c("green","blue","red"))

plot(rownames(tables), tables[,2], type = 'l', col ='black', xlab = "Document", lty = 2, ylab = "Receiver F1Score", ylim = c(min(tables[,2]), max(tables[,2])), main = "Dare")
#lines(colnames(Dare_topic), Dare_topic[2,], lty = 3,col = "blue")
#lines(colnames(Dare_topic), Dare_topic[3,], lty = 4, col = "green")
#legend(65, -405, pch = 21, legend = c("IPTM with C = 3","IPTM with C = 2","LDA"), col = c("green","blue","red"))


plot(rownames(tables), tables[,3], type = 'l', col ='black', xlab = "Document", lty = 2, ylab = "Time MSE", ylim = c(min(tables[,3]), max(tables[,3])), main = "Dare")
#lines(colnames(Dare_topic), Dare_topic[2,], lty = 3,col = "blue")
#lines(colnames(Dare_topic), Dare_topic[3,], lty = 4, col = "green")
#legend(65, -405, pch = 21, legend = c("IPTM with C = 3","IPTM with C = 2","LDA"), col = c("green","blue","red"))


table_sum = matrix(rep(colMeans(tables), 7), 3)
rownames(table_sum) =  c("sender", "receiver", "time")
colnames(table_sum) = c("2", "5", "10", "25", "50", "75", "100")


setwd('/Users/bomin8319/Desktop/IPTM/Dare_PPE')
tables = list()
for (d in 1:3) {
  D = selectD[d]
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
    MSE = 0
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
}

tables = tables[order(selectD[1:3])]
all = Reduce('+', tables) / length(tables)

par(mfrow = c(1,3))
plot(colnames(table_sum), table_sum[1,], type = 'b', pch = 19, col = 'red', xlab = "Number of Topics", lty = 2, ylab = "Sender AUC", ylim = c(table_sum[1,1] - 0.04, table_sum[1,1] + 0.02))
lines(colnames(table_sum), all[1:7,1], lty = 3,col = "blue", type = 'b', pch = 19)
lines(colnames(table_sum), all[8:14,1], lty = 4,col = "green",type = 'b', pch = 19)
lines(colnames(table_sum), all[15:21,1], lty = 5,col = "purple",type = 'b', pch = 19)
legend(60, 0.53, pch = 19, bty = "n", lty = 2:5, legend = c("IPTM with C = 3","IPTM with C = 2","IPTM with C = 1", "Regressions"), col = c("purple", "green", "blue", "red"))

plot(colnames(table_sum), table_sum[2,], type = 'b', pch = 19,  col = 'red', xlab = "Number of Topics", lty = 2, ylab = "Receiver F1Score", ylim = c(table_sum[2,1] - 0.35, table_sum[2,1] + 0.15))
lines(colnames(table_sum), all[1:7,2], lty = 3,col = "blue",type = 'b', pch = 19)
lines(colnames(table_sum), all[8:14,2], lty = 4,col = "green",type = 'b', pch = 19)
lines(colnames(table_sum), all[15:21,2], lty = 5,col = "purple",type = 'b', pch = 19)
legend(60, 1.05,  pch = 19,bty = "n", lty = 2:5, legend = c("IPTM with C = 3","IPTM with C = 2","IPTM with C = 1", "Regressions"), col = c("purple", "green", "blue", "red"))

plot(colnames(table_sum), table_sum[3,], type = 'b', pch = 19, col = 'red', xlab = "Number of Topics", lty = 2, ylab = "Time MSE", ylim = c(0, 0.1))
lines(colnames(table_sum), all[1:7,3], lty = 3,col = "blue",type = 'b', pch = 19)
lines(colnames(table_sum), all[8:14,3], lty = 4,col = "green",type = 'b', pch = 19)
lines(colnames(table_sum), all[15:21,3], lty = 5,col = "purple",type = 'b', pch = 19)
legend(60, 0.1,  pch = 19, bty = "n",lty = 2:5, legend = c("IPTM with C = 3","IPTM with C = 2","IPTM with C = 1", "Regressions"), col = c("purple", "green", "blue", "red"))


