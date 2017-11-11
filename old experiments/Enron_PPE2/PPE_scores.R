library(MLmetrics)
library(pROC)
load('/Users/bomin8319/Desktop/IPTM/enron/Enron.RData')
set.seed(1)
selectD = sample(1963:3925, 200, replace = FALSE)
setwd('/Users/bomin8319/Desktop/IPTM/Enron_PPE2')
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
  observeds = tabulate(Enron$edge[[D]]$sender, 33)
  observedr = tabulate(Enron$edge[[D]]$receiver, 33)
  observedt = Enron$edge[[D]]$timestamp - Enron$edge[[D-1]]$timestamp
   filename = paste0("EnronSimulatedDocumentsD",D,".RData")
  load(filename)
    MSE = c()
    F1score = 0
    AUC = 0
    for (o in 1:100) {
        AUC = AUC + as.numeric(multiclass.roc(observeds, tabulate(simulated_document_matrix[o,1], 33))$auc)
        F1score = F1score + F1_Score(observedr, tabulate(which(simulated_document_matrix[o,2:34] == 1), 33))
        MSE[o] = abs(simulated_document_matrix[o,35]- observedt)
    }
    tables[d, 1] = AUC / 100
    tables[d, 2] = F1score / 100
    tables[d, 3] = median(MSE)
 d = d + 1
  }

# #tables = tables[-which(tables[,3] > 50),]
# par(mfrow = c(1,3))
# plot(rownames(tables), tables[,1], type = 'l', col = 'black', xlab = "Document", lty = 2, ylab = "Sender AUC", ylim = c(min(tables[,1]), max(tables[,1])), main = "Enron")
# #lines(colnames(Dare_topic), Dare_topic[2,], lty = 3,col = "blue")
# #lines(colnames(Dare_topic), Dare_topic[3,], lty = 4, col = "green")
# #legend(65, -405, pch = 21, legend = c("IPTM with C = 3","IPTM with C = 2","LDA"), col = c("green","blue","red"))

# plot(rownames(tables), tables[,2], type = 'l', col = 'red', xlab = "Document", lty = 2, ylab = "Receiver F1Score", ylim = c(min(tables[,2]), max(tables[,2])), main = "Enron")
# #lines(colnames(Dare_topic), Dare_topic[2,], lty = 3,col = "blue")
# #lines(colnames(Dare_topic), Dare_topic[3,], lty = 4, col = "green")
# #legend(65, -405, pch = 21, legend = c("IPTM with C = 3","IPTM with C = 2","LDA"), col = c("green","blue","red"))


# plot(rownames(tables), tables[,3], type = 'l', col = 'red', xlab = "Document", lty = 2, ylab = "Time MSE", ylim = c(min(tables[,3]), max(tables[,3])), main = "Enron")
# #lines(colnames(Dare_topic), Dare_topic[2,], lty = 3,col = "blue")
# #lines(colnames(Dare_topic), Dare_topic[3,], lty = 4, col = "green")
# #legend(65, -405, pch = 21, legend = c("IPTM with C = 3","IPTM with C = 2","LDA"), col = c("green","blue","red"))



table_sum2 = matrix(rep(colMeans(tables), 7), 3)
rownames(table_sum2) =  c("sender", "receiver", "time")
colnames(table_sum2) = c("2", "5", "10", "25", "50", "75", "100")



tables2 = list()
d = 1
for (D in commondocs) {
  observeds = tabulate(Enron$edge[[D]]$sender, 33)
  observedr = tabulate(Enron$edge[[D]]$receiver, 33)
  observedt = Enron$edge[[D]]$timestamp
  tables2[[d]] = matrix(0, nrow = 21, ncol = 3)
  rownames(tables2[[d]]) = 1:21
  iter = 1
  for (nIP in 1:3) {
    for (K in c(2, 5, 10, 25, 50, 75, 100)){
    filename = paste0("PPE_",D,"_",nIP,"_",K,".RData")
    load(filename)
    rownames(tables2[[d]])[iter] = filename
    MSE = c()
    F1score = 0
    AUC = 0
    it = 1
    for (o in 1:10) {
      for (r in 1:50) {
          if (r > 10 && (r-10) %% 4 == 0) {
        AUC = AUC + as.numeric(multiclass.roc(observeds, tabulate(PPE[[o]][[r]]$edge$sender, 33))$auc)
        F1score = F1score + F1_Score(observedr, tabulate(PPE[[o]][[r]]$edge$receiver, 33))
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
  d = d+1
}
all2 = Reduce('+', tables2) / length(tables2)

par(mfrow = c(1,3))
plot(colnames(table_sum2), table_sum2[1,], type = 'b', pch = 19, col = 'red', xlab = "Number of Topics", lty = 2, ylab = "Sender AUC", ylim = c(0.51, 0.56), main = "Enron")
lines(colnames(table_sum2), all2[1:7,1], lty = 3,col = "blue", type = 'b', pch = 19)
lines(colnames(table_sum2), all2[8:14,1], lty = 4,col = "green",type = 'b', pch = 19)
lines(colnames(table_sum2), all2[15:21,1], lty = 5,col = "purple",type = 'b', pch = 19)
legend(60, 0.56, pch = 19, bty = "n", lty = 2:5, legend = c("IPTM with C = 3","IPTM with C = 2","IPTM with C = 1", "Regressions"), col = c("purple", "green", "blue", "red"))

plot(colnames(table_sum2), table_sum2[2,], type = 'b', pch = 19,  col = 'red', xlab = "Number of Topics", lty = 2, ylab = "Receiver F1Score", ylim = c(0.935, 0.955))
lines(colnames(table_sum2), all2[1:7,2], lty = 3,col = "blue",type = 'b', pch = 19)
lines(colnames(table_sum2), all2[8:14,2], lty = 4,col = "green",type = 'b', pch = 19)
lines(colnames(table_sum2), all2[15:21,2], lty = 5,col = "purple",type = 'b', pch = 19)
legend(60, 0.955,  pch = 19,bty = "n", lty = 2:5, legend = c("IPTM with C = 3","IPTM with C = 2","IPTM with C = 1", "Regressions"), col = c("purple", "green", "blue", "red"))

plot(colnames(table_sum2), table_sum2[3,], type = 'b', pch = 19, col = 'red', xlab = "Number of Topics", lty = 2, ylab = "Time MAD", ylim = c(9.5, 10.0))
lines(colnames(table_sum2), all2[1:7,3], lty = 3,col = "blue",type = 'b', pch = 19)
lines(colnames(table_sum2), all2[8:14,3], lty = 4,col = "green",type = 'b', pch = 19)
lines(colnames(table_sum2), all2[15:21,3], lty = 5,col = "purple",type = 'b', pch = 19)
legend(60, 10.0,  pch = 19, bty = "n",lty = 2:5, legend = c("IPTM with C = 3","IPTM with C = 2","IPTM with C = 1", "Regressions"), col = c("purple", "green", "blue", "red"))
