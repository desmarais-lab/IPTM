load('/Users/bomin8319/Desktop/IPTM/enron/Enron.RData')
Enron$edge[[2484]]

load('/Users/bomin8319/Desktop/Enron_PPE/PPE_2484_1_5.RData')

sender = c()
for (o in 1:10) {
  for (r in 1:50) {
    sender = c(sender, PPE[[o]][[r]]$edge$sender)
  }
}
sum(sender == Enron$edge[[2484]]$sender) / 500

receiver = c()
for (o in 1:10) {
  for (r in 1:50) {
    receiver = c(receiver, PPE[[o]][[r]]$edge$receiver)
  }
}
sum(receiver %in% Enron$edge[[2484]]$receiver) / 500

lreceiver = c()
for (o in 1:10) {
  for (r in 1:50) {
    lreceiver = c(lreceiver, which(PPE[[o]][[r]]$iJi[Enron$edge[[2484]]$sender,] ==1))
  }
}
sum(lreceiver %in% Enron$edge[[2484]]$receiver) / 500

multicast = c()
for (o in 1:10) {
  for (r in 1:50) {
    multicast = c(multicast, length(PPE[[o]][[r]]$edge$receiver))
  }
}

time = c()
for (o in 1:10) {
  for (r in 1:50) {
    time = c(time, PPE[[o]][[r]]$edge$timestamp)
  }
}

hist(time)
abline(v = Enron$edge[[2484]]$timestamp, col = 'red')


MSE = sum((time-Enron$edge[[2484]]$timestamp)^2) /500
print(MSE)