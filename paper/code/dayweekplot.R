setwd('/Users/bomin8319/Desktop/IPTM/full')
library(fields)
library(ergm)
library(anytime)
load('Darenew.RData')
# 762 - 
attach(Dare)
Dare$edge = Dare$edge[-which(sapply(Dare$text, function(d){length(d)})==0)]
Dare$text = Dare$text[-which(sapply(Dare$text, function(d){length(d)})==0)]

wday = function(x) {
	as.POSIXlt(x)$wday
}

weeks = sapply(Dare$edge, function(x) wday(anytime(x$unixtime,  tz = "America/New_York")) )

timeinc = c(sapply(1:2209, function(x) Dare$edge[[x+1]]$unixtime - Dare$edge[[x]]$unixtime), 0)
boxplot(timeinc~weeks, names = c("Sun", "Mon", "Tue", "Wed", "Thur", "Fri", "Sat"))

rm(list=ls())
setwd('/Users/bomin8319/Desktop/IPTM/full')
library(fields)
library(ergm)
library(lubridate)
load('Darenew.RData')
# 762 - 
attach(Dare)
Dare$edge = Dare$edge[-which(sapply(Dare$text, function(d){length(d)})==0)]
Dare$text = Dare$text[-which(sapply(Dare$text, function(d){length(d)})==0)]
timeinc = c(sapply(1:2209, function(x) Dare$edge[[x+1]]$unixtime - Dare$edge[[x]]$unixtime), 0)

day = sapply(Dare$edge, function(x) hour(ymd_hms(anytime(x$unixtime,  tz = "America/New_York"))))

day = cut(day, c(0,6,9,14,Inf), c("Overnight", "Morning", "Afternoon", "Prime"))
day[is.na(day)] = "Overnight" 


boxplot(timeinc~day)

