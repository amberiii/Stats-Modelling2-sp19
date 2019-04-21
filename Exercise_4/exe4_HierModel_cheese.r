# Exercise 04
# Cheese Case - Price Elasticity Modeling (Hierarchical Models)

rm(list=ls())
library(mosaic)
library(lme4)
library(dplyr)
library(ggplot2)

data = read.csv('/Users/ippjune/Documents/2019_Spring/StatsModeling/cheese.csv', head = T)
attach(data)

# display vs price -- location
#Read in data.

#Read in functions.
source('/Users/ippjune/Documents/2019_Spring/StatsModeling/exe4_FUNCTIONS_cheese.r')

#Add log price and log volume variables.
data$logP = log(data$price)
data$logQ = log(data$vol)

#Add a week variable.
data$week = ave(rep(1:nrow(data)), data$store, FUN=seq_along)

#Number of stores.
s = length(unique(data$storenum))

#Add store numbers, then convert to list.
store.names.all = data$store
store.names = levels(store.names.all)

levels(data$store) = 1:length(levels(data$store))
data$store = as.numeric(data$store)

#Order data by store number, then week.
data = data[with(data,order(data$store,data$week)), ]
rownames(data) = 1:nrow(data)

#================================================================
# Exploratory ===================================================
#================================================================

### Inspect grouping structure.
xtabs( ~ factor(week) + factor(store), data=data)
#Shows that each week has one observation per store.
#Notice there are 68 weeks, and some stores are missing a few obs.
xtabs( ~ factor(week), data=data)
#Shows how many observations each week has. 
#Can see there are some weeks near the end with 87, and
#a few weeks with only 28.  Not a balanced model.	

#A few stores have almost all ad, or almost no ad weeks.  (Add pictures for this.)
#Is a good reason to use pooling in this case.
#ADD MORE EDA PLOTS.	

### Does avg price appear to vary by store?

vol.avg.by.store = aggregate(vol,list(store),mean)
#par(mfrow=c(1,2))
#boxplot(vol ~ store)
plot(vol.avg.by.store$x,pch=19,main='Avg Volume by Store Over All Weeks',xlab='Store',ylab='Vol')


#Exploratory Data Analysis
#boxplots of log vol ~ disp
boxplot(logQ ~ disp, data=data)

#Obvious confounders? (Corr with both predictor and response)
plot(logQ ~ logP, data=data) #Globally  
boxplot(logP ~ disp, data=data)

#Interactions.
xyplot(logQ ~ logP | store, data=data)

#================================================================
# Naive Approach ================================================
#================================================================

# This approach ignores the pseudoreplication due to weekly 
# measurements of same stores.

# Main effects
lm0 = lm(logQ ~ logP + disp, data = data)
lm1 = lm(logQ ~ logP + disp + store, data=data)
lm2 = lm(logQ ~ logP + store, data=data)

summary(lm0)
summary(lm1)
anova(lm1)

anova(lm1,lm0)	#F-test, concludes store is significant.
anova(lm2,lm0)	#F-test, concludes disp is significant.

#================================================================
# Hierarchical Model: Empirical Bayes ===========================
#================================================================

#----------------------------------------------------------------
#	Model: logQ_{it} = log(alpha_{i}) + b_{i}*log(P_{it}) + gamma_{i} * ad_{it} + v_{i} + e_{it}
#		for i = 1...88 indexing stores, and t indexing time points.
#		logQ_{it}	= Volume at each store at each time.
#
#	This model uses empirical Bayes; errors estimated from data.
#----------------------------------------------------------------

#Hierarchical linear model; allows intercept, logQ and disp to change among stores.
#	Includes ad*price interaction.
hlm = lmer(logQ ~ (logP + disp + disp:logP | store), data=data)
summary(hlm)

#Plot coefficients.
coef(hlm)

#Plot conditional modes (same as means, for linear mixed models) of random effects.
resid = ranef(hlm, condVar = T)
dotplot(resid,scales=list(cex=c(.5,.5)),layout=c(3,1),main=T, main.title='Random Effects by by Store',)
#dev.off()

plot(hlm)
#dev.off()


#================================================================
# Hierarchical Model: Fully Bayesian w Gibbs Sampler ============
#================================================================

#Data setup for function.
y = data$logQ																	#Responses for all i,t
X = as.data.frame(model.matrix(logQ ~ logP + disp + logP*disp,data=data))		#Covariates for all i,t
idx = data[,c("store","week")]	

#OLS for muB estimates.
ols = lm(logQ ~ logP + disp + logP*disp,data=data)

#Noninformative hyperpriors for muB,v,C,d.
p = ncol(X)
muB = ols$coef
V = diag(p)
C = diag(p)
d = p+1

#Call gibbs sampler function.
output = gibbs.cheese(y,X,idx,muB,V,C,d,iter=5100,burn=100,thin=2)

#Set column names for posterior bi vector.
rownames(output$bi.pm) = c("Intercept","logP","disp","logP:disp")

#----------------------------------------------------------------
#Plot results.

#Trace plots.
par(mfrow=c(2,2))
plot(output$bi[1,])
plot(output$Beta[1])
plot(output$bi[2,])
plot(output$Beta[2])
dev.off()

#Demand curves for each store.

# Add ad and non-ad demand curves for this specific store.
# Formula for curve:  
#		AD = NO group: exp(log(Intercept)) * x ^ (log(Price)))
#		AD = YES group: 	exp(log(Intercept + disp) * x ^ (log(Price) - log(Price):disp))
#                             log(price)     disp log(price):disp (Intercept)

# Under my model formulation:
# AD = NO:	exp(log(B + Intercept_i))

par(mfrow=c(4,5),oma=c(1,1,0,0) + .1, mar=c(0,0,1,1)+1)

#Extract values for demand curves for store i.
int 	= output$Beta.pm[1] + output$bi.pm[1,]
logp 	= output$Beta.pm[2] + output$bi.pm[2,]
dsp 	= output$Beta.pm[3] + output$bi.pm[3,]
logp.d 	= output$Beta.pm[4] + output$bi.pm[4,]

sort.avg.price = order(aggregate(data$price, list(data$store), mean)$x,decreasing=T)

for (i in sort.avg.price){
  
  #Plot 'AD = 1' points in red.
  plot(vol ~ price, data=subset(data, data$store == i),
       xlab=paste('Store',i),ylab='',col='red',xaxt='n',yaxt='n')
  
  #Add the 'AD = 0' points in blue.
  points(vol ~ price, col='blue',
         data=subset(data, data$store == i & data$disp == 0))
  
  #Fitted demand curve for disp=0 group at store i.
  curve(exp(int[i])*x^(logp[i]), add=TRUE, col='blue')
  
  #Fitted demand curve for disp=1 group at store i.
  curve(exp(int[i] + dsp[i])*x^(logp[i] + logp.d[i]), add=TRUE, col='red')
}

dev.off()