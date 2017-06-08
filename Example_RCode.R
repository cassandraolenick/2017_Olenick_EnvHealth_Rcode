
# Example R code for Stage 2 models described by O'Lenick et al. 
# in "Ozone and childhood respiratory disease in three US cities: 
# evaluation of effect measure modification by neighborhood socioeconomic status 
# using a Bayesian hierarchical approach"  
# (Environmental Health, 2017, DOI: 10.1186/s12940-017-0244-2)

#Requires example dataset, ENHE-D-16-00179-R1_Additional_Files_Example.Data.csv

library(stats)
library(tlnise)
options()
options(digits=6, scipen=6)


#call dataset with ZCTA-specific parameters from simulated Stage 1 analyses
data <- read.csv("~/Desktop/ENHE-D-16-00179-R1_Additional_Files_Example.Data.csv")

###################
#### Estimating overall associations between ozone and respiratory disease for each city  ####
###################

### extract parameter estimates for meta-regression (assocations between ozone and respiratory disease)
logOR = data$logOR
stderr = data$SE

### make indicator variable for city

city1= as.numeric (data$city=="ATL")
city2= as.numeric (data$city=="DFW")
city3= as.numeric (data$city=="STL")
X = cbind(city1, city2, city3) # in overall models, there is an indicator for each city. There is not an intercept. 

library(tlnise)

### fit model to estimate overall associations between ozone and respiratory disease for each city
overall_fit=tlnise (logOR, stderr^2, X, intercept=F, prnt = F) # no intercept model
overall_fit$gamma #prints betas and standard errors

###################
#### Examining continous ZCTA-level SES as an effect modifier of ozone-related respiratory disease
###################

#### Create the design matrix with an intercept, 
#       an indicator variable for 2 of the 3 cities, and continuous ZCTA-level % below poverty.

#   The intercept will represent a ZCTA in Atlanta with 0% poverty
poverty=data$pov_percent    
Xses=cbind (1, city2, city3, poverty, poverty^2, poverty^3) #design matrix

### fit model to examine non-linear (cubic) ZCTA-level % below poverty as an effect modifier of associations between ozone and respiratory disease for each city
effect.mod_fit=tlnise (logOR, stderr^2, Xses, intercept=T, prnt = F) # model fit with an intercept (intercept=T)

### extracting parameter estimates and estimating lower and upper PIs (posterior intervals) 
### est/se z-score for contribution of each variable to the model
effect.mod_fit$gamma #effect estimates and standard errors
effect.mod_fit$Dgamma # covariance matrix

covar.matrix<-effect.mod_fit$Dgamma #extract covariance matrix
var<-Xses%*%covar.matrix%*%t(Xses) #estimate variance for each zcta log OR
diagonal.var<-diag(var) #extract diagonals
std.err=sqrt(diagonal.var) #standard error for each ZCTA log OR

effect.est<-effect.mod_fit$gamma[,1] #extract log ORs
odds.ratio=exp((Xses%*%effect.est)*25) # estimate Odds Ratio: exponentiate and scale to 25 ppb
odds.ratio.lowerCI=exp(((Xses%*%effect.est)*25) - (1.96*(std.err*25))) #estimate 95% lower CI with standard error:exponentiate and scale to 25 ppb
odds.ratio.upperCI=exp(((Xses%*%effect.est)*25) + (1.96*(std.err*25)))#estimate 95% upper CI with standard error:exponentiate and scale to 25 ppb

#Merge output into a new dataset that includes estimated ORs and 95% CIs for each ZCTA
data.merge<-(cbind(data$ZCTA,data$city,data$pov_percent,odds.ratio,odds.ratio.lowerCI,odds.ratio.upperCI))
colnames(data.merge)<-c("ZCTA","city", "per_pov" ,"OR","LowerCI","UpperCI")
data.results=as.data.frame(data.merge)
#order by %poverty (not by ZCTA)
data.results<-data.results[with(data.results,order(per_pov)),]
View(data.results)

##########################
#### EXAMPLE PLOTS - Atlanta
###########################

# Example R code for results plots similar to those presented in Figure 3a in
# "Ozone and childhood respiratory disease in three US cities: 
# evaluation of effect measure modification by neighborhood socioeconomic status 
# using a Bayesian hierarchical approach"  
# (O'Lenick, et al. Environmental Health, 2017, DOI: 10.1186/s12940-017-0244-2)
# plots may need adjusting depending on plot margins

atlanta<- subset(data.results, city==1) 
View(atlanta)

par(mar=c(4,4.1,3,2.4))
par(bty="n")

#plot estimated mean OR (from analysis above) for each ZCTA by ZCTA-level poverty
plot(atlanta$OR~atlanta$per_pov,type="l",
     yaxt="n",lab=c(6,5,7), 
     xlab= "",
     ylab="",
     lty="dashed",
     ylim=c(.7,1.3),
     xlim=c(0,52),
     cex.lab=1,las=1,cex.axis=1,main="",cex.main=1)

# add polygon and overlplot the lines

# gray polygons represent 95% PI
polygon(c(atlanta$per_pov,rev(atlanta$per_pov)),c(atlanta$LowerCI,rev(atlanta$UpperCI)),col="gray",border=NA)
# Black line represents mean OR for each ZCTA
lines(atlanta$per_pov,atlanta$OR,lty="solid",col="black",lwd=1.5)
# add line to show where OR = 1.0
abline(h=1, lwd=1)
#Adjust y-axis
axis(2,at=c(0.9,1,1.1,1.2,1.3),las=2, cex.axis=1, font=2)

#extract 2.5 and 97.5 percentile of ZCTA-level % below poverty
atl_qts <- quantile(atlanta$per_pov,probs=c(.025,.975))
atl_qts #3.73, 35.6
abline(v=atl_qts[1],lty=2)
abline(v=atl_qts[2],lty=2)

#labeling the y-axis
mtext("OR per 25 ppb ozone", side=2, line=2.7, las=0, cex=1,adj=0, font=2)
#labeling the x-axis
mtext("% Below Poverty", side=1, line=3, las=0, cex=1.2, font=2)

#overlay new plot
par(new=T)
par(mar=c(4,4.1,8,2.4))
par(bty="n")

# adding histrogram below graph
hist <- hist(atlanta$per_pov,breaks=seq(0,52,by=2),ann=F,col=grey(0.95),axes=F,plot=T)
# add axis with frequency-count for histogram
axis(4,at=c(0,5,10,15,20),pos=53,las =2, cex.axis=0.7)

#adding verticle dotted lines to indicate the 2.5 and 97.5 percentiles of ZCTA % below poverty
abline(v=atl_qts[1],lty=2)
abline(v=atl_qts[2],lty=2)


