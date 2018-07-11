#Library
library(vars)
library(ucminf)
library(rgenoud)
library(matlib)
##library(MTS)

#----------------------------------------------------------------------------------------

########################################################################
##Data Pre-processing
##
##T-bill:       Statistic Canada 1-month and 3 month rate
##Bond yield:   Statistic Canada 2 year, 3 year, 5 year, 7 year, 10 year, Long term
##Inflation:    Statistic Canada CPI
##Stock return: Statistic Canada
##              CANADA: S&P/TSX composite index
##              
##
##
##Time span:    Data available from Jan 1986 - June 2016
##
##All the datas are the continuously compounded annual rate, with sample taken in monthly frequency
##Data used to estimate VAR model is 1-month T-bill, 10 year bond, inflation, US and Canada Stock return
##########################################################################

#bonds
CDNbonds <- read.table("bond yield partial.csv", header=TRUE, sep=",", na.strings=" ")

plot(CDNbonds$X2.year, type="l",xaxt='n', col="green", main="Bond Yields",ylab="annualized bond yields",lwd=1,xlab="time")
axis(side=1, at=seq(1, 306, by=60),labels=CDNbonds$Date[seq(1, 306, by=60)])
legend(x=10,y=0.09,legend=c("3 month","1 year", "2 year","5 year", "10 year", "15 year"), col=c("red", "green", "purple", "orange", "black"),lty=1,lwd=1, horiz = T)
lines(CDNbonds$X3.month, type="l",lwd=1,col="red")
lines(CDNbonds$X5.year, type="l",lwd=1,col="purple")
lines(CDNbonds$X10.year, type="l",lwd=1,col="orange")
lines(CDNbonds$X15.year, type = "l", lwd = 1, col = "black")
lines(CDNbonds$X1.year, type = "l", lwd = 1, col = "blue")

BondsMonthly = as.matrix(CDNbonds[c(2:10)])
BondsMonthly = BondsMonthly / 12

plot(BondsMonthly[, 2], type = "l", xaxt = 'n', col = "green", main = "Bond Yields", ylab = "Monthly Rate", lwd = 1.5, xlab = "Date")
axis(side = 1, at = seq(1, 306, by = 60), labels = CDNbonds$Date[seq(1, 306, by = 60)])
legend(x = 95, y = 0.0076,legend = c("3-month", "1-year", "2-year", "5-year", "10-year", "15-year"), col = c("red", "green", "purple", "orange", "black"), lty = 1, lwd = 1, horiz = T, cex =0.6)
ylim = -0.001:0.009
lines(BondsMonthly[, 3], type = "l", lwd = 1, col = "red")
lines(BondsMonthly[, 4], type = "l", lwd = 1, col = "purple")
lines(BondsMonthly[, 6], type = "l", lwd = 1, col = "orange")
lines(BondsMonthly[, 8], type = "l", lwd = 1, col = "black")
lines(BondsMonthly[, 9], type = "l", lwd = 1, col = "blue")


#inflation
CDNinflation <- read.table("canada cpi.csv", header=TRUE, sep=",", na.strings=" ")
CDNinfrate<-CDNinflation$continuous

plot(CDNinfrate, type="l", main = "Inflation Rate", xaxt='n', xlab="Date", ylab="Monthly Rate")
axis(side=1, at=seq(1, 306, by=60),labels=CDNbonds$Date[seq(1, 306, by=60)])

#Canadian stock TSX 
TSXindex <- read.table("TSX with dividen.csv", header=TRUE, sep=",", na.strings=" ")
CDNstock <- TSXindex$conti
#CDNstock <- TSXindex$conti + CDNdividend

plot(12*CDNstock, type="l", main = "TSX stock return", xaxt='n', xlab="time",ylab="continuously compounded annual rate")
axis(side=1, at=seq(1, 306, by=60),labels=CDNbonds$Date[seq(1, 306, by=60)])

#Excess Stock Resturn
ExCDNstock <- CDNstock - BondsMonthly[,2]

plot(ExCDNstock, type="l", main = "Excess Stock Return", xaxt='n', xlab="Date",ylab="Monthly Rate")
axis(side=1, at=seq(1, 306, by=60),labels=CDNbonds$Date[seq(1, 306, by=60)])

#Canadian stock dividend yield
CDNdividend <-TSXindex$DividenMonth

plot(CDNdividend, type="l", main = "Dividend Yield", xaxt='n', xlab="Date",ylab="Monthly Rate")
axis(side = 1, at = seq(1, 306, by = 60), labels = CDNbonds$Date[seq(1, 306, by = 60)])

#make z4 as stock return + dividend yield - bondsMonthly
ExCDNstock <- ExCDNstock + CDNdividend

#-----------------------------------------------------------------------------------------



#----------------------------------------------------------------------------------------

##VAR Model

##VAR input data

##Data Without Dividend Yield
VARinput <- data.frame(tbill = BondsMonthly[, 2], bonds15 = BondsMonthly[, 9], infrate = CDNinfrate, TSXexcess = ExCDNstock)
MeanData <- apply(X = VARinput, MARGIN = 2, FUN = mean)
StdData <- apply(X = VARinput, MARGIN = 2, FUN = sd)

##Combine the input of VAR data into one data frame
VARinput <- data.frame(tbill = BondsMonthly[,2], bonds15 = BondsMonthly[,9], infrate = CDNinfrate, TSXexcess = ExCDNstock)
MeanData <- apply(X=VARinput, MARGIN=2, FUN=mean)
StdData  <- apply(X=VARinput, MARGIN=2, FUN=sd)

##Center the data at mean=0
ZeroMeanData <- VARinput
ZeroMeanData$tbill <- VARinput$tbill-MeanData[1]
ZeroMeanData$bonds15 <- VARinput$bonds15-MeanData[2]
ZeroMeanData$infrate <- VARinput$infrate-MeanData[3]
ZeroMeanData$TSXexcess <- VARinput$TSXexcess - MeanData[4]


##Estimate the VAR(1) model
VAR1model <- VAR(ZeroMeanData,p = 1, type="none")
VAR1result<- summary(VAR1model)

#VAR2model<-VAR(VARinput,p=2,type="const")

##Coefficient of the VAR(1) model
B=Bcoef(VAR1model)
Nu=(diag(4)-B)%*%MeanData
Sigma = VAR1result$covres

##choleski decomposition
##Knowledge about Cholesky Transfermation
##http://blogs.sas.com/content/iml/2012/02/08/use-the-cholesky-transformation-to-correlate-and-uncorrelate-variables.html
CDSigma = chol(Sigma)

##Stationary Check
abseigen=abs(eigen(B)$value)
#----------------------------------------------------------------------------------------



#----------------------------------------------------------------------------------------

##Affine Term Structure model



CDNbonds_full <- read.table("bond yield monthly full.csv", header=TRUE, sep=",", na.strings=" ")



#data used to fit the affine term model
rankedCDNbonds=CDNbonds_full[,c(1:8)]
histdata=as.matrix(rankedCDNbonds)
VARaffine = as.matrix(VARinput)
Nhist=nrow(histdata)

#delta in pricing kernel
delta0=0
delta1=c(1,0,0,0)

#Estimating the parameter by minimize the value of sum of squares
Optimizefunc = function(sl){
  sl1matrix=matrix(data = c(sl[1:4],sl[5:8],0,0,0,0,B[4,]),4,4,byrow=TRUE)
  sl0vector=c(sl[9],sl[10],0,Nu[4])
  
  AnVector=numeric(181)
  BnMatrix=matrix(data=0,nrow=4,ncol=181)
  
  for(i in 2:181){
    BnMatrix[,i] = -delta1+t(B-sl1matrix)%*%BnMatrix[,i-1]
    AnVector[i]  = -delta0+AnVector[i-1]+t(BnMatrix[,i-1])%*%(Nu-sl0vector)+0.5*t(BnMatrix[,i-1])%*%t(CDSigma)%*%CDSigma%*%BnMatrix[,i-1]
  }
  
  yfit=matrix(data=0,nrow=Nhist,ncol=8)
  
  for (i in 1:Nhist){
	yfit[i, 1] = -1 / 12 * (AnVector[13] - sum(BnMatrix[, 13] * (VARaffine[i,])))
	yfit[i, 2] = -1 / 36 * (AnVector[37] + sum(BnMatrix[, 37] * (VARaffine[i,])))
	yfit[i, 3] = -1 / 60 * (AnVector[61] + sum(BnMatrix[, 61] * (VARaffine[i,])))
	yfit[i, 4] = -1 / 84 * (AnVector[85] + sum(BnMatrix[, 85] * (VARaffine[i,])))
	yfit[i, 5] = -1 / 120 * (AnVector[121] + sum(BnMatrix[, 121] * (VARaffine[i,])))
	yfit[i, 6] = -1 / 144 * (AnVector[145] + sum(BnMatrix[, 145] * (VARaffine[i,])))
	yfit[i, 7] = -1 / 168 * (AnVector[169] + sum(BnMatrix[, 169] * (VARaffine[i,])))
	yfit[i, 8] = -1 / 180 * (AnVector[181] + sum(BnMatrix[, 181] * (VARaffine[i,])))
  }
  
  sumsquared = sum((yfit-histdata)^2)
  return(sumsquared)

}

slbest=vector("double",10)
minfunc=100

for(i in 1:1){
  sl<-ucminf(par=runif(10,min=-0.0000,max=0.0000),fn=Optimizefunc,control = list(xtol = 1e-30))$par
  vfunc<-Optimizefunc(sl)
  if(vfunc<minfunc){
    slbest=sl
    B=vfunc
  }
}


L1 = matrix(data=c(slbest[1:4],slbest[5:8],0,0,0,0,B[4,]),4,4,byrow=TRUE)
L0 = c(slbest[9],slbest[10],0,Nu[4])





#slsaved<-sl
sl <- slbest
slsaved11 <- sl

slbest = vector("double", 12)
minfunc = 100

for (i in 1:100) {
	sl <- ucminf(par = runif(12, min = -0.02, max = 0.02), fn = Optimizefunc, control = list(xtol = 1e-10))$par
	vfunc <- Optimizefunc(sl)
	if (vfunc < minfunc) {
		slbest = sl
		minfunc = vfunc
	}
}


#slsaved<-sl
sl <- slbest
slsaved22 <- sl


slbest = vector("double", 12)
minfunc = 100


for (i in 1:100) {
	sl <- ucminf(par = runif(12, min = -0.02, max = 0.02), fn = Optimizefunc, control = list(xtol = 1e-10))$par
	vfunc <- Optimizefunc(sl)
	if (vfunc < minfunc) {
		slbest = sl
		minfunc = vfunc
	}
}


#slsaved<-sl
sl <- slbest
slsaved33 <- sl


slbest = vector("double", 12)
minfunc = 100
for (i in 1:100) {
	sl <- ucminf(par = runif(12, min = -0.035, max = 0.035), fn = Optimizefunc, control = list(xtol = 1e-10))$par
	vfunc <- Optimizefunc(sl)
	if (vfunc < minfunc) {
		slbest = sl
		minfunc = vfunc
	}
}


#slsaved<-sl
sl <- slbest
slsaved44 <- sl


slbest = vector("double", 12)
minfunc = 100
for (i in 1:100) {
	sl <- ucminf(par = runif(12, min = -0.035, max = 0.035), fn = Optimizefunc, control = list(xtol = 1e-10))$par
	vfunc <- Optimizefunc(sl)
	if (vfunc < minfunc) {
		slbest = sl
		minfunc = vfunc
	}
}


#slsaved<-sl
sl <- slbest
slsaved55 <- sl

slbest = vector("double", 12)
minfunc = 100
for (i in 1:100) {
	sl <- ucminf(par = runif(12, min = -0.035, max = 0.035), fn = Optimizefunc, control = list(xtol = 1e-10))$par
	vfunc <- Optimizefunc(sl)
	if (vfunc < minfunc) {
		slbest = sl
		minfunc = vfunc
	}
}


#slsaved<-sl
sl <- slbest
slsaved6 <- sl


slbest = vector("double", 12)
minfunc = 100
for (i in 1:100) {
	sl <- ucminf(par = runif(12, min = -0.035, max = 0.035), fn = Optimizefunc, control = list(xtol = 1e-10))$par
	vfunc <- Optimizefunc(sl)
	if (vfunc < minfunc) {
		slbest = sl
		minfunc = vfunc
	}
}


#slsaved<-sl
sl <- slbest
slsaved7 <- sl


sl <- slsaved1
#sl <-genoud(Optimizefunc,nvars=12,max=FALSE,pop.size=500,max.generations=200)
#sl <-sl$par

#slOriginal<-Original

##zero starting value optimization


SigmaLambda1 = matrix(data = c(sl[1:5],sl[6:10],0,0,0,0,0,B[4,],0,0,0,0,0),5,5,byrow=TRUE)
SigmaLambda0 = c(sl[11],sl[12],0,Nu[4],0)

An=numeric(181)
Bn=matrix(data=0,nrow=5,ncol=181)

for(i in 2:181){
  Bn[,i] = -delta1+t(B-SigmaLambda1)%*%Bn[,i-1]
  An[i]  = -delta0+An[i-1]+t(Bn[,i-1])%*%(Nu-SigmaLambda0)+0.5*t(Bn[,i-1])%*%t(CDSigma)%*%CDSigma%*%Bn[,i-1]
}

fitvalue=matrix(data=0,nrow=Nhist,ncol=9)

for (i in 1:Nhist){
	fitvalue[i, 1] = -(An[2] + sum(Bn[, 2] * (VARaffine[i,])))
	fitvalue[i, 2] = -1 / 3 * (An[4] + sum(Bn[, 4] * (VARaffine[i,])))
	fitvalue[i, 3] = -1 / 12 * (An[13] + sum(Bn[, 13] * (VARaffine[i,])))
	fitvalue[i, 4] = -1 / 24 * (An[25] + sum(Bn[, 25] * (VARaffine[i,])))
	fitvalue[i, 5] = -1 / 36 * (An[37] + sum(Bn[, 37] * (VARaffine[i,])))
	fitvalue[i, 6] = -1 / 60 * (An[61] + sum(Bn[, 61] * (VARaffine[i,])))
	fitvalue[i, 7] = -1 / 84 * (An[85] + sum(Bn[, 85] * (VARaffine[i,])))
	fitvalue[i, 8] = -1 / 120 * (An[121] + sum(Bn[, 121] * (VARaffine[i,])))
	fitvalue[i, 9] = -1 / 180 * (An[181] + sum(Bn[, 181] * (VARaffine[i,])))

}

x11(h=16,w=14,pointsize=12)
par(mfrow=c(4,2))
plot(fitvalue[,1],type="l",col="Red",xaxt='n',main="1 month maturity")
legend(x=300,y=0.12,legend=c("fitted value","historical"), col=c("red",  "black"),lty=1,lwd=1)
axis(side=1, at=seq(1, 302, by=60),labels=CDNbonds$Date[seq(1, 302, by=60)])
lines(BondsMonthly[, 1])


plot(fitvalue[, 3], type = "l", col = "Red", xaxt = 'n', main = "1 year maturity")
legend(x = 300, y = 0.12, legend = c("fitted value", "historical"), col = c("red", "black"), lty = 1, lwd = 1)
axis(side = 1, at = seq(1, 302, by = 60), labels = CDNbonds$Date[seq(1, 302, by = 60)])
lines(BondsMonthly[, 3])

plot(fitvalue[,4],type="l",col="Red",xaxt='n',main="2 year maturity")
legend(x=300,y=0.12,legend=c("fitted value","historical"), col=c("red",  "black"),lty=1,lwd=1)
axis(side = 1, at = seq(1, 302, by = 60), labels = CDNbonds$Date[seq(1, 302, by = 60)])
lines(BondsMonthly[, 4])

plot(fitvalue[,5],type="l",col="Red",xaxt='n',main="3 year maturity")
legend(x=300,y=0.12,legend=c("fitted value","historical"), col=c("red",  "black"),lty=1,lwd=1)
axis(side = 1, at = seq(1, 302, by = 60), labels = CDNbonds$Date[seq(1, 302, by = 60)])
lines(BondsMonthly[, 5])

plot(fitvalue[,6],type="l",col="Red",xaxt='n',main="5 year maturity")
legend(x=300,y=0.12,legend=c("fitted value","historical"), col=c("red",  "black"),lty=1,lwd=1)
axis(side = 1, at = seq(1, 302, by = 60), labels = CDNbonds$Date[seq(1, 302, by = 60)])
lines(BondsMonthly[, 6])

plot(fitvalue[,7],type="l",col="Red",xaxt='n',main="7 year maturity")
legend(x=300,y=0.12,legend=c("fitted value","historical"), col=c("red",  "black"),lty=1,lwd=1)
axis(side = 1, at = seq(1, 302, by = 60), labels = CDNbonds$Date[seq(1, 302, by = 60)])
lines(BondsMonthly[, 7])

plot(fitvalue[,8],type="l",col="Red",xaxt='n',main="10 year maturity")
legend(x=300,y=0.12,legend=c("fitted value","historical"), col=c("red",  "black"),lty=1,lwd=1)
axis(side = 1, at = seq(1, 302, by = 60), labels = CDNbonds$Date[seq(1, 302, by = 60)])
lines(BondsMonthly[, 8])

plot(fitvalue[,9],type="l",col="Red",xaxt='n',main="15 year maturity")
legend(x=300,y=0.12,legend=c("fitted value","historical"), col=c("red",  "black"),lty=1,lwd=1)
axis(side = 1, at = seq(1, 302, by = 60), labels = CDNbonds$Date[seq(1, 302, by = 60)])
lines(BondsMonthly[, 9])
#----------------------------------------------------------------------------------------



#----------------------------------------------------------------------------------------

#Forecast Plot
z0=as.numeric(VARinput[302,])
z0[1] = z0[1] + 0.0025 / 12
z0[2] = z0[2] + 0.0025 / 12


VAR1pred = predict(VAR1model, n.ahead = 480, ci = 0.95, dumvar = NULL)

x11(h=7,w=12,pointsize=12)
par(mfrow=c(1,1))
plot(VAR1pred)

#----------------------------------------------------------------------------------------



#----------------------------------------------------------------------------------------
##Simulation Result

ProjYear = 50
nscen=10000
projmonth= ProjYear*12


scenarios1 = SimulationFunc(ProjYear, nscen)
scenarios2 = SimulationFunc(ProjYear, nscen)
scenarios3 = SimulationFunc(ProjYear, nscen)
scenarios4 = SimulationFunc(ProjYear, nscen)
scenarios5 = SimulationFunc(ProjYear, nscen)
scenariosL = LambdaSimulationFunc (ProjYear,nscen)

scenarios = scenarios2[[1]]

StateScen1 = StateVariableScen(scenarios1)
StateScen2 = StateVariableScen(scenarios2)
StateScen3 = StateVariableScen(scenarios3)
StateScen4 = StateVariableScen(scenarios4)
StateScen5 = StateVariableScen(scenarios5)
StateScenL = StateVariableScen(scenariosL)

MValue1 = MValueFunc(scenarios1, StateScen1)
MValue2 = MValueFunc(scenarios2, StateScen2)
MValue3 = MValueFunc(scenarios3, StateScen3)
MValue4 = MValueFunc(scenarios4, StateScen4)
MValue5 = MValueFunc(scenarios5, StateScen5)
MValueL = MValueFunc(scenariosL, StateScenL)

YieldCurveScn2 = YieldCurveFunc(scenarios2)
ContiAnnualReturnScn2 = AnnualReturnFunc(StateScen2, YieldCurveScn2)
EffAnnualReturnScn2 = EffAnnualReturnFunc(ContiAnnualReturnScn2)
BondYield = BondYieldAnnual(YieldCurveScn2)

InfAnnual = EffAnnualReturnScn2$EffInfAnnual
BondAnnual = EffAnnualReturnScn2$EffBdAnnual
SrAnnual = EffAnnualReturnScn2$EffSrAnnual
DivAnnual = EffAnnualReturnScn2$EffDivAnnual
BondYield15 = BondYield$BondYield15
BondYield1m = BondYield$BondYield1m

#YieldCurveScnL = YieldCurveFunc(scenariosL)
#ContiAnnualReturnScnL = AnnualReturnFunc(StateScenL, YieldCurveScnL)
#EffAnnualReturnScnL = EffAnnualReturnFunc(ContiAnnualReturnScnL)
#BondYieldL = BondYieldAnnual(YieldCurveScnL)

#InfAnnual = EffAnnualReturnScnL$EffInfAnnual
#BondAnnual = EffAnnualReturnScnL$EffBdAnnual
#SrAnnual = EffAnnualReturnScnL$EffSrAnnual
#DivAnnual = EffAnnualReturnScnL$EffDivAnnual
#BondYield15 = BondYieldL$BondYield15
#BondYield1m = BondYieldL$BondYield1m



plot(SrAnnual[, 2], col = 'red', type = 'l', ylim = c(-0.2, 1))
lines(MValue2$MtAnnual[, 2], type = 'l')
lines(BondYield1m[,2],type = 'l')
#-------------------------------------------------------------------------------------------------



#---------------Some Plot-------------------------------------------------------------------------


x11(h=7,w=12,pointsize=12)
par(mfrow=c(1,1))
plot(ContiAnnualReturnScn2$BondAnnual[, 1], type = "l")

qqplot(MValue1[[3]][26,], MValue2[[3]][26,])
qqplot(MValue1[[3]][26,], MValue3[[3]][26,])
qqplot(MValue1[[3]][26,], MValue4[[3]][26,])
qqplot(MValue1[[3]][26,], MValue5[[3]][26,])
qqplot(MValue2[[3]][26,], MValue3[[3]][26,])
qqplot(MValue2[[3]][26,], MValue4[[3]][26,])
qqplot(MValue2[[3]][26,], MValue5[[3]][26,])
qqplot(MValue3[[3]][26,], MValue4[[3]][26,])
qqplot(MValue3[[3]][26,], MValue5[[3]][26,])
qqplot(MValue4[[3]][26,], MValue5[[3]][26,])



hist(MValue1[[3]][26,])
hist(MValue2[[3]][26,])
hist(MValue3[[3]][26,])
hist(MValue4[[3]][26,])
hist(MValue5[[3]][26,])
hist(MValueL[[3]][26,])
#-------------------------------------------------------------------------------------------------

