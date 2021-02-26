maledata <- c(180,  278,
  186,  277,
  206, 308,
  184,  290,
  177,  273,
  177,  284,
  176,  267,
  200,  281,
  191,  287,
  193,  271,
  212,  302,
  181,  254,
  195,  297,
  187,  281,
  190,  284,
  185,  282,
  195,  285,
  183,  276,
  202,  308,
  177,  254,
  177,  268,
  170,  260,
  186,  274,
  177,  272,
  178,  266,
  192,  281,
  204,  276,
  191,  290,
  178,  265,
  177,  275,
  284,  277,
  176,  281,
  185,  287,
  191,  295,
  177,  267,
  197,  310,
  199,  299,
  190,  273,
  180,  278,
  189,  280,
  194,  290,
  186,  287,
  191,  286,
  187,  288,
  186,  275)
T6.11 <- matrix(data=maledata,ncol=2,byrow=TRUE)

femaledata <- c(191,  284,
  197,  285,
  208,  288,
  180,  273,
  180,  275,
  188,  280,
  210,  283,
  196,  288,
  191,  271,
  179,  257,
  208,  289,
  202,  285,
  200,  272,
  192,  282,
  199,  280,
  186,  266,
  197,  285,
  201,  295,
  190,  282,
  209,  305,
  187,  285,
  207,  297,
  178,  268,
  202,  271,
  205,  285,
  190,  280,
  189,  277,
  211,  310,
  216,  305,
  189,  274,
  173,  271,
  194,  280,
  198,  300,
  180,  272,
  190,  292,
  191,  286,
  196,  285,
  207,  286,
  209,  303,
  179,  261,
  186,  262,
  174,  245,
  181,  250,
  189,  262,
  188,  258)
T5.12 <- matrix(data=femaledata,ncol=2,byrow=TRUE)

#Create Males data tables and scatterplot
males=data.frame(TailLength=T6.11[,1],WingLength=T6.11[,2])
plot(males,main="Scatterplot Males Wing Length vs. Males Tail Length")

#Male Data Summary
summary(males)
Mcov=cov(males)
malesbar<-apply(males,2,mean)

#Create Females data tables and scatterplot
females=data.frame(TailLength=T5.12[,1],WingLength=T5.12[,2])
plot(females,main="Scatterplot Female Wing Length vs. Female Tail Length")

#Female Data Summary
summary(females)
femalesbar<-apply(females,2,mean)
Fcov=cov(females)

#Create pooled data
n1=dim(males)[1]
n2=dim(females)[1]
p=dim(males)[2]
Spooled=((n1-1)*Mcov+(n2-1)*Fcov)/(n1+n2-2)

#Determine T2, obs. T2 and p-values
library(MASS)
T2=t(malesbar-femalesbar)%*%ginv((1/n1+1/n2)*Spooled)%*%(malesbar-femalesbar)
print(T2)
((n1+n2-2)*p)/(n1+n2-p-1)*qf(1-.05,df1=p,df2=n1+n2-p-1)
Pvalue<-1-pf(T2*(n1+n2-p-1)/(p*(n1+n2-2)),df1=p,df2=n1+n2-p-1)
print(Pvalue)

#Change outlier
newmales=males
newmales[31,1]=184

#NewMales Data Summary
summary(newmales)
NMcov=cov(newmales)
Nmalesbar<-apply(newmales,2,mean)

#Create New Pooled Data
NSpooled=((n1-1)*NMcov+(n2-1)*Fcov)/(n1+n2-2)

#Determine New T2, New obs. T2 and New p-values
library(MASS)
NT2=t(Nmalesbar-femalesbar)%*%ginv((1/n1+1/n2)*NSpooled)%*%(Nmalesbar-femalesbar)
print(NT2)
(((n1+n2-2)*p)/(n1+n2-p-1))*qf(1-.05,df1=p,df2=n1+n2-p-1)
NPvalue<-1-pf(NT2*(n1+n2-p-1)/(p*(n1+n2-2)),df1=p,df2=n1+n2-p-1)
print(NPvalue)

#Determine 95% confidence ellipse
Xbar=Nmalesbar-femalesbar
eigen(NSpooled)

y1=sqrt(259.86310)*(sqrt((1/n1+1/n2)*((n1+n2-2)*p)/(n1+n2-p-1)*qf(p=1-.05,df1=p,df2=n1+n2-p-1)))
y2=sqrt(32.67327)*(sqrt((1/n1+1/n2)*((n1+n2-2)*p)/(n1+n2-p-1)*qf(p=1-.05,df1=p,df2=n1+n2-p-1)))
print(y1)
print(y2)

#95% Confidence Interval
mu1.L=Xbar[1]-sqrt((((n1+n2-2)*p)/(n1+n2-p-1))*qf(1-.05,df1=p,df2=n1+n2-p-1))*sqrt((1/n1+1/n2)*NSpooled[1,1])
mu1.U=Xbar[1]+sqrt((((n1+n2-2)*p)/(n1+n2-p-1))*qf(1-.05,df1=p,df2=n1+n2-p-1))*sqrt((1/n1+1/n2)*NSpooled[1,1])

mu2.L=Xbar[2]-sqrt((((n1+n2-2)*p)/(n1+n2-p-1))*qf(1-.05,df1=p,df2=n1+n2-p-1))*sqrt((1/n1+1/n2)*NSpooled[2,2])
mu2.U=Xbar[2]+sqrt((((n1+n2-2)*p)/(n1+n2-p-1))*qf(1-.05,df1=p,df2=n1+n2-p-1))*sqrt((1/n1+1/n2)*NSpooled[2,2])

c(mu1.L,mu1.U)
c(mu2.L,mu2.U)


#Plot Confidence ellipse
library(ellipse)

plot(ellipse(NSpooled,centre=Xbar,t=sqrt((1/n1+1/n2)*((n1+n2-2)*p)/(n1+n2-p-1)*qf(1-.05,df1=p,df2=n1+n2-p-1))), type="l", main="Hook-billed Kites Data 95% Confidence Regions for
 Differenced Mean Vectors", xlim=c(-20,10),ylim=c(-10,10))
points(Xbar[1],Xbar[2])

#plot(ellipse(NSpooled,centre=Xbar,t=sqrt(((n1+n2-2)*p)/(n1+n2-p-1))*qf(.05,df1=p,df2=n1+n2-p-1)), type="l", main="Hook-billed Kites Data \n95% Density 
#Region for Normal Population \n and 95% Confidence Regions for Differenced Mean Vectors", xlim=c(-20,20),ylim=c(-10,10))
#points(Xbar[1],Xbar[2])

#library(car)
#rad <- sqrt(((n1+n2-2)*p)/(n1+n2-p-1)*qf(1-.05,df1=p,df2=n1+n2-p-1))
#ellipse(center=c(Xbar[1],Xbar[2]),shape=(1/n1+1/n2)*NSpooled,radius=rad,draw=TRUE)

#ellipse(center=c((xbar1.2[1,1]-xbar2[1,1]),(xbar1.2[2,1]-xbar2[2,1])),shape=(1/n1+1/n2)*spooled.2,radius=sqrt(c2.2),draw=TRUE)

lines(c(mu1.L,mu1.L),c(-20,mu2.U),lty=2)
lines(c(mu1.U,mu1.U),c(-20,mu2.U),lty=2)
lines(c(-25,mu1.U),c(mu2.L,mu2.L),lty=2)
lines(c(-25,mu1.U),c(mu2.U,mu2.U),lty=2)
