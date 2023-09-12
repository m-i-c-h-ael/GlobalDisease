#Find the best way to model a timecourse which is made up of two linear effects
 #the second effect starts at a given time after the first

rm(list=ls())
cat('\014')
graphics.off()

library('tidyverse')

#an alternative model would be that the effect of a2 depends on level at the
 #given time point

set.seed(01092023)

a0=2
a1=1
a2= -1.5

t= 0:100
t0= 40
#t0= rep(40,length(t))  #time when effect 2 starts

y_noNoise= a0 + a1*t + a2*ifelse(t>t0,t-t0,0)
#noise1= rnorm(length(t),mean=0,sd=1)  #constant noise level
noise2= rnorm(length(t),mean=0,sd=0.1*abs(y_noNoise))  #noise relative to data
y= y_noNoise+noise2
y[y<0]= 0   #make sure no negative values
DF= cbind.data.frame(t=t,y=y)

plot(t,y,ylim=c(0,1.1*max(y)))

twoLine_fun= function(t, pms, obs){  #ssq minimization between observed data
                                       #and two-line model
  a0= pms[1]; a1=pms[2]; a2=pms[3]; t0= pms[4]
  y= a0 + a1*t + a2*ifelse(t>t0,t-t0,0)
  
  ssq= sum( (obs-y)^2 )
  return(ssq)
}

pms= c(a0,a1,a2,t0)
twoLine_fun(t,pms,obs=y)

nlmOUT= nlm(f=twoLine_fun,p=c(a0=0,a1=1,a2=1,t0=max(t)/2),t=t,obs=y)
est= nlmOUT$estimate

lines(t, est[1] + est[2]*t + est[3]*ifelse(t>est[4],t-est[4],0),col='blue')


###******************************************************####
#extimate parameters based on Bayesian model
library('rjags')
source('H:/My Drive/20230815_GlobBurdDisease/DBDA2E-utilities_mod.R')  #modified to fit markdown


modStr= 
'
data{
  #Ntotal= length(y)
   #transform data by rescaling y and shifting t to t1/2
    t_half= (max(t)-min(t))/2
    mean_y= mean(y)
    sd_y= sd(y)
    for(i in 1:Ntotal){
      t_z[i]= t[i]-t_half
      y_z[i]= (y[i]-mean_y)/sd_y
    }
}

model{
  for(i in 1:Ntotal){
    t2[i]= ifelse(t_z[i]>t0_z,t_z[i]-t0_z,0)  #only subtracted if t_z[i]>t0_z
    y_z[i] ~ dnorm( a0 + a1*t_z[i] + a2*t2[i], 1/sig_y^2 ) #homog. variance
  }
   #priors
   a0 ~ dnorm(0,1)      #orig. scale mean= mean_y
   a1 ~ dgamma(.05,2)     #orig. scale mean= 0.025*sd_y, ESS= 2
   a2 = -a2_neg         #ensure a2 is negative on same scale as a1
   a2_neg~ dgamma(.05,2)
   t0_z= t0-t_half
   t0 ~ dunif(0,2*t_half)  #t0 anywhere between 0 and max(t)
   sig_y ~ dgamma(1/10,1) #on transformed scale
   
   #re-transform
   b0= a0*sd_y - a1*sd_y*t_half + mean_y
   b1= a1*sd_y
   b2= a2*sd_y
}'

#y=c(); t=c()  #if you want priors
dataList= list(
  t= t,
  y= y,
  Ntotal= length(y)
)

jagMod= jags.model(file=textConnection(modStr),data=dataList,n.chains=3)  
cS= coda.samples(model=jagMod,variable.names = c('b0','b1','b2','t0'),n.iter=10000)

#diagnostics
for(i in 1:length(varnames(cS))){
  diagMCMC(codaObject= cS,parName=varnames(cS)[i]  )
}

# MCMC results
DF1= data.frame(cS[[1]])
head(DF1)
par(mfrow=c(2,2))
hist(DF1$b0,col='skyblue',breaks=50)
lines(x=HDIofMCMC(DF1$b0,credMass=.95),y=c(0,0),lwd=5)

hist(DF1$b1,col='skyblue',breaks=50)
lines(x=HDIofMCMC(DF1$b1,credMass=.95),y=c(0,0),lwd=5)

hist(DF1$b2,col='skyblue',breaks=50)
lines(x=HDIofMCMC(DF1$b2,credMass=.95),y=c(0,0),lwd=5)

hist(DF1$t0,col='skyblue',breaks=50)
lines(x=HDIofMCMC(DF1$t0,credMass=.95),y=c(0,0),lwd=5)


par(mfrow=c(1,1))
plot(DF1$b1,DF1$b2,col='skyblue')


##***********************************************###
# Test prediction power for some data

a0=2
t= 0:100
t0= 40
coefDF= expand.grid(a1=c(0.1,1,10), a2=c(-.01,-.1,-0.2,-0.5))   #a1>0, a2<0


dev.new(width=9,height=9,noRStudioGD = TRUE)
par(mfrow=c(nrow=ceiling(dim(coefDF)[1]/3),ncol=3))
DF_result= data.frame(matrix(NA,nrow=dim(coefDF)[1],ncol=6))
colnames(DF_result)= c('b1','b2','b1_correct','b2_correct','b1_CV', 'b2_CV')
for(i in 1:dim(coefDF)[1]){
  a1= coefDF$a1[i]; a2= coefDF$a2[i]
  DF_result[i,c(1,2)]= c(a1,a2)
  y_noNoise= a0 + a1*t + a2*ifelse(t>t0,t-t0,0)
  noise= rnorm(length(t),mean=0,sd=0.1*abs(y_noNoise))  #noise relative to data
  y= y_noNoise+noise
  y[y<0]= 0   #make sure no negative values
  plot(t,y,ylim=c(0,1.1*max(y)),main=paste('b1:',a1,', b2:',a2)) #change naming: 'a' are simulated, 'b' estimated params
  
  dataList= list( t= t, y= y, Ntotal= length(y) )
  
  jagMod= jags.model(file=textConnection(modStr),data=dataList,n.chains=1)  
  cS= coda.samples(model=jagMod,variable.names = c('b0','b1','b2','t0'),n.iter=10000)
  
  DF1= data.frame(cS[[1]])
  HDIlower_b1= HDIofMCMC(DF1$b1,credMass=.95)[1]  #lower bound of HDI -> diff. from 0?
  DF_result$b1_correct[i]= ifelse(a1>0 & HDIlower_b1>0 | a1==0 & HDIlower_b1<=0,1,0) #b1 correctly classified
  DF_result$b1_CV[i]= abs(mean(DF1$b1)-a1)/abs(a1)  #error of mean estimate
  HDIupper_b2= HDIofMCMC(DF1$b2,credMass=.95)[2]  #lower bound of HDI -> diff. from 0?
  DF_result$b2_correct[i]= ifelse(a2<0 & HDIupper_b2<0 | a2==0 & HDIupper_b2>=0,1,0) #b2 correctly classified
  DF_result$b2_CV[i]= abs(mean(DF1$b2)-a2)/abs(a2)  #error of mean estimate
  
  lines(t, mean(DF1$b0) + mean(DF1$b1)*t + mean(DF1$b2)*ifelse(t>mean(DF1$t0),t-mean(DF1$t0),0),col='blue')
}
#dev.off()
dev.new(width=4,height=7,noRStudioGD = TRUE)

#manually reformat
DF_result_long= cbind.data.frame(
  coef= c(DF_result$b1,DF_result$b2),
  x= factor(c(rep(1,dim(DF_result)[1]),rep(2,dim(DF_result)[1])),levels=c(1,2)),
  y= rep(1:dim(DF_result)[1],times=2),
  correct= factor(c(DF_result$b1_correct,DF_result$b2_correct),levels=c(0,1)),
  CV= c(DF_result$b1_CV,DF_result$b2_correct_CV),
  log10_recipCV= log(1/c(DF_result$b1_CV,DF_result$b2_correct_CV),10) #log10 of recip. CV
)

ggplot()+
  geom_point(data=DF_result_long,aes(x=x,y=y,color=correct,size=log10_recipCV))+  
  scale_color_manual(limits=as.character(c(0,1)),values=c('red','green'),labels=c('no','yes'))+
  scale_size_continuous(breaks=log(1/c(0.01,0.1,0.3),10),labels=c(0.01,0.1,0.3),name='CV')+
  geom_text(data=DF_result_long,aes(x=as.numeric(x)-0.2,y=y, label=coef))+
  theme_minimal()+
  scale_x_discrete(breaks=c(1,2),labels=c('b1','b2'))+
  theme(axis.title= element_blank(),axis.text.y=element_blank())
