#Find the best way to model a timecourse which is made up of two linear effects
 #the second effect starts at a given time after the first

#an alternative model would be that the effect of a2 depends on level at the
 #given time point

set.seed(01092023)

a0=2
a1=1
a2=-1.5

t= 0:100
t0= 40
#t0= rep(40,length(t))  #time when effect 2 starts

y_noNoise= a0 + a1*t + a2*ifelse(t>t0,t-t0,0)
#noise1= rnorm(length(t),mean=0,sd=1)  #constant noise level
noise2= rnorm(length(t),mean=0,sd=0.5*y_noNoise)  #noise relative to data
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
source('G:/My Drive/20230815_GlobBurdDisease/DBDA2E-utilities_mod.R')  #modified to fit markdown


modStr= 
'data{
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
    y_z[i] ~ dnorm( a0 + a1*t_z[i] + a2*(t_z[i]-t0_z), 1/sig_y^2 ) #homog. variance
  }
   #priors
   a0 ~ dgamma(1,1)      #orig. scale mean= mean_y, ESS= 1
   a1 ~ dgamma(1,10)     #orig. scale mean= sd_y/10, ESS= 10
   a2 = -a2_neg         #ensure a2 is negative on same scale as a1
   a2_neg~ dgamma(1,10)
   t0_z= t0-t_half
   t0 ~ dunif(0,2*t_half)  #t0 anywhere between 0 and max(t)
   sig_y ~ dgamma(1/10,1)
   
   #re-transform
   b0= a0*sd_y+mean_y
   b1= a1*sd_y
   b2= a2*sd_y
}'

dataList= list(
  t=t,
  y=y,
  Ntotal= length(y)
)

jagMod= jags.model(file=textConnection(modStr),data=dataList,n.chains=3)  
cS= coda.samples(model=jagMod,variable.names = c('b0','b1','b2','t0','sig_y'),n.iter=10000)

#diagnostics
for(i in 1:length(varnames(cS))){
  diagMCMC(codaObject= cS,parName=varnames(cS)[i]  )
}

# MCMC results
DF1= data.frame(cS[[1]])
head(DF1)
par(mfrow=c(2,2))
hist(DF1$b0,col='skyblue',breaks=50)
lines(x= HDIofMCMC(DF1$b0,credMass=.95), y=c(0,0),lwd=5)
hist(DF1$b1,col='skyblue',breaks=50)
hist(DF1$b2,col='skyblue',breaks=50)
hist(DF1$t0,col='skyblue',breaks=50)

par(mfrow=c(1,1))
plot(DF1$b1,DF1$b2,col='skyblue')
