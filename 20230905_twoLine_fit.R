twoLine_fit= function(t, y){
  #fit combination of lines with slope>=0 starting at t=0 and line with slope<=0
   #starting at t>=0 to y
  
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
     
     #Prior for a1 has tail in the right and for a2 in the left direction ->
      #it is clear which is the positive and which the negative slope
     #To allow parameter estimates of 0, shift them so that the mean of each
      #prior is at 0
     #(cases with only small slopes will have high autocorrelation)
     a1_shift ~ dgamma(.05,2)     #orig. scale mean= 0.025*sd_y, ESS= 2
     a1= a1_shift-0.025  #mean at 0
     a2 = -a2_neg_shift+0.025         #ensure a2 is negative on same scale as a1
     a2_neg_shift~ dgamma(.05,2)
     
     t0_z= t0-t_half
     t0 ~ dunif(0,2*t_half)  #t0 anywhere between 0 and max(t)
     sig_y ~ dgamma(1/10,1) #on transformed scale
     
     #re-transform
     b0= a0*sd_y - a1*sd_y*t_half + mean_y
     b1= a1*sd_y
     b2= a2*sd_y
  }'

  dataList= list( t= t, y= y, Ntotal= length(y) )
  
  jagMod= jags.model(file=textConnection(modStr),data=dataList,n.chains=1,quiet=TRUE)  
  cS= coda.samples(model=jagMod,variable.names = c('b0','b1','b2','t0'),n.iter=10000)
  
  DF1= data.frame(cS[[1]])
  HDIlower_b1= HDIofMCMC(DF1$b1,credMass=.95)[1]  #lower bound of HDI -> diff. from 0?
  b1_correct= ifelse(a1>0 & HDIlower_b1>0 | a1==0 & HDIlower_b1<=0,1,0) #b1 correctly classified
  b1_relErr= abs(mean(DF1$b1)-a1)/abs(a1+1E3)  #relative error of mean estimate; w. flooring
  HDIupper_b2= HDIofMCMC(DF1$b2,credMass=.95)[2]  #lower bound of HDI -> diff. from 0?
  b2_correct= ifelse(a2<0 & HDIupper_b2<0 | a2==0 & HDIupper_b2>=0,1,0) #b2 correctly classified
  b2_relErr= abs(mean(DF1$b2)-a2)/abs(a2+1E3)  #relative error of mean estimate
  
  return( list(MCMC= DF1,
               eval= c(b1_correct,b2_correct,b1_relErr,b2_relErr)) )
}
