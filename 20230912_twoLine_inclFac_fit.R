twoLine_fit= function(t, y){
  #fit combination of lines with slope>=0 starting at t=0 and line with slope<=0
   #starting at t>=0 to y
  # source('G:/My Drive/20230815_GlobBurdDisease/DBDA2E-utilities_mod.R')  #modified to fit markdown
  # today= gsub('-','',Sys.Date())
  
  #to do: its probably better to use a model appropriate for bounding at lower end; e.g.exponential + Poisson
  
  t_orig= t
  t= t_orig-min(t_orig)  #make sure t starts at 0
  
  if(all(y== mean(y))){   #no variation
    df= data.frame(matrix(data= c(mean(y),rep(0,6)),nrow=1,ncol=7))
    colnames(df)= c('b0','b1','b2','t0','iF1','iF2','t0')
    return(df)
  }
  
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
      y_z[i] ~ dnorm( a0 + iF1*a1*t_z[i] + iF2*a2*t2[i], 1/sig_y^2 ) #homog. variance
    }
     #priors
     a0 ~ dnorm(0,1)      #orig. scale mean= mean_y
     iF1 ~ dbern(0.1)   #should a1 be included?; bias against it
     iF2 ~ dbern(0.1)   #should a2 be included?; bias against it
     
     # a1 can only be positive, a2 only negative
     a1 ~ dgamma(.05,2)     #orig. scale mean= 0.025*sd_y, ESS= 2
     a2 = -a2_neg         #ensure a2 is negative on same scale as a1
     a2_neg~ dgamma(.05,2)
     
     t0_z= t0-t_half
     t0 ~ dunif(0,2*t_half)  #t0 anywhere between 0 and max(t)
     sig_y ~ dgamma(1/10,1) #on transformed scale
     
     #re-transform
     b0= a0*sd_y - a1*sd_y*t_half + mean_y
     b1= iF1*a1*sd_y  #if indicator is 0, parameter is 0
     b2= iF2*a2*sd_y
  }'

  dataList= list( t= t, y= y, Ntotal= length(y) )
  
  jagMod= jags.model(file=textConnection(modStr),data=dataList,n.chains=2,quiet=TRUE)  
  cS= coda.samples(model=jagMod,variable.names = c('b0','b1','b2','t0','iF1','iF2'),n.iter=20000)
  
  # print("Diagnostic plots saved to file.")
  # for(i in 1:length(varnames(cS))){
  #   diagMCMC(codaObject= cS, parName=varnames(cS)[i],
  #            saveName = paste(getwd(),'/plots/',today,'_diagPlots_',sep='') )
  #   dev.off()
  # }
  
  DF1= data.frame(cS[[1]])
  DF1$t0 = DF1$t0+min(t_orig)
  #DF1$iF1= factor(DF1$iF1, levels=c(0,1)) #to make sure both levels are listed even if only one observed
  #DF1$iF2= factor(DF1$iF2, levels=c(0,1))
  
  return( DF1 )
}
