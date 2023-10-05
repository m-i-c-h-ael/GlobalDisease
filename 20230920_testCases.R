# Simulated Test Cases
rm(list=ls())
set.seed(03102023)


# 1,
 #Deaths: Three causes: Linearly increasing, decreasing and constant (from 1990)
  #No subgroups (location, age, sex)
  #No noise
 #Prevalence/Incidence: Constant
parDF= rbind.data.frame(
  c(200,  1,  0,  0,  0),
  c(200,  0, -10,  0,  0),
  c(200,  0,  0,  0,  0)
)
colnames(parDF)= c('a0','a1','a2','t0','noiseSD')

causes=     c('HIV/AIDS','Tuberculosis','Stroke')
locations=  c('lA','lA','lA')
ages=       c('25-30 years','25-30 years','25-30 years')   #string, not number
sexes=      c('M','M','M')

population= c(100000,100000,100000)  #population per location
names(population)= locations

t= 0:29
#years= 1990+t

DF= data.frame(matrix(NA,nrow=0,ncol=6))
for(measure in c('Incidence','Prevalence','Deaths')){
  if(measure == 'Deaths'){
    for (i in 1:dim(parDF)[1]){
      a0= parDF$a0[i]; a1= parDF$a1[i]; a2= parDF$a2[i]; t0= parDF$t0[i]; noiseSD= parDF$noiseSD[i]
      y_noNoise= a0 + a1*t + a2*ifelse(t>t0,t-t0,0)
      #noise1= rnorm(length(t),mean=0,sd=1)  #constant noise level
      noise2= rnorm(length(t),mean=0,sd=noiseSD*abs(y_noNoise))  #noise relative to data
      y= y_noNoise+noise2
      y[y<0]= 0   #make sure no negative values
      y[y>population[locations[i]]]= population[locations[i]]
    
      DF= rbind.data.frame( DF, 
        cbind.data.frame(
          matrix( rep(c('Deaths',locations[i],sexes[i],ages[i],causes[i],'Number'), each=length(t)), 
                  ncol=6,byrow=FALSE), t+1990,y,rep(NA,length(t)),rep(NA,length(t)))
      )
      #plot(t,y,ylim=c(0,1.1*max(y)))
    }
  } else {
    for (i in 1:dim(parDF)[1]){
      y= rep(1,length(t))
      DF= rbind.data.frame( DF, 
              cbind.data.frame(
                matrix( rep(c(measure,locations[i],sexes[i],ages[i],causes[i],'Number'), each=length(t)), 
                      ncol=6,byrow=FALSE), t+1990,y,rep(NA,length(t)),rep(NA,length(t)))
      )
    }
  }
}
colnames(DF)= c('measure','location','sex','age','cause','metric','year','val','upper','lower')

DF_percent= DF
pop_vec= population[DF$location]
DF_percent$val= DF$val / pop_vec * 100
DF_percent$metric= 'Percent'
DF= rbind.data.frame(DF,DF_percent)

write.csv(DF,file=paste('./testData/test1.csv'),row.names = FALSE)

ggplot(data=DF[DF$measure=='Deaths' & DF$metric=='Percent',], aes(x=year,y=val,color=age))+
  geom_line()+
  ggtitle('Deaths')+
  facet_wrap(~location+cause)
