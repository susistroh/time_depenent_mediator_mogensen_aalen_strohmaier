load("~/yourpath/csl_CP.Rda")

library(survival)
library(timereg)
library(mvna)


#####################################################################  
# Functions for the stochastic integral  in equation B.5
# input parameter: Data in counting process format, treatment indicator, start -, stop intervals, event indicator
# output: stochastic integral and increments of the stochastic integral, unique event times

f_stoch_int<-function(data_CP, treat_ind, start_t, stop_t, event_ind) {
  # group a=1
  treat<-data_CP[data_CP[[treat_ind]]==1,]
  # group a=0
  placebo<-data_CP[data_CP[[treat_ind]]==0,]
  # Cox model with time-depend. mediator in group a=0 (equation B.1)
  cox.pla <- coxph(Surv(placebo[[start_t]], placebo[[stop_t]], placebo[[event_ind]]==1)~age+sex+trans_prot,
                   data = placebo, method = "breslow")
  # storing the linear predictor
  pla.lp<-cox.pla$linear.predictors
  data.pla=cbind(placebo,pla.lp=pla.lp)
  # event times
  obestime_x= data.pla[[stop_t]][data.pla[[event_ind]]==1]
  obstimes = sort(unique(obestime_x))
  
  K<-length(obstimes)
  h<-rep(0,K)
  
  # predict and store the linear predictor for the group a=1 based on the the model for a=0;
  treat.pred.lp<-predict(cox.pla,newdata=treat,type="lp")
  data.treat=cbind(treat,treat.pred.lp=treat.pred.lp)
  # loop over event times 
  for(i in 1:K){
    # risk set in group a=0
    risk.set.pla<-subset(data.pla,data.pla[[start_t]]<obstimes[i] & data.pla[[stop_t]]>=obstimes[i])
    # linear predictors for individuals in the risk set (a=0)
    lp_test<-risk.set.pla$pla.lp
    # risk
    risk<-exp(lp_test)
    # risk set in group a=1 for event times of group a=0
    risk.set.treat<-subset(data.treat,data.treat[[start_t]]<obstimes[i] & data.treat[[stop_t]]>=obstimes[i])
    # liner predictor for individuals in that risk set
    lp_treat_at_risk<-risk.set.treat$treat.pred.lp
    # corresponding risk 
    risk_treat<-exp(lp_treat_at_risk)
    # components of the increments in the stochastic integral 
    numerator_treat=sum(risk_treat)*sum(data.pla[[event_ind]][data.pla[[stop_t]]==obstimes[i]])
    denumerator=dim( risk.set.treat)[1]*sum(risk)
    # increment 
    h[i]=numerator_treat/denumerator
    #h
    
  }
  int=cumsum(h)
  list(incr_h=h, int=int, obstimes= obstimes)
}

# evaluate the stochastic integral 
pla.stoch.int<-f_stoch_int(data_CP=csl_CP, treat_ind="treat", start_t="lt", stop_t="rt", event_ind="e_ind")
#####################################################################
# step -function to evaluate the stochastic integral at any (sensible) time point
eval_stoch_int<-stepfun(pla.stoch.int$obstimes, c(0,pla.stoch.int$int), f=0)
# evaluate at event times
event_times<-sort(unique(csl_CP$rt[csl_CP$e_ind==1]))
event.eval.stoch.int<-eval_stoch_int(event_times)


#####################################################################
# function to evaluate the Nelson Aalen estimator at any time point 
# input: data with one line per subject, event indicator, event time, desired evalution times
# output: evaluation time points, value of the NA estimator at these times 

# take data in counting process format, create variable with one line per subject 
#head(csl_CP)
csl_oo<-csl_CP[csl_CP$numbering==csl_CP$max_ob,]


eval_Nel_Aal<-function(data_oo, id, event_ind, event_t, predict_times){
  # use data set with one observation per subject
  # create a dataframe needed for the mvna package 
  mvna.data<-data.frame(id=data_oo[[id]],
                        from=data_oo[[event_ind]],
                        to=data_oo[[event_ind]],
                        time= data_oo[[event_t]])
  
  mvna.data$from<-c(rep(0,dim(data_oo)[1]))
  
  tra2<-matrix(ncol=2,nrow=2,FALSE)
  tra2[1,2]<-TRUE
  
  # recode 0 in to to cens
  mvna.data$to[mvna.data$to==0]<-"cens"
  
  # obtain Nelson Aalen estimate  
  na.estimate<-mvna(mvna.data,c("0","1"),tra2,"cens")
 # predict NA value at the specified times
  predict.ob<-predict(na.estimate,times=predict_times)
  times<-predict.ob$`0 1`$time
  na.predictions<-predict.ob$`0 1`$na
  
  list(pred.times=times,na.predictions=na.predictions)
}

# evaluate the NA based on group a=1 at all event times 
treat_oo<-csl_oo[csl_oo$treat==1,]
event.eval.Nel.Aal<-eval_Nel_Aal(data_oo=treat_oo, id="id", event_ind = "e_ind", event_t="rt", predict_times = event_times)

##################################################
# obtain R(t)
roh=event.eval.Nel.Aal$na.predictions-event.eval.stoch.int


##################################################
# Kaplan Meier group a=1; evaluated at all event times 
KM_treat.fit <- survfit(Surv(rt,e_ind)~1, data=treat_oo)
KM_treat<-summary(KM_treat.fit, times = event_times)
est_KM_treat<-KM_treat$surv

# Kaplan Meier group a=0; evaluated at all event times 
placebo_oo<-csl_oo[csl_oo$treat==0,]
KM_placebo.fit <- survfit(Surv(rt,e_ind)~1, data=placebo_oo)
KM_placebo<-summary(KM_placebo.fit, times = event_times)
est_KM_placebo<-KM_placebo$surv

# KM ratio
KM_ratio<-est_KM_treat/est_KM_placebo

#Effects on the relative survival scale 
# direct effect on the relative survival scale 
SDE<-exp(-roh)

#indirect effect 
SIE=exp(roh)*KM_ratio

# total effect 
STE=SDE*SIE


####################################################################
####################################################################
# Bootstrap confindence intervall
####################################################################
####################################################################
start.time.bootstrap <- Sys.time()
set.seed(123)
N=100 # number of boostrap samples
effects=c() # object ot store the effects
#Bootstrap loop:
#for(i in 1:N){
for(i in 1:N){
  print(i)

  # sample individuals 
  boot.id<-unique(csl_CP$id)
  boot.sample<-sample(boot.id,replace=T)
  bootset<-csl_CP[csl_CP$id %in% boot.sample,]
 
  #event times 
  boot_event_times<-sort(unique(bootset$rt[bootset$e_ind==1]))
  #apply stoch int 
  boot_stoch_int<-f_stoch_int(data_CP=bootset, treat_ind="treat", start_t="lt", stop_t="rt", event_ind="e_ind")
  # function to evaluate the stochastic integral at any time 
  boot.eval_stoch_int<-stepfun(boot_stoch_int$obstimes, c(0,boot_stoch_int$int), f=0)
  # evaluate at event times
  boot.event.eval.stoch.int<-boot.eval_stoch_int(boot_event_times)
  
  # bootset with one observation per subject 
  bootset_oo<- bootset[bootset$numbering== bootset$max_ob,]
  #  group a=1
  boot.treat_oo<-bootset_oo[bootset_oo$treat==1,]
  #  group a=0
  boot.placebo_oo<-bootset_oo[bootset_oo$treat==0,]
  
  
  # function to evaluate Nelson Aalen estimator at any time point 
  boot.event.eval.Nel.Aal<-eval_Nel_Aal(data_oo=boot.treat_oo, id="id", event_ind = "e_ind", event_t="rt", 
                                        predict_times = boot_event_times)
  
  # evaluate R(t)
  boot.roh= boot.event.eval.Nel.Aal$na.predictions-boot.event.eval.stoch.int
  
  # direct effect on the relative survival scale 
  boot.SDE<-exp(-boot.roh)
  
  # indirect effect on the relative survival scale 
  # Kaplan Meier group a=1; evaluated at all event times 
  boot.KM_treat.fit <- survfit(Surv(rt,e_ind)~1, data=boot.treat_oo)
  boot.KM_treat<-summary(boot.KM_treat.fit, times =boot_event_times)
  boot.est_KM_treat<-boot.KM_treat$surv
  
  # Kaplan Meier group a=0; evaluated at all event times 
  boot.KM_placebo.fit <- survfit(Surv(rt,e_ind)~1, data=boot.placebo_oo)
  boot.KM_placebo<-summary(boot.KM_placebo.fit, times = boot_event_times)
  boot.est_KM_placebo<-boot.KM_placebo$surv
  
  boot.KM_ratio<-boot.est_KM_treat/boot.est_KM_placebo
  #indirect effect 
  boot.SIE=exp(boot.roh)*boot.KM_ratio
  
  # total effect on the relative survival scale (product SDE*SIE)
  boot.STE=boot.SDE*boot.SIE
  
  #estimates
  boot.estimates<-cbind(boot_event_times, boot.SDE, boot.SIE, boot.STE)
  
  effects<-rbind(effects, boot.estimates)
  #output
  list(effects=effects)
}

end.time.bootstrap <- Sys.time()
time.taken.bootstrap <- end.time.bootstrap - start.time.bootstrap
time.taken.bootstrap


#####################################################################

##############################################################
obstimes.imput = sort(unique(csl_CP$rt[csl_CP$e_ind==1] ))

# objects to store the estimates:
direct.SDE.l=c()
direct.SDE.u=c()

indirect.SIE.l=c()
indirect.SIE.u=c()

total.l=c()
total.u=c()

for (j in 1:length(obstimes.imput)){
  # select the corresponding columns form the 'lm.effect' object,
  #which is the output from the bootstrap loop
  subseti=subset(effects,effects[,1]==obstimes.imput[j])
  
  temp.direct.SDE.l=quantile(subseti[,2],0.025,names=FALSE)
  temp.direct.SDE.u=quantile(subseti[,2],0.975,names=FALSE)
  
  # pointwise CI for the indirect effectss
  temp.indirect.SIE.l=quantile(subseti[,3],0.025,names=FALSE)
  temp.indirect.SIE.u=quantile(subseti[,3],0.975,names=FALSE)
  
  # pointwise CI for the total effect
  temp.total.l=quantile(subseti[,4],0.025,names=FALSE)
  temp.total.u=quantile(subseti[,4],0.975,names=FALSE)
  
  # combine pointswise values to one vector:
  direct.SDE.l=rbind(direct.SDE.l,temp.direct.SDE.l)
  direct.SDE.u=rbind(direct.SDE.u,temp.direct.SDE.u)
  
  
  indirect.SIE.l=rbind(indirect.SIE.l, temp.indirect.SIE.l)
  indirect.SIE.u=rbind(indirect.SIE.u, temp.indirect.SIE.u)
  
  total.l=rbind(total.l,temp.total.l)
  total.u=rbind(total.u,temp.total.u)
}


event_times = sort(unique(csl_CP$rt[csl_CP$e_ind==1] ))


##################################################
# basic plots 
##################################################
par(mfrow=c(1,3))
# direct effect of treatment on outcome
plot(obstimes.imput,SDE,type="l", 
     main='Direct effect of treatment',
     xlab='Years since randomisation', xlim=c(-0.001,8),
     ylab='Relative Survival', ylim=c(0.7,1.5),cex.lab = 1.5)
lines(obstimes.imput,rep(1,length(obstimes.imput)))
lines(obstimes.imput,direct.SDE.l,type='l',col="grey")
lines(obstimes.imput,direct.SDE.u,type='l',col="grey")

plot(obstimes.imput,SIE,type="l",
     main='Indirect effect through prothrombin',
     xlab='Years since randomisation', xlim=c(-0.001,8),
     ylab='Relative Survival', ylim=c(0.9,1.5),cex.lab = 1.5)
lines(obstimes.imput,rep(1,length(obstimes.imput)))
lines(obstimes.imput,indirect.SIE.l,type='l',col="grey")
lines(obstimes.imput,indirect.SIE.u,type='l',col="grey")

plot(obstimes.imput,STE,type="l", 
     main='Total effect of treatment',
     xlab='Years since randomisation', xlim=c(-0.001,8),
     ylab='Relative Survival', ylim=c(0.8,1.5),cex.lab = 1.5)
lines(obstimes.imput,rep(1,length(obstimes.imput)))
lines(obstimes.imput,total.l,type='l',col="grey")
lines(obstimes.imput,total.u,type='l',col="grey")

###################################################

