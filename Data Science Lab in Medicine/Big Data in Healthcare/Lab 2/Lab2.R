### -----------------------------
### LAB 2: Cox regression model
### -----------------------------

### Exercise 1

#   0. Read the data and the call the "survival" package

ova<-read.table(file.choose(), na.strings=".",header=T,row.names=NULL)
names(ova)
head(ova)

ova$survtime<-ova$survtime/365.25
ova$turesc<-ifelse(ova$tures<=4, 0, ifelse(ova$tures>=7, 2, 1)) 
ova$karc<-ifelse(ova$kar<=70, 0, 1) 
ova$histoc<-ifelse(ova$histo==1 | ova$histo==7, 0,1)

# 1. Estimate the KM survival function for each treatment

library(survival)
fit<-survfit(Surv(survtime, surv) ~ treat, data=ova)
plot(fit[treat=1]$time, fit[treat=1]$surv, type='s',col=1, ylim=c(0,1) , xlim=c(0,13), xlab='time', ylab='survival probability')
lines(fit[treat=2]$time, fit[treat=2]$surv, type='s', lwd=2,col=2)
lines(fit[treat=3]$time, fit[treat=3]$surv, type='s', lwd=3,col=3)
legend("top",legend=c("CAP","CP","P"),col=c(1,2,3),lwd=c(1,2,3),lty=1)


# 2. Test the null hypothesis of equality of the survival curves, by using nonparametric tests 

survdiff(Surv(survtime, surv) ~ treat,data=ova)

# 3. Test the null hypothesis of equality of the survival curves (specific for treatment), by using nonparametric hypothesis tests stratifying for the prognostic factors

survdiff(Surv(survtime, surv) ~ treat + strata(karc), data=ova)
survdiff(Surv(survtime, surv) ~ treat + strata(turesc), data=ova)

# 4. Fit the Cox model including the treat covariate

model<-coxph( Surv(survtime, surv) ~ factor(treat), data = ova)
summary(model)

treat.r<-relevel(factor(ova$treat),ref=3)
model<-coxph(formula = Surv(survtime, surv) ~ treat.r, data = ova)
summary(model)

# 5. Fit the Cox model including the treat covariate and adjusting for the other prognostic factors

model5<-coxph(Surv(survtime, surv) ~ factor(treat)+factor(figo) +factor(karc)+factor(turesc)+factor(histoc),data=ova)
summary(model5)
anova(model,model5)

# 6. From the Cox model in 5., estimate and plot:
#1. the baseline survival curve 
#2. the survival curve for a patient treated with CP having the following characteristics: FIGO= III, Karnofsky Index =80, residual tumor size=3 cm, histotype=3

model<-coxph(Surv(survtime, surv) ~ factor(treat)+factor(figo) +factor(karc)+factor(turesc)+factor(histoc),data=ova)
bas<-basehaz(model,centered=FALSE)
bas.surv<- exp( -bas[,1] )
plot(bas$time, bas.surv, type='s',col=1, ylim=c(0,1) , xlim=c(0,13),lty=2, xlab='time', ylab='survival probability')

# CP:treat =2, FIGO= III, Karnofsky Index =80 -> karc=1 & residual tumor size=3 cm ->turesc=1, histotype=3 -> histoc=1
surv_patientCP<-bas.surv^(exp(0.1858+ 0.4119 -0.3931 + 0.7019 +0.2385 ))
lines(bas$time, surv_patientCP, type='s', lwd=2,col=2)

# same result, using survfit:
fit<-survfit(model,newdata=data.frame(treat=2,figo=3,karc=1,turesc=1,histoc=1))
surv_patientCP<-fit$surv
lines(bas$time, surv_patientCP, type='s', lwd=2,col=2)


# 7. Compare the survival curve in 6. with the survival curves of two patients with the same characteristics, but treated with P and CAP, respectively

plot(bas$time, bas.surv, type='s',col=1, ylim=c(0,1) , xlim=c(0,13),lty=2, xlab='time', ylab='survival probability')
lines(bas$time, surv_patientCP, type='s', lwd=2,col=2)

fit<-survfit(model,newdata=data.frame(treat=1,figo=3,karc=1,turesc=1,histoc=1))
surv_patientCAP<-fit$surv
lines(bas$time, surv_patientCAP, type='s', lwd=2,col=1)

fit<-survfit(model,newdata=data.frame(treat=3,figo=3,karc=1,turesc=1,histoc=1))
surv_patientP  <-fit$surv
lines(bas$time, surv_patientP, type='s', lwd=2,col=3)
legend("top",legend=c("baseline","CAP","CP","P"),col=c(1,1,2,3),lwd=c(1,1,2,3),lty=c(2,1,1,1))

# 8. Fit the model of point 5. and check the PH assumption for karc and for histoc by the schoenfeld residuals plot
model5<-coxph(Surv(survtime, surv) ~ factor(treat)+factor(figo)+factor(karc)+factor(turesc)+factor(histoc),data=ova)
summary(model5)
par(mfrow=c(1,2),mar=c(4,4,2,2))
checkPH.kar<-cox.zph(model5)[3]
plot(checkPH.kar,main="Check PH assumption of karc")
points(checkPH.kar$x,checkPH.kar$y,pch=16,col="lightgray")
abline(h=0,lty=2,col=2)

library(survminer)
ggcoxzph(cox.zph(model5))

checkPH.hist<-cox.zph(model5)[5]
plot(checkPH.hist,main="Check PH assumption of histoc")
points(checkPH.hist$x,checkPH.hist$y,pch=16,col="lightgray")
abline(h=0,lty=2,col=2)

# 9. Check the PH assumption for karc and for histoc by the log(-log(Survival)) plot

km.karc<-survfit(Surv(survtime, surv) ~ factor(karc),data=ova)
plot(km.karc, col=c("black", "red"), fun="cloglog",ylab="log(-log(Survival))",xlab="log(time)",main="Check PH assumption of karc")

km.histoc<-survfit(Surv(survtime, surv) ~ factor(histoc),data=ova)
plot(km.histoc, col=c("black", "red"), fun="cloglog",ylab="log(-log(Survival))",xlab="log(time)",main="Check PH assumption of histoc")


# 10. Fit the model of point 5. but stratifying on the histoc variable

model.s<-coxph(Surv(survtime, surv) ~ factor(turesc)+strata(histoc) + factor(treat)+factor(figo) +factor(karc),data=ova)
summary(model.s)


# 11. Fit the model of point 5. but with kar (numerical) instead of karc. Check the linearity assumption using the martingale residuals plot and splines

model11<-coxph(Surv(survtime, surv) ~ factor(treat)+factor(figo) +kar+factor(turesc)+factor(histoc),data=ova)
summary(model11)

par(mfrow=c(1,2),mar=c(4,4,2,2))
mar.res<-resid(model11,type='martingale')
plot(ova$kar, mar.res,
     xlab="Time", ylab="Martingale Residuals",
     main="Check functional form of kar")
lines(lowess(ova$kar, mar.res),col='red')

library(splines)
model.kar.bs <- coxph(Surv(survtime, surv>0) ~ bs(kar,4), data = ova)
library(Greg)
par(mar=c(4,4,1,1))
plotHR(model.kar.bs, term="kar", plot.bty="o", ylog=T, xlim = c(30, 100), rug="density", xlab="Karnofsky",polygon_ci=T)

# 12. Fit the model of point 11. but with a time-dependent effect for histoc: include an interaction term with the function of time g(t)=ln(t)-ln(3)

par(mfrow=c(1,1),mar=c(4,4,2,2))
cut.points<- unique(ova$survtime[ova$surv==1])
ova$entry<-0
ova2<-survSplit(data=ova, cut=cut.points, end="survtime", start="entry", event="surv")
ova3<-ova2[order(ova2$npaz),]
ova3$td1<-ifelse(ova3$histoc==1, log(ova3$survtime/3), 0)

modeltd<-coxph(Surv(entry, survtime, surv) ~ factor(treat)+factor(figo) +kar+factor(turesc)+factor(histoc)+td1,data=ova3)
summary(modeltd)
coefftd<-modeltd$coeff

t<-seq(0.01,12,0.01)
HRtd<-exp( coefftd[8]+coefftd[9]*log(t/3))
plot(t, HRtd, type='l', xlab="time", ylab="HR(t)")

head(cbind(t, HRtd))
