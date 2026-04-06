### -----------------------------
### LAB 4: competing risks and multistate models
### -----------------------------


### Exercise 1

# 0. Read the data remdesivir.txt and call the "survival", "prodlim" and "cmprsk" packages

library(survival)
library(prodlim)
library(cmprsk)

rdv  <- read.table(file.choose(), na.strings=".",header=T,row.names=NULL)
str(rdv)
head(rdv)


# 1. Estimate the crude incidence functions of estubation and death in ICU for the two treatment groups. 
# Estimate also the proportion of patients still remaining in ICU over time.

par(mar=c(4,2,3,1),mfrow=c(1,2))
crFit <- prodlim(Hist(time,event)~treatment,data=rdv)
summary(crFit)
plot(crFit,cause=1,xlab="Time since intubation (days)",xlim=c(0,30),
     legend.x="topleft",legend.legend=c("no RDV","RDV"),atrisk.title="",atrisk.labels=c("no RDV","RDV"))
mtext("Estubation",side=3,at=14,cex=2,line=1)

plot(crFit,cause=2,xlab="Time since intubation (days)",xlim=c(0,30),
     legend.x="topleft",legend.legend=c("no RDV","RDV"),atrisk.title="",atrisk.labels=c("no RDV","RDV"))
mtext("Death in ICU",side=3,at=14,cex=2,line=1)

par(mar=c(4,2,3,1),mfrow=c(1,1))
fit <- prodlim(Hist(time,event>0)~treatment,data=rdv)
summary(fit)
plot(fit,xlab="Time since intubation (days)",xlim=c(0,30),
     legend.x="topright",legend.legend=c("no RDV","RDV"),atrisk.title="",atrisk.labels=c("no RDV","RDV"))
mtext("Still in ICU",side=3,at=14,cex=2,line=1)


# 2. Compare the crude incidence curves under the two treatments by the Gray test
ci<-with(rdv,cuminc(time,event,treatment))
ci$Tests


# 3. Estimate the cumulative cause-specific hazards by treatment
# Estimate the cumulative hazard of the composite end-point (Estubation or death in ICU)

# Estubation
fit.est<-survfit(Surv(time,event==1) ~ treatment, data=rdv)
summary(fit.est)

na.haz.est0<-cumsum(fit.est["treatment=0"]$n.event/fit.est["treatment=0"]$n.risk)
na.haz.est1<-cumsum(fit.est["treatment=1"]$n.event/fit.est["treatment=1"]$n.risk)

par(mfrow=c(1,2),mar=c(4,4,3,1))
plot( c(0,fit.est["treatment=0"]$time), c(0,na.haz.est0), xlab='time', ylab="cumulative cause-specific hazard", type="s", lwd=3, xlim=c(0,30), ylim=c(0,1.5))
lines( c(0,fit.est["treatment=1"]$time), c(0,na.haz.est1), type="s", lwd=3, col="goldenrod1")
legend("top",c("No RDV","RDV"),lwd=3,col=c("black","goldenrod1"),bty="n")
mtext("Estubation",side=3,at=14,cex=2,line=1)

# Death in ICU
fit.d<-survfit(Surv(time,event==2) ~ treatment, data=rdv)
summary(fit.d)

na.haz.d0<-cumsum(fit.d["treatment=0"]$n.event/fit.d["treatment=0"]$n.risk)
na.haz.d1<-cumsum(fit.d["treatment=1"]$n.event/fit.d["treatment=1"]$n.risk)

plot( c(0,fit.d["treatment=0"]$time), c(0,na.haz.d0), xlab='time', ylab="cumulative cause-specific hazard", type="s", lwd=3, xlim=c(0,30), ylim=c(0,1.5))
lines( c(0,fit.d["treatment=1"]$time), c(0,na.haz.d1), type="s", lwd=3, col="goldenrod1")
legend("top",c("No RDV","RDV"),lwd=3,col=c("black","goldenrod1"),bty="n")
mtext("Death in ICU",side=3,at=14,cex=2,line=1)


# Composite endpoint
fit.tot<-survfit(Surv(time,event>0) ~ treatment, data=rdv)
summary(fit.tot)

na.haz.tot0<-cumsum(fit.tot["treatment=0"]$n.event/fit.tot["treatment=0"]$n.risk)
na.haz.tot1<-cumsum(fit.tot["treatment=1"]$n.event/fit.tot["treatment=1"]$n.risk)

par(mfrow=c(1,1),mar=c(4,4,3,1))
plot( c(0,fit.tot["treatment=0"]$time), c(0,na.haz.tot0), xlab='time', ylab="cumulative hazard", type="s", lwd=3, xlim=c(0,30), ylim=c(0,2))
lines( c(0,fit.tot["treatment=1"]$time), c(0,na.haz.tot1), type="s", lwd=3, col="goldenrod1")
legend("top",c("No RDV","RDV"),lwd=3,col=c("black","goldenrod1"),bty="n")
mtext("Composite end-point",side=3,at=14,cex=2,line=1)


# 4. Compare the Aalen-Johansen estimate of estubation and death in ICU with the Kaplan-Meier estimate (where additional censoring on competing event is added)
par(mar=c(4,2,3,1),mfrow=c(1,2))
crFit <- prodlim(Hist(time,event)~treatment,data=rdv)
plot(crFit,cause=1,xlab="Time since intubation (days)",xlim=c(0,30),confint = F,
     legend.x="topleft",legend.legend=c("no RDV AJ","RDV AJ"),atrisk.title="",atrisk.labels=c("no RDV","RDV"))

kmFit1 <- prodlim(Hist(time,event==1)~treatment,data=rdv)
plot(kmFit1,type="cuminc",add=T,lty=3,confint=F)
mtext("AJ (solid) vs KM (dashed) estubation",side=3,at=12,cex=1,line=1)


plot(crFit,cause=2,xlab="Time since intubation (days)",xlim=c(0,30),confint = F,
     legend.x="topleft",legend.legend=c("no RDV AJ","RDV AJ"),atrisk.title="",atrisk.labels=c("no RDV","RDV"))

kmFit2 <- prodlim(Hist(time,event==2)~treatment,data=rdv)
plot(kmFit2,type="cuminc",add=T,lty=3,confint=F)
mtext("AJ (solid) vs KM (dashed) death in ICU",side=3,at=12,cex=1,line=1)


# this is a cause-specific Cox model for the cause-specific hazard of event 1
m<-coxph(Surv(time,event==2)~treatment,data=rdv)
summary(m)

# this is a cause-specific Cox model for the cause-specific hazard of event 2
m<-coxph(Surv(time,event==1)~treatment,data=rdv)
summary(m)


### Exercise 2

# 0. Call the "mstate" package and load the ebmt3 dataset

library(mstate)
data(ebmt3)

str(ebmt3)
head(ebmt3)

# 1. Build the transition matrix for the illness-death model with Transplant, PR and RelDeath states
tmat <- transMat(x = list(c(2, 3), c(3), c()), names = c("Tx", "PR","RelDeath"))
tmat

# 2. Transform time from days to years and reshape data from wide to long format
ebmt3$prtime <- ebmt3$prtime/365.25
ebmt3$rfstime <- ebmt3$rfstime/365.25
msbmt <- msprep(time = c(NA, "prtime", "rfstime"), status = c(NA,"prstat", "rfsstat"), 
                data = ebmt3, trans = tmat)
head(msbmt)

# 3. Estimate and plot the cumulative transition hazards
c1 <- coxph(Surv(Tstart, Tstop, status) ~  strata(trans), data = msbmt, method = "breslow")
msf1 <- msfit(c1, newdata = msbmt, trans = tmat)
par(mfrow=c(1,1),mar=c(4,4,2,2))
plot(msf1, cols = rep(1, 3), lwd = 2, lty = 1:3, xlab = "Years since transplant", ylab = "Cumulative hazard", legend.pos = c(3, 0.4),ylim=c(0,0.9))

# 4. Estimate and plot the state probabilities in a "stacked" plot
pt <- probtrans(msf1, predt = 0)
plot(pt, ord = c(3, 2, 1), lwd = 2, xlab = "Years since transplant",
     ylab = "Prediction probabilities", cex = 0.75, type="filled", 
     legend = c("Alive", "PR", "Relapse or death"))
