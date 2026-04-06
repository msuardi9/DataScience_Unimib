### -----------------------------
### LAB 3 - Assessment of the performance of a risk predictor
### -----------------------------

# 0. Load data and R packages

htrap<-read.table(file.choose(), na.strings=".",header=T,row.names=NULL)
names(htrap)
head(htrap)

library(survival)
library(pROC)
library(riskRegression)

# 1. Estimate the overall Kaplan-Meier survival function and compute the estimated risk of death within 1 year after transplant

fit<-survfit(Surv(time.y, event) ~ 1, data=htrap)
plot(fit, xlab='Years since transplant', ylab='Survival probability')
summary(fit,times=1)
1-summary(fit,times=1)$surv


# 2. Estimate the overall censoring function and check the number of subject censored within 1 year after transplant

fit.rev<-survfit(Surv(time.y, event==0) ~ 1, data=htrap)
plot(fit.rev, xlab='Years since transplant', ylab='Censoring probability')
summary(fit.rev,times=1)
summary(fit.rev,times=1)$n.event


# 3. Exclude patients with a potential follow-up < 1 year (i.e. patients transplanted from 01/01/2017 onwards). This will unbiasedly produce a sample where all alive patients have been observed for at least 1 year. 

htrap.1y<-htrap[grepl(2017, htrap$date.htx)==F,]
table(htrap.1y$time.y>1,htrap.1y$event)


# 4. Working on the sample created in 3., estimate the overall Kaplan-Meier survival function and compute the estimated risk of death within 1 year after transplant

fit<-survfit(Surv(time.y, event) ~ 1, data=htrap.1y)
plot(fit, xlab='Years since transplant', ylab='Survival probability')
summary(fit,times=1)
1-summary(fit,times=1)$surv


# 5. Working on the sample created in 3., fit two Cox PH models: the first including as predictors "gender", "age", "bmi.cat" and "status"; the second including also "crea" and "azot".
# For each patient, calculate the predicted risk of death within one year according to both models.

model1 <- coxph(Surv(time.y, event==1)~gender+age+bmi.cat+status,data=htrap.1y)
summary(model1)
fit1<-survfit(model1,newdata=htrap.1y)
htrap.1y$riskdeath1<-1-as.numeric(summary(fit1,times=1)$surv)

model2 <- coxph(Surv(time.y, event==1)~gender+age+bmi.cat+status+crea+azot,data=htrap.1y)
summary(model2)
fit2<-survfit(model2,newdata=htrap.1y)
htrap.1y$riskdeath2<-1-as.numeric(summary(fit2,times=1)$surv)


# 6. To check the model calibration, create the calibration plot for predicted risks of both models.

#Risk quantiles of 1
q1 <- quantile(htrap.1y$riskdeath1,  probs = seq(0.1, 0.9, 0.1))

# Risk classes of 1:
riskclass1 <- ifelse(htrap.1y$riskdeath1 < q1[1], 1,
                     ifelse(htrap.1y$riskdeath1 < q1[2], 2,
                            ifelse(htrap.1y$riskdeath1 < q1[3], 3,
                                   ifelse(htrap.1y$riskdeath1 < q1[4], 4,
                                          ifelse(htrap.1y$riskdeath1 < q1[5], 5,
                                                 ifelse(htrap.1y$riskdeath1 < q1[6], 6,
                                                        ifelse(htrap.1y$riskdeath1 < q1[7], 7,
                                                               ifelse(htrap.1y$riskdeath1 < q1[8], 8,
                                                                      ifelse(htrap.1y$riskdeath1 < q1[9], 9, 
                                                                             10)))))))))
# Observed proportion of events in risk classes of 1:
obsclass1 <- c()
for (i in 1:10) {
  group <- i
  ratio <- sum(htrap.1y$event[riskclass1 == group]) / sum(riskclass1 == group)
  obsclass1 <- c(obsclass1, ratio)
}

# Average predicted risk in classes of 1:
meanriskclass1 <- c()
for (i in 1:10) {
  group <- i
  meanrisk <- mean(htrap.1y$riskdeath1[riskclass1 == group])
  meanriskclass1 <- c(meanriskclass1, meanrisk)
}

#Risk quantiles of 2
q2 <- quantile(htrap.1y$riskdeath2,  probs = seq(0.1, 0.9, 0.1))

# Risk classes of 1:
riskclass2 <- ifelse(htrap.1y$riskdeath2 < q2[1], 1,
                     ifelse(htrap.1y$riskdeath2 < q2[2], 2,
                            ifelse(htrap.1y$riskdeath2 < q2[3], 3,
                                   ifelse(htrap.1y$riskdeath2 < q2[4], 4,
                                          ifelse(htrap.1y$riskdeath2 < q2[5], 5,
                                                 ifelse(htrap.1y$riskdeath2 < q2[6], 6,
                                                        ifelse(htrap.1y$riskdeath2 < q2[7], 7,
                                                               ifelse(htrap.1y$riskdeath2 < q2[8], 8,
                                                                      ifelse(htrap.1y$riskdeath2 < q2[9], 9, 
                                                                             10)))))))))
# Observed proportion of events in risk classes of 2:
obsclass2 <- c()
for (i in 1:10) {
  group <- i
  ratio <- sum(htrap.1y$event[riskclass2 == group]) / sum(riskclass2 == group)
  obsclass2 <- c(obsclass2, ratio)
}

# Average predicted risk in classes of 2:
meanriskclass2 <- c()
for (i in 1:10) {
  group <- i
  meanrisk <- mean(htrap.1y$riskdeath2[riskclass2 == group])
  meanriskclass2 <- c(meanriskclass2, meanrisk)
}


# Calibration plot

par(mar=c(5,5,1,1))
plot(meanriskclass1, obsclass1, main = '', pch = 1, xlim = c(0,1), cex = 1.5,
     ylim = c(0,1), lwd = 3, ylab = 'Observed event rate', cex.lab = 1.7, cex.axis = 1.7,
     xlab = 'Average risk within decile', xaxt = "n", yaxt = "n", frame = F) 
axis(2, at = c(0,0.2,0.4,0.6,0.8,1), labels = NA, pos = 0)
axis(2, at = c(0,0.2,0.4,0.6,0.8,1), labels = c(0,0.2,0.4,0.6,0.8,1), cex.axis = 1.7, pos = 0)
axis(1, at = c(0,0.2,0.4,0.6,0.8,1), labels = c(0,0.2,0.4,0.6,0.8,1), cex.axis = 1.7, pos = 0)

points(meanriskclass2, obsclass2, lwd = 3, pch = 2, cex = 1.5)

lines(c(0, 1), c(0, 1), lty = 2, lwd = 2) 

lines(c(0,1), c(1,1), lty = 1)
lines(c(1,1), c(0,1), lty = 1)

legend(x = 0, y = 1, c("riskdeath1","riskdeath2"), lwd = c(3, 3), pch = c(1, 2), 
       bty = 'n', cex = 1.8, lty = c(NA,NA)) 


# 7. Compare the calibration of the two risk predictors also by calculating the Brier Score.

# indicator of death within 1 year
htrap.1y$event.1y<-ifelse(htrap.1y$time.y<=1,1,0)

# Brier Score
(BS_1 <- mean((htrap.1y$event.1y - htrap.1y$riskdeath1) ^ 2)) 
(BS_2 <- mean((htrap.1y$event.1y - htrap.1y$riskdeath2) ^ 2))

# Brier Score under strong calibration
BS_1sc  <- mean(htrap.1y$riskdeath1*(1-htrap.1y$riskdeath1))
BS_2sc <- mean(htrap.1y$riskdeath2*(1-htrap.1y$riskdeath2))

BS_1-BS_1sc
BS_2-BS_2sc


# 8. Check the discrimination performance of the two risk predictors using ROC curve and AUC index. Find also the optimal cut-off value according to the Youden index.

roc1<-roc(htrap.1y$event.1y, htrap.1y$riskdeath1)
plot(1 - roc1$specificities, roc1$sensitivities, 
     type = 'l', ylab = 'TPF', xlab = 'FPF', lwd = 3, xaxt = "n", yaxt = "n", 
     xlim = c(0,1), cex.lab = 1.7, frame = F)
axis(1, at = c(0,0.25,0.5,0.75,1), labels = NA, pos = 0)
axis(1, at = c(0,0.25,0.5,0.75,1), labels = c(0,0.25,0.5,0.75,1), cex.axis = 1.7, pos = 0)
axis(2, at = c(0,0.25,0.5,0.75,1), labels = c(0,0.25,0.5,0.75,1), cex.axis = 1.7, pos = 0)
Youden1<-roc1$sensitivities+roc1$specificities-1
optimal.cut.off1<-roc1$thresholds[Youden1==max(Youden1)]
cbind(optimal.cut.off1,Youden=max(Youden1))
points(1-roc1$specificities[roc1$thresholds==optimal.cut.off1],roc1$sensitivities[roc1$thresholds==optimal.cut.off1],pch=0,cex=1.7)

roc2<-roc(htrap.1y$event.1y, htrap.1y$riskdeath2)
lines(1 - roc2$specificities, roc2$sensitivities, 
      lwd = 3, lty = 3)
Youden2<-roc2$sensitivities+roc2$specificities-1
optimal.cut.off2<-roc2$thresholds[Youden2==max(Youden2)]
cbind(optimal.cut.off2,Youden=max(Youden2))
points(1-roc2$specificities[roc2$thresholds==optimal.cut.off2],roc2$sensitivities[roc2$thresholds==optimal.cut.off2],pch=0,cex=1.7)

# Chance line:
abline(a=0, b=1, lty = 2, lwd = 2)
lines(c(0,1), c(1,1), lty = 1)
lines(c(1,1), c(0,1), lty = 1)

legend(x = 0, y = 1, c("model1","model2"), lwd = c(3,3), lty = c(1,3), bty = 'n', cex = 1.7)

# AUC
(AUC1 <- roc1$auc)
(AUC2 <- roc2$auc)   
AUC2 - AUC1


# 9. Working on the original dataset (N=805 patients) compute the Brier score (adjusted for censoring) and the AUC index. Perform internal validation using bootstrap to correct for the over-optimism.

#remeber to set x=T in order to use the output for the "Score" function in "riskRegression" package
model1 <- coxph(Surv(time.y, event==1)~gender+age+bmi.cat+status,data=htrap,x=T)
summary(model1)
fit1<-survfit(model1,newdata=htrap)
htrap$riskdeath1<-1-as.numeric(summary(fit1,times=1)$surv)


model2 <- coxph(Surv(time.y, event==1)~gender+age+bmi.cat+status+crea+azot,data=htrap,x=T)
summary(model2)
fit2<-survfit(model2,newdata=htrap)
htrap$riskdeath2<-1-as.numeric(summary(fit2,times=1)$surv)


score<- Score(list("model1"=model1,"model2"=model2),
              formula=Surv(time.y, event==1)~1,
              data=htrap,conf.int=T,
              times=seq(1,5,1),
              plots=c("calibration","ROC"))

plotROC(score,times=1,cens.method="local")
title(main="time-dependent ROC at 1 years")
plotROC(score,times=3,cens.method="local")
title(main="time-dependent ROC at 3 years")
plotCalibration(score,times=1,cens.method="local",method="quantile",q=10)
title(main="calibration at 1 years")
plotCalibration(score,times=3,cens.method="local",method="quantile",q=10)
title(main="calibration at 3 years")

# internal validation using bootstrap
CVscore<- Score(list("model1"=model1,"model2"=model2),
                formula=Surv(time.y, event==1)~1,
                data=htrap,conf.int=T,
                times=seq(1,5,1),
                split.method="loob",B=100,seed=1803)
CVscore


# 10. Calculate and plot the expected Net Benefit according to the risk threshold for both risk predictors.

library(dcurves)
dca(Surv(time.y, event==1) ~ riskdeath1 + riskdeath2, data = htrap, time = 1)