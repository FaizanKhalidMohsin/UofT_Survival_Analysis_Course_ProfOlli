
rm(list=ls())
library(survival)
library(Epi)

head(jasa)

# Create time variables:

jasa$acceptyr <- as.numeric(cal.yr(as.Date(jasa$accept.dt, format="%Y-%m-%d")))
jasa$deathyr <- as.numeric(cal.yr(as.Date(jasa$fu.date, format="%Y-%m-%d")))
jasa$birthyr <- as.numeric(cal.yr(as.Date(jasa$birth.dt, format="%Y-%m-%d")))
jasa$txyr <- as.numeric(cal.yr(as.Date(jasa$tx.date, format="%Y-%m-%d")))
jasa$timetodeath <- round(as.numeric(jasa$deathyr - jasa$acceptyr), 3)
jasa$timetotransplant <- round(as.numeric(jasa$txyr - jasa$acceptyr), 3)

# Check that the event times are positive:

eps <- 0.00001
jasa$timetotransplant <- ifelse(jasa$timetotransplant==0, eps, jasa$timetotransplant)
jasa$timetodeath <- ifelse(jasa$timetodeath==0, eps, jasa$timetodeath)
sel <- jasa$transplant==1
jasa$timetodeath[sel] <- ifelse(jasa$timetodeath[sel]==jasa$timetotransplant[sel], jasa$timetodeath[sel]+eps, jasa$timetodeath[sel])
head(jasa)

# Time and event indicator variables for competing risks analysis:

jasa$timetofirstevent <- ifelse(is.na(jasa$timetotransplant), jasa$timetodeath, jasa$timetotransplant)
head(jasa)

jasa$firstevent <- ifelse(is.na(jasa$timetotransplant), 2 * jasa$fustat, 1)
table(jasa$firstevent)
table(jasa$transplant)
table(jasa$fustat==1 & jasa$transplant==0)
head(jasa)

# K-M curve for probability of receiving transplantation:

fit <- survfit(Surv(timetofirstevent,firstevent==1)~1, data=jasa)
summary(fit)

plot(fit, col=2, xlab="Years since entering the waiting list", ylim=c(0,1),
     xlim=c(0,1), ylab="Probability", lwd=2, conf.int=T, mark.time=TRUE,
     main="Waiting time for transplantation", fun="event")
legend("bottomright", col=2, legend=c("1-KM"), lwd=2, bty='n')

# But less than 70% actually received the transplant during the study:

sum(jasa$transplant)/nrow(jasa)

# Non-parametric cumulative incidence analysis:

library(cmprsk)
cifit <- cuminc(jasa$timetofirstevent, jasa$firstevent, cencode=0)

plot(cifit, curvlab=c('Transplant','Death'), xlab="Years since entering the waiting list", lwd=2, xlim=c(0,1))
lines(fit, lwd=2, col=2, conf.int=F, fun="event")
legend("bottomright", col=2, legend=c("1-KM"), lwd=2, bty='n')

# Cox models for transplantation and death before transplantation:

jasa$acceptyr <- jasa$acceptyr - mean(jasa$acceptyr)

model1 <- coxph(Surv(jasa$timetofirstevent,jasa$firstevent==1)~age + surgery + acceptyr, data=jasa)
summary(model1)

model2 <- coxph(Surv(jasa$timetofirstevent,jasa$firstevent==2)~age + surgery + acceptyr, data=jasa)
summary(model2)

# Calculate the cumulative incidence of receiving transplant at average covariate values:

L1 <- basehaz(model1, centered=FALSE) 
L1$hazard <- L1$hazard * as.numeric(exp(t(colMeans(model.matrix(model1))) %*% coef(model1)))
names(L1) <- c('ch1','time')
L2 <- basehaz(model2, centered=FALSE) 
L2$hazard <- L2$hazard * as.numeric(exp(t(colMeans(model.matrix(model2))) %*% coef(model2)))
names(L2) <- c('ch2','time')
intersect(names(L1),names(L2))
ci <- merge(L1, L2)

# Cumulative hazard increments:

ci$h1 <- c(ci$ch1[1],ci$ch1[2:nrow(ci)]-ci$ch1[1:(nrow(ci)-1)])
ci$h2 <- c(ci$ch2[1],ci$ch2[2:nrow(ci)]-ci$ch2[1:(nrow(ci)-1)])
head(ci)
tail(ci)

# Overall survival probability (two alternative estimators):

ci$surv <- exp(-(c(0.0,ci$ch1[1:(nrow(ci)-1)]) + c(0.0,ci$ch2[1:(nrow(ci)-1)])))
ci$kmsurv <- cumprod(1.0 - (c(0.0,ci$h1[1:(nrow(ci)-1)]) + c(0.0,ci$h2[1:(nrow(ci)-1)])))

head(ci)
tail(ci)

# plot(ci$time, ci$surv, type='s')
# lines(ci$time, ci$kmsurv, type='s', col='red')

# Calculate cumulative incidence:

ci$ci1 <- cumsum(ci$h1 * ci$kmsurv)
ci$ci2 <- cumsum(ci$h2 * ci$kmsurv)
head(ci)
tail(ci)

# Plot:

plot(ci$time, ci$ci1, type='s', ylim=c(0,1), xlim=c(0,1), lwd=2, lty=c('solid'), xlab="Years since entering the waiting list", 
    ylab="Cumulative incidence", main='Cumulative incidence functions')
lines(ci$time, ci$ci2, type='s', lwd=2, lty=c('dashed'))
legend('bottomright', legend=c('Transplant','Death'), lwd=2, lty=c('solid','dashed'), bty='n')

# Fine & Gray models:

mm <- cbind(jasa$age, jasa$surgery, jasa$acceptyr)
colnames(mm) <- c('age','surgery','acceptyr')

fgmodel1 <- crr(jasa$timetofirstevent, jasa$firstevent, cov1=mm, failcode=1, cencode=0)
summary(fgmodel1)

fgmodel2 <- crr(jasa$timetofirstevent, jasa$firstevent, cov1=mm, failcode=2, cencode=0)
summary(fgmodel2)

# Cumulative incidence at average covariate values:

fgci1 <- predict(fgmodel1, cov1=t(colMeans(mm)))
fgci2 <- predict(fgmodel2, cov1=t(colMeans(mm)))

plot(ci$time, ci$ci1, type='s', ylim=c(0,1), xlim=c(0,1), lwd=2, lty=c('solid'), xlab="Years since entering the waiting list", ylab="Cumulative incidence",
     main='Cumulative incidence functions')
lines(ci$time, ci$ci2, type='s', lwd=2, lty=c('dashed'))
lines(fgci1[,1], fgci1[,2], type='s', lwd=2, lty=c('solid'), col='red')
lines(fgci2[,1], fgci2[,2], type='s', lwd=2, lty=c('dashed'), col='red')
legend('bottomright', legend=c('Transplant (F-G)','Death (F-G)','Transplant (Cox)','Death (Cox)'), lwd=2, 
       lty=c('solid','dashed','solid','dashed'), col=c('red','red','black','black'), bty='n')

# Comparison of predictions: calculate 2-month cumulative incidences of receiving transplant for everyone in the dataset:

ciall <- rep(NA, nrow(jasa))
fgciall <- rep(NA, nrow(jasa))
s <- 2/12
     
for (i in 1:nrow(jasa)) {    
    L1 <- basehaz(model1, centered=FALSE) 
    L1$hazard <- L1$hazard * as.numeric(exp(t(model.matrix(model1)[i,]) %*% coef(model1)))
    names(L1) <- c('ch1','time')
    L2 <- basehaz(model2, centered=FALSE) 
    L2$hazard <- L2$hazard * as.numeric(exp(t(model.matrix(model2)[i,]) %*% coef(model2)))
    names(L2) <- c('ch2','time')
    ci <- merge(L1, L2)
    
    ci$h1 <- c(ci$ch1[1],ci$ch1[2:nrow(ci)]-ci$ch1[1:(nrow(ci)-1)])
    ci$h2 <- c(ci$ch2[1],ci$ch2[2:nrow(ci)]-ci$ch2[1:(nrow(ci)-1)])

    ci$surv <- exp(-(c(0.0,ci$ch1[1:(nrow(ci)-1)]) + c(0.0,ci$ch2[1:(nrow(ci)-1)])))
    ci$kmsurv <- cumprod(1.0 - (c(0.0,ci$h1[1:(nrow(ci)-1)]) + c(0.0,ci$h2[1:(nrow(ci)-1)])))
        
    ci$ci1 <- cumsum(ci$h1 * ci$kmsurv)
    ci$ci2 <- cumsum(ci$h2 * ci$kmsurv)
    ciall[i] <- ci$ci1[findInterval(s, ci$time)]
    if (i %% 10 == 0)
        print(i)
}
sum(ciall)

fgci <- predict(fgmodel1, cov1=mm)
fgciall <- fgci[findInterval(s, fgci[,1]),2:ncol(fgci)]
sum(fgciall)

plot(ciall, fgciall, xlab='Cox model based CI', ylab='F-G model based CI', xlim=c(0,1), ylim=c(0,1),
     main='Comparison of predictions from the two models')
abline(a=0, b=1, lty='dotted', col='blue')




