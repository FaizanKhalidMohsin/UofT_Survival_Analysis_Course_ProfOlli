
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

# Wrong analysis:

# K-M comparison of survival in 'transplanted' and 'not transplanted' groups:

fit <- survfit(Surv(timetodeath,fustat)~transplant, data=jasa)
summary(fit)

plot(fit, col=1:2, xlab="Years since entering the waiting list",
     ylab="Survival", lwd=2, conf.int=T, mark.time=TRUE,
     main="Kaplan-Meier curves for mortality")
legend("topright",col=1:2, legend=c("Not transplanted","Transplanted"), lwd=2)

# log-rank test:

diff <- survdiff(Surv(timetodeath,fustat)~transplant, data=jasa)
diff

# Cox model with a fixed exposure variable:

model <- coxph(Surv(futime,fustat)~transplant, data=jasa)
summary(model)

# Add baseline covariates:

model <- coxph(Surv(timetodeath,fustat)~transplant + age + surgery + acceptyr, data=jasa)
summary(model)

# Correct analysis:

# Create long format dataset:

nrow(jasa)
jasa2 <- NULL
for (i in 1:nrow(jasa)) {
    if (jasa$transplant[i] == 0) {
        jasa2 <- rbind(jasa2, data.frame(jasa[i,], id=i, start=0, stop=jasa$timetodeath[i], event=jasa$fustat[i], tr=0)) 
    }
    else if (jasa$transplant[i] == 1) {
        jasa2 <- rbind(jasa2, data.frame(jasa[i,], id=i, start=0, stop=jasa$timetotransplant[i], event=0, tr=0))
        jasa2 <- rbind(jasa2, data.frame(jasa[i,], id=i, start=jasa$timetotransplant[i], stop=jasa$timetodeath[i], event=jasa$fustat[i], tr=1))
    }
}
nrow(jasa2)
head(jasa2)
head(jasa2[(ncol(jasa2)-4):ncol(jasa2)])

# Both starting and ending time of the follow-up interval have to be specified now in the survival outcome:

model <- coxph(Surv(start,stop,event)~tr + age + surgery + acceptyr, data=jasa2)
summary(model)

# Use tt to create the time-dependent exposure variable (the results should be the same):

jasa$timetotransplant <- ifelse(is.na(jasa$timetotransplant), Inf, jasa$timetotransplant)
head(jasa)

model <- coxph(Surv(timetodeath,fustat)~tt(timetotransplant) + age + surgery + acceptyr, tt=function(x,t, ...) I(x<t), data=jasa)
summary(model)

# Allow interaction between transplant and time since transplant:

model <- coxph(Surv(timetodeath,fustat)~tt(timetotransplant) + age + surgery + acceptyr, 
tt=function(x,t, ...) cbind(I(x<t), I(x<t) * (t-ifelse(is.finite(x), x, 0))), data=jasa)
summary(model)

# Allow interaction between transplant and wait time:

model <- coxph(Surv(timetodeath,fustat)~tt(timetotransplant) + age + surgery + acceptyr, 
tt=function(x,t, ...) cbind(I(x<t), I(x<t) * ifelse(is.finite(x), x, 0)), data=jasa)
summary(model)

# Investigate whether the mismatch scores and rejection predict subsequent survival through interactions with the transplant indicator:

table(jasa$mismatch)
jasa2$mismatchia <- ifelse(jasa2$tr == 1, jasa2$mismatch-1, 0)
model <- coxph(Surv(start,stop,event)~tr + factor(mismatchia) + age + surgery + acceptyr, data=jasa2)
summary(model)

table(jasa$hla.a2)
jasa2$hla.a2ia <- ifelse(jasa2$tr == 1, jasa2$hla.a2, 0)
model <- coxph(Surv(start,stop,event)~tr + factor(hla.a2ia) + age + surgery + acceptyr, data=jasa2)
summary(model)

table(jasa$mscore)
jasa2$mscoreia <- ifelse(jasa2$tr == 1, jasa2$mscore, 0)
model <- coxph(Surv(start,stop,event)~tr + mscoreia + age + surgery + acceptyr, data=jasa2)
summary(model)

table(jasa$reject)
jasa2$rejectia <- ifelse(jasa2$tr == 1, jasa2$reject, 0)
model <- coxph(Surv(start,stop,event)~tr + rejectia + age + surgery + acceptyr, data=jasa2)
summary(model)

# Investigate interactions with the baseline covariates:

jasa2$ageia <- ifelse(jasa2$tr == 1, jasa2$age, 0)
jasa2$surgeryia <- ifelse(jasa2$tr == 1, jasa2$surgery, 0)
jasa2$yria <- ifelse(jasa2$tr == 1, jasa2$acceptyr, 0)
model <- coxph(Surv(start,stop,event)~tr + rejectia + ageia + surgeryia + yria + age + surgery + acceptyr, data=jasa2)
summary(model)

# Trouble, big correlations between the parameter estimates:

cov2cor(vcov(model))

plot(jasa2$tr, jasa2$yria)

# Center the age and year variables and try again:

jasa2$age <- jasa2$age - mean(jasa$age)
jasa2$acceptyr <- jasa2$acceptyr - mean(jasa$acceptyr)
head(jasa2)

jasa2$ageia <- ifelse(jasa2$tr == 1, jasa2$age, 0)
jasa2$surgeryia <- ifelse(jasa2$tr == 1, jasa2$surgery, 0)
jasa2$yria <- ifelse(jasa2$tr == 1, jasa2$acceptyr, 0)
model <- coxph(Surv(start,stop,event)~tr + rejectia + ageia + surgeryia + yria + age + surgery + acceptyr, data=jasa2)
summary(model)

cov2cor(vcov(model))

# Choose this as the final model and interpret the results:

model <- coxph(Surv(start,stop,event)~tr + rejectia + yria + age + surgery + acceptyr, data=jasa2)
summary(model)

# We could calculate robust/sandwich variance estimates which allow for dependencies
# between observations by adding a cluster term into the model. However, this is not necessary here. (Why?)

model <- coxph(Surv(start,stop,event)~tr + rejectia + yria + age + surgery + acceptyr + cluster(id), data=jasa2)
summary(model)


