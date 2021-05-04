
setwd('~/Dropbox/work/CHL5209H_2018/data')
getwd()

library(survival)
library(Epi)
library(Hmisc)
head(transplant)
brain <- read.csv('brain.csv')
head(brain)

table(brain$treat, brain$event)

brain <- brain[!is.na(brain$local),]
model <- coxph(Surv(weeks, event) ~ treat + resect75 + age + interval + karn + race + local + male + nitro + factor(path) + grade, data=brain)
summary(model)

martingaleres <- residuals(model, type=c('martingale'))
devianceres <- residuals(model, type=c('deviance'))
dfbeta <- residuals(model, type=c('dfbeta'))
dfbetas <- residuals(model, type=c('dfbetas'))

outpath <- '~/Dropbox/work/CHL5209H_2018/slides'

# Check martingale and deviance residuals for continuous covariates:

pdf(file.path(outpath, paste('martingale_age.pdf', sep='')), height=7, width=7, paper='special')
op <- par(mar=c(4.5,4.5,1,1))
plot(brain$age, martingaleres, xlab='Age', ylab='Martingale residual')
lines(lowess(brain$age, martingaleres), lwd=2, col='blue')
abline(h=0, lty='dotted')
par(op)
dev.off()

pdf(file.path(outpath, paste('martingale_interval.pdf', sep='')), height=7, width=7, paper='special')
op <- par(mar=c(4.5,4.5,1,1))
plot(brain$interval, martingaleres, xlab='Interval', ylab='Martingale residual')
lines(lowess(brain$interval, martingaleres), lwd=2, col='blue')
abline(h=0, lty='dotted')
par(op)
dev.off()

pdf(file.path(outpath, paste('deviance_age.pdf', sep='')), height=7, width=7, paper='special')
op <- par(mar=c(4.5,4.5,1,1))
plot(brain$age, devianceres, xlab='Age', ylab='Deviance residual')
lines(lowess(brain$age, devianceres), lwd=2, col='blue')
abline(h=0, lty='dotted')
par(op)
dev.off()

pdf(file.path(outpath, paste('deviance_interval.pdf', sep='')), height=7, width=7, paper='special')
op <- par(mar=c(4.5,4.5,1,1))
plot(brain$interval, devianceres, xlab='Interval', ylab='Deviance residual')
lines(lowess(brain$interval, devianceres), lwd=2, col='blue')
abline(h=0, lty='dotted')
par(op)
dev.off()

# Unscaled dfbeta influence measures:

pdf(file.path(outpath, paste('dfbeta_age.pdf', sep='')), height=7, width=7, paper='special')
op <- par(mar=c(4.5,4.5,1,1))
plot(brain$age, dfbeta[,3], xlab='Age', ylab='dfbeta', ylim=c(-2/sqrt(nrow(dfbeta)), 2/sqrt(nrow(dfbeta))))
lines(lowess(brain$age, dfbeta[,3]), lwd=2, col='blue')
abline(h=0, lty='dotted')
par(op)
dev.off()

pdf(file.path(outpath, paste('dfbeta_interval.pdf', sep='')), height=7, width=7, paper='special')
op <- par(mar=c(4.5,4.5,1,1))
plot(brain$interval, dfbeta[,4], xlab='Interval', ylab='dfbeta')
lines(lowess(brain$interval, dfbeta[,4]), lwd=2, col='blue')
abline(h=0, lty='dotted')
par(op)
dev.off()

# Scaled dfbeta influence measures (compare to the threshold of 2/sqrt(n)):

pdf(file.path(outpath, paste('dfbetas_age.pdf', sep='')), height=7, width=7, paper='special')
op <- par(mar=c(4.5,4.5,1,1))
plot(brain$age, dfbetas[,3], xlab='Age', ylab='dfbetas')
lines(lowess(brain$age, dfbetas[,3]), lwd=2, col='blue')
abline(h=c(-2/sqrt(nrow(dfbetas)), 0, 2/sqrt(nrow(dfbetas))), lty='dotted')
par(op)
dev.off()

pdf(file.path(outpath, paste('dfbetas_interval.pdf', sep='')), height=7, width=7, paper='special')
op <- par(mar=c(4.5,4.5,1,1))
plot(brain$interval, dfbetas[,4], xlab='Interval', ylab='dfbeta')
lines(lowess(brain$interval, dfbetas[,4]), lwd=2, col='blue')
abline(h=0, lty='dotted')
abline(h=c(-2/sqrt(nrow(dfbetas)), 0, 2/sqrt(nrow(dfbetas))), lty='dotted')
par(op)
dev.off()

# Checks for proportionality:

cox.zph(model, global=FALSE)

pdf(file.path(outpath, paste('schoenfeld_age.pdf', sep='')), height=7, width=7, paper='special')
op <- par(mar=c(4.5,4.5,1,1))
plot(cox.zph(model, global=FALSE), var=3)
abline(h=0, lty='dotted')
par(op)
dev.off()

pdf(file.path(outpath, paste('schoenfeld_interval.pdf', sep='')), height=7, width=7, paper='special')
op <- par(mar=c(4.5,4.5,1,1))
plot(cox.zph(model, global=FALSE), var=4)
abline(h=0, lty='dotted')
par(op)
dev.off()

pdf(file.path(outpath, paste('schoenfeld_male.pdf', sep='')), height=7, width=7, paper='special')
op <- par(mar=c(4.5,4.5,1,1))
plot(cox.zph(model, global=FALSE), var=8)
abline(h=0, lty='dotted')
par(op)
dev.off()

# Examples of tests for covariate-time interactions:

model <- coxph(Surv(weeks, event) ~ treat + tt(treat) + resect75 + age + interval + 
                                    karn + race + local + male + tt(male) + nitro + 
                                    factor(path) + grade, 
                                    tt=function(x,t, ...) x * t, data=brain)
summary(model)

model <- coxph(Surv(weeks, event) ~ treat + resect75 + age + interval + karn + 
                                    race + local + male + tt(male) + nitro + 
                                    factor(path) + grade, 
                                    tt=function(x,t, ...) x * t, data=brain)
summary(model)
