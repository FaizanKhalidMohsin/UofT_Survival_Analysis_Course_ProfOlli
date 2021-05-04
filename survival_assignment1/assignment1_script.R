age = rep(22:29, 2)
death = c(433, 412, 337, 331, 287, 242, 215, 192, 
      24, 36, 66, 102, 138, 171, 185, 200)
personyears = c(91444, 86835, 75892, 63241, 52023, 42123, 36915, 32215,
                8556, 12708, 23203, 35415, 46207, 55675, 60470, 64770) 
z = c(rep(0, 8), rep(1, 8))

model1 = glm(death ~ z + age + offset(log(personyears)), family = poisson(link="log") )
summary(model1)

model2 = glm(death ~ z + factor(age) + offset(log(personyears)), family = poisson(link="log") )
summary(model2)

model3 = glm(death ~ z + factor(age) + z:age offset(log(personyears)), family = poisson(link="log") )
summary(model3)

# Question 4
#http://individual.utoronto.ca/osaarela/finrisk82.csv
data=read.csv("http://individual.utoronto.ca/osaarela/finrisk82.csv", sep = ";")
str(data)
summary(data)


###### Function to transform the data. 

rm(list=ls())

outpath <- '~/Dropbox/work/CHL5209H_2016/data'
dataset <- read.csv2(file.path(outpath, 'finrisk82.csv'), colClasses='character')
dataset = data
dataset$events <- as.numeric(dataset$events)
dataset$followupyears <- as.numeric(dataset$followupyears)
dataset$year <- as.numeric(dataset$year)
ncol(dataset)
nrow(dataset)
by(as.numeric(dataset$events), dataset$endpoint, sum)
by(as.numeric(dataset$followupyears), dataset$endpoint, sum)

# Plot Lexis diagram:

str(dataset)
table(dataset$endpoint)
frdeaths <- dataset[dataset$endpoint == 'DEATH',]
chddeaths <- dataset[dataset$endpoint == 'CHD2',!(names(dataset) %in% c('endpoint','followupyears'))]
names(chddeaths)[names(chddeaths) == 'events'] <- 'chd'
frdeaths <- merge(frdeaths, chddeaths)
frdeaths <- frdeaths[,!(names(frdeaths) %in% 'endpoint')]
table(frdeaths$events >= frdeaths$chd)

frdeaths$yearmid <- 1984.5 * (frdeaths$year >= 1982 & frdeaths$year < 1987) +
  1989.5 * (frdeaths$year >= 1987 & frdeaths$year < 1992) +
  1994.5 * (frdeaths$year >= 1992 & frdeaths$year < 1997) +
  1999.5 * (frdeaths$year >= 1997 & frdeaths$year < 2002) +
  2004.5 * (frdeaths$year >= 2002 & frdeaths$year < 2007) +
  2009 * (frdeaths$year >= 2007 & frdeaths$year < 2011)
frdeaths$agemid <- 30 * (frdeaths$agegr == '<35') +
  40 * (frdeaths$agegr == '35-44') +
  50 * (frdeaths$agegr == '45-54') +
  60 * (frdeaths$agegr == '55-64') +
  67.5 * (frdeaths$agegr == '65-69') +
  72.5 * (frdeaths$agegr == '70-74') +
  77.5 * (frdeaths$agegr == '75-79') +
  82.5 * (frdeaths$agegr == '80-84') +
  89.0 * (frdeaths$agegr == '>85')
fragg <- aggregate(frdeaths[,c('events','followupyears')], by=list(frdeaths$yearmid, frdeaths$agemid), FUN=sum)

# postscript(file.path(outpath, 'frlexis.eps'), width=6, height=6, paper='special', horizontal=FALSE)
op <- par(mar=c(4,4,0,0), mgp=c(2,1,0))
minyr <- 1982
maxyr <- 2011
minage <- 25
maxage <- 94

plot(NULL, NULL, type='n', xlim=c(minyr, maxyr), ylim=c(minage, maxage), axes=FALSE,
     xlab='calendar year', ylab='age', main='')
axis(1, at=c(seq(minyr, maxyr, by=5), 2011), las=1, pos=minage)
axis(2, at=c(25,35,45,55,65,70,75,80,85,94), las=1, pos=minyr)

ygrid <- c(seq(minyr, 2007, by=5), 2011)
agrid <- c(25,35,45,55,65)
segments(rep(minyr, length(agrid)), agrid, pmin(minyr + (maxage - agrid), maxyr), pmin(agrid + (maxyr - minyr), maxage), col='gray80')

ygrid <- c(seq(minyr, 2007, by=5), 2011)
agrid <- c(25,35,45,55,65,70,75,80,85,94)
segments(rep(minyr, length(agrid)), agrid, rep(maxyr, length(agrid)), agrid, col='gray50')
segments(ygrid, rep(minage, length(ygrid)), ygrid, rep(maxage, length(ygrid)), col='gray50')

lines(c(maxyr, maxyr), c(minage, maxage))

ygrid <- seq(1840, 1950, by=10)
agrid <- seq(20, 90, by=5)
for(i in 1:nrow(fragg)) {
  text(fragg[i,'Group.1'], fragg[i,'Group.2'], fragg[i,'events'], pos=3, offset=0.15, cex=0.8, col='red')
  text(fragg[i,'Group.1'], fragg[i,'Group.2'], round(fragg[i,'followupyears']), pos=1, offset=0.15, cex=0.8, col='blue')    
}
par(op)
# dev.off()

# Analysis variables:

frdeaths$ageg <- 1 * (frdeaths$agegr == '<35') +
  2 * (frdeaths$agegr == '35-44') +
  3 * (frdeaths$agegr == '45-54') +
  4 * (frdeaths$agegr == '55-64') +
  5 * (frdeaths$agegr == '65-69') +
  6 * (frdeaths$agegr == '70-74') +
  7 * (frdeaths$agegr == '75-79') +
  8 * (frdeaths$agegr == '80-84') +
  9 * (frdeaths$agegr == '>85')
frdeaths$yearg <- frdeaths$year - min(frdeaths$year) + 1
frdeaths$sexg <- 0 * (frdeaths$sex == 'men') +
  1 * (frdeaths$sex == 'women')
frdeaths$area <- 0 * (frdeaths$rua == 'FIN-EASa') +
  1 * (frdeaths$rua == 'FIN-WESa')
nyears <- length(unique(frdeaths$yearg))
nagegroups <- length(unique(frdeaths$ageg))

# Poisson regression for total mortality:
model <- glm(events ~ as.factor(yearg) + as.factor(ageg) + sexg + area, offset=log(followupyears), data=frdeaths, family=poisson(link='log'))
summary(model)





