# install the 'remotes' package
install.packages("remotes")

# install the 'esmpack' package (version 0.1-17)
remotes::install_github("wviechtb/esmpack@v0.1-17")

# install the 'latticeExtra' package
install.packages("latticeExtra")


library(tidyverse)
library(lme4)
library(lmerTest)

# load the 'nlme' package
library(nlme)

# load the 'esmpack' package
library(esmpack)

# load the 'lattice' and 'latticeExtra' packages
library(lattice)
library(latticeExtra)


#load data
dat <- read.table("example_1.txt", header=TRUE, sep="\t", na.strings="", as.is=TRUE)


# check that there are no missing values in design and time-related variables
check.nomiss(id, id, data=dat)
check.nomiss(day, id, data=dat)
check.nomiss(beep, id, data=dat)
check.nomiss(obs, id, data=dat)
check.nomiss(beeptime, id, data=dat)

# check that age, sex, and status are really time-invariant within subjects
# (and that they are either completely missing or complete within subjects)
check.timeinvar(age, id, data=dat, na.rm=FALSE)
check.timeinvar(sex, id, data=dat, na.rm=FALSE)
check.timeinvar(status, id, data=dat, na.rm=FALSE)

# age, sex, and status are actually never missing in this dataset
check.nomiss(age, id, data=dat)
check.nomiss(sex, id, data=dat)

# number of subjects and number of rows of data
nsub(dat$id)
nrow(dat)


# frequency table for each item (to check for out-of-range values)
table(dat$mood_cheerf)
table(dat$mood_relaxed)
table(dat$mood_satisfi)
table(dat$mood_irritat)
table(dat$mood_anxious)
table(dat$mood_down)
table(dat$mood_guilty)
table(dat$mood_insecur)
table(dat$mood_lonely)
table(dat$eventpl)
table(dat$soc_alone)
table(dat$soc_pleasant)
table(dat$act_else)
table(dat$use_coffee)
table(dat$use_alcohol)
table(dat$location)

# frequency table for day, beep, and obs counters
table(dat$day)
table(dat$beep)
table(dat$obs)


# check that there are no duplicated day-beep combinations within each subject
check.nodup(interaction(day, beep), id, data=dat)
# check that there are no duplicated obs values within each subject
check.nodup(obs, id, data=dat)
# check that there are no duplicated beeptime values within each day within each subject
check.nodup(beeptime, interaction(id, day), data=dat)


# check that the beep times are consistent with the blocks (beep numbers)
all((dat$beeptime >=  450 & dat$beeptime <  540)[dat$beep ==  1])
all((dat$beeptime >=  540 & dat$beeptime <  630)[dat$beep ==  2])
all((dat$beeptime >=  630 & dat$beeptime <  720)[dat$beep ==  3])
all((dat$beeptime >=  720 & dat$beeptime <  810)[dat$beep ==  4])
all((dat$beeptime >=  810 & dat$beeptime <  900)[dat$beep ==  5])
all((dat$beeptime >=  900 & dat$beeptime <  990)[dat$beep ==  6])
all((dat$beeptime >=  990 & dat$beeptime < 1080)[dat$beep ==  7])
all((dat$beeptime >= 1080 & dat$beeptime < 1170)[dat$beep ==  8])
all((dat$beeptime >= 1170 & dat$beeptime < 1260)[dat$beep ==  9])
all((dat$beeptime >= 1260 & dat$beeptime < 1350)[dat$beep == 10])


# cross-tabulation of soc_alone and soc_pleasant; there should not be any
# observations for soc_pleasant when soc_alone is 1; can also check if there
# are any observations for soc_pleasant when soc_alone is missing
table(dat$soc_alone, dat$soc_pleasant, useNA="ifany")

# compute response delays
dat$delay <- dat$resptime - dat$beeptime

# barplot of the delay values
barplot(table(dat$delay), xlab="Delay in Minutes")

# mean, SD, and range of the delay values
round(mean(dat$delay, na.rm=TRUE), digits=1)
round(sd(dat$delay, na.rm=TRUE), digits=1)
range(dat$delay, na.rm=TRUE)


# compute the response frequencies (the number of non-missing response times within subjects)
respfreq <- calc.nomiss(resptime, id, data=dat)

# histogram of the response frequencies
hist(respfreq, col="lightgray", xlab="Response Frequency", right=FALSE)

# mean and SD of the compliance rates
round(mean(respfreq) / 60 * 100, digits=1)
round(sd(respfreq) / 60 * 100, digits=1)

# range of the compliance rates
round(range(respfreq / 60 * 100), digits=1)


############################################################################

# mean, SD, and range of the age of the subjects
round(mean(get.timeinvar(age, id, data=dat)), digits=1)
round(sd(get.timeinvar(age, id, data=dat)), digits=1)
round(range(get.timeinvar(age, id, data=dat)), digits=1)

# number of female/male subjects
table(get.timeinvar(sex, id, data=dat))

# as percentages
round(prop.table(table(get.timeinvar(sex, id, data=dat))) * 100, digits=1)

# number of control, depressed, and psychotic subjects
table(get.timeinvar(status, id, data=dat))

# compute variables posaff and negaff (average of the 3 PA and 6 NA items)
dat$posaff <- combitems(c("mood_cheerf", "mood_relaxed", "mood_satisfi"), data=dat)
dat$negaff <- combitems(c("mood_irritat", "mood_anxious", "mood_down", "mood_guilty", "mood_insecur", "mood_lonely"), data=dat)

# add the person-level means of the eventpl variable to the dataset
dat$meventpl <- calc.mean(eventpl, id, data=dat, expand=TRUE)

# create a new variable for the within-person mean centered eventpl variable
dat$ceventpl <- calc.mcent(eventpl, id, data=dat)

# create a lagged version of the eventpl variable
dat$leventpl <- lagvar(eventpl, id=id, day=day, data=dat)

# create a lagged version of the posaff variable
dat$lposaff <- lagvar(posaff, id=id, day=day, data=dat)

# examine part of the data for subjects c100 and c101
dat[dat$id == "c100",c("age", "sex", "day", "beep", "eventpl", "meventpl", "ceventpl", "leventpl", "posaff")][c(1:4,60),]
dat[dat$id == "c101",c("age", "sex", "day", "beep", "eventpl", "meventpl", "ceventpl", "leventpl", "posaff")][c(1:4,60),]

############################################################################

# a random sample of 6 subjects from each group

ids <- c("c104", "c115", "c119", "c154", "c199", "c209", "d215", "d242", "d272", "d282", "d304", "d309", "p328", "p352", "p367", "p373", "p385", "p395")
sub <- dat[dat$id %in% ids,]

# plot posaff as a function of tothours

col3 <- c("blue1","orange2","green3")
col3.light <- c(rgb(130,130,205,maxColorValue=255), rgb(255,204,50,maxColorValue=255), rgb(130,205,130,maxColorValue=255))
col3.alpha <- apply(col2rgb(col3), MARGIN=2, FUN=function(x) rgb(x[1], x[2], x[3], alpha=50, maxColorValue=255))

xyplot(posaff ~ tothours | id, data=sub, groups=id, subset = !is.na(posaff),
       xlab="Time (in Hours)", ylab="Positive Affect",
       scales=list(x=list(at=seq(0,144,24))),
       col=rep(col3, each=6), fill=rep(col3, each=6), col.line=rep(col3.light, each=6),
       layout=c(6,3), aspect=1, pch=21,
       panel = function(x,y,...){
         xval <- 1:400
         panel.xblocks(xval, xval >=  0 & xval <=  24, col = "lightgray", alpha=0.4)
         panel.xblocks(xval, xval >= 48 & xval <=  72, col = "lightgray", alpha=0.4)
         panel.xblocks(xval, xval >= 96 & xval <= 120, col = "lightgray", alpha=0.4)
         panel.xyplot(x, y, type=c("l","p"), ...)
         panel.xyplot(x, y, type="p", ...)
       },
       par.settings = list(strip.background=list(col="lightgrey")))

# plot posaff as a function of eventpl

xyplot(posaff ~ eventpl | id, data=sub, groups=id, type=c("p","r"),
       xlab="Pleasantness of Event", ylab="Positive Affect",
       col=rep(col3.light, each=6), col.line=rep(col3, each=6),
       layout=c(6,3), aspect=1, pch=19, lwd=2,
       par.settings = list(strip.background=list(col="lightgrey")))

############################################################################


### analysis / model fitting

# estimate the average PA and the between- and within-person variability
res <- lme(posaff ~ 1, random = ~ 1 | id, data=dat, na.action=na.omit)
summary(res)

# intercept variance
round(getVarCov(res)[1,1], digits=3)

# error variance
round(sigma(res)^2, digits=3)

# approximate 95% CI for mu
tab <- coef(summary(res))
round(tab[1] + c(-1,1)*qnorm(.975) * tab[2], digits=2)

# approximate 95% interval for the subject-specific averages
round(fixef(res) + c(-1,1)*qnorm(.975) * sqrt(getVarCov(res)[1,1]), digits=2)

# interval into which 95% of the residuals should fall
round(c(-1,1)*qnorm(.975) * sigma(res), digits=2)

# histogram of the BLUPs
hist(coef(res)[[1]], breaks=seq(0-.125,8+.125,by=0.25), freq=FALSE, main="", xlab="Average Positive Affect")
curve(dnorm(x, mean=fixef(res), sd=sqrt(getVarCov(res)[1,1])), add=TRUE, lwd=3)

# histogram of the residuals
hist(resid(res), breaks=seq(-6-.1,6+.1,by=0.2), xlim=c(-4,4), freq=FALSE, main="", xlab="Residual")
curve(dnorm(x, mean=0, sd=sigma(res)), add=TRUE, lwd=3)

# ICC (intraclass correlation coefficient)
round(getVarCov(res)[1,1] / (getVarCov(res)[1,1] + sigma(res)^2), digits=2)

############################################################################

# examine differences in average level of positive affect as a function of sex,
# age, and mental health status
res <- lme(posaff ~ sex + I(age-35) + status, random = ~ 1 | id, data=dat, na.action=na.omit)
summary(res)
round(coef(summary(res)), digits=3)

# intercept and error variance
round(getVarCov(res)[1,1], digits=3)
round(sigma(res)^2, digits=3)

# proportion of intercept variance accounted for
res0 <- lme(posaff ~ 1, random = ~ 1 | id, data=dat, na.action=na.omit)
res1 <- lme(posaff ~ sex + I(age-35) + status, random = ~ 1 | id, data=dat, na.action=na.omit)
round((getVarCov(res0)[1,1] - getVarCov(res1)[1,1]) / getVarCov(res0)[1,1], digits=2)

############################################################################

# model to examine the association between positive affect and event pleasantness
res <- lme(posaff ~ eventpl, random = ~ eventpl | id, data=dat, na.action=na.omit)
summary(res)
round(coef(summary(res)), digits=3)

# intercept variance, slope variance, intercept-slope correlation, and error variance
round(getVarCov(res)[1,1], digits=3)
round(getVarCov(res)[2,2], digits=4)
round(getVarCov(res)[1,2] / sqrt(getVarCov(res)[1,1] * getVarCov(res)[2,2]), digits=2)
round(sigma(res)^2, digits=3)

# proportion of error variance accounted for by eventpl
res0 <- lme(posaff ~ 1,       random = ~ 1       | id, data=dat, na.action=na.omit)
res1 <- lme(posaff ~ eventpl, random = ~ eventpl | id, data=dat, na.action=na.omit)
round((sigma(res0)^2 - sigma(res1)^2) / sigma(res0)^2, digits=2)

############################################################################

# histograms of the BLUPs, residuals, and scatterplots of the BLUPs

par(mfrow=c(2,2), mar=c(5,5,4,0.5))

hist(coef(res)[[1]], breaks=seq(0-.125,8+.125,by=0.25), freq=FALSE, main="", xlab="Intercept")
curve(dnorm(x, mean=fixef(res)[1], sd=sqrt(getVarCov(res)[1,1])), add=TRUE, lwd=3)
title("(a) Distribution of Estimated Intercepts", font.main=1)

hist(coef(res)[[2]], breaks=seq(-0.1-0.0125,0.5+0.0125,by=0.025), freq=FALSE, main="", xlab="Slope")
curve(dnorm(x, mean=fixef(res)[2], sd=sqrt(getVarCov(res)[2,2])), add=TRUE, lwd=3)
title("(b) Distribution of Estimated Slopes", font.main=1)

hist(resid(res), breaks=seq(-6-.1,6+.1,by=0.2), xlim=c(-4,4), freq=FALSE, main="", xlab="Residual")
curve(dnorm(x, mean=0, sd=sigma(res)), add=TRUE, lwd=3)
title("(c) Distribution of the Residuals", font.main=1)

plot(coef(res)[[1]], coef(res)[[2]], pch=19, xlab="Intercept", ylab="Slope", cex=1, col="gray40")
abline(v=fixef(res)[1], lty="dotted")
abline(h=fixef(res)[2], lty="dotted")
title("(d) Scatterplot of Estimated Intercepts versus Slopes", font.main=1)

############################################################################

# model to examine the differences between the mental health status groups in
# the association between positive affect and event pleasantness
res <- lme(posaff ~ status + eventpl + status:eventpl, random = ~ eventpl | id, data=dat, na.action=na.omit)
summary(res)
round(coef(summary(res)), digits=3)

# intercept variance, slope variance, intercept-slope correlation, and error variance
round(getVarCov(res)[1,1], digits=3)
round(getVarCov(res)[2,2], digits=4)
round(getVarCov(res)[1,2] / sqrt(getVarCov(res)[1,1] * getVarCov(res)[2,2]), digits=2)
round(sigma(res)^2, digits=3)

# test difference in intercepts between the depression and psychosis groups
anova(res, L=c(0,1,-1,0,0,0))

# test difference in slopes between the depression and psychosis groups
anova(res, L=c(0,0,0,0,-1,1))

# proportion of slope variance accounted for by group
res0 <- lme(posaff ~ eventpl, random = ~ eventpl | id, data=dat, na.action=na.omit)
res1 <- lme(posaff ~ status + eventpl + status:eventpl, random = ~ eventpl | id, data=dat, na.action=na.omit)
round((getVarCov(res0)[2,2] - getVarCov(res1)[2,2]) / getVarCov(res0)[2,2], digits=2)

############################################################################

# plot intercepts and slopes per group

res <- lme(posaff ~ 0 + status + status:eventpl, random = ~ eventpl | id, data=dat, na.action=na.omit)
round(coef(summary(res)), digits=3)

coefs <- coef(res, augFrame=TRUE, which=which(names(dat) == "status"))

par(mfrow=c(1,3))

for (j in unique(dat$status)) {
  plot(NA, NA, xlim=c(-3,3), ylim=c(1,7), xlab="Pleasantness of Event", ylab="Positive Affect")
  title(paste("Group:", j), cex.main=1.4)
  for (i in 1:nrow(coefs)) {
    if (coefs$status[i] == j)
      abline(coefs[i,which(unique(dat$status) == j)] + coefs[i,7], coefs[i,3+which(unique(dat$status) == j)] + coefs[i,8], col=col3.alpha[which(unique(dat$status) == j)])
  }
  abline(fixef(res)[which(unique(dat$status) == j)], fixef(res)[3+which(unique(dat$status) == j)], lwd=4, col=col3[which(unique(dat$status) == j)])
  box()
}

############################################################################

# disentangling the within- and between-person associations between event pleasantness and positive affect
res <- lme(posaff ~ meventpl + ceventpl, random = ~ ceventpl | id, data=dat, na.action=na.omit)
summary(res)
round(coef(summary(res)), digits=3)

# intercept variance, slope variance, intercept-slope correlation, and error variance
round(getVarCov(res)[1,1], digits=3)
round(getVarCov(res)[2,2], digits=4)
round(getVarCov(res)[1,2] / sqrt(getVarCov(res)[1,1] * getVarCov(res)[2,2]), digits=2)
round(sigma(res)^2, digits=3)

# test if the slopes of the within- and between-person relationships differ significantly
fixef(res)[2] / fixef(res)[3]
z <- (fixef(res)[2] - fixef(res)[3]) / sqrt(vcov(res)[2,2] + vcov(res)[3,3] - 2*vcov(res)[2,3])
round(c(zval=z, pval=2*pnorm(abs(z), lower.tail=FALSE)), digits=3)

# proportion of error variance accounted for by ceventpl
res0 <- lme(posaff ~ 1, random = ~ 1 | id, data=dat, na.action=na.omit)
res1 <- lme(posaff ~ meventpl + ceventpl, random = ~ ceventpl | id, data=dat, na.action=na.omit)
round((sigma(res0)^2 - sigma(res1)^2) / sigma(res0)^2, digits=2)

############################################################################

# allow the within- and between-person associations to differ across groups
res <- lme(posaff ~ status + meventpl + ceventpl + status:meventpl + status:ceventpl,
           random = ~ ceventpl | id, data=dat, na.action=na.omit)
summary(res)
round(coef(summary(res)), digits=3)

############################################################################

# fit model with lagged event pleasantness as predictor
res <- lme(posaff ~ leventpl, random = ~ leventpl | id, data=dat, na.action=na.omit)
summary(res)
round(coef(summary(res)), digits=3)

# intercept variance, slope variance, intercept-slope correlation, and error variance
round(getVarCov(res)[1,1], digits=3)
round(getVarCov(res)[2,2], digits=4)
round(getVarCov(res)[1,2] / sqrt(getVarCov(res)[1,1] * getVarCov(res)[2,2]), digits=2)
round(sigma(res)^2, digits=3)

# proportion of error variance accounted for by leventpl
res0 <- lme(posaff ~ 1, random = ~ 1 | id, data=dat, na.action=na.omit)
res1 <- lme(posaff ~ leventpl, random = ~ leventpl | id, data=dat, na.action=na.omit)
round((sigma(res0)^2 - sigma(res1)^2) / sigma(res0)^2, digits=2)

# amount of data used for model (3)
sum(!is.na(dat$posaff) & !is.na(dat$eventpl))

# amount of data used for model (6)
sum(!is.na(dat$posaff) & !is.na(dat$leventpl))

############################################################################

# examine group differences in the lagged association
res <- lme(posaff ~ status + leventpl + status:leventpl,
           random = ~ leventpl | id, data=dat, na.action=na.omit)
summary(res)

# test difference in slopes between the depression and psychosis groups
anova(res, L=c(0,0,0,0,1,-1))

############################################################################

# model to estimate the autocorrelation and the lagged relationship
res <- lme(posaff ~ lposaff + leventpl,
           random = ~ lposaff + leventpl | id, data=dat, na.action=na.omit)
summary(res)

# examine group differences in the autocorrelation and the lagged relationship
res <- lme(posaff ~ status + lposaff + leventpl + status:lposaff + status:leventpl,
           random = ~ lposaff + leventpl | id, data=dat, na.action=na.omit)
summary(res)

# also add resphour number within a day as predictor (and allow this to interact with status)
res <- lme(posaff ~ status + lposaff + leventpl + resphour + status:lposaff + status:leventpl + status:resphour,
           random = ~ lposaff + leventpl + resphour | id, data=dat, na.action=na.omit)
summary(res)
round(coef(summary(res)), digits=3)

# intercept variance, slope variances, and error variance
round(getVarCov(res)[1,1], digits=3)
round(getVarCov(res)[2,2], digits=4)
round(getVarCov(res)[3,3], digits=4)
round(getVarCov(res)[4,4], digits=4)
round(sigma(res)^2, digits=3)

############################################################################
### analysis / model fitting

# estimate the average PA and the between- and within-person variability
res <- lme(posaff ~ 1, random = ~ 1 | id, data=dat, na.action=na.omit)
summary(res)

# intercept variance
round(getVarCov(res)[1,1], digits=3)

# error variance
round(sigma(res)^2, digits=3)

# approximate 95% CI for mu
tab <- coef(summary(res))
round(tab[1] + c(-1,1)*qnorm(.975) * tab[2], digits=2)

# approximate 95% interval for the subject-specific averages
round(fixef(res) + c(-1,1)*qnorm(.975) * sqrt(getVarCov(res)[1,1]), digits=2)

# interval into which 95% of the residuals should fall
round(c(-1,1)*qnorm(.975) * sigma(res), digits=2)

# histogram of the BLUPs
hist(coef(res)[[1]], breaks=seq(0-.125,8+.125,by=0.25), freq=FALSE, main="", xlab="Average Positive Affect")
curve(dnorm(x, mean=fixef(res), sd=sqrt(getVarCov(res)[1,1])), add=TRUE, lwd=3)

# histogram of the residuals
hist(resid(res), breaks=seq(-6-.1,6+.1,by=0.2), xlim=c(-4,4), freq=FALSE, main="", xlab="Residual")
curve(dnorm(x, mean=0, sd=sigma(res)), add=TRUE, lwd=3)

# ICC (intraclass correlation coefficient)
round(getVarCov(res)[1,1] / (getVarCov(res)[1,1] + sigma(res)^2), digits=2)

############################################################################

# examine differences in average level of positive affect as a function of sex,
# age, and mental health status
res <- lme(posaff ~ sex + I(age-35) + status, random = ~ 1 | id, data=dat, na.action=na.omit)
summary(res)
round(coef(summary(res)), digits=3)

# intercept and error variance
round(getVarCov(res)[1,1], digits=3)
round(sigma(res)^2, digits=3)

# proportion of intercept variance accounted for
res0 <- lme(posaff ~ 1, random = ~ 1 | id, data=dat, na.action=na.omit)
res1 <- lme(posaff ~ sex + I(age-35) + status, random = ~ 1 | id, data=dat, na.action=na.omit)
round((getVarCov(res0)[1,1] - getVarCov(res1)[1,1]) / getVarCov(res0)[1,1], digits=2)

############################################################################

# model to examine the association between positive affect and event pleasantness
res <- lme(posaff ~ eventpl, random = ~ eventpl | id, data=dat, na.action=na.omit)
summary(res)
round(coef(summary(res)), digits=3)

# intercept variance, slope variance, intercept-slope correlation, and error variance
round(getVarCov(res)[1,1], digits=3)
round(getVarCov(res)[2,2], digits=4)
round(getVarCov(res)[1,2] / sqrt(getVarCov(res)[1,1] * getVarCov(res)[2,2]), digits=2)
round(sigma(res)^2, digits=3)

# proportion of error variance accounted for by eventpl
res0 <- lme(posaff ~ 1,       random = ~ 1       | id, data=dat, na.action=na.omit)
res1 <- lme(posaff ~ eventpl, random = ~ eventpl | id, data=dat, na.action=na.omit)
round((sigma(res0)^2 - sigma(res1)^2) / sigma(res0)^2, digits=2)

############################################################################

# histograms of the BLUPs, residuals, and scatterplots of the BLUPs

par(mfrow=c(2,2), mar=c(5,5,4,0.5))

hist(coef(res)[[1]], breaks=seq(0-.125,8+.125,by=0.25), freq=FALSE, main="", xlab="Intercept")
curve(dnorm(x, mean=fixef(res)[1], sd=sqrt(getVarCov(res)[1,1])), add=TRUE, lwd=3)
title("(a) Distribution of Estimated Intercepts", font.main=1)

hist(coef(res)[[2]], breaks=seq(-0.1-0.0125,0.5+0.0125,by=0.025), freq=FALSE, main="", xlab="Slope")
curve(dnorm(x, mean=fixef(res)[2], sd=sqrt(getVarCov(res)[2,2])), add=TRUE, lwd=3)
title("(b) Distribution of Estimated Slopes", font.main=1)

hist(resid(res), breaks=seq(-6-.1,6+.1,by=0.2), xlim=c(-4,4), freq=FALSE, main="", xlab="Residual")
curve(dnorm(x, mean=0, sd=sigma(res)), add=TRUE, lwd=3)
title("(c) Distribution of the Residuals", font.main=1)

plot(coef(res)[[1]], coef(res)[[2]], pch=19, xlab="Intercept", ylab="Slope", cex=1, col="gray40")
abline(v=fixef(res)[1], lty="dotted")
abline(h=fixef(res)[2], lty="dotted")
title("(d) Scatterplot of Estimated Intercepts versus Slopes", font.main=1)

############################################################################

# model to examine the differences between the mental health status groups in
# the association between positive affect and event pleasantness
res <- lme(posaff ~ status + eventpl + status:eventpl, random = ~ eventpl | id, data=dat, na.action=na.omit)
summary(res)
round(coef(summary(res)), digits=3)

# intercept variance, slope variance, intercept-slope correlation, and error variance
round(getVarCov(res)[1,1], digits=3)
round(getVarCov(res)[2,2], digits=4)
round(getVarCov(res)[1,2] / sqrt(getVarCov(res)[1,1] * getVarCov(res)[2,2]), digits=2)
round(sigma(res)^2, digits=3)

# test difference in intercepts between the depression and psychosis groups
anova(res, L=c(0,1,-1,0,0,0))

# test difference in slopes between the depression and psychosis groups
anova(res, L=c(0,0,0,0,-1,1))

# proportion of slope variance accounted for by group
res0 <- lme(posaff ~ eventpl, random = ~ eventpl | id, data=dat, na.action=na.omit)
res1 <- lme(posaff ~ status + eventpl + status:eventpl, random = ~ eventpl | id, data=dat, na.action=na.omit)
round((getVarCov(res0)[2,2] - getVarCov(res1)[2,2]) / getVarCov(res0)[2,2], digits=2)

############################################################################

# plot intercepts and slopes per group

res <- lme(posaff ~ 0 + status + status:eventpl, random = ~ eventpl | id, data=dat, na.action=na.omit)
round(coef(summary(res)), digits=3)

coefs <- coef(res, augFrame=TRUE, which=which(names(dat) == "status"))

par(mfrow=c(1,3))

for (j in unique(dat$status)) {
   plot(NA, NA, xlim=c(-3,3), ylim=c(1,7), xlab="Pleasantness of Event", ylab="Positive Affect")
   title(paste("Group:", j), cex.main=1.4)
   for (i in 1:nrow(coefs)) {
      if (coefs$status[i] == j)
         abline(coefs[i,which(unique(dat$status) == j)] + coefs[i,7], coefs[i,3+which(unique(dat$status) == j)] + coefs[i,8], col=col3.alpha[which(unique(dat$status) == j)])
   }
   abline(fixef(res)[which(unique(dat$status) == j)], fixef(res)[3+which(unique(dat$status) == j)], lwd=4, col=col3[which(unique(dat$status) == j)])
   box()
}

############################################################################

# disentangling the within- and between-person associations between event pleasantness and positive affect
res <- lme(posaff ~ meventpl + ceventpl, random = ~ ceventpl | id, data=dat, na.action=na.omit)
summary(res)
round(coef(summary(res)), digits=3)

# intercept variance, slope variance, intercept-slope correlation, and error variance
round(getVarCov(res)[1,1], digits=3)
round(getVarCov(res)[2,2], digits=4)
round(getVarCov(res)[1,2] / sqrt(getVarCov(res)[1,1] * getVarCov(res)[2,2]), digits=2)
round(sigma(res)^2, digits=3)

# test if the slopes of the within- and between-person relationships differ significantly
fixef(res)[2] / fixef(res)[3]
z <- (fixef(res)[2] - fixef(res)[3]) / sqrt(vcov(res)[2,2] + vcov(res)[3,3] - 2*vcov(res)[2,3])
round(c(zval=z, pval=2*pnorm(abs(z), lower.tail=FALSE)), digits=3)

# proportion of error variance accounted for by ceventpl
res0 <- lme(posaff ~ 1, random = ~ 1 | id, data=dat, na.action=na.omit)
res1 <- lme(posaff ~ meventpl + ceventpl, random = ~ ceventpl | id, data=dat, na.action=na.omit)
round((sigma(res0)^2 - sigma(res1)^2) / sigma(res0)^2, digits=2)

############################################################################

# allow the within- and between-person associations to differ across groups
res <- lme(posaff ~ status + meventpl + ceventpl + status:meventpl + status:ceventpl,
           random = ~ ceventpl | id, data=dat, na.action=na.omit)
summary(res)
round(coef(summary(res)), digits=3)

############################################################################

# fit model with lagged event pleasantness as predictor
res <- lme(posaff ~ leventpl, random = ~ leventpl | id, data=dat, na.action=na.omit)
summary(res)
round(coef(summary(res)), digits=3)

# intercept variance, slope variance, intercept-slope correlation, and error variance
round(getVarCov(res)[1,1], digits=3)
round(getVarCov(res)[2,2], digits=4)
round(getVarCov(res)[1,2] / sqrt(getVarCov(res)[1,1] * getVarCov(res)[2,2]), digits=2)
round(sigma(res)^2, digits=3)

# proportion of error variance accounted for by leventpl
res0 <- lme(posaff ~ 1, random = ~ 1 | id, data=dat, na.action=na.omit)
res1 <- lme(posaff ~ leventpl, random = ~ leventpl | id, data=dat, na.action=na.omit)
round((sigma(res0)^2 - sigma(res1)^2) / sigma(res0)^2, digits=2)

# amount of data used for model (3)
sum(!is.na(dat$posaff) & !is.na(dat$eventpl))

# amount of data used for model (6)
sum(!is.na(dat$posaff) & !is.na(dat$leventpl))

############################################################################

# examine group differences in the lagged association
res <- lme(posaff ~ status + leventpl + status:leventpl,
           random = ~ leventpl | id, data=dat, na.action=na.omit)
summary(res)

# test difference in slopes between the depression and psychosis groups
anova(res, L=c(0,0,0,0,1,-1))

############################################################################

# model to estimate the autocorrelation and the lagged relationship
res <- lme(posaff ~ lposaff + leventpl,
           random = ~ lposaff + leventpl | id, data=dat, na.action=na.omit)
summary(res)

# examine group differences in the autocorrelation and the lagged relationship
res <- lme(posaff ~ status + lposaff + leventpl + status:lposaff + status:leventpl,
           random = ~ lposaff + leventpl | id, data=dat, na.action=na.omit)
summary(res)

# also add resphour number within a day as predictor (and allow this to interact with status)
res <- lme(posaff ~ status + lposaff + leventpl + resphour + status:lposaff + status:leventpl + status:resphour,
           random = ~ lposaff + leventpl + resphour | id, data=dat, na.action=na.omit)
summary(res)
round(coef(summary(res)), digits=3)

# intercept variance, slope variances, and error variance
round(getVarCov(res)[1,1], digits=3)
round(getVarCov(res)[2,2], digits=4)
round(getVarCov(res)[3,3], digits=4)
round(getVarCov(res)[4,4], digits=4)
round(sigma(res)^2, digits=3)

############################################################################

#MODEL 1: CONTROL + Positive Affect
dat.control <- dat %>% 
  filter(status == "control")

#step one look at intercept variance in positive affect
res <- lme(posaff ~ 1, random =~ 1 | id, data = dat, na.action = na.omit)
summary(res)

# ICC (intraclass correlation coefficient)
round(getVarCov(res)[1,1] / (getVarCov(res)[1,1] + sigma(res)^2), digits=2)

#step two; include negative affect in positive affect
res <- lme(posaff ~ 1 + negaff, random =~ 1 | id, data = dat.control, na.action = na.omit)
summary(res)

#step three; disentangle within- and between-person negative affect

# add the person-level means of the eventpl variable to the dataset
dat$mnegaff <- calc.mean(negaff, id, data=dat, expand=TRUE)

# create a new variable for the within-person mean centered eventpl variable
dat$cnegaff <- calc.mcent(negaff, id, data=dat)

dat.control <- dat %>% 
  filter(status == "control")

res <- lme(posaff ~ 1 + mnegaff + cnegaff, random =~ 1 | id, data = dat.control, na.action = na.omit)
summary(res)

#MODEL 2: Psychosis & Negative affect


#step three; disentangle within- and between-person negative affect

dat.psychosis<- dat %>% 
  filter(status == "psychotic")

res <- lme(posaff ~ 1 + mnegaff + cnegaff, random =~ 1 | id, data = dat.psychosis, na.action = na.omit)
summary(res)

#MODEL 3: Depression & Negative affect

dat.depression <- dat %>% 
  filter(status == "depressed")

res <- lme(posaff ~ 1 + mnegaff + cnegaff, random =~ 1 | id, data = dat.depression, na.action = na.omit)
summary(res)


#repeat the same ordeal, but using positive affect as predictor of negative affect
#MODEL 1: CONTROL + Negative Affect
dat.control <- dat %>% 
  filter(status == "control")

#step one look at intercept variance in negative affect
res <- lme(negaff ~ 1, random =~ 1 | id, data = dat.control, na.action = na.omit)
summary(res)

# ICC (intraclass correlation coefficient)
round(getVarCov(res)[1,1] / (getVarCov(res)[1,1] + sigma(res)^2), digits=2)

#step two; include positive affect as predictor
res <- lme(negaff ~ 1 + posaff, random =~ 1 | id, data = dat.control, na.action = na.omit)
summary(res)

#step three; disentangle within- and between-person negative affect

# add the person-level means of the eventpl variable to the dataset
dat$mposaff <- calc.mean(posaff, id, data=dat, expand=TRUE)

# create a new variable for the within-person mean centered eventpl variable
dat$cposaff <- calc.mcent(posaff, id, data=dat)

dat.control <- dat %>% 
  filter(status == "control")

res <- lme(negaff ~ 1 + mposaff + cposaff, random =~ 1 | id, data = dat.control, na.action = na.omit)
summary(res)

#MODEL 2: Psychosis & Negative affect


#step three; disentangle within- and between-person negative affect

dat.psychosis<- dat %>% 
  filter(status == "psychotic")

res <- lme(negaff ~ 1 + mposaff + cposaff, random =~ 1 | id, data = dat.psychosis, na.action = na.omit)
summary(res)

#MODEL 3: Depression & Negative affect

dat.depression <- dat %>% 
  filter(status == "depressed")

res <- lme(negaff ~ 1 + mposaff + cposaff, random =~ 1 | id, data = dat.depression, na.action = na.omit)
summary(res)



############# MULTIVARIATE MODEL #######################

#transform data into a 'long' format with a single column for the different variables

dat.long <- pivot_longer(data = dat, 
                         cols = c("posaff", "negaff"), 
                         names_to = "var", values_to = "val") %>% 
            mutate(
              #add the dummy columns
              d_posaff = ifelse(var == "posaff", 1, 0),
              d_negaff = ifelse(var == "negaff", 1, 0)
            )

#check that the transformation worked
table(filter(dat.long, var == "posaff")$val == dat$posaff)
table(filter(dat.long, var == "negaff")$val == dat$negaff)

res <- lme(fixed = val ~ 0 + var, random = ~ 0 + var | id, data=dat.long,
           correlation = corSymm(form = ~ 1 | id/obs),
           weights = varIdent(form = ~ 1 | var), na.action = na.omit)
summary(res)

intervals(res)
coef(res)
