###############################################################################
# R script file to fit the Bayesian GLM to PM data using Markov Chain
# Monte Carlo (MCMC).
#
# September, 2022
# Rob Laurie

library(ggplot2)
library(rstan)
library(TeachingDemos)
library(ggpubr)
library(bayesplot)
library(ggbreak)

mildew <- read.csv("Data.csv") 
head(mildew)

###############################################################################
# Model 1.

mod1 <- glm(FS_Binary ~ Second_Year, data = mildew, family=binomial)
summary(mod1)

idx1 <- 1:dim(mildew)[1]
na.idx1 <- idx1[with(mildew, is.na(FS_Binary) | is.na(Second_Year))]
length(idx1[-na.idx1])

mildew1 <- with(mildew, data.frame(FS_Binary[-na.idx1], Second_Year[-na.idx1]))
names(mildew1) <- c("FS_Binary", "Second_Year")

mildewData1 <- list(N = nrow(mildew1),
                    p = 2,
                    FS_Binary = mildew1$FS_Binary,
                    Second_Year = mildew1$Second_Year)

fit1a <- stan(file="mildew1.stan",data=mildewData1, iter = 1000, chains = 4)

traceplot(fit1a, "beta")
print(fit1a, "beta")

fit1b <- stan(fit=fit1a, data=mildewData1, iter=10000, chains=4)
traceplot(fit1b, par = "beta", inc_warmup=F)
results1 <- extract(fit1b, pars="beta",permuted = F, inc_warmup = FALSE) 
str(results1) 
print(fit1b, par = "beta")

Beta1 <- data.frame(Intercept = as.vector(results1[,1:4,1]),
                    X1 = as.vector(results1[,1:4,2]))

OutputMatrix1 <- data.frame(Risk = as.factor(c("EST","2nd")),
                            ProbEst = numeric(2),
                            Lower = numeric(2),
                            Upper = numeric(2))

### 0
est = exp(Beta1[,1])/(1+exp(Beta1[,1]))
OutputMatrix1$ProbEst[1] = round(median(est),2)
OutputMatrix1[1,3:4] = round(emp.hpd(est),2)

### 1
est = exp(Beta1[,1] + Beta1[,2])/(1+exp(Beta1[,1] + Beta1[,2]))
OutputMatrix1$ProbEst[2] = round(median(est),2)
OutputMatrix1[2,3:4] = round(emp.hpd(est),2)

Name_vec1 = c("EST","2nd")

OutputMatrix1$Risk = factor(OutputMatrix1$Risk,levels=Name_vec1, ordered = T)

###############################################################################
# Model 2.

mod2 <- glm(FS_Binary ~ Second_Year + Suscept_Binary, data = mildew, 
            family=binomial)
summary(mod2)

idx2 <- 1:dim(mildew)[1]
na.idx2 <- idx2[with(mildew, is.na(FS_Binary) | is.na(Second_Year) | 
                       is.na(Suscept_Binary))]
length(idx2[-na.idx2])

mildew2 <- with(mildew, data.frame(FS_Binary[-na.idx2], Second_Year[-na.idx2], 
                                   Suscept_Binary[-na.idx2]))
names(mildew2) <- c("FS_Binary", "Second_Year", "Suscept_Binary")

mildewData2 <- list(N = nrow(mildew2),
                    p = 3,
                    FS_Binary = mildew2$FS_Binary,
                    Second_Year = mildew2$Second_Year,
                    Suscept_Binary = mildew2$Suscept_Binary)

fit2a <- stan(file="mildew2.stan",data=mildewData2, iter = 1000, chains = 4)

traceplot(fit2a, "beta")
print(fit2a, "beta")

fit2b <- stan(fit=fit2a, data=mildewData2, iter=10000, chains=4)
traceplot(fit2b, par = "beta", inc_warmup=F)
results2 <- extract(fit2b, pars="beta",permuted = F, inc_warmup = FALSE) 
str(results2) 
print(fit2b, par = "beta")

Beta2 <- data.frame(Intercept = as.vector(results2[,1:4,1]),
                    X1 = as.vector(results2[,1:4,2]),
                    X2 = as.vector(results2[,1:4,3]))

OutputMatrix2 <- data.frame(Risk = as.factor(c("Est","2nd","Est Susceptible",
                                               "2nd Susceptible")),
                            ProbEst = numeric(4),
                            Lower = numeric(4),
                            Upper = numeric(4))

### 00
est = exp(Beta2[,1])/(1+exp(Beta2[,1]))
OutputMatrix2$ProbEst[1] = round(median(est),3)
OutputMatrix2[1,3:4] = round(emp.hpd(est),3)

### 10
est = exp(Beta2[,1] + Beta2[,2])/(1+exp(Beta2[,1] + Beta2[,2]))
OutputMatrix2$ProbEst[2] = round(median(est),3)
OutputMatrix2[2,3:4] = round(emp.hpd(est),3)

### 01
est = exp(Beta2[,1] + Beta2[,3])/(1+exp(Beta2[,1] + Beta2[,3]))
OutputMatrix2$ProbEst[3] = round(median(est),3)
OutputMatrix2[3,3:4] = round(emp.hpd(est),3)

### 11
est = exp(Beta2[,1] + Beta2[,3] + Beta2[,2])/(1+exp(Beta2[,1] + Beta2[,3] + 
                                                      Beta2[,2]))
OutputMatrix2$ProbEst[4] = round(median(est),3)
OutputMatrix2[4,3:4] = round(emp.hpd(est),3)

Name_vec2 = c("Est","2nd","Est Susceptible","2nd Susceptible")

OutputMatrix2$Risk = factor(OutputMatrix2$Risk,levels=Name_vec2, ordered = T)

###############################################################################
# Model 3.

mod3 <- glm(FS_Binary ~ Second_Year + Suscept_Binary + Poor_Prune, 
            data = mildew, family=binomial)
summary(mod3)

idx3 <- 1:dim(mildew)[1]
na.idx3 <- idx3[with(mildew, is.na(FS_Binary) | is.na(Second_Year) | 
                       is.na(Suscept_Binary) | is.na(Poor_Prune))]
length(idx3[-na.idx3])

mildew3 <- with(mildew, data.frame(FS_Binary[-na.idx3], Second_Year[-na.idx3], 
                              Suscept_Binary[-na.idx3], Poor_Prune[-na.idx3]))
names(mildew3) <- c("FS_Binary", "Second_Year", "Suscept_Binary", "Poor_Prune")

mildewData3 <- list(N = nrow(mildew3),
                    p = 4,
                    FS_Binary = mildew3$FS_Binary,
                    Second_Year = mildew3$Second_Year,
                    Suscept_Binary = mildew3$Suscept_Binary,
                    Poor_Prune = mildew3$Poor_Prune)

fit3a <- stan(file="mildew3.stan",data=mildewData3, iter = 1000, chains = 4)

traceplot(fit3a, "beta")
print(fit3a, "beta")

fit3b <- stan(fit=fit3a, data=mildewData3, iter=10000, chains=4)
traceplot(fit3b, par = "beta", inc_warmup=F)
results3 <- extract(fit3b, pars="beta",permuted = F, inc_warmup = FALSE) 
str(results3) 
print(fit3b, par = "beta")

Beta3 <- data.frame(Intercept = as.vector(results3[,1:4,1]),
                    X1 = as.vector(results3[,1:4,2]),
                    X2 = as.vector(results3[,1:4,3]),
                    X3 = as.vector(results3[,1:4,4]))

OutputMatrix3 <- data.frame(Risk = as.factor(c("EST","2nd","EST P","2nd P",
                                               "EST SUS","2nd SUS","EST SUS+P",
                                               "2nd SUS+P")),
                            ProbEst = numeric(8),
                            Lower = numeric(8),
                            Upper = numeric(8))

### 000
est = exp(Beta3[,1])/(1+exp(Beta3[,1]))
OutputMatrix3$ProbEst[1] = round(median(est),4)
OutputMatrix3[1,3:4] = round(emp.hpd(est),4)

### 100
est = exp(Beta3[,1] + Beta3[,2])/(1+exp(Beta3[,1] + Beta3[,2]))
OutputMatrix3$ProbEst[2] = round(median(est),4)
OutputMatrix3[2,3:4] = round(emp.hpd(est),4)

### 001
est = exp(Beta3[,1] + Beta3[,4])/(1+exp(Beta3[,1] + Beta3[,4]))
OutputMatrix3$ProbEst[3] = round(median(est),4)
OutputMatrix3[3,3:4] = round(emp.hpd(est),4)

### 101
est = exp(Beta3[,1] + Beta3[,2] + Beta3[,4])/(1+exp(Beta3[,1] + Beta3[,2] + 
                                                      Beta3[,4]))
OutputMatrix3$ProbEst[4] = round(median(est),4)
OutputMatrix3[4,3:4] = round(emp.hpd(est),4)

### 010
est = exp(Beta3[,1] + Beta3[,3])/(1+exp(Beta3[,1] + Beta3[,3]))
OutputMatrix3$ProbEst[5] = round(median(est),4)
OutputMatrix3[5,3:4] = round(emp.hpd(est),4)

### 110
est = exp(Beta3[,1] + Beta3[,2] + Beta3[,3])/(1+exp(Beta3[,1] + Beta3[,2] + 
                                                      Beta3[,3]))
OutputMatrix3$ProbEst[6] = round(median(est),4)
OutputMatrix3[6,3:4] = round(emp.hpd(est),4)

### 011
est = exp(Beta3[,1] + Beta3[,3] + Beta3[,4])/(1+exp(Beta3[,1] + Beta3[,3] + 
                                                      Beta3[,4]))
OutputMatrix3$ProbEst[7] = round(median(est),4)
OutputMatrix3[7,3:4] = round(emp.hpd(est),4)

### 111
est = exp(Beta3[,1] + Beta3[,2] + Beta3[,3] + Beta3[,4])/(1+exp(Beta3[,1] + 
                                            Beta3[,2] + Beta3[,3] + Beta3[,4]))
OutputMatrix3$ProbEst[8] = round(median(est),4)
OutputMatrix3[8,3:4] = round(emp.hpd(est),4)

Name_vec3 = c("EST","2nd","EST P","2nd P","EST SUS","2nd SUS","EST SUS+P",
              "2nd SUS+P")

OutputMatrix3$Risk = factor(OutputMatrix3$Risk,levels=Name_vec3, ordered = T)

###############################################################################
# Model 4.

mod4 <- glm(FS_Binary ~ Second_Year + Suscept_Binary + Mech_Prune, 
            data = mildew, family=binomial)
summary(mod4)

idx4 <- 1:dim(mildew)[1]
na.idx4 <- idx4[with(mildew, is.na(FS_Binary) | is.na(Second_Year) | 
                       is.na(Suscept_Binary) | is.na(Mech_Prune))]
length(idx4[-na.idx4])

mildew4 <- with(mildew, data.frame(FS_Binary[-na.idx4], Second_Year[-na.idx4], 
                                Suscept_Binary[-na.idx4], Mech_Prune[-na.idx4]))
names(mildew4) <- c("FS_Binary", "Second_Year", "Suscept_Binary", "Mech_Prune")

mildewData4 <- list(N = nrow(mildew4),
                    p = 4,
                    FS_Binary = mildew4$FS_Binary,
                    Second_Year = mildew4$Second_Year,
                    Suscept_Binary = mildew4$Suscept_Binary,
                    Mech_Prune = mildew4$Mech_Prune)

fit4a <- stan(file="mildew4.stan",data=mildewData4, iter = 1000, chains = 4)

traceplot(fit4a, "beta")
print(fit4a, "beta")

fit4b <- stan(fit=fit4a, data=mildewData4, iter=10000, chains=4)
traceplot(fit4b, par = "beta", inc_warmup=F)
results4 <- extract(fit4b, pars="beta",permuted = F, inc_warmup = FALSE) 
str(results4) 
print(fit4b, par = "beta")

Beta4 <- data.frame(Intercept = as.vector(results4[,1:4,1]),
                    X1 = as.vector(results4[,1:4,2]),
                    X2 = as.vector(results4[,1:4,3]),
                    X3 = as.vector(results4[,1:4,4]))

OutputMatrix4 <- data.frame(Risk = as.factor(c("EST","2nd","EST NM","2nd NM",
                                               "EST SUS","2nd SUS","EST SUS+NM",
                                               "2nd SUS+NM")),
                            ProbEst = numeric(8),
                            Lower = numeric(8),
                            Upper = numeric(8))

### 000
est = exp(Beta4[,1])/(1+exp(Beta4[,1]))
OutputMatrix4$ProbEst[1] = round(median(est),4)
OutputMatrix4[1,3:4] = round(emp.hpd(est),4)

### 100
est = exp(Beta4[,1] + Beta4[,2])/(1+exp(Beta4[,1] + Beta4[,2]))
OutputMatrix4$ProbEst[2] = round(median(est),4)
OutputMatrix4[2,3:4] = round(emp.hpd(est),4)

### 001
est = exp(Beta4[,1] + Beta4[,4])/(1+exp(Beta4[,1] + Beta4[,4]))
OutputMatrix4$ProbEst[3] = round(median(est),4)
OutputMatrix4[3,3:4] = round(emp.hpd(est),4)

### 101
est = exp(Beta4[,1] + Beta4[,2] + Beta4[,4])/(1+exp(Beta4[,1] + Beta4[,2] + 
                                                      Beta4[,4]))
OutputMatrix4$ProbEst[4] = round(median(est),4)
OutputMatrix4[4,3:4] = round(emp.hpd(est),4)

### 010
est = exp(Beta4[,1] + Beta4[,3])/(1+exp(Beta4[,1] + Beta4[,3]))
OutputMatrix4$ProbEst[5] = round(median(est),4)
OutputMatrix4[5,3:4] = round(emp.hpd(est),4)

### 110
est = exp(Beta4[,1] + Beta4[,2] + Beta4[,3])/(1+exp(Beta4[,1] + Beta4[,2] + 
                                                      Beta4[,3]))
OutputMatrix4$ProbEst[6] = round(median(est),4)
OutputMatrix4[6,3:4] = round(emp.hpd(est),4)

### 011
est = exp(Beta4[,1] + Beta4[,3] + Beta4[,4])/(1+exp(Beta4[,1] + Beta4[,3] + 
                                                      Beta4[,4]))
OutputMatrix4$ProbEst[7] = round(median(est),4)
OutputMatrix4[7,3:4] = round(emp.hpd(est),4)

### 111
est = exp(Beta4[,1] + Beta4[,2] + Beta4[,3] + Beta4[,4])/(1+exp(Beta4[,1] + 
                                            Beta4[,2] + Beta4[,3] + Beta4[,4]))
OutputMatrix4$ProbEst[8] = round(median(est),4)
OutputMatrix4[8,3:4] = round(emp.hpd(est),4)

Name_vec4 = c("EST","2nd","EST NM","2nd NM","EST SUS","2nd SUS","EST SUS+NM",
              "2nd SUS+NM")

OutputMatrix4$Risk = factor(OutputMatrix4$Risk,levels=Name_vec4, ordered = T)

###############################################################################
# Model 5.

mod5 <- glm(FS_Binary ~ Second_Year + Suscept_Binary + Poor_Prune + Mech_Prune, 
            data = mildew, family=binomial)
summary(mod5)

idx5 <- 1:dim(mildew)[1]
na.idx5 <- idx5[with(mildew, is.na(FS_Binary) | is.na(Second_Year) | 
                is.na(Suscept_Binary) | is.na(Poor_Prune) | is.na(Mech_Prune))]
length(idx5[-na.idx5])

mildew5 <- with(mildew, data.frame(FS_Binary[-na.idx5], Second_Year[-na.idx5], 
                                   Suscept_Binary[-na.idx5], 
                                   Poor_Prune[-na.idx5], Mech_Prune[-na.idx5]))
names(mildew5) <- c("FS_Binary", "Second_Year", "Suscept_Binary", "Poor_Prune",
                    "Mech_Prune")

mildewData5 <- list(N = nrow(mildew5),
                    p = 5,
                    FS_Binary = mildew5$FS_Binary,
                    Second_Year = mildew5$Second_Year,
                    Suscept_Binary = mildew5$Suscept_Binary,
                    Poor_Prune = mildew5$Poor_Prune,
                    Mech_Prune = mildew5$Mech_Prune)

fit5a <- stan(file="mildew5.stan",data=mildewData5, iter = 1000, chains = 4)

traceplot(fit5a, "beta")
print(fit5a, "beta")

fit5b <- stan(fit=fit5a, data=mildewData5, iter=10000, chains=4)
traceplot(fit5b, par = "beta", inc_warmup=F)
results5 <- extract(fit5b, pars="beta",permuted = F, inc_warmup = FALSE) 
str(results5) 
print(fit5b, par = "beta")

Beta5 <- data.frame(Intercept = as.vector(results5[,1:4,1]),
                    X1 = as.vector(results5[,1:4,2]),
                    X2 = as.vector(results5[,1:4,3]),
                    X3 = as.vector(results5[,1:4,4]),
                    X4 = as.vector(results5[,1:4,5]))

OutputMatrix5 <- data.frame(Risk = as.factor(c("1 Est","2 2yo","3 Est","4 2yo",
                                               "5 Est","6 2yo","7 Est","8 2yo",
                                               "9 Est","10 2yo","11 Est",
                                               "12 2yo","13 Est","14 2yo",
                                               "15 Est","16 2yo")),
                            ProbEst = numeric(16),
                            Lower = numeric(16),
                            Upper = numeric(16))

### 0000
est = exp(Beta5[,1])/(1+exp(Beta5[,1]))
OutputMatrix5$ProbEst[1] = round(median(est),5)
OutputMatrix5[1,3:4] = round(emp.hpd(est),5)

### 1000
est = exp(Beta5[,1] + Beta5[,2])/(1+exp(Beta5[,1] + Beta5[,2]))
OutputMatrix5$ProbEst[2] = round(median(est),5)
OutputMatrix5[2,3:4] = round(emp.hpd(est),5)

### 0010
est = exp(Beta5[,1] + Beta5[,4])/(1+exp(Beta5[,1] + Beta5[,4]))
OutputMatrix5$ProbEst[3] = round(median(est),5)
OutputMatrix5[3,3:4] = round(emp.hpd(est),5)

### 1010
est = exp(Beta5[,1] + Beta5[,2] + Beta5[,4])/(1+exp(Beta5[,1] + Beta5[,2] + 
                                                      Beta5[,4]))
OutputMatrix5$ProbEst[4] = round(median(est),5)
OutputMatrix5[4,3:4] = round(emp.hpd(est),5)

### 0001
est = exp(Beta5[,1] + Beta5[,5])/(1+exp(Beta5[,1] + Beta5[,5]))
OutputMatrix5$ProbEst[5] = round(median(est),5)
OutputMatrix5[5,3:4] = round(emp.hpd(est),5)

### 1001
est = exp(Beta5[,1] + Beta5[,2] + Beta5[,5])/(1+exp(Beta5[,1] + Beta5[,2] + 
                                                      Beta5[,5]))
OutputMatrix5$ProbEst[6] = round(median(est),5)
OutputMatrix5[6,3:4] = round(emp.hpd(est),5)

### 0011
est = exp(Beta5[,1] + Beta5[,4] + Beta5[,5])/(1+exp(Beta5[,1] + Beta5[,4] + 
                                                      Beta5[,5]))
OutputMatrix5$ProbEst[7] = round(median(est),5)
OutputMatrix5[7,3:4] = round(emp.hpd(est),5)

### 1011
est = exp(Beta5[,1] + Beta5[,2] + Beta5[,4] + Beta5[,5])/(1+exp(Beta5[,1] + 
                                            Beta5[,2] + Beta5[,4] + Beta5[,5]))
OutputMatrix5$ProbEst[8] = round(median(est),5)
OutputMatrix5[8,3:4] = round(emp.hpd(est),5)

### 0100
est = exp(Beta5[,1] + Beta5[,3])/(1+exp(Beta5[,1] + Beta5[,3]))
OutputMatrix5$ProbEst[9] = round(median(est),5)
OutputMatrix5[9,3:4] = round(emp.hpd(est),5)

### 1100
est = exp(Beta5[,1] + Beta5[,2] + Beta5[,3])/(1+exp(Beta5[,1] + Beta5[,2] + 
                                                      Beta5[,3]))
OutputMatrix5$ProbEst[10] = round(median(est),5)
OutputMatrix5[10,3:4] = round(emp.hpd(est),5)

### 0110
est = exp(Beta5[,1] + Beta5[,3] + Beta5[,4])/(1+exp(Beta5[,1] + Beta5[,3] + 
                                                      Beta5[,4]))
OutputMatrix5$ProbEst[11] = round(median(est),5)
OutputMatrix5[11,3:4] = round(emp.hpd(est),5)

### 1110
est = exp(Beta5[,1] + Beta5[,2] + Beta5[,3] + Beta5[,4])/(1+exp(Beta5[,1] + 
                                            Beta5[,2] + Beta5[,3] + Beta5[,4]))
OutputMatrix5$ProbEst[12] = round(median(est),5)
OutputMatrix5[12,3:4] = round(emp.hpd(est),5)

### 0101
est = exp(Beta5[,1] + Beta5[,3] + Beta5[,5])/(1+exp(Beta5[,1] + Beta5[,3] + 
                                                      Beta5[,5]))
OutputMatrix5$ProbEst[13] = round(median(est),5)
OutputMatrix5[13,3:4] = round(emp.hpd(est),5)

### 1101
est = exp(Beta5[,1] + Beta5[,2] + Beta5[,3] + Beta5[,5])/(1+exp(Beta5[,1] + 
                                            Beta5[,2] + Beta5[,3] + Beta5[,5]))
OutputMatrix5$ProbEst[14] = round(median(est),5)
OutputMatrix5[14,3:4] = round(emp.hpd(est),5)

### 0111
est = exp(Beta5[,1] + Beta5[,3] + Beta5[,4] + Beta5[,5])/(1+exp(Beta5[,1] + 
                                            Beta5[,3] + Beta5[,4] + Beta5[,5]))
OutputMatrix5$ProbEst[15] = round(median(est),5)
OutputMatrix5[15,3:4] = round(emp.hpd(est),5)

### 1111
est = exp(Beta5[,1] + Beta5[,2] + Beta5[,3] + Beta5[,4] + Beta5[,5])/
  (1+exp(Beta5[,1] + Beta5[,2] + Beta5[,3] + Beta5[,4] + Beta5[,5]))
OutputMatrix5$ProbEst[16] = round(median(est),5)
OutputMatrix5[16,3:4] = round(emp.hpd(est),5)

Name_vec0 = c("1 Est","2 2yo","5 Est","6 2yo","3 Est","4 2yo","7 Est","8 2yo",
              "9 Est","10 2yo","13 Est","14 2yo","11 Est","12 2yo","15 Est",
              "16 2yo")

OutputMatrix5$Risk = factor(OutputMatrix5$Risk,levels=Name_vec0, ordered = T)

###############################################################################
# Plots, analysis, etc.
# Manual data input from 'BetaX' and 'OutputMatrixX' data frames for some plots.

###############################################################################
# Histogram with number of flag shoots/transect.
# New data set created to include actual number of flag shoots in results.

FS <- read.csv("Flagshoots.csv")

gghistogram(FS, "Flag1", bins = 90, xlab = "Number of flag shoots", ylab = 
              "Number of transects") + 
  geom_density(adjust = 5, color="cadetblue3", fill="cadetblue3", alpha = 0.4) +
  scale_y_break(c(40, 2800)) + 
  scale_y_continuous(breaks = seq(0, 2820, by = 10), limits = c(0, 2820))

###############################################################################
# Box plot describing yard age in each of the sampling years and flag shoot-
# -presence or absence.
# New data set created to remove NA's.

Yardage <- read.csv("Yardage.csv")
YEAR <- as.factor(Yardage$Year)
FS_BINARY <- as.factor(Yardage$FS_Binary)

ggplot(Yardage, aes(x = YEAR, y = Yard_Age, fill = FS_BINARY)) +
  geom_boxplot(outlier.size = 0, outlier.colour = "white") +
  scale_fill_manual(values=c("cadetblue3", "cadetblue1")) +
  guides(fill=guide_legend(title="Flag shoot occurrence")) +
  geom_jitter(shape = 16, size = 0.4, position=position_jitter(0.2)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_x_discrete(name ="Year") +
  scale_y_continuous(name ="Yard age") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

###############################################################################
# Histograms showing 95% HPD intervals for each parameter.
# One-factor.

onehist <- ggplot(Beta1) + xlim(-7, -2) +
  geom_histogram(bins = 60, aes(x = Intercept), color = "dark blue", fill = 
                   "cadetblue3") + 
  geom_vline(aes(xintercept = -3.07), linetype = "dashed") +
  geom_vline(aes(xintercept = -2.43), linetype = "dashed") +
  geom_vline(aes(xintercept = median(x = Beta1$Intercept)), linetype = 
               "longdash") +
  labs(y = "Count", x = "") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) 

onehist1 <- ggplot(Beta1) +  xlim(-2, 1.5) +
  geom_histogram(bins = 60, aes(x = X1), color = "dark blue", fill = 
                   "cadetblue3") + 
  geom_vline(aes(xintercept = -0.90), linetype = "dashed") +
  geom_vline(aes(xintercept = 0.93), linetype = "dashed") +
  geom_vline(aes(xintercept = median(x = Beta1$X1)), linetype = "longdash") +
  labs(y = "", x = "") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) 

onehist
onehist1

###############################################################################
# Histograms showing 95% HPD intervals for each parameter.
# Two-factor.

twohist <- ggplot(Beta2) + xlim(-7, -2) + 
  geom_histogram(bins = 60, aes(x = Intercept), color = "dark blue", fill = 
                   "cadetblue3") + 
  geom_vline(aes(xintercept = -4.72), linetype = "dashed") +
  geom_vline(aes(xintercept = -2.82), linetype = "dashed") +
  geom_vline(aes(xintercept = median(x = Beta2$Intercept)), linetype = 
               "longdash") +
  labs(y = "Count", x = "") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) 

twohist1 <- ggplot(Beta2) + xlim(-2, 1.5) +
  geom_histogram(bins = 60, aes(x = X1), color = "dark blue", fill = 
                   "cadetblue3") + 
  geom_vline(aes(xintercept = -0.74), linetype = "dashed") +
  geom_vline(aes(xintercept = 1.10), linetype = "dashed") +
  geom_vline(aes(xintercept = median(x = Beta2$X1)), linetype = "longdash") +
  labs(y = "", x = "") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) 

twohist2 <- ggplot(Beta2) + xlim(-0.5, 3.5) +
  geom_histogram(bins = 60, aes(x = X2), color = "dark blue", fill = 
                   "cadetblue3") + 
  geom_vline(aes(xintercept = 0.20), linetype = "dashed") +
  geom_vline(aes(xintercept = 2.17), linetype = "dashed") +
  geom_vline(aes(xintercept = median(x = Beta2$X2)), linetype = "longdash") +
  labs(y = "", x = "") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) 

twohist
twohist1
twohist2

###############################################################################
# Histograms showing 95% HPD intervals for each parameter.
# Incomplete pruning.

threehist <- ggplot(Beta3) + xlim(-7, -2) + 
  geom_histogram(bins = 60, aes(x = Intercept), color = "dark blue", fill = 
                   "cadetblue3") + 
  geom_vline(aes(xintercept = -5.58), linetype = "dashed") +
  geom_vline(aes(xintercept = -3.43), linetype = "dashed") +
  geom_vline(aes(xintercept = median(x = Beta3$Intercept)), linetype = 
               "longdash") +
  labs(y = "Count", x = "") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) 

threehist1 <- ggplot(Beta3) + xlim(-2, 1.5) +
  geom_histogram(bins = 60, aes(x = X1), color = "dark blue", fill = 
                   "cadetblue3") + 
  geom_vline(aes(xintercept = -0.81), linetype = "dashed") +
  geom_vline(aes(xintercept = 1.08), linetype = "dashed") +
  geom_vline(aes(xintercept = median(x = Beta3$X1)), linetype = "longdash") +
  labs(y = "", x = "") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) 

threehist2 <- ggplot(Beta3) + xlim(-0.5, 3.5) +
  geom_histogram(bins = 60, aes(x = X2), color = "dark blue", fill = 
                   "cadetblue3") + 
  geom_vline(aes(xintercept = 0.28), linetype = "dashed") +
  geom_vline(aes(xintercept = 2.24), linetype = "dashed") +
  geom_vline(aes(xintercept = median(x = Beta3$X2)), linetype = "longdash") +
  labs(y = "", x = "") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) 

threehist3 <- ggplot(Beta3) + xlim(-0.5, 2.5) +
  geom_histogram(bins = 60, aes(x = X3), color = "dark blue", fill = 
                   "cadetblue3") + 
  geom_vline(aes(xintercept = 0.33), linetype = "dashed") +
  geom_vline(aes(xintercept = 1.74), linetype = "dashed") +
  geom_vline(aes(xintercept = median(x = Beta3$X3)), linetype = "longdash") +
  labs(y = "", x = "") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) 

threehist
threehist1
threehist2
threehist3

###############################################################################
# Histograms showing 95% HPD intervals for each parameter.
# Non-mechanical model.

fourhist <- ggplot(Beta4) + xlim(-7, -2) + 
  geom_histogram(bins = 60, aes(x = Intercept), color = "dark blue", fill = 
                   "cadetblue3") + 
  geom_vline(aes(xintercept = -5.80), linetype = "dashed") +
  geom_vline(aes(xintercept = -3.41), linetype = "dashed") +
  geom_vline(aes(xintercept = median(x = Beta4$Intercept)), linetype = 
               "longdash") +
  labs(y = "Count", x = "") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) 

fourhist1 <- ggplot(Beta4) + xlim(-2, 1.5) +
  geom_histogram(bins = 60, aes(x = X1), color = "dark blue", fill = 
                   "cadetblue3") + 
  geom_vline(aes(xintercept = -0.86), linetype = "dashed") +
  geom_vline(aes(xintercept = 1.00), linetype = "dashed") +
  geom_vline(aes(xintercept = median(x = Beta4$X1)), linetype = "longdash") +
  labs(y = "", x = "") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) 

fourhist2 <- ggplot(Beta4) + xlim(-0.5, 3.5) +
  geom_histogram(bins = 60, aes(x = X2), color = "dark blue", fill = 
                   "cadetblue3") + 
  geom_vline(aes(xintercept = 0.31), linetype = "dashed") +
  geom_vline(aes(xintercept = 2.30), linetype = "dashed") +
  geom_vline(aes(xintercept = median(x = Beta4$X2)), linetype = "longdash") +
  labs(y = "", x = "") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) 

fourhist3 <- ggplot(Beta4) + xlim(-1, 2.5) +
  geom_histogram(bins = 60, aes(x = X3), color = "dark blue", fill = 
                   "cadetblue3") + 
  geom_vline(aes(xintercept = 0.21), linetype = "dashed") +
  geom_vline(aes(xintercept = 1.83), linetype = "dashed") +
  geom_vline(aes(xintercept = median(x = Beta4$X3)), linetype = "longdash") +
  labs(y = "", x = "") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) 

fourhist
fourhist1
fourhist2
fourhist3

###############################################################################
# Histograms showing 95% HPD intervals for each parameter.
# Full model.

fivehist <- ggplot(Beta5) + xlim(-7, -2) +
  geom_histogram(bins = 60, aes(x = Intercept), color = "dark blue", fill = 
                   "cadetblue3") + 
  geom_vline(aes(xintercept = -6.18), linetype = "dashed") +
  geom_vline(aes(xintercept = -3.69), linetype = "dashed") +
  geom_vline(aes(xintercept = median(x = Beta5$Intercept)), linetype = 
               "longdash") +
  labs(y = "Count", x = "Parameter value") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) 

fivehist1 <- ggplot(Beta5) + xlim(-2, 1.5) +
  geom_histogram(bins = 60, aes(x = X1), color = "dark blue", fill = 
                   "cadetblue3") + 
  geom_vline(aes(xintercept = -0.88), linetype = "dashed") +
  geom_vline(aes(xintercept = 1.00), linetype = "dashed") +
  geom_vline(aes(xintercept = median(x = Beta5$X1)), linetype = "longdash") +
  labs(y = "", x = "Parameter value") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) 

fivehist2 <- ggplot(Beta5) + xlim(-0.5, 3.5) +
  geom_histogram(bins = 60, aes(x = X2), color = "dark blue", fill = 
                   "cadetblue3") + 
  geom_vline(aes(xintercept = 0.31), linetype = "dashed") +
  geom_vline(aes(xintercept = 2.32), linetype = "dashed") +
  geom_vline(aes(xintercept = median(x = Beta5$X2)), linetype = "longdash") +
  labs(y = "", x = "Parameter value") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) 

fivehist3 <- ggplot(Beta5) + xlim(-0.5, 2.5) +
  geom_histogram(bins = 60, aes(x = X3), color = "dark blue", fill = 
                   "cadetblue3") + 
  geom_vline(aes(xintercept = 0.09), linetype = "dashed") +
  geom_vline(aes(xintercept = 1.57), linetype = "dashed") +
  geom_vline(aes(xintercept = median(x = Beta5$X3)), linetype = "longdash") +
  labs(y = "", x = "Parameter value") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) 

fivehist4 <- ggplot(Beta5) + xlim(-1, 2.5) +
  geom_histogram(bins = 60, aes(x = X4), color = "dark blue", fill = 
                   "cadetblue3") + 
  geom_vline(aes(xintercept = -0.13), linetype = "dashed") +
  geom_vline(aes(xintercept = 1.59), linetype = "dashed") +
  geom_vline(aes(xintercept = median(x = Beta5$X4)), linetype = "longdash")+
  labs(y = "", x = "Parameter value") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) 

fivehist
fivehist1
fivehist2
fivehist3
fivehist4

###############################################################################
# Histograms showing 95% HPD intervals for each parameter.
# Arrangement of all models (Figure 1).

figure1 <- ggarrange(onehist, onehist1,
                     labels = c("Intercept", "Yard age"),
                     font.label = list(size = 15),
                     ncol = 2, nrow = 1)

figure2 <- ggarrange(twohist, twohist1, twohist2,
                     labels = c("Intercept", "Yard age", "Susceptibility"),
                     font.label = list(size = 15),
                     ncol = 3, nrow = 1)

figure3 <- ggarrange(threehist, threehist1, threehist2, threehist3,
                     labels = c("Intercept", "Yard age", "Susceptibility", 
                                "Pruning quality"),
                     font.label = list(size = 15),
                     ncol = 4, nrow = 1)

figure4 <- ggarrange(fourhist, fourhist1, fourhist2, fourhist3,
                     labels = c("Intercept", "Yard age", "Susceptibility", 
                                "Pruning type"),
                     font.label = list(size = 15),
                     ncol = 4, nrow = 1)

figure5 <- ggarrange(fivehist, fivehist1, fivehist2, fivehist4, fivehist3,
                     labels = c("Intercept", "Yard age", "Susceptibility", 
                                "Pruning type", "Pruning quality"),
                     font.label = list(size = 15),
                     ncol = 5, nrow = 1)

figure1
figure2
figure3
figure4
figure5

###############################################################################
# Markov Chain Monte Carlo visualizations for full model.

t1 <- mcmc_trace(fit5b, c("beta[1]"), iter1 = 5001) + ggtitle("Trace plots")
t2 <- mcmc_trace(fit5b, c("beta[2]"), iter1 = 5001)
t3 <- mcmc_trace(fit5b, c("beta[3]"), iter1 = 5001) 
t4 <- mcmc_trace(fit5b, c("beta[4]"), iter1 = 5001) 
t5 <- mcmc_trace(fit5b, c("beta[5]"), iter1 = 5001) 

trace1 <- ggarrange(t1, t2, t3, t4, t5, 
                    ncol = 1, nrow = 5) 

rankoverlay <- mcmc_rank_overlay(fit5b, c("beta[1]", "beta[2]", "beta[3]", 
                                          "beta[4]", "beta[5]")) + 
                                          ggtitle("Rank overlay")

rankhist <- mcmc_rank_hist(fit5b, c("beta[1]", "beta[2]", "beta[3]", "beta[4]", 
                        "beta[5]"), ref_line = TRUE) + ggtitle("Rank histogram")

areas <- mcmc_areas(fit5b, vars(starts_with("beta")), area_method = 
                      "equal area") +  ggtitle("Density curves")

acf <- mcmc_acf(fit5b, c("beta[1]", "beta[2]", "beta[3]", "beta[4]", "beta[5]"))

trace1
rankoverlay
rankhist
areas
acf

###############################################################################
# Full model posterior probability point estimate plot (Figure 2).

ggplot(OutputMatrix5, aes(x = Risk, y = ProbEst)) +
  ylab("Probability estimate") +
  xlab("Yard age") +
  geom_point(aes(shape = Risk, fill = Risk), size=3) +
  scale_shape_manual(values=c(21,21,21,21,21,21,21,21,21,21,21,21,21,21,21,21))+
  scale_fill_manual(values=c('lightblue','lightblue','deepskyblue3','deepskyblue3',
                             'blueviolet','blueviolet','red','red',
                             'lightblue','lightblue','deepskyblue3','deepskyblue3',
                             'blueviolet','blueviolet','red','red')) +
  theme(legend.position='none') +
  geom_linerange(aes(ymin = Lower, ymax = Upper)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_segment(aes(x="1 Est",xend="8 2yo",y = 0.224, yend = 0.224)) +
  geom_segment(aes(x="9 Est",xend="16 2yo",y = 0.224, yend = 0.224)) +
  geom_text(x = "6 2yo", y = 0.231, label = "Low susceptibility cultivars") +
  geom_text(x = "14 2yo", y = 0.231, label = "Susceptible cultivars") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank()) 

###############################################################################
# Chi-squared tests (Table 3).

# Yard age * thoroughness.

chi1 <- matrix(c(35, 49, 280, 373),nrow=2,ncol=2)
chisq.test(chi1,correct=FALSE)

# Yard age * pruning type.

chi2 <- matrix(c(14, 71, 218, 435),nrow=2,ncol=2)
chisq.test(chi2,correct=FALSE)

# Pruning type * thoroughness.

chi3 <- matrix(c(174, 75, 156, 361),nrow=2,ncol=2)
chisq.test(chi3,correct=FALSE)

# Tables.
chi1
chi2
chi3

###############################################################################