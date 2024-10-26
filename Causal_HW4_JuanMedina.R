library(causaldata)
#1
data("nhefs_complete")
#a)
#For this mediation analysis I will be exploring the effect of having diabetes in 1971 on weight change,
#mediated through income.

#b)
#Make household income dichotomous (above median or below median) at $10,000
table(nhefs_complete$income)
nhefs_complete$income_binary <- ifelse(nhefs_complete$income %in% c("19","20","21","22"), 1, 0)
table(nhefs_complete$income_binary)

#Select subset of observations where exposure is not missing (diabetes!=2)
nhefs_subset <- subset(nhefs_complete,diabetes != "2")

# Outcome model (continuous outcome):
ymodel <- glm(wt82_71 ~ diabetes*income_binary + sex + age + race, family=gaussian(link='identity'), data=nhefs_subset)
summary(ymodel)

# Mediator model:
zmodel <- glm(income_binary ~ diabetes + sex + age + race, family=binomial(link='logit'), data=nhefs_subset)
summary(zmodel)

# Expected potential outcomes:
newdat11 <- data.frame(diabetes=1,income_binary=1,sex=nhefs_subset$sex, age=nhefs_subset$age, race=nhefs_subset$race)
predy11 <- predict(ymodel, newdata=newdat11)

newdat10 <- data.frame(diabetes=1,income_binary=0,sex=nhefs_subset$sex, age=nhefs_subset$age, race=nhefs_subset$race)
predy10 <- predict(ymodel, newdata=newdat10)

newdat01 <- data.frame(diabetes=0,income_binary=1,sex=nhefs_subset$sex, age=nhefs_subset$age, race=nhefs_subset$race)
predy01 <- predict(ymodel, newdata=newdat01)

newdat00 <- data.frame(diabetes=0,income_binary=0,sex=nhefs_subset$sex, age=nhefs_subset$age, race=nhefs_subset$race)
predy00 <- predict(ymodel, newdata=newdat00)

# Expected potential mediators:
predz1 <- predict(zmodel, newdata=newdat11, type='response')
predz0 <- predict(zmodel, newdata=newdat00, type='response')

# Calculate the natural direct effect estimate:
nde <- mean(predy11 * predz0 + predy10 * (1.0 - predz0)) - mean(predy01 * predz0 + predy00 * (1.0 - predz0))
nde

# Calculate the natural indirect effect estimate:
nie <- mean(predy11 * predz1 + predy10 * (1.0 - predz1)) - mean(predy11 * predz0 + predy10 * (1.0 - predz0))
nie

# Total effect:
te <- nde + nie
te

#Proportion mediated
p_med <- nie/te
p_med

#Bootstrap 95% CI 
nobs<- nrow(nhefs_subset)
nb <- 1000
nie<- rep(NA, nb)
nde<- rep(NA, nb)
te<- rep(NA, nb)
p_med<- rep(NA, nb)
set.seed(1)
for (i in 1:nb) {
  bootidx <- sample(c(1:nobs), size=nobs, replace=TRUE)
  nhefs_boot <- nhefs_subset[bootidx,]
  # Outcome model (continuous outcome):
  ymodel <- glm(wt82_71 ~ diabetes*income_binary + sex + age + race, family=gaussian(link='identity'), data=nhefs_boot)
  # Mediator model:
  zmodel <- glm(income_binary ~ diabetes + sex + age + race, family=binomial(link='logit'), data=nhefs_boot)
  summary(zmodel)
  
  # Expected potential outcomes:
  newdat11 <- data.frame(diabetes=1,income_binary=1,sex=nhefs_boot$sex, age=nhefs_boot$age, race=nhefs_boot$race)
  predy11 <- predict(ymodel, newdata=newdat11)
  
  newdat10 <- data.frame(diabetes=1,income_binary=0,sex=nhefs_boot$sex, age=nhefs_boot$age, race=nhefs_boot$race)
  predy10 <- predict(ymodel, newdata=newdat10)
  
  newdat01 <- data.frame(diabetes=0,income_binary=1,sex=nhefs_boot$sex, age=nhefs_boot$age, race=nhefs_boot$race)
  predy01 <- predict(ymodel, newdata=newdat01)
  
  newdat00 <- data.frame(diabetes=0,income_binary=0,sex=nhefs_boot$sex, age=nhefs_boot$age, race=nhefs_boot$race)
  predy00 <- predict(ymodel, newdata=newdat00)
  
  # Expected potential mediators:
  predz1 <- predict(zmodel, newdata=newdat11, type='response')
  predz0 <- predict(zmodel, newdata=newdat00, type='response')
  
  # Calculate the natural direct effect estimate:
  nde_boot <- mean(predy11 * predz0 + predy10 * (1.0 - predz0)) - mean(predy01 * predz0 + predy00 * (1.0 - predz0))
  nde[i] <- nde_boot
  
  # Calculate the natural indirect effect estimate:
  nie_boot <- mean(predy11 * predz1 + predy10 * (1.0 - predz1)) - mean(predy11 * predz0 + predy10 * (1.0 - predz0))
  nie[i] <- nie_boot
  # Total effect:
  te_boot <- nde_boot + nie_boot
  te[i] <- te_boot
  #Proportion mediated
  p_med_boot <- nie_boot/te_boot
  p_med[i] <- p_med_boot
  
}
quantile(nde, prob=c(0.025,0.975))
quantile(nie, prob=c(0.025,0.975))
quantile(te, prob=c(0.025,0.975))
quantile(p_med, prob=c(0.025,0.975))





#2
#b
sim.r <- function(samplesize = 500)
{
  set.seed(123)
  expit <- function(x){exp(x)/(1+exp(x))}
  #covariates;
  L1 <- runif(n=samplesize,0,1)
  L2 <- runif(n=samplesize,0,1)
  U <- rbinom(n=samplesize, size = 1, prob = 0.2)
  
  #treatment;
  Aprob <- expit(3*L1-3*L2+4*U)
  A <- rbinom(n=samplesize, size=1, prob=Aprob)
  
  #outcome;
  Yprob <- expit(3*A+3*L1-3*L2+4*U)
  Y <- rbinom(n = samplesize, size = 1, prob = Yprob)
  dat <- cbind(L1, L2, U, A, Y)
  dat <- data.frame(dat)
  return(dat)
}

simdat <- sim.r(samplesize = 500)
summary(simdat)

#i
library(survey)
#Method 1
ps_model<- glm(A ~ L1 + L2, data=simdat, family=binomial(link="logit"))
simdat$ps <- predict(ps_model, type="response")
simdat$iptw <- ifelse(simdat$A == 1, 1 / simdat$ps, 1 / (1 - simdat$ps))
simdat$iptw
 #Model
design <- svydesign(ids = ~1, weights = ~iptw, data = simdat)
model <- svyglm(Y ~ A, family = binomial(link = "logit"), design = design)
#Predict Risk Ratio
rr <- mean(predict(model, newdata = data.frame('A'=1), type = "response"))/mean(predict(model, newdata = data.frame('A'=0), type = "response"))
rr

#Bootstrap CI
nobs<- nrow(simdat)
nb <- 1000
rr_b<- rep(NA, nb)
set.seed(1)
for (i in 1:nb) {
  bootidx <- sample(c(1:nobs), size=nobs, replace=TRUE)
  simdat_boot <- simdat[bootidx,]
  ps_model<- glm(A ~ L1 + L2, data=simdat_boot, family=binomial(link="logit"))
  simdat_boot$ps <- predict(ps_model, type="response")
  simdat_boot$iptw <- ifelse(simdat_boot$A == 1, 1 / simdat_boot$ps, 1 / (1 - simdat_boot$ps))
  #Model
  design_boot <- svydesign(ids = ~1, weights = ~iptw, data = simdat_boot)
  model_boot <- svyglm(Y ~ A, family = binomial(link = "logit"), design = design_boot)
  #Predict Risk Ratio
  rr_b[i] <- mean(predict(model_boot, newdata = data.frame('A'=1), type = "response"))/mean(predict(model_boot, newdata = data.frame('A'=0), type = "response"))
}
quantile(rr_b, prob=c(0.025,0.975))
warnings()

#Method II
ps_model<- glm(A ~ L1 + L2, data=simdat, family=binomial(link="logit"))
ps <- predict(ps_model, type="response")
risk_ratio <- mean((simdat$A==1)*(simdat$Y/ps))/mean((simdat$A==0)*(simdat$Y/(1-ps)))
#risk_ratio <- mean((simdat$A == 1) * simdat$Y * mean(simdat$A) / ps) / mean((simdat$A == 0) * simdat$Y * (1 - mean(simdat$A)) / (1 - ps))
risk_ratio
#Bootstap 95% CI
nobs<- nrow(simdat)
nb <- 1000
risk_ratio_b<- rep(NA, nb)
set.seed(1)
for (i in 1:nb) {
  bootidx <- sample(c(1:nobs), size=nobs, replace=TRUE)
  simdat_boot <- simdat[bootidx,]
  ps_model_b<- glm(A ~ L1 + L2, data=simdat_boot, family=binomial(link="logit"))
  psb <- predict(ps_model_b, type="response")
  risk_ratio_b[i] <- mean((simdat$A==1)*(simdat$Y/psb))/mean((simdat$A==0)*(simdat$Y/(1-psb)))
}
quantile(risk_ratio_b, prob=c(0.025,0.975))
#Method III
library(WeightIt)
IPTW <- weightit(A ~ L1 + L2,
                 data = simdat,
                 method = "glm",
                 stabilize = TRUE)

summary(IPTW)
risk_ratio <- mean((simdat$A==1)*(simdat$Y*IPTW$weights))/mean((simdat$A==0)*(simdat$Y*IPTW$weights))
risk_ratio
#Bootstap 95% CI
nobs<- nrow(simdat)
nb <- 1000
risk_ratio_b<- rep(NA, nb)
set.seed(1)
for (i in 1:nb) {
  bootidx <- sample(c(1:nobs), size=nobs, replace=TRUE)
  simdat_boot <- simdat[bootidx,]
  IPTW_boot <- weightit(A ~ L1 + L2,
                   data = simdat_boot,
                   method = "glm",
                   stabilize = TRUE)
  risk_ratio_b[i] <- risk_ratio <- mean((simdat_boot$A==1)*(simdat_boot$Y*IPTW_boot$weights))/mean((simdat_boot$A==0)*(simdat_boot$Y*IPTW_boot$weights))
}
quantile(risk_ratio_b, prob=c(0.025,0.975))
#ii
library(EValue)
E_val <- evalues.RR(3.029558, lo=2.430075, hi=3.844333)
E_val
#iii
ps_model<- glm(A ~ L1 + L2 + U, data=simdat, family=binomial(link="logit"))
ps <- predict(ps_model, type="response")
risk_ratio <- mean((simdat$A==1)*(simdat$Y/ps))/mean((simdat$A==0)*(simdat$Y/(1-ps)))
risk_ratio
#Bootstap 95% CI
nobs<- nrow(simdat)
nb <- 1000
risk_ratio_b<- rep(NA, nb)
set.seed(1)
for (i in 1:nb) {
  bootidx <- sample(c(1:nobs), size=nobs, replace=TRUE)
  simdat_boot <- simdat[bootidx,]
  ps_model_b<- glm(A ~ L1 + L2 + U, data=simdat_boot, family=binomial(link="logit"))
  psb <- predict(ps_model_b, type="response")
  risk_ratio_b[i] <- mean((simdat$A==1)*(simdat$Y/psb))/mean((simdat$A==0)*(simdat$Y/(1-psb)))
}
quantile(risk_ratio_b, prob=c(0.025,0.975))
#Method II
library(WeightIt)
IPTW <- weightit(A ~ L1 + L2 + U,
                 data = simdat,
                 method = "glm",
                 stabilize = TRUE)
#Notice some extreme weights
summary(IPTW)
#Trim weights
IPTW.trim <- trim(IPTW, at = .99)
summary(IPTW.trim)
risk_ratio <- mean((simdat$A==1)*(simdat$Y*IPTW.trim$weights))/mean((simdat$A==0)*(simdat$Y*IPTW.trim$weights))
risk_ratio
#Bootstap 95% CI
nobs<- nrow(simdat)
nb <- 1000
risk_ratio_b<- rep(NA, nb)
set.seed(1)
for (i in 1:nb) {
  bootidx <- sample(c(1:nobs), size=nobs, replace=TRUE)
  simdat_boot <- simdat[bootidx,]
  IPTW_boot <- weightit(A ~ L1 + L2 + U,
                        data = simdat_boot,
                        method = "glm",
                        stabilize = TRUE)
  #Trim weights
  IPTW_boot.trim <- trim(IPTW_boot, at = .99)
  risk_ratio_b[i] <- risk_ratio <- mean((simdat_boot$A==1)*(simdat_boot$Y*IPTW_boot.trim$weights))/mean((simdat_boot$A==0)*(simdat_boot$Y*IPTW_boot.trim$weights))
}
quantile(risk_ratio_b, prob=c(0.025,0.975))