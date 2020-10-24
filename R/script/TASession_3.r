## ----data---------------------------------------------------------------------
dt <- read.csv(
  file = "./data/titanic.csv", 
  header = TRUE,  sep = ",", row.names = NULL,  stringsAsFactors = FALSE)

dt$female <- ifelse(dt$sex == "female", 1, 0)
dt <- subset(dt, !is.na(survived)&!is.na(age)&!is.na(fare)&!is.na(female))

dt <- dt[,c("survived", "age", "fare", "female")]
head(dt)


## ----ols----------------------------------------------------------------------
model <- survived ~ factor(female) + age + fare
LPM <- lm(model, data = dt)


## ----RobustSE-----------------------------------------------------------------
# heteroskedasticity-robust standard errors
dt$"(Intercept)" <- 1
X <- as.matrix(dt[,c("(Intercept)", "female", "age", "fare")])
u <- diag(LPM$residuals^2)

XX <- t(X) %*% X
avgXX <- XX * nrow(X)^{-1}
inv_avgXX <- solve(avgXX)

uXX <- t(X) %*% u %*% X
avguXX <- uXX * nrow(X)^{-1} 

vcov_b <- (inv_avgXX %*% avguXX %*% inv_avgXX) * nrow(X)^{-1}
rse_b <- sqrt(diag(vcov_b))

label <- c("(Intercept)", "factor(female)1", "age", "fare")
names(rse_b) <- label

# homoskedasticity-based standard errors
se_b <- sqrt(diag(vcov(LPM)))

print("The Variance of OLS"); vcov(LPM)
print("The Robust variance of OLS"); vcov_b
print("The Robust se using White method"); rse_b
print("The Robust t-value using White method"); coef(LPM)/rse_b


## ----lmtest-------------------------------------------------------------------
library(lmtest) #use function `coeftest`
library(sandwich) #use function `vcovHC`
coeftest(LPM, vcov = vcovHC(LPM, type = "HC0"))[, "Std. Error"]
coeftest(LPM, vcov = vcovHC(LPM, type = "HC0"))[, "t value"]


## ----LPM_result, results = "asis"---------------------------------------------
# t-stats
t_b <- coef(LPM)/se_b 
rt_b <- coef(LPM)/rse_b
# p-value Pr( > |t|)
p_b <- pt(abs(t_b), df = nrow(X)-ncol(X), lower = FALSE)*2
rp_b <- pt(abs(rt_b), df = nrow(X)-ncol(X), lower = FALSE)*2

library(stargazer)
stargazer(
  LPM, LPM,
  se = list(se_b, rse_b), t = list(t_b, rt_b), p = list(p_b, rp_b),
  t.auto = FALSE, p.auto = FALSE,
  report = "vcstp", keep.stat = c("n"),
  covariate.labels = c("Female = 1"),
  add.lines = list(
    c("Standard errors", "Homoskedasticity-based", "Heteroskedasticity-robust")),
  title = "Results of Linear Probability Model", label = "LPM",
  type = "latex", header = FALSE, font.size = "small",
  omit.table.layout = "n", table.placement = "h"
)


## ----probit-------------------------------------------------------------------
Y <- dt$survived
female <- dt$female; age <- dt$age; fare <- dt$fare

# log-likelihood
LnLik <- function(b, model = c("probit", "logit")) {

  xb <- b[1]+ b[2]*female + b[3]*age + b[4]*fare

  if (model == "probit") {
    L <- pnorm(xb)
  } else {
    L <- 1/(1 + exp(-xb))
  }

  LL_i <- Y * log(L) + (1 - Y) * log(1 - L)
  LL <- -sum(LL_i)

  return(LL)
}

#Newton-Raphson
init <- c(0.218, 0.512, -0.002, 0.001)
probit <- nlm(LnLik, init, model = "probit", hessian = TRUE)

label <- c("(Intercept)", "factor(female)1", "age", "fare")
names(probit$estimate) <- label
colnames(probit$hessian) <- label; rownames(probit$hessian) <- label

b_probit <- probit$estimate
vcov_probit <- solve(probit$hessian); se_probit <- sqrt(diag(vcov_probit))
LL_probit <- -probit$minimum

#glm function
model <- survived ~ factor(female) + age + fare
probit_glm <- glm(model, data = dt, family = binomial("probit"))

#result
print("The MLE of probit model using nlm"); b_probit
print("The Variance of probit model using nlm"); vcov_probit
print("The se of probit model using nlm"); se_probit
print("The coefficients of probit using glm"); coef(probit_glm)
print("The se of probit using glm"); sqrt(diag(vcov(probit_glm)))


## ----logit--------------------------------------------------------------------
#Newton-Raphson
logit <- nlm(LnLik, init, model = "logit", hessian = TRUE)

label <- c("(Intercept)", "factor(female)1", "age", "fare")
names(logit$estimate) <- label
colnames(logit$hessian) <- label; rownames(logit$hessian) <- label

b_logit <- logit$estimate
vcov_logit <- solve(logit$hessian); se_logit <- sqrt(diag(vcov_logit))
LL_logit <- -logit$minimum

#glm function
logit_glm <- glm(model, data = dt, family = binomial("logit"))

#result
print("The MLE of logit model"); b_logit
print("The Variance of logit model"); vcov_logit
print("The se of logit model"); se_logit
print("The coefficients of logit using glm"); coef(logit_glm)
print("The se of logit using glm"); sqrt(diag(vcov(logit_glm)))


## ----summary_probit_logit, results = "asis"-----------------------------------
# z-value
z_probit <- b_probit/se_probit
z_logit <- b_logit/se_logit

# Pr(>|z|)
p_probit <- pnorm(abs(z_probit), lower = FALSE)*2
p_logit <- pnorm(abs(z_logit), lower = FALSE)*2

stargazer(
  probit_glm, logit_glm,
  coef = list(b_probit, b_logit), se = list(se_probit, se_logit),
  t = list(z_probit, z_logit), p = list(p_probit, p_logit),
  t.auto = FALSE, p.auto = FALSE,
  report = "vcstp", keep.stat = c("n"),
  covariate.labels = c("Female = 1"),
  add.lines = list(
    c("Log-Likelihood", round(LL_probit, 3), round(LL_logit, 3))),
  title = "Results of Probit and Logit model",
  label = "probit_logit",
  type = "latex", header = FALSE, font.size = "small",
  table.placement = "h", omit.table.layout = "n"
)


## ----AME----------------------------------------------------------------------
DeltaAME <- function(b, X, vcov, jbin = NULL, model = c("probit", "logit")) {
  Xb <- numeric(nrow(X))
  for (i in 1:length(b)) {
     Xb <- Xb + b[i] * X[,i]
  }

  if (model == "probit") {
    dens <- dnorm(Xb)
    grad <- -Xb * dens
  } else {
    dens <- exp(-Xb)/(1 + exp(-Xb))^2
    grad <- dens * (-1+2*exp(-Xb)/(1+exp(-Xb)))
  }

  ame <- mean(dens) * b
  if (!is.null(jbin)) {
    for (i in jbin) {
      val1 <- X[,-i] %*% matrix(b[-i], ncol = 1) + b[i]
      val0 <- X[,-i] %*% matrix(b[-i], ncol = 1)
      if (model == "probit") {
        amed <- mean(pnorm(val1) - pnorm(val0))
      } else { 
        amed <- mean((1/(1 + exp(-val1))) - (1/(1 + exp(-val0))))
      }
      ame[i] <- amed
    }
  }

  e <- NULL
  for (i in 1:length(b)) {
    e <- c(e, rep(mean(X[,i] * grad), length(b)))
  }

  Jacob <- matrix(e, nrow = length(b), ncol = length(b))

  for (i in 1:nrow(Jacob)) {
    Jacob[i,] <- b[i] * Jacob[i,]
  }
  diag(Jacob) <- diag(Jacob) + rep(mean(dens), length(b))

  if (!is.null(jbin)) {
    for (i in jbin) {
      val1 <- X[,-i] %*% matrix(b[-i], ncol = 1) + b[i]
      val0 <- X[,-i] %*% matrix(b[-i], ncol = 1)
      de <- NULL
      if (model == "probit") {
        for (j in 1:length(b)) {
           if (j != i) {
            dep <- X[,j] * (dnorm(val1) - dnorm(val0))
            de <- c(de, mean(dep))
           } else {
            dep <- dnorm(val1)
            de <- c(de, mean(dep))
           }
        }        
      } else {
        for (j in 1:length(b)) {
           if (j != i) {
            dep <- X[,j] * 
              ((exp(-val1)/(1 + exp(-val1))^2) - (exp(-val0)/(1 + exp(-val0))^2))
            de <- c(de, mean(dep))
           } else {
            dep <- exp(-val1)/(1 + exp(-val1))^2
            de <- c(de, mean(dep))
           }
        } 
      }
      Jacob[i,] <- de
    }
  }

  label <- names(b)
  colnames(Jacob) <- label; rownames(Jacob) <- label

  vcov_ame <- Jacob %*% vcov %*% t(Jacob)
  se_ame <- sqrt(diag(vcov_ame))
  z_ame <- ame/se_ame
  p_ame <- pnorm(abs(z_ame), lower = FALSE)*2

  return(list(AME = ame[-1], SE = se_ame[-1], zval = z_ame[-1], pval = p_ame[-1]))
}

X <- as.matrix(dt[,c("(Intercept)", "female", "age", "fare")])
ame_probit <- DeltaAME(b_probit, X, vcov_probit, jbin = 2, model = "probit")
ame_logit <- DeltaAME(b_logit, X, vcov_logit, jbin = 2, model = "logit")

print("AME of probit estimates"); ame_probit$AME
print("AME of logit estimates"); ame_logit$AME 
print("SE of AME of probit estimates"); ame_probit$SE 
print("SE of AME of logit estimates"); ame_logit$SE


## ----Margin-------------------------------------------------------------------
library(margins)
summary(margins(probit_glm))
summary(margins(logit_glm))


## ----pcp----------------------------------------------------------------------
Y <- dt$survived
X <- as.matrix(dt[,c("(Intercept)", "female", "age", "fare")])

Xb_lpm <- X %*% matrix(coef(LPM), ncol = 1)
Xb_probit <- X %*% matrix(b_probit, ncol = 1)
Xb_logit <- X %*% matrix(b_logit, ncol = 1)

hatY_lpm <- ifelse(Xb_lpm > 0.5, 1, 0)
hatY_probit <- ifelse(pnorm(Xb_probit) > 0.5, 1, 0)
hatY_logit <- ifelse(1/(1 + exp(-Xb_logit)) > 0.5, 1, 0)

pcp_lpm <- round(sum(Y == hatY_lpm)/nrow(X), 4)
pcp_probit <- round(sum(Y == hatY_probit)/nrow(X), 4)
pcp_logit <- round(sum(Y == hatY_logit)/nrow(X), 4)


## ----pr2----------------------------------------------------------------------
Y2 <- Y^2

hatu_lpm <- (Y - Xb_lpm)^2
hatu_probit <- (Y - pnorm(Xb_probit))^2
hatu_logit <- (Y - 1/(1 + exp(-Xb_logit)))^2

pr2_lpm <- round(1 - sum(hatu_lpm)/sum(Y2), 4)
pr2_probit <- round(1 - sum(hatu_probit)/sum(Y2), 4)
pr2_logit <- round(1 - sum(hatu_logit)/sum(Y2), 4)


## ----BinaryModelResult, results = "asis"--------------------------------------
stargazer(
  LPM, probit_glm, logit_glm,
  coef = list(coef(LPM), ame_probit$AME, ame_logit$AME),
  se = list(rse_b, ame_probit$SE, ame_logit$SE),
  t = list(rt_b, ame_probit$zval, ame_logit$zval),
  p = list(rp_b, ame_probit$pval, ame_logit$pval),
  t.auto = FALSE, p.auto = FALSE,
  omit = c("Constant"), covariate.labels = c("Female = 1"),
  report = "vcstp", keep.stat = c("n"),
  add.lines = list(
    c("Percent correctly predicted", pcp_lpm, pcp_probit, pcp_logit),
    c("Pseudo R-squared", pr2_lpm, pr2_probit, pr2_logit)
  ),
  omit.table.layout = "n", table.placement = "t",
  title = "Titanic Survivors: LPM, Probit (AME), and Logit (AME)",
  label = "titanic",
  type = "latex", header = FALSE
)


## ----data2--------------------------------------------------------------------
house <- read.csv(file = "./data/housing.csv", header = TRUE,  sep = ",")
house <- house[,c("Level", "lnPrice", "Top25")]
house$Levelnum <- ifelse(
  house$Level == 1, 25, 
  ifelse(house$Level == 2, 50, 
  ifelse(house$Level == 3, 75, 100)))
head(house)


## ----orderModel---------------------------------------------------------------
library(MASS)
library(tidyverse) #use case_when()

ols <- lm(Levelnum ~ lnPrice + Top25, data = house)

model <- factor(Level) ~ lnPrice + Top25
oprobit <- polr(model, data = house, method = "probit")
ologit <- polr(model, data = house, method = "logistic")

a_oprobit <- round(oprobit$zeta, 3)
a_ologit <- round(ologit$zeta, 3)

xb_oprobit <- oprobit$lp 
xb_ologit <- ologit$lp

hatY_oprobit <- case_when(
  xb_oprobit <= oprobit$zeta[1] ~ 1,
  xb_oprobit <= oprobit$zeta[2] ~ 2,
  xb_oprobit <= oprobit$zeta[3] ~ 3,
  TRUE ~ 4
)
hatY_ologit <- case_when(
  xb_ologit <= ologit$zeta[1] ~ 1,
  xb_ologit <= ologit$zeta[2] ~ 2,
  xb_ologit <= ologit$zeta[3] ~ 3,
  TRUE ~ 4
)

pred_oprobit <- round(sum(house$Level == hatY_oprobit)/nrow(house), 3)
pred_ologit <- round(sum(house$Level == hatY_ologit)/nrow(house), 3)


## ----predict------------------------------------------------------------------
quantef <- function(model) {
  b <- coef(model)
  val1 <- mean(house$lnPrice)*b[1] + b[2]
  val0 <- mean(house$lnPrice)*b[1]

  prob <- matrix(c(rep(val1, 3), rep(val0, 3)), ncol = 2, nrow = 3)
  for (i in 1:3) {
    for (j in 1:2) {
      prob[i,j] <- pnorm(model$zeta[i] - prob[i,j])
    }
  }
  Ey1 <- 25*prob[1,1] + 50*(prob[2,1]-prob[1,1]) + 
    75*(prob[3,1]-prob[2,1]) + 100*(1-prob[3,1])
  Ey0 <- 25*prob[1,2] + 50*(prob[2,2]-prob[1,2]) + 
    75*(prob[3,2]-prob[2,2]) + 100*(1-prob[3,2])
  
  return(Ey1 - Ey0)
} 

ef_oprobit <- round(quantef(oprobit), 3)
ef_ologit <- round(quantef(ologit), 3)


## ----tab_orderModel, results = "asis"-----------------------------------------
stargazer(
  ols, oprobit, ologit,
  report = "vcstp", keep.stat = c("n"),
  omit = c("Constant"),
  add.lines = list(
    c("Cutoff value at 1|2", "", a_oprobit[1], a_ologit[1]),
    c("Cutoff value at 2|3", "", a_oprobit[2], a_ologit[2]),
    c("Cutoff value at 3|4", "", a_oprobit[3], a_ologit[3]),
    c("Quantitative Effect of Top25", "", ef_oprobit, ef_ologit),
    c("Percent correctly predicted", "", pred_oprobit, pred_ologit)
  ),
  omit.table.layout = "n", table.placement = "t",
  title = "Floor Level of House: Ordered Probit and Logit Model",
  label = "housing",
  type = "latex", header = FALSE
)


## ----supermarket--------------------------------------------------------------
library(AER)
data(BankWages)
dt <- BankWages
dt$job <- as.character(dt$job)
dt$job <- factor(dt$job, levels = c("admin", "custodial", "manage"))
head(BankWages, 5)


## ----multinomial, results = "hide"--------------------------------------------
library(nnet)
est_mlogit <- multinom(job ~ education + gender, data = dt)

# observations and percent correctly predicted
pred <- est_mlogit$fitted.value
pred <- colnames(pred)[apply(pred, 1, which.max)]
n <- length(pred)
pcp <- round(sum(pred == dt$job)/n, 3)

# Log-likelihood and pseudo R-sq
loglik1 <- as.numeric(nnet:::logLik.multinom(est_mlogit))
est_mlogit0 <- multinom(job ~ 1, data = dt)
loglik0 <- as.numeric(nnet:::logLik.multinom(est_mlogit0))
pr2 <- round(1 - loglik1/loglik0, 3)


## ----tab_multinomial, results = "asis"----------------------------------------
stargazer(
  est_mlogit,
  covariate.labels = c("Education", "Female = 1"),
  report = "vcstp", omit.stat = c("aic"),
  add.lines = list(
    c("Observations", n, ""),
    c("Percent correctly predicted", pcp, ""),
    c("Log-likelihood", round(loglik1, 3), ""),
    c("Pseudo R-sq", pr2, "")
  ),
  omit.table.layout = "n", table.placement = "t",
  title = "Multinomial Logit Model of Job Position",
  label = "job",
  type = "latex", header = FALSE  
)
