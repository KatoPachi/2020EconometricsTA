library(AER) # Including lmtest, sandwich and ivreg(2SLS)
library(gmm) # Including gmm function
library(MASS) # Using ginv: generalized inverse operation
library(stargazer)

# Preparations 
setwd("C:/Users/Administrator/Desktop/D1/D1 AW/Econometrics 2/TA/TA session/TA session #10")

dt = read.csv(
  file = "boston.csv",
  header = TRUE, sep = ",", row.names = NULL, stringsAsFactors = FALSE)
dt = dt[complete.cases(dt),] # complete.cases returns TURE if one row doesn't contain NA

n = nrow(dt)
head(dt[, c("value", "crime", "industrial", "distance", "black",  "ptratio")])

# Extract varibales and instruments
name_x = c("constant", "crime", "industrial", "distance")
name_z = c("constant", "black",  "ptratio", "industrial", "distance") # means L = d + 1 = 5
constant = rep(1, n)

X = cbind(constant, dt[, c("crime", "industrial", "distance")])
colnames(X) = name_x

Z = cbind(constant, dt[, c("black", "ptratio", "industrial", "distance")])
colnames(Z) = name_z

Y = dt[, "value"]




# Numerical 2SGMM Method to A Simple Linear Over-identified case

# Since we tested out the heteroskedasticity, 
# this dataset provides an appropriate example for studying 2SGMM.


# Firstly, recall HA-robust ols estimates for later comparison
model = value ~ crime + industrial + distance

ols = lm(model, data = dt)
#summary(ols)
coeftest(ols, vcov. = NeweyWest)

cov_hac = NeweyWest(ols)
se_hac = sqrt(diag(cov_hac))

t_hac = coef(ols)/se_hac
p_hac = pt(abs(t_hac), df = nrow(dt) - ncol(X), lower.tail = FALSE)*2


# Secondly, perform 2SGMM and 
# the 1st step is like follows (looks like 2SLS).

Z = as.matrix(Z)
X = as.matrix(X)

ZX = t(Z) %*% X/n # or ZX = crossprod(Z, X)/n
ZY = t(Z) %*% Y/n

var_g_1 = t(Z) %*% Z/n

theta_1 = solve(t(ZX) %*% solve(var_g_1) %*% ZX) %*% 
            t(ZX) %*% solve(var_g_1) %*% ZY

cov_1   = solve(t(ZX) %*% solve(var_g_1) %*% ZX)


# Caculate the estimates and check
se_1  = sqrt(diag(cov_1))
t_1   = theta_1/se_1
p_1   = pt(abs(t_1), df = nrow(dt) - ncol(X), lower.tail = FALSE)*2

print("1st step estimates:"); theta_1

print("1st step se estimates:"); se_1

print("T statistics :"); t_1

print("P values:"); p_1


# Do 2nd step 2SGMM

# Theoretically, we should calculate the estimate of the variance with theta_1 by

Zu_1    = t(Z) %*% (Y - X %*% theta_1)
var_g_2 = Zu_1 %*% t(Zu_1)/n
# or var_g_1 = tcrossprod(Zu_1, Zu_1)/n

# The previous var_g_2 is derived under the iid assumption. 
# It is equal to set vcov = "iid" in gmm function, which would slow down the calculation.

# In practice, let's build a HAC covariance estimate from 2sls residuals and X.
# For simplicity, set L = 10(<= n). 
u_1 = Y - X %*% theta_1
i = 1
l = 1
L = 10
weight = rep(0, L)

Q = matrix(0, nrow = ncol(X), ncol = ncol(X))
Q_1 = matrix(0, nrow = ncol(X), ncol = ncol(X))
Q_2 = matrix(0, nrow = ncol(X), ncol = ncol(X))

for (l in 1:L) {
  weight[l] = 1 - l/(L + 1)
}

for (i in 1:n) {
  Q_1 = Q_1 + u_1[i]^2 * (X[i,] %*% t(X[i,])) # X[i,] here returns column vector
}
Q_1 = Q_1/n

for (l in 1:L) {
  for (i in (l+1):n) {
    Q_2 = Q_2 + weight[l]*u_1[i]*u_1[i-l]*
      (X[i,] %*% t(X[i-l,]) + X[i-l,] %*% t(X[i,]))
  }
}
Q_2 = Q_2/n
Q = Q_1 + Q_2


XXT = X %*% t(X)
omega = ginv(XXT) %*% (X %*% Q %*% t(X)) %*% ginv(XXT)
# Use ginv here to suppress error warnings from traditional inverse.
# omega = ginv(t(X)) %*% Q %*% ginv(X) returns almostly same values.
var_g_2 = t(Z) %*% omega %*% Z /n

theta_2 = ginv(t(ZX) %*% ginv(var_g_2) %*% ZX) %*% 
          t(ZX) %*% ginv(var_g_2) %*% ZY

cov_2   = ginv(t(ZX) %*% ginv(var_g_2) %*% ZX) 


# Calculate the estimates and check
se_2  = sqrt(diag(cov_2))
se_2  = as.vector(se_2)
t_2   = theta_2/se_2
p_2   = pt(abs(t_2), df = nrow(dt) - ncol(X), lower.tail = FALSE)*2

print("2nd step estimates:"); theta_2

print("2nd step se estimates:"); se_2

print("T statistics :"); t_2

print("P values:"); p_2

# Finally, theta_1 is not optimal
se_1 >= se_2

# Display the results
stargazer(
  ols, ols, ols, 
  coef = list(coef(ols), theta_1, theta_2), se = list(se_hac, se_1, se_2),
  t = list(t_hac, t_1, t_2), p = list(p_hac, p_1, p_2),
  t.auto = FALSE, p.auto = FALSE,
  report = "vcstp", keep.stat = c("n"),
  add.lines = list(
    c("Type", "HA-Roubusted OLS", "1st step GMM", "2nd step GMM")),
  title = "Results of rols and numerical 2 step GMM",
  label = "Numeric",
  type = "latex", header = FALSE, font.size = "small",
  table.placement = "htb", omit.table.layout = "n"
)




# Built-in R Functions


# 2SLS and Hausman Wu Test for Endogeneity (By ivreg).
model_iv = value ~ crime + industrial + distance | black + ptratio + industrial + distance

twoSLS = ivreg(model_iv, data = dt)

# Durbin-Wu-Hausman Test
dwh.test = function(model.iv, model.ols){
  cf_diff = coef(model.iv) - coef(model.ols)
  vc_diff = vcovHC(model.iv, "HC0") - vcovHC(model.ols, "HC0") 
  # NeweyWest() doesn't fit well to model.iv so we use white estimator.
  x2_diff = as.vector(t(cf_diff) %*% solve(vc_diff) %*% cf_diff)
  pvalue = pchisq(q = x2_diff, df = dim(vc_diff)[1], lower.tail = F)
  
  result = list(x2_diff, dim(vc_diff)[1], pvalue)
  names(result) = c("DWH Statistic", "df", "P value")
  return(result)
}

dwh.test(twoSLS, ols) # "http://klein.uk/R/myfunctions.R" from professor Thilo Klein.

# Use "diagnostics = TRUE" calls a Wu-Hausman(2 step) F test(available in Wooldridge(2003)).
summary(twoSLS, vcov. = NeweyWest, df = Inf, diagnostics = TRUE)


# 2SGMM and Sargan's J Test (By gmm)
value = dt[, "value"]
crime = dt[, "crime"]
industrial = dt[, "industrial"]
distance = dt[, "distance"]

instruments = dt[, c("black", "ptratio", "industrial", "distance")]

twoSGMM_iid = gmm(g = model, x = instruments, vcov = "iid")
twoSGMM_HAC = gmm(g = model, x = instruments, vcov = "HAC")

summary(twoSGMM_iid)
summary(twoSGMM_HAC)

# Display
stargazer(
  twoSLS, twoSGMM_iid, twoSGMM_HAC,
  report = "vcstp", keep.stat = c("n"),
  add.lines = list(
    c("Type", "2SLS", "2SGMM(iid)", "2SGMM(HAC)")),
  title = "Results of 2SLS and 2SGMM by R",
  label = "byR",
  type = "latex", header = FALSE, font.size = "small",
  table.placement = "htb", omit.table.layout = "n"
)

setwd("~")
