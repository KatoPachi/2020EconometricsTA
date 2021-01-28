library(AER) # Including lmtest and sandwich
library(nlme) # Do gls

# Preparations and preliminary tests

dt = read.csv(
  file = "C:/Users/Administrator/Desktop/D1/D1 AW/Econometrics 2/TA/TA session/TA session #4/boston.csv",
  header = TRUE, sep = ",", row.names = NULL, stringsAsFactors = FALSE)
dt = dt[complete.cases(dt),] # complete.cases returns TURE if one row doesn't contain NA
dt = dt[,c("value", "crime", "industrial", "distance")]

n = nrow(dt)
head(dt)

model = value ~ crime + industrial + distance
ols = lm(model, data = dt)
ols_summary = summary(ols)

# Checking heteroskedasticity and autocorrelation

# Breusch=Pagan=Godfrey test against heteroskedasticity.

(bpgtest = bptest(ols, data = dt, studentize = F)) # from lmtest
# If studentize is set to TRUE Koenker's studentized version of the test statistic will be used.

# Durbin-Watson test for autocorrelation of disturbances.

(dwtest = dwtest(ols, data = dt, alternative = "two.sided")) #alternative = c("greater", "two.sided", "less")


# OLS method with HAC covariance matrix estirmator

cov_hac = NeweyWest(ols)
se_hac = sqrt(diag(cov_hac))

t_hac = coef(ols)/se_hac
p_hac = pt(abs(t_hac), df = nrow(dt) - 4, lower.tail = FALSE)*2


print("NeweyWest covariance matrix estimate:"); cov_hac

print("NeweyWest se estimates:"); se_hac

print("T statistics by NeweyWest covariance:"); t_hac

print("P values:"); p_hac

# Comparing with lmtest::coeftest 
print("Test using NeweyWest estimator:"); coeftest(ols, vcov. = NeweyWest)


# FGLS method by numerical calculation

X = as.matrix(cbind(rep(1, n), dt$crime, dt$industrial, dt$distance))
label = c("(Intercept)", "crime", "industrial", "distance")
colnames(X) = label

Y = as.vector(dt[, "value"])

# Using HCCME(HC_0)

cov_hat = diag(resid(ols)^2)


# Deriving estimates by closed-form estimators 
b_fgls =  solve(t(X) %*% solve(cov_hat) %*% X) %*% (t(X) %*% solve(cov_hat) %*% Y)
cov_fgls = solve(t(X) %*% solve(cov_hat) %*% X)

se_fgls = sqrt(diag(cov_fgls))
t_fgls = b_fgls/se_fgls
p_fgls = pt(abs(t_fgls), df = n - ncol(X), lower.tail = FALSE)*2

# Displaying results 

print("FGLS covariance matrix estimate:"); cov_fgls

print("FGLS se estimates:"); se_fgls

print("T statistics by FGLS covariance:"); t_fgls

print("P values:"); p_fgls


# Iterative numerical method

termination = function(b, omega){
  grad = -2 * t(X) %*% solve(omega) %*% (Y - X %*% b)
  Inner = t(grad) %*% grad
  Inner_sqrt = sqrt(Inner)
  return(Inner_sqrt)
}

b_loop = b_fgls # step1
cov_u = cov_hat
n_loop = 0

while (termination(b_loop, cov_u) > 10^-12) { # judge
  u_hat = as.vector(Y - X %*% b_loop) # step2
  cov_u = diag(u_hat^2)   # step3
  cov_b_loop = solve(t(X) %*% solve(cov_u) %*% X)
  b_loop =  solve(t(X) %*% solve(cov_u) %*% X) %*% (t(X) %*% solve(cov_u) %*% Y) # step4
  n_loop = n_loop + 1
  }


se_loop = sqrt(diag(cov_b_loop))
t_loop = b_loop/se_loop
p_loop = pt(abs(t_loop), df = n - ncol(X), lower.tail = FALSE)*2


# Checking number of iterations and displaying results

print("Number of iterations:"); n_loop

print("FGLS covariance matrix estimate:"); cov_b_loop

print("FGLS se estimates:"); se_loop

print("T statistics by FGLS covariance:"); t_loop

print("P values:"); p_loop



# Estimating by nlme::gls 

# Caculating AR(1) estimate from dwtest
rho = 1 - 0.5 * as.numeric(dwtest$statistic)

# Estimating
gls_esti = gls(model, data = dt, correlation = corAR1(rho), weights = varPower())


(gls_summary = summary(gls_esti))



# Making a table to show each model
library(stargazer)

R2_fgls = 1 - as.numeric(crossprod(Y - X %*% b_fgls)) / sum((Y - mean(Y))^2)
R2_loop = 1 - as.numeric(crossprod(Y - X %*% b_loop)) / sum((Y - mean(Y))^2)
R2_byR  = 1 - as.numeric(crossprod(Y - X %*% coef(gls_esti))) / sum((Y - mean(Y))^2)

se_byR  = coef(gls_summary)[, "Std.Error"]
t_byR   = coef(gls_summary)[, "t-value"]
p_byR   = coef(gls_summary)[, "p-value"]

stargazer(
  ols, gls_esti, gls_esti, gls_esti, 
  coef = list(coef(ols), b_fgls, b_loop, coef(gls_esti)), se = list(se_hac, se_fgls, se_loop, se_byR),
  t = list(t_hac, t_fgls, t_loop, t_byR), p = list(p_hac, p_fgls, p_loop, p_byR),
  t.auto = FALSE, p.auto = FALSE,
  report = "vcstp", keep.stat = c("n"),
  add.lines = list(
    c("Type", "HA-Roubusted OLS", "fgls", "fgls_loop", "Built-in R gls"),
    c("R-Squared", round(ols_summary$r.squared, 4), round(R2_fgls, 4), round(R2_loop, 4), round(R2_byR, 4))),
  title = "Results of linear model estimations",
  label = "LS",
  type = "latex", header = FALSE, font.size = "small",
  table.placement = "htb", omit.table.layout = "n"
)