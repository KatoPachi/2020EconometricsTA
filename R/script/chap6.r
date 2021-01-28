## ----data---------------------------------------------------------------------
dt <- read.csv("./data/wages.csv")
head(dt, 14)


## ----summary------------------------------------------------------------------
dt <- dt[,c("id", "time", "exp", "lwage")]
dt$sqexp <- dt$exp^2
summary(dt)


## ----model--------------------------------------------------------------------
model <- lwage ~ -1 + exp + sqexp


## ----pooledOLSE---------------------------------------------------------------
bols1 <- lm(model, data = dt)

library(plm)
bols2 <- plm(model, data = dt, model = "pooling", index = c("id", "time"))


## ----Inference of Pooled OLSE-------------------------------------------------
# Setup
N <- length(unique(dt$id)); T <- length(unique(dt$time))
X <- model.matrix(bols1); k <- ncol(X)

# Inference
uhat <- bols1$residuals
uhatset <- matrix(0, nrow = nrow(X), ncol = nrow(X))

i_from <- 1; j_from <- 1
for (i in 1:max(dt$id)) {
  x <- as.numeric(rownames(dt))[dt$id == i]
  usq <- uhat[x] %*% t(uhat[x])
  i_to <- i_from + nrow(usq) - 1
  j_to <- j_from + ncol(usq) - 1
  uhatset[i_from:i_to, j_from:j_to] <- usq
  i_from <- i_to + 1; j_from <- j_to + 1
}

Ahat <- t(X) %*% X
Bhat <- t(X) %*% uhatset %*% X
vcovols <- solve(Ahat) %*% Bhat %*% solve(Ahat)
seols <- sqrt(diag(vcovols))

# Easy way
library(lmtest)
library(sandwich)
easy_cluster <- coeftest(
  bols2, vcov = vcovHC(bols2, type = "HC0", cluster = "group"))


## ----FGLS estimator-----------------------------------------------------------
# Setup
X <- model.matrix(model, dt); k <- ncol(X)
y <- dt$lwage
N <- length(unique(dt$id)); T <- length(unique(dt$time))

# Estimator of Omega
uhat <- bols1$residuals

Omega_sum <- matrix(0, ncol = T, nrow = T)
for (i in 1:N) {
  x <- as.numeric(rownames(dt))[dt$id == i]
  Omega_sum <- uhat[x] %*% t(uhat[x]) + Omega_sum
}
Omega <- Omega_sum/N

# FGLS estimator
kroOmega <- diag(N) %x% solve(Omega)
bfgls <- solve(t(X) %*% kroOmega %*% X) %*% (t(X) %*% kroOmega %*% y)

# Easy way!!!
easy_fgls <- pggls(
  model, data = dt, index = c("id", "time"), model = "pooling")


## ----Inference of FGLS estimator----------------------------------------------
ufgls <- y - X %*% bfgls
uhatset <- matrix(0, nrow = nrow(X), ncol = nrow(X))
i_from <- 1; j_from <- 1
for (i in 1:max(dt$id)) {
  x <- as.numeric(rownames(dt))[dt$id == i]
  usq <- uhat[x] %*% t(uhat[x])
  i_to <- i_from + nrow(usq) - 1
  j_to <- j_from + ncol(usq) - 1
  uhatset[i_from:i_to, j_from:j_to] <- usq
  i_from <- i_to + 1; j_from <- j_to + 1
}

Ahat <- t(X) %*% kroOmega %*% X
Bhat <- t(X) %*% kroOmega %*% uhatset %*% kroOmega %*% X
vcovfgls <- solve(Ahat) %*% Bhat %*% solve(Ahat)
sefgls <- sqrt(diag(vcovfgls))


## ----time-demeaning matrix----------------------------------------------------
# extract outcome variables for i = 1
i <- as.numeric(rownames(dt))[dt$id == 1]
y1 <- dt$lwage[i]

# deviation from mean
Ydev1 <- y1 - mean(y1)
print("Deviation from mean across time"); Ydev1

# time demean-matrix
T <- length(y1)
vec1 <- rep(1, T)
Qt <- diag(T) - vec1 %*% solve(t(vec1) %*% vec1) %*% t(vec1)
Ydev2 <- Qt %*% y1
print("Time-demeaning matrix"); Ydev2


## ----FE estimator-------------------------------------------------------------
# Setup
X <- model.matrix(model, dt); k <- ncol(X)
y <- dt$lwage
N <- length(unique(dt$id)); T <- length(unique(dt$time))

# FE estimator
i <- rep(1, T)
Qt <- diag(T) - i %*% solve(t(i) %*% i) %*% t(i)
Ydev <- diag(N) %x% Qt %*% y
Xdev <- diag(N) %x% Qt %*% X
bfe <- solve(t(Xdev) %*% Xdev) %*% t(Xdev) %*% Ydev

# Awesome way !!!
plmfe <- plm(model, data = dt, index = c("id", "time"), model = "within")


## ----Inference of FE estimator------------------------------------------------
uhat <- Ydev - Xdev %*% bfe
sigmahat <- sum(uhat^2)/(N*(T-1)-k)
vcovfe <- sigmahat * solve(t(Xdev) %*% Xdev)
sefe <- sqrt(diag(vcovfe))


## ----FGLS-type RE estimator---------------------------------------------------
# Setup
X <- model.matrix(model, dt)
y <- dt$lwage
k <- ncol(X)
N <- length(unique(dt$id))
T <- length(unique(dt$time))

# estimator of Omega
pols <- lm(model, dt)
vhat <- pols$residuals
sigmav <- sum(vhat^2)/(N*T - k)

sumuc <- matrix(0, nrow = N, ncol = T-1)
for (i in 1:N) {
  for (t in 1:T-1) {
    it <- as.numeric(rownames(dt))[dt$id == i & dt$time == t]
    is <- as.numeric(rownames(dt))[dt$id == i & dt$time > t]
    sumuc[i,t] <- vhat[it] * sum(vhat[is])
  }
}
sigmac <- sum(colSums(sumuc))/((N*T*(T-1))/2-k)
sigmau <- sigmav - sigmac

i <- rep(1, T)
Omega <- sigmau * diag(T) + sigmac * i %*% t(i)
kroOmega <- diag(N) %x% solve(Omega)

# Random effect
bre <- solve(t(X) %*% kroOmega %*% X) %*% t(X) %*% kroOmega %*% y


## ----Inference of RE estimator------------------------------------------------
vcovre <- solve(t(X) %*% kroOmega %*% X)
sere <- sqrt(diag(vcovre))




## ----table, results = "asis", echo = FALSE------------------------------------
library(stargazer)
stargazer(
  bols1, bols1, bols1, bols1, 
  column.labels = c("Pooled OLS", "FGLS", "Fixed Effect", "Random Effect"),
  coef = list(coef(bols1), bfgls, bfe, bre),
  se = list(seols, sefgls, sefe, sere),
  keep.stat = c("n"), report = "vcs",
  omit.table.layout = "n", table.placement = "t",
  header = FALSE,
  type = "latex",
  title = "Effect of Experience on Wages (Standard errors are in parentheses)",
  label = "pdm"
)


## ----HausmanTest--------------------------------------------------------------
delta <- bre - bfe
diffv <- vcovre - vcovfe
H <- t(delta) %*% solve(diffv) %*% delta
qtchi <- qchisq(0.99, nrow(delta))
paste("The test statistics of Hausman test is ", round(H, 3))
paste("The 1% quantile value of chi-sq dist is", round(qtchi, 3))

