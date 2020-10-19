## ----setup, include = FALSE, echo = FALSE-------------------------------------
# This is options to make pdf file. You should ignore
knitr::opts_chunk$set(
  message = FALSE, 
  warning = FALSE, 
  echo = TRUE, 
  cache = FALSE,
  fig.pos = "h")


## ----addition_subtraction, eval = FALSE, echo = FALSE-------------------------
## A <- matrix(c(1, 2, 3, 4), nrow = 2); B <- matrix(c(5, 6, 7, 8), nrow = 2)
## print("The matrix A is"); A
## print("The matrix B is"); B
## print("A+B is"); A + B
## print("A-B is"); A - B


## ----multiplication, eval = FALSE, echo = FALSE-------------------------------
## A <- matrix(c(1, 2, 3, 4), nrow = 2); B <- matrix(c(5, 6, 7, 8), nrow = 2)
## print("AB is"); A %*% B # A*B returns Hadamard product (element-wise product)
## print("BA is"); B %*% A # Non-commutativity
## print("The Kronecker prodct of two matricies is"); A %x% B


## ----transposed_matrix, eval = FALSE, echo = FALSE----------------------------
## A <- matrix(c(1, 2, 3, 4), nrow = 2); B <- matrix(c(5, 6, 7, 8), nrow = 2)
## print("The transposed matrix of A is"); t(A)
## # first property
## print("(AB)' is"); t(A %*% B)
## print("B'A' is"); t(B) %*% t(A)
## # third property
## print("Using the Kronecker product, (AB)' is"); t(A %x% B)
## print("Using the Kronecker product, A'B' is"); t(A) %x% t(B)


## ----diagonal_matrix, eval = FALSE, echo = FALSE------------------------------
## print("The diagonal matrix is"); diag(c(2, 4, 6))
## print("The 3-by-3 identity matrix is"); diag(3)
## 
## A <- matrix(c(1, 2, 3, 4), nrow = 2)
## print("The diagonal element of A is"); diag(A)
## sprintf("The trace of A is %1d", sum(diag(A)))


## ----inverse_matrix, eval = FALSE, echo = FALSE-------------------------------
## A <- matrix(c(1, 2, 3, 4), nrow = 2); B <- matrix(c(5, 6, 7, 8), nrow = 2)
## print("The inverse matrix of A is"); solve(A)
## # First property
## print("The transposed matrix of the inverse matrix of A is"); t(solve(A))
## print("The inverse matrix of transposed matrix of A is"); solve(t(A))
## # Third property
## print("The example of third property is"); (A %x% B) %*% (solve(A) %x% solve(B))


## ----matrix 3, eval = FALSE, echo = FALSE-------------------------------------
## A <- matrix(1, nrow = 2, ncol = 2)
## print("The matrix A is"); A
## print("The characteristic roots are"); eigen(A)$values


## ----WLLN, fig.cap = "Simulation Result of WLLN", fig.width = 6, fig.height = 4----
set.seed(120504)

data <- data.frame(
  trial = 1:20000,
  success = rbinom(n = 20000, size = 1, prob = .5)
)
data$sum_success <- cumsum(data$success)
data$prob <- data$sum_success/data$trial

plot(
  log(data$trial), data$prob, type = "l", col = "blue", 
  ylim = c(0.3, 0.7), xlab = "logged trials", ylab = "Pr(head)")
lines(c(0, 10), c(0.5, 0.5), lwd = 1, col = "red")


## ----CLT, fig.cap = "Simulation Result of CLT", fig.width = 6, fig.height = 4----
set.seed(120504)
m <- 10000; n <- c(3, 100, 1000); p <- 0.5
a <- seq(-4, 4, .01); b <- dnorm(a)

dt <- list("n = 3"=numeric(m), "n = 100"=numeric(m), "n = 1000"=numeric(m))
for (i in 1:3) {
  dt[[i]] <- rbinom(n = m, size = n[i], prob = p)
  dt[[i]] <- sqrt(n[i])*(dt[[i]]/n[i] - p)/sqrt(p*(1-p))
}

par(mfrow=c(2,2), mai = c(0.5, 0.5, 0.35, 0.35)) 
for (i in 1:3) {
  hist(dt[[i]], col = "grey", freq = FALSE, 
    xlab = "", main = names(dt)[i], xlim = c(-4, 4))
  par(new = TRUE)
  plot(a, b, type = "l", col = "red", axes = FALSE, 
    xlab = "", ylab = "", main = "")
}

