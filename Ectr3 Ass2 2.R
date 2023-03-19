rm(list=ls())
data <- read.csv2("~/Documents/GitHub/Econometrics3-2/data_assignment2.csv") # Round the big values to integers

library(vars) # Only used to find lag order
library(tsDyn) # Only used to check 

library(geigen)

# a

date = data$X
OR = data$OR

# Transform variables
m = log(data$M)
p = log(data$P)
y = log(data$Y)
whh = log(data$WHH)

rms = m-p
dwhh = whh[2:92]-whh[1:91]

# Plot timeseries
plot(ts(rms, start = c(1990,1), frequency = 4), ylab="(m-p)")
plot(ts(y, start = c(1990,1), frequency = 4), ylab="y")
plot(ts(dwhh, start = c(1990,4), frequency = 4), ylab="dwhh")
plot(ts(OR, start = c(1990,1), frequency = 4), ylab="OR")



# b

# Generate y matrix
dikke_y = matrix(nrow= 92,ncol = 3)
dikke_y[,1] = rms 
dikke_y[,2] = y
dikke_y[,3] = OR
colnames(dikke_y) = list("rms", "y", "OR")

# # Find best p
VAR_1 = VAR(dikke_y,p=1,type='none',season=NULL,ic="AIC")
VAR_2 = VAR(dikke_y,p=2,type='none',season=NULL,ic="AIC")
VAR_3 = VAR(dikke_y,p=3,type='none',season=NULL,ic="AIC")
AIC(VAR_1)
AIC(VAR_2)
AIC(VAR_3)

# Lowest AIC for p=2 in VAR, so p=1 in VECM

X = dikke_y




# Kopie van tut

## Build matrices
build_dX <- function(dataset, p){
  return(dataset[(p+2):nrow(dataset),] - dataset[(p+1):(nrow(dataset)-1),])
}
build_Xlag <- function(dataset, p){
  return(dataset[(p+1):(nrow(dataset)-p),])
}
build_W <- function(dataset, p){
  dX_full <- t(dataset[2:(nrow(dataset)-1),] - dataset[1:(nrow(dataset)-2),])
  W <- matrix(c(dX_full[, p:1]), 1)
  for (t in (p+1):ncol(dX_full)){
    W <- rbind(W, matrix(c(dX_full[,(t-p+1):t]), 1))
  }
  return(W)
}
build_Mw <- function(W){
  Mw <- (W %*% solve(t(W) %*% W) %*% t(W))
  return(diag(1, nrow(Mw)) - Mw)
}
build_S00 <- function(Mw, dX){
  return((t(dX) %*% Mw %*% dX) / nrow(dX))
}
build_S11 <- function(Mw, Xlag){
  return((t(Xlag) %*% Mw %*% Xlag) / nrow(Xlag))
}
build_S01 <- function(Mw, Xlag, dX){
  return((t(dX) %*% Mw %*% Xlag) / nrow(Xlag))
}

# Define estimators as functions of beta
est_alpha <- function(beta, S01, S11){
  return(S01 %*% beta %*% solve(t(beta) %*% S11 %*% beta))
}
est_phi <- function(beta, dX, Xlag, W, alpha){
  return(t(dX - (Xlag %*% beta %*% t(alpha))) %*% W %*% solve(t(W) %*% W))
}
est_omega <- function(beta, S00, S01, S11){
  S00 - (S01 %*% beta %*% solve(t(beta) %*% S11 %*% beta) %*% t(beta) %*% t(S01))
}

### Finally, estimate beta by the eigenvectors
est_beta <- function(S11, S00, S01, r){
  A <- (t(S01) %*% solve(S00)) %*% (S01)
  B <- S11
  eigenvalues <- geigen(A, B, symmetric=FALSE)
  beta <- eigenvalues$vectors[,order(-eigenvalues$values)]
  beta <- matrix(beta[, 1:r], ncol=r)
  lambdas <- eigenvalues$values[order(-eigenvalues$values)]
  return(list(beta, lambdas))
}

# Put it in a function
est_VECM <- function(dataset, p, r){
  dX <- build_dX(dataset, p)
  Xlag <- build_Xlag(dataset, p)
  W <- build_W(dataset, p)
  Mw <- build_Mw(W)
  S00 <- build_S00(Mw, dX)
  S11 <- build_S11(Mw, Xlag)
  S01 <- build_S01(Mw, Xlag, dX)
  eigenvalues_results <- est_beta(S11, S00, S01, r=r)
  beta <- eigenvalues_results[[1]]
  lambdas <- eigenvalues_results[[2]]
  alpha <- est_alpha(beta, S01, S11)
  phi <- est_phi(beta, dX, Xlag, W, alpha)
  omega <- est_omega(beta, S00, S01, S11)
  return(list(beta=beta, alpha=alpha, phi=phi, omega=omega,
              lambdas=lambdas, S00=S00, S11=S11))
}

p <- 1
r <- 1

VECM_results <- est_VECM(X, p, r)
vecm_check <- VECM(X, lag=p, r=r, include=c('none'), estim="ML")#, LRinclude="none")
summary(vecm_check)

VECM_results$phi # Our parameter estimates match

VECM_results$beta/VECM_results$beta[1]

# Rank test
iT = 92

LRstat <- function(lambdas, n=iT-2){
  -n*sum(log(1 - lambdas))
}
lambdas <- est_VECM(X, p=1, r=3)$lambdas
LRstat(lambdas)
LRstat(lambdas[1])
LRstat(lambdas[2])
LRstat(lambdas[3])

# So, p=1 and r=1   HIER HEBBEN WE DUS GEEN BEWIJS VOOR MOETEN WE NOG FF UITZOEKEN

# Check 
summary(rank.test(VECM(X, lag=p, r=2, include='none', estim="ML",
                       LRinclude="none"), cval=0.01))


# c

data2 <- data[data$X<"2007-3",]


date2 = data2$X
OR2 = data2$OR
m2 = log(data2$M)
p2 = log(data2$P)
y2 = log(data2$Y)
whh2 = log(data2$WHH)

rms2 = m2-p2
dwhh2 = whh2[2:70]-whh2[1:69]

dikke_y2 = matrix(nrow= 69,ncol = 4)
dikke_y2[,1] = rms2[2:70] 
dikke_y2[,2] = y2[2:70] 
dikke_y2[,3] = OR2[2:70] 
dikke_y2[,4] = dwhh2
colnames(dikke_y2) = list("rms", "y", "OR", "dwhh")

build_dX <- function(dataset, p){
  return(dataset[(p+2):nrow(dataset),] - dataset[(p+1):(nrow(dataset)-1),])
}


est_VECM2 <- function(dataset, p, r){
  constant <- rep(1,67)
  trend <- seq(1, 67)
  dX <- cbind(build_dX(dataset, p), constant)
  Xlag <- cbind(build_Xlag(dataset, p), trend)
  W <- cbind(build_W(dataset, p), constant)
  Mw <- build_Mw(W)
  S00 <- build_S00(Mw, dX)
  S11 <- build_S11(Mw, Xlag)
  S01 <- build_S01(Mw, Xlag, dX)
  eigenvalues_results <- est_beta(S11, S00, S01, r=r)
  beta <- eigenvalues_results[[1]]
  lambdas <- eigenvalues_results[[2]]
  alpha <- est_alpha(beta, S01, S11)
  phi <- est_phi(beta, dX, Xlag, W, alpha)
  omega <- est_omega(beta, S00, S01, S11)
  return(list(beta=beta, alpha=alpha, phi=phi, omega=omega,
              lambdas=lambdas, S00=S00, S11=S11))
}


summary(VECM(dikke_y2, lag=p, r=r, include=c('const'), estim=c("ML"), LRinclude="trend")) # Add our test
VECM_results2<-est_VECM2(dikke_y2, p, r) # Just storing the results
VECM_results2$beta/VECM_results2$beta[1] # Relation
