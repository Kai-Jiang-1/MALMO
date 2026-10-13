rm(list = ls())

# set working directory 
setwd("~/Desktop/MALMO/simulation")

# package 
library(R.matlab)

# parameters 
n <- 50 # number of samples
q <- 3 # number of outcomes
p <- 100 # number of features
snr <- 1 # signal to noise ratio (SNR = mean(abs(c_ij))/\sigma)
gen_miss <- "true"  
missing_type <- "mnar0.33" # missing type, we have 2 combinations:"1/3 mnar + 2/3 mar", "2/3 mnar + 1/3 mar"
proportion <- 0.2 

# setting
set.seed(123)
q_true_num <- q # number of true outcomes
p_true_num <- 5 # number of true features

# print error if q_true_num is greater than q: 
if (q_true_num > q) {
  stop("Error: q_true_num cannot be greater than q.")
}

q_true <- sort(sample(1:q, q_true_num)) # randomly sample two selected outcomes for outcome true
p_true <- sort(sample(1:p, p_true_num)) # randomly sample two selected features for feature true
q_true
p_true

# C matrix (coefficient matrix)
C <- matrix(0, nrow = p, ncol = q)
# the true c_ij is between -2 to -1, and between 1 to 2
C[p_true, q_true] <- matrix(round(sample(c(1,-1), length(p_true) * length(q_true), replace = T) * runif(length(p_true) * length(q_true), 1, 2), 2),
                            nrow = length(p_true), ncol = length(q_true))
C_idx <- C
C_idx[which(C_idx != 0)] <- 1 # create a binary matrix for C matrix, 1 for true coefficients, 0 for false coefficients

# set seed for genearte X matrix
set.seed(123)

# X matrix
# step 1: create a covariance matrix for W matrix
Wcor <- diag(0, p) # covariance matrix of W
# for selected features, we set the covariance to 0.75 - |i - j|
# for unselected features, we set the covariance to 0.2
for (i in 1:p) {
  for (j in 1:p) {
    if (i != j) {
      Wcor[i, j] <- 0.15
    }
  }
}

for (i in 1:p_true_num) {
  for (j in 1:p_true_num) {
    Wcor[p_true[i], p_true[j]] <- 0.75 - 0.0015*abs(i - j)
  }
}

Wcor <- (Wcor + t(Wcor)) / 2 
diag(Wcor) <- 1 # set the diagonal to 1

# step 2: mean vector for W matrix
# for selected features, we set the theta as log(0.5p), for non-selected features, we set the theta as 0
theta <- rep(0, p)
theta[p_true] <- log(0.5 * p)

# step 3: generate W matrix
W <- MASS::mvrnorm(n, mu = theta, Sigma = Wcor)

X <- W 
X_taro <- X

# X <- sweep(X, 1, geo_mean, FUN = "-")
X <- scale(X, center = TRUE, scale = TRUE)
X_taro <- scale(X_taro, center = TRUE, scale = TRUE)

# Error matrix 
sigma <- mean(abs(C[p_true, q_true])) / snr 
set.seed(123)
E <- matrix(rnorm(n * q, mean = 0, sd = sigma), nrow = n, ncol = q)

# Outcome matrix 
Y_true <- X %*% C + E # Y = X * C + E

# missing value 
if (gen_miss == "false") {
  Y_miss <- Y_true
  data_type <- "complete"
  cat("No missing values generated.\n")
} else if (missing_type == "mnar0.33") {
  Y_miss <- Y_true 
  
  # 1/3 mnar 
  lod <- quantile(Y_true, probs = c(proportion * (1/3))) # calculate the lod
  Y_miss[Y_miss < lod] <- 0 # set values below the lod to 0
  
  # 2/3 mar 
  remain_idx <- which(Y_miss != 0) # remaining values
  miss_idx <- sample(remain_idx, size = ceiling((n*q) * proportion * (2/3))) # randomly sample 2/3 of the remaining values
  miss_idx <- arrayInd(miss_idx, .dim = c(n, q))
  Y_miss[miss_idx] <- 0 # set missing values to 0
  
  data_type <- "mnar0.33"
  cat("Proportion of missing values:", sum(Y_miss == 0)/(n * q), "\n") # proportion of missing values
} else if (missing_type == "mnar0.67") {
  Y_miss <- Y_true 
  
  # 2/3 mnar 
  lod <- quantile(Y_true, probs = c(proportion * (2/3))) # calculate the lod
  Y_miss[Y_miss < lod] <- 0 # set values below the lod to 0
  
  # 1/3 mar 
  remain_idx <- which(Y_miss != 0) # remaining values
  miss_idx <- sample(remain_idx, size = ceiling((n*q) * proportion * (1/3))) # randomly sample 1/3 of the remaining values
  miss_idx <- arrayInd(miss_idx, .dim = c(n, q))
  Y_miss[miss_idx] <- 0 # set missing values to 0
  
  data_type <- "mnar0.67"
  cat("Proportion of missing values:", sum(Y_miss == 0)/(n * q), "\n") # proportion of missing values
} 

mydata <- list(
  Y_miss = Y_miss,
  Y_true = Y_true,
  C = C,
  C_idx = C_idx,
  X = X,
  X_taro = X_taro,
  data_type = data_type
)

if (gen_miss == "true") {
  cat("Missing data saved.\n")
  
  writeMat(paste0(mydata$data_type, "_n", n, "_q", q, "_p", p, "_snr", snr, "_prop", proportion, ".mat"),
           X = mydata$X,
           X_taro = mydata$X_taro,
           Y_true = mydata$Y_true,
           Y_miss = mydata$Y_miss,
           C = mydata$C,
           C_idx = mydata$C_idx,
           data_type = mydata$data_type)
  
  sum(mydata$Y_miss != mydata$Y_true)/(n*q) # check the number of missing values
  cat("Missing proportion:", sum(mydata$Y_miss == 0)/(n * q), "\n")
  cat("Dimensions: X =", dim(X), ", Y_true =", dim(Y_true), ", C =", dim(C), "\n")
  cat("SNR used:", snr, ", Implied sigma:", sigma, "\n")
  
  prop_missing <- colMeans(Y_miss ==0)
  prop_missing
  
} else if (gen_miss == "false") {
  cat("Complete data saved.\n")
  
  # save the matlab file
  writeMat(paste0("complete_n", n, "_q", q, "_p", p, "_snr", snr, ".mat"),
           X = mydata$X,
           X_taro = mydata$X_taro,
           Y_true = mydata$Y_true,
           Y_miss = mydata$Y_miss,
           C = mydata$C,
           C_idx = mydata$C_idx,
           data_type = mydata$data_type)
}

names(mydata)