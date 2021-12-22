# Data 1: Friedman Example
Friedman = function(n, p, sigma, binary) {
  
  X = matrix(runif(n*p), nrow = n, ncol = p)
  f0 = 10 * sin(pi * X[, 1] * X[,2]) + 20 * (X[,3] - 0.5)^2 + 10 * X[,4] + 5 * X[,5]
  X = as.data.frame(X)
  dimnames(X)[[2]] = sapply(1:p, function(s) paste("X", s, sep=""))
  
  if (binary) {
    f0 = scale(f0, center = T, scale = F)
    p = pnorm(f0, mean = 0, sd = 1)
    Y = c()
    for (i in 1:n) {
      Y[i] = rbinom(1, 1, prob = p[i])
    }
    return(list(X=X, Y=Y, f0=f0, p=p))
    
  } else {
    Y = rnorm(n, mean = f0, sd = sigma)
    return(list(X=X, Y=Y, f0=f0, sigma=sigma))
  }
  
}

# Data 2: Checkerboard
Checkerboard = function(n, p, sigma, binary) {
  
  cov_mtx = matrix(NA, nrow = p, ncol = p)
  for (i in 1:p) {
    for (j in 1:p) {
      cov_mtx[i, j] = 0.3^abs(i-j)
    }
  }
  X = rmvnorm(n, sigma = cov_mtx)
  
  f0 = 2*X[,1]*X[,4] + 2*X[,7]*X[,10]
  X = as.data.frame(X)
  dimnames(X)[[2]] = sapply(1:p, function(s) paste("X", s, sep=""))
  
  if (binary) {
    f0 = scale(f0, center = T, scale = F)
    p = pnorm(f0, mean = 0, sd = 1)
    Y = c()
    for (i in 1:n) {
      Y[i] = rbinom(1, 1, prob = p[i])
    }
    return(list(X=X, Y=Y, f0=f0, p=p))
    
  } else {
    Y = rnorm(n, mean = f0, sd = sigma)
    return(list(X=X, Y=Y, f0=f0, sigma=sigma))
  }
  
}

# Data 3: Linear
Linear = function(n, p, sigma, binary) {
  
  cov_mtx = matrix(rep(0.6, p^2), nrow = p, ncol = p) + diag(x = 0.4, nrow = p)
  X = rmvnorm(n, sigma = cov_mtx)
  
  f0 = X[,1] + 2*X[,2] + 3*X[,3] - 2*X[,4] - X[,5]
  X = as.data.frame(X)
  dimnames(X)[[2]] = sapply(1:p, function(s) paste("X", s, sep=""))
  
  if (binary) {
    f0 = scale(f0, center = T, scale = F)
    p = pnorm(f0, mean = 0, sd = 1)
    Y = c()
    for (i in 1:n) {
      Y[i] = rbinom(1, 1, prob = p[i])
    }
    return(list(X=X, Y=Y, f0=f0, p=p))
    
  } else {
    Y = rnorm(n, mean = f0, sd = sigma)
    return(list(X=X, Y=Y, f0=f0, sigma=sigma))
  }
  
}

# Data 4: Tree (not finished yet)
Tree = function(n, p, sigma, binary) {
  cov_mtx = matrix(NA, nrow = p, ncol = p)
  for (i in 1:p) {
    for (j in 1:p) {
      cov_mtx[i, j] = 0.3^abs(i-j)
    }
  }
  X = rmvnorm(n, sigma = cov_mtx)
  X = as.data.frame(X)
  dimnames(X)[[2]] = sapply(1:p, function(s) paste("X", s, sep=""))
  
  f0 = prior_trees(x.train = X[,1:5], sigmaf = sqrt(3))$fp
  
  if (binary) {
    f0 = scale(f0, center = T, scale = F)
    p = pnorm(f0, mean = 0, sd = 1)
    Y = c()
    for (i in 1:n) {
      Y[i] = rbinom(1, 1, prob = p[i])
    }
    return(list(X=X, Y=Y, f0=f0, p=p))
    
  } else {
    Y = rnorm(n, mean = f0, sd = sigma)
    return(list(X=X, Y=Y, f0=f0, sigma=sigma))
  }
}

# Data 5: Mixed-type I
MixOne = function(n, p, sigma, binary) {
  k = p/2
  X = matrix(NA, nrow = n, ncol = p)
  X[, 1:k] = matrix(rbinom(n*k, size = 1, prob = 0.5), nrow = n, ncol = k)
  X[, (k+1):p] = matrix(runif(n*k), nrow = n, ncol = k)
  
  f0 = 10 * sin(pi * X[, k+1] * X[, k+2]) + 20 * (X[, k+3] - 0.5)^2 + 10 * X[, 1] + 5 * X[, 2]
  #f0 = 10 * sin(pi * X[, 1] * X[, k+1]) + 20 * (X[, k+3] - 0.5)^2 + 10 * X[, 2] + 5 * X[, k+2]
  X = as.data.frame(X)
  for (j in 1:k) {
    X[, j] = factor(X[, j], levels = c("0", "1"), labels = c("0", "1"))
  }
  dimnames(X)[[2]] = sapply(1:p, function(s) paste("X", s, sep=""))
  
  if (binary) {
    f0 = scale(f0, center = T, scale = F)
    p = pnorm(f0, mean = 0, sd = 1)
    Y = c()
    for (i in 1:n) {
      Y[i] = rbinom(1, 1, prob = p[i])
    }
    return(list(X=X, Y=Y, f0=f0, p=p))
    
  } else {
    Y = rnorm(n, mean = f0, sd = sigma)
    return(list(X=X, Y=Y, f0=f0, sigma=sigma))
  }
  
}

# Data 6: Mixed-type II
MixTwo = function(n, sigma, binary) {
  X = matrix(NA, nrow = n, ncol = 84)
  X[, 1:20] = matrix(rbinom(n*20, size = 1, prob = 0.2), nrow = n, ncol = 20)
  X[, 21:40] = matrix(rbinom(n*20, size = 1, prob = 0.5), nrow = n, ncol = 20)
  covmat = matrix(rep(0.3, 1936), nrow = 44, ncol = 44) + diag(x = 0.7, nrow = 44)
  X[, 41:84] = rmvnorm(n, mean = rep(0, 44), sigma = covmat)
  
  f0 = -4 + X[, 1] + sin(pi * X[, 1] * X[, 44]) - X[, 21] + 
    0.6 * X[, 41] * X[, 42] - exp(-2 * (X[, 42] + 1)^2) - X[, 43]^2 + 
    0.5 * X[, 44]
  X = as.data.frame(X)
  for (j in 1:40) {
    X[, j] = factor(X[, j], levels = c("0", "1"), labels = c("0", "1"))
  }
  dimnames(X)[[2]] = sapply(1:84, function(s) paste("X", s, sep=""))
  
  if (binary) {
    f0 = scale(f0, center = T, scale = F)
    p = pnorm(f0, mean = 0, sd = 1)
    Y = c()
    for (i in 1:n) {
      Y[i] = rbinom(1, 1, prob = p[i])
    }
    return(list(X=X, Y=Y, f0=f0, p=p))
    
  } else {
    Y = rnorm(n, mean = f0, sd = sigma)
    return(list(X=X, Y=Y, f0=f0, sigma=sigma))
  }
  
}