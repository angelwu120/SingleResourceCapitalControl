#You are asked to write codes to implement some of the methods we talked about in class. 
#The inputs to your codes are:
#  capacity level C,
#  number of classes n,
#  fair price vector f,
#  means mu and standard deviations sigma of different demand classes (assuming normal distributed).
#(a) Compute the optimal protection levels.
#(b) Compute the EMSR-a protection levels.
#(c) Compute the EMSR-b protection levels.
#(d) Implement the procedure that combines censored forecasting with EMSR-b protection levels.
#(e) Implement the adaptive algorithm. Feel free to adjust the step
#    sizes to see the difference.
#    Please test your codes on the fare and demand data in Table 2.5 in the textbook 
#    (The Theory and Practice of Revenue Management).

##############################
### INPUTS for (a) (b) (c) ###
##############################
C <- 1000 #capacity
n <- 4  #number of classes
f <- matrix(c(1050, 567, 534, 520), nrow = n, ncol = 1, byrow = TRUE)  #fair price
m <- matrix(c(17.3, 45.1, 39.6, 34.0), nrow = n, ncol = 1, byrow = TRUE)  #mean
s <- matrix(c(5.8, 15.0, 13.2, 11.3), nrow = n, ncol = 1, byrow = TRUE)  #deviation

###########################################
#(a) Compute the optimal protection levels#
###########################################
# Initial Setup
y <- array(data = NA, c(n,1)) #protection level
D <- array(data = 0, c(n,C)) #demand
u <- array(data = 0, c(n,C)) #aggregated future demand

z <- C
K <- array(data = 0, c(n,C))
K[1,] <- (1:z)

# Generate random demands from normal distribution N(m,s)
for (i in 1:n) {
    D[i,] <- rnorm(n = C, m = m[i], sd = s[i])
}

# Calculate aggregated future demands
for (i in 1:C) {
  for (j in 1:n) {
    for (k in 1:j) {
      u[j,i] <- u[j,i] + D[k,i]
    }
  }
}

# Calculate optimal protection level
for (i in 1:(n-1)) {
  U <- as.matrix(sort(u[i,K[i,(1:z)]],decreasing = FALSE))
  l <- floor((f[i+1] / f[i]) * z)
  if (l > 0) {
    y[i] <- (1/2) * (U[l] + U[l+1])
    k <- 0
  } else {
    y[i] <- U[l+1]
  }
  for (j in K[i,(1:z)]) {
    if (u[i,j] > y[i]) {
      k <- k + 1
      K[i+1,k] <- j
    }
  }
  z <- k
}

##########################################
#(b) Compute the EMSR-a protection levels#
##########################################
EMSRa <- function(C, n, f, m, s) {
  # Initial Setup
  ya <- array(data = NA, c(n,1)) #protection level for EMSR-a
  mya <- array(data = 0, c(n-1,n-1))
  
  # Calculate EMSR-a
  for (i in n:2) {
    for (j in (i-1):1) {
      mya[i-1,j] <- qnorm(1-f[i]/f[j],mean = m[j], sd = s[j])
    }
    ya[i-1] <- sum(mya[i-1,]) #EMSR-a
  }
  return(ya)
}
ya <- EMSRa(C, n, f, m, s)

##########################################
#(c) Compute the EMSR-b protection levels#
##########################################
EMSRb <- function(C, n, f, m, s) {
  # Initial setup
  yb <- array(data = NA, c(n,1)) #protection level for EMSR-b
  ff <- array(data = NA, c(n-1,1)) #average f
  ss <- array(data = NA, c(n-1,1)) #average s
  
  # Calsulate EMSR-b
  for (i in 1:(n-1)) {
    a <- 0 
    d <- 0 #mean
    v <- 0 #variance
    for (j in 1:i) {
      a <- a + f[j]*m[j]
      d <- d + m[j]
      v <- v + s[j]^2
    }
    ff[i] <- a/d
    ss[i] <- v^(1/2)
    yb[i] <- qnorm(1-f[i+1]/ff[i],mean = d, sd = ss[i]) #EMSR-b
  }
  return(yb)
}
yb <- EMSRb(C, n, f, m, s)

##########################
### INPUTS for (d) (e) ###
##########################
f2 <- matrix(c(1050, 567, 527, 350), nrow = n, ncol = 1, byrow = TRUE)  #fair price
m2 <- matrix(c(17.3, 45.1, 73.6, 19.8), nrow = n, ncol = 1, byrow = TRUE)  #mean
s2 <- matrix(c(5.8, 15.0, 17.4, 6.6), nrow = n, ncol = 1, byrow = TRUE)  #deviation

################################################################
#(d) Implement the procedure that combines censored forecasting#
################################################################
yd <- EMSRb(C, n, f2, m2, s2)

######################################
#(e) Implement the adaptive algorithm#
######################################
adapted <- function(C, n, f, m, s) {
  yy <- array(data = 0, c(n,1))
  myy <- array(data = NA, c(n,1))
  D2 <- array(data = 0, c(n,1))
  mD2 <- array(data = 0, c(n,1))
  
  for (i in 1:n) {
    yy[i] <- C/(2^(n-i))
  }
  
  for (i in 1:C) {
    for (l in 1:n) {
      D2[l] <- rnorm(n = 1, m = m2[l], sd = s2[l])
    }
    for (j in 1:(n-1)) {
      mD2[j] <- sum(D2[(1:j)])
      if (sum(mD2[(1:j),]>yy[(1:j)]) == j) {
        myy[j] <- yy[j] - (C/(n*i)) * (f2[j+1]/f2[1]-1)
      } else {
        myy[j] <- yy[j] - (C/(n*i)) * (f2[j+1]/f2[1])
      }
    }
    yy <- myy
  }
  return(yy)
}
yy <- adapted(C, n, f2, m2, s2)


##############
### OUTPUT ###
##############
# install.packages("abind") #!!!Run this line if the package has not been installed
library("abind")
abind(f,m,s,y,ya,yb,new.names=c("Price","Mean","StandardDeviation","y*","EMSR-a","EMSR-b"))
abind(f2,m2,s2,yd,yy,new.names=c("Price","Mean","StandardDeviation","y-EMSR","adapted-y*"))

