#Introduction to Software - R
#Assignment 1 - Out-Of-Sample Model Comparison 
#Lucas Hoff, i6203030
#Toby Pfeiffer, i6203090
library(matrixStats);  library(reshape2); library(ggplot2); library(mvtnorm)

#First we load our data, real GDP of the US and Japan respectively
#raw.df  <- read.csv("USA_EMP_txt.csv",TRUE, sep =",") #Lucas
raw.df <- read.csv("JPN_EMP_txt.csv",TRUE, sep =",") # Tobi

#next we extract all relevant data for the regression (13 corresponds to investment and 17 to the real GDP of the countries)
Y <- matrix(log(raw.df [,17][-1]))
X <- cbind(rep(1), log(lag(raw.df [,17][-65])), log(raw.df [,13][-1])) 
nvar <- ncol(X)

# Initialise Priors
#Prior Nr.1
# B0 <- as.matrix(c(rep(0, nvar))) #beata vector filled with means of our prior(s)
# sigma0 <- diag(1,nvar) #variance matrix
# T1 = nrow(Y) #number of observations n
# T0 = 1 #prior degrees of freedom
# w = 0.1 #prior scale (theta0)
# sigma2 = 1 #value for variance

#Prior Nr.2 USA
# B0 <- as.matrix(c(1,0.9,0.02))
# sigma0 <- matrix(c(s1,0,0,0,s2,0,0,0,s3), nrow = 3, ncol = 3)
# T1 = nrow(Y)
# T0 = tbd
# sigma2 = tbd

#Prior Nr.2 JPN
B0 <- as.matrix(c(1,0.9,0.03))
sigma0 <- matrix(c(0.64,0,0,0,0.005625,0,0,0,0.0225), nrow = 3, ncol = 3)
T1 = nrow(Y)
T0 = 1
sigma2 = 0.09

reps = 10000 #number of iterations
forcast_length = 3 #length of our forecast/out-of-sample obseravtions --> watch for logGDP in line 120

set.seed(123)

gibbs_sampler <- function(X,Y,sigma0,sigma2){
  
  nvar <- ncol(X)
  draws <- matrix(0, nrow = reps, ncol = nvar + 1) 
  
  forc <- matrix(0, nrow = reps, ncol = forcast_length)
  
  for(i in 1:reps){
    
    Var <- solve(sigma0 + (1/sigma2) * t(X) %*% X) #variance of posterior 
    
    Mean <- Var %*% (sigma0 %*% B0 + (1/sigma2) * t(X) %*% Y)  #mean of posterior
    
    B_hat <- t(rmvnorm(1,Mean,Var))  #draw for B_hat conditional on posterior distribution
    
    resids <- Y- X%*%B_hat
    
    #draw sigma Inverse Gamma
    
    sigma2 = 1/rgamma(1,shape = (T1-1)/2, scale = ((T1-1)*t(resids) %*% resids) / 2 )
    
    #we store the draws in the draws-matrix
    draws[i,] <- t(matrix(c(t(B_hat),sigma2)))
    
    #next, we perform a static forecast of our time frame 
    cfactor = sqrt(sigma2)
    #we initiate a vector of approprite length to store the forecasted values
    y_forecast <- vector("numeric", forcast_length) 
    
    #the comp-vector will contain regressor data: [1, y(t-1), x(t)]
    comp = c(1, rep(0,2))
    
    #we run through the last #time.frame elements of the data to forecast them
    for(h in (length(Y) - forcast_length+1):length(Y)){ 
      comp[2] = Y[h-1]
      comp[3] = X[,3][h]
      y_forecast[h-(length(Y)-forcast_length)] = comp %*% (B_hat)
    }
    forc[i,] <- y_forecast 
  }
  
  return = list(draws,forc)
}

gibbs <- gibbs_sampler(X,Y,sigma0,sigma2)

coef <- gibbs[[1]]
forecasts <- gibbs[[2]]



##################################################


b_0 <- mean(coef[,1])
b_1 <- mean(coef[,2])
b_2 <- mean(coef[,3])
sigma <- mean(coef[,4])

qplot(coef[,1], geom = "histogram", bins = 50, main = 'Dist of the Constant',colour="red")+ geom_vline(xintercept = b_0,colour = "green")

qplot(coef[,2], geom = "histogram", bins = 50,main = 'Dist of  Beta1',colour="red") + geom_vline(xintercept = b_1,colour = "green")

qplot(coef[,3], geom = "histogram", bins = 50,main = 'Dist of  Beta2',colour="red") + geom_vline(xintercept = b_2,colour = "green")

qplot(coef[,4], geom = "histogram", bins = 60,main = 'Dist of  Sigma',colour="red")+ geom_vline(xintercept = sigma,colour = "green")

quantiles <- colQuantiles(forecasts,prob = c(0.05,0.95))
Y_int = cbind(Y[1:(length(Y)-forcast_length)],Y[1:(length(Y)-forcast_length)])

HC_int <- rbind(Y_int, quantiles)

logGDP <- as.matrix(c(Y[1:(length(Y)-forcast_length)],mean(forecasts[,1]),mean(forecasts[,2]),mean(forecasts[,3])))

forecastsComb <- cbind.data.frame(HC_int[,1],logGDP, HC_int[,2])
names(forecastsComb) <- c('lower', 'logGDP', 'upper')

Date <-raw.df[,4][-1]

data.plot <- cbind.data.frame(Date, forecastsComb)

forcast.frame <- 64 - forcast_length

ggplot(data.plot[55:64,], aes(x = Date, y = logGDP)) + geom_line(colour = "blue", lwd = 1.2) + geom_ribbon(data = data.plot[forcast.frame:64,],aes(ymin = lower, ymax = upper ,  alpha = 0.2))

#plot the ACF functions to check convergence
acf(coef[,1])
acf(coef[,2])
acf(coef[,3])
acf(coef[,4])