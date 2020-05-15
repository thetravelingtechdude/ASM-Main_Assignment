install.packages("magrittr")
install.packages("dplyr")    
install.packages("ggplot2")
install.packages("plotly")
install.packages("ggridges")
install.packages("ggpubr")
install.packages("viridis") 

library(viridis)
library(ggpubr)
library(magrittr) 
library(dplyr)
library(ggridges)
library(ggplot2)
library(plotly)

#Read the data
data <-read.csv("/Users/abishekvaithylingam/Desktop/Projects/AppliedStatisticalModeling/Assignment2/assignment2/winemag-data-130k-v2.csv")

#Filter dataset to contain reviews about Italian wines < 20 euros
data <- data[(data$country == "Italy") & (data$price < 20),]

#Size of dataset after filtering
dim(data)

#Handling missing values
data <- data[!is.na(data$points),]

table1 <- table(data$region_1)

#Filter  out countries that have less than 4 reviews
data1 <- subset(data, region_1 %in% names(table1[table1 > 3]))

data1 <- data1 %>% group_by(region_1)
nlevels(data1$region_1)
data1 <- droplevels(data1)
region_levels <- nlevels(data1$region_1)


ggplot(data1) + geom_boxplot(aes(x = reorder(region_1, points, median), points, fill = reorder(region_1, points, median)), show.legend=FALSE)
ggplot(data1, aes(x = reorder(region_1, region_1, length))) + stat_count()							   
ggplot(data1, aes(points)) + stat_bin()

ggplot(data.frame(size = tapply(data1$points, data1$region_1, length),
                  mean_score = tapply(data1$points, data1$region_1, mean)), aes(size, mean_score)) +
  geom_point() + xlab("Region 1 sample size") + ylab("Mean Points") +
  ggtitle("Effect size versus sample size")

compare_m_gibbs <- function(y, ind, mu0 = 50, tau0 = 1/400,
                            a0 = 1, b0 = 50, alpha0 = 1, beta0 = 50, maxiter = 5000)
{
  ### weakly informative priors
  
  a0 <- 1/2 ; b0 <- 50 ## tau_w hyperparameters
  alpha0 <-1/2 ; beta0 <- 50 ## tau_b hyperparameters
  mu0<-50 ; tau0 <- 1/25
  
  y <- y + rnorm(length(y),1,1)/10000
  
  
  ### starting values
  
  m <- nlevels(ind)
  ybar <- theta <- tapply(y, ind, mean)
  tau_w <- mean(1 / tapply(y, ind, var)) ##within group precision
  mu <- mean(theta)
  tau_b <-var(theta) ##between group precision
  n_m <- tapply(y, ind, length)
  alphan <- alpha0 + sum(n_m)/2
  
  ### setup MCMC
  
  theta_mat <- matrix(0, nrow=maxiter, ncol=m)
  mat_store <- matrix(0, nrow=maxiter, ncol=3)
  
  ### MCMC algorithm
  for(s in 1:maxiter)
  {
    # sample new values of the thetas
    
    for(j in 1:m)
    {
      taun <- n_m[j] * tau_w + tau_b
      thetan <- (ybar[j] * n_m[j] * tau_w + mu * tau_b) / taun
      theta[j]<-rnorm(1, thetan, 1/sqrt(taun))
    }
    
    #sample new value of tau_w
    ss <- 0
    for(j in 1:m){
      ss <- ss + sum((y[ind == j] - theta[j])^2)
    }
    betan <- beta0 + ss/2
    tau_w <- rgamma(1, alphan, betan)
    
    #sample a new value of mu
    taum <- m * tau_b + tau0
    mum <- (mean(theta) * m * tau_b + mu0 * tau0) / taum
    mu <- rnorm(1, mum, 1/ sqrt(taum))
    
    # sample a new value of tau_b
    am <- a0 + m/2
    bm <- b0 + sum((theta - mu)^2) / 2
    tau_b <- rgamma(1, am, bm)
    
    #store results
    theta_mat[s,] <- theta
    mat_store[s, ] <- c(mu, tau_w, tau_b)
  }
  
  colnames(mat_store) <- c("mu", "tau_w", "tau_b")
  return(list(params = mat_store, theta = theta_mat))
}

fit <- compare_m_gibbs(data1$points, data1$region_1)

apply(fit$params, 2, mean)

theta_hat <- apply(fit$theta, 2, mean)

names(theta_hat) <- 1:region_levels

sort(theta_hat, decreasing = TRUE)

theta_ci <- apply(fit$theta, 2, quantile, prob = c(0.025, .975))

df_error <- data.frame(lower = theta_ci[1, ], upper = theta_ci[2, ], mean = theta_hat, region_1 = factor(1:168))

ggplot(df_error, aes(x = reorder(region_1, mean), mean)) + geom_errorbar(aes(ymin = lower, ymax = upper))					   

theta_df <- data.frame(samples = as.numeric(fit$theta), region_1 = rep(1:ncol(fit$theta), each = nrow(fit$theta)))

ggplot(theta_df) + geom_boxplot(aes(x = reorder(region_1, samples, median), samples, fill = reorder(region_1, samples, median)), show.legend=FALSE)

ggplot(data.frame(size = tapply(data1$points, data1$region_1, length), theta_hat = theta_hat), aes(size, theta_hat)) + geom_point()

ggplot(data.frame(ybar = tapply(data1$points, data1$region_1, mean), theta_hat = theta_hat), aes(ybar, theta_hat)) + geom_point()

mu <- sort(fit$params[1:region_levels,1], decreasing = TRUE)

regions <- sum(theta_hat > mu)
regions