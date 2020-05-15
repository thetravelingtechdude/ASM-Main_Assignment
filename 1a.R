library(dplyr)
library(ggplot2)
library(MCMCpack)

data <-read.csv("/Users/abishekvaithylingam/Desktop/Projects/AppliedStatisticalModeling/Assignment2/assignment2/winemag-data-130k-v2.csv")

df1 <- data[(data$country == 'South Africa') & (data$variety == 'Sauvignon Blanc') & (data$price == 15),]
df2 <- data[(data$country == 'Chile') & (data$variety == 'Chardonnay') & (data$price == 15),]
df3 <- bind_rows(df1, df2)
df3 <- df3[!is.na(df3$points),]

ggplot(df3) + geom_boxplot(aes(variety, points, fill = variety)) + geom_jitter(aes(variety, points, shape = variety))

t.test(points ~ variety, data=df3, var.equal = TRUE)

compare_2_gibbs <- function(df3, variety, mu0 = 50, tau0 = 1/400, del0 = 0, gamma0 = 1/400,
                            a0 = 1, b0 = 50, maxiter = 5000)
{
  y1 <- df3[variety == 'Sauvignon Blanc']
  y2 <- df3[variety == 'Chardonnay']
  n1 <- length(y1)
  n2 <- length(y2)
  ##### starting values
  mu <- (mean(y1) + mean(y2)) / 2
  del <- (mean(y1) - mean(y2)) / 2
  mat_store <- matrix(0, nrow = maxiter, ncol = 3)
  #####
  ##### Gibbs sampler
  an <- a0 + (n1 + n2)/2
  for(s in 1 : maxiter)
  {
    ##update tau
    bn <- b0 + 0.5 * (sum((y1 - mu - del) ^ 2) + sum((y2 - mu + del) ^ 2))
    tau <- rgamma(1, an, bn)
    ##
    ##update mu
    taun <- tau0 + tau * (n1 + n2)
    mun <- (tau0 * mu0 + tau * (sum(y1 - del) + sum(y2 + del))) / taun
    mu <- rnorm(1, mun, sqrt(1/taun))
    ##
    ##update del
    gamman <- tau0 + tau*(n1 + n2)
    deln <- ( del0 * tau0 + tau * (sum(y1 - mu) - sum(y2 - mu))) / gamman
    del<-rnorm(1, deln, sqrt(1/gamman))
    ##
    ## store parameter values
    mat_store[s, ] <- c(mu, del, tau)
  }
  colnames(mat_store) <- c("mu", "del", "tau")
  return(mat_store)
}

fit <- compare_2_gibbs(df3$points, as.factor(df3$variety))
plot(as.mcmc(fit))
raftery.diag(as.mcmc(fit))
better_value <- apply(fit, 2, mean)
better <- better_value[2] * 2
better

y1_sim <- rnorm(5000, fit[, 1] + fit[, 2], sd = 1/sqrt(fit[, 3]))
y2_sim <- rnorm(5000, fit[, 1] - fit[, 2], sd = 1/sqrt(fit[, 3]))
ggplot(data.frame(y_sim_diff = y1_sim - y2_sim)) + stat_bin(aes(y_sim_diff))
mean(y1_sim > y2_sim)
ggplot(data.frame(y1_sim, y2_sim)) + geom_point(aes(y1_sim, y2_sim), alpha = 0.3) +
  geom_abline(slope = 1, intercept = 0)