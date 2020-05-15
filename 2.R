library(mclust)
library(ggplot2)


data <-read.csv("/Users/abishekvaithylingam/Desktop/Projects/AppliedStatisticalModeling/Assignment2/assignment2/winemag-data-130k-v2.csv")

data <- data[(data$country == "US") & (!is.na(data$price)),]

us_wine_price <- data$price
us_wine_points <- data$points

wine_data <- do.call(rbind, Map(data.frame, A=us_wine_points, B=us_wine_price))

colnames(wine_data) <- c("Points", "Price")

plot(wine_data, log='y')
fit <- Mclust(wine_data)


print(fit)
summary(fit)

fit$parameters$pro
fit$parameters$mean
plot(fit, what = "classification")
plot(fit, what = "uncertainty")
head(fit$classification)
head(fit$uncertainty)
plot(fit, what = "BIC")
fit$BIC


fit <- Mclust(wine_data, G = 9, modelNames = "VVV")

plot(fit, what = "classification", log='x')
plot(fit, what = "uncertainty", log='x')

summary(fit2)

region_mean <- aggregate(wine_italy_filtered$points, list(wine_italy_filtered$region_1), mean)
overall_mean <- mean(region_mean$x)
sum(region_mean$x > overall_mean)
region_mean$Group.1[region_mean$x > overall_mean]

