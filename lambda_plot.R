
library(glmnet)
URL <- "https://archive.ics.uci.edu/ml/machine-learning-databases/wine-quality/winequality-red.csv"
wine <- read.csv2(URL, na.strings="N/A", dec=".")


x <- model.matrix(~., wine[1:(dim(wine)[2]-1)])
y <- wine$quality
lambda_grid <- 10^seq(from = 2.5, to = -2.5, length = 100)
ridge_model <- glmnet(x, y, alpha = 0, albmda = lambda_grid)
lasso_model <- glmnet(x, y, alpha = 1, lambda = lambda_grid)
png(filename="lambda_plot.png", width=800, height=400)
par(mfrow=c(1,2))
plot(lasso_model, xvar='lambda', xlim=c(-5,5), ylim=c(-1,1))
plot(ridge_model, xvar='lambda', xlim=c(-5,5), ylim=c(-1,1))
dev.off()

