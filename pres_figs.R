rm(list=ls())

source("model_functions.R")

library(ggplot2)
library(MASS)

n <- 50
p <- 2
beta <- rep(1, p)
intercepts <- c(0, 3)
K <- 3
dev_prob <- 1/3
dev_size <- .5
cat_colors <- c("purple", "red")

set.seed(1254721)

# Generate uniformly distributed X
if(p > 1){
    X <- matrix(runif(n*p, min=-1, max=1), nrow=n, ncol=p)
} else{
    X <- runif(n, min=-1, max=1)
}

# Generate beta_mat

beta_mat <- matrix(0, nrow=p, ncol=K - 1)
beta_mat[, 1] <- beta

# beta_mat[, 2] <- beta

# Generate sparse random deviations
dev_mat <- matrix(0, nrow=p, ncol=K - 2)

for(k in 1:(K - 2)){
    # Select random indices for deviations
    inds <- which(as.logical(rbinom(p, size=1, prob=dev_prob)))
    dev_mat[inds, k] <- dev_size
}

# Create random sign flips to multiply deviations by
sign_mat_content <- sample(c(-1, 1), size=p*(K - 2), replace=TRUE)
sign_mat <- matrix(sign_mat_content, nrow=p, ncol=K - 2)

# Multiply deviations by random sign flips
dev_mat <- dev_mat * sign_mat

# Add sparse random deviations to beta_mat
for(k in 2:(K-1)){
    beta_mat[, k] <- beta_mat[, k - 1] + dev_mat[, k - 1]
}

# Generate X, y, probabilities
ret <- sim_data(n, p, K, intercepts, beta_mat, dev_num, dev_size)

X <- ret$X
y <- ret$y
beta_final <- ret$beta_final
ycoarse <- y
ycoarse[ycoarse==2] <- 1
ycoarse[ycoarse==3] <- 2
ycoarse <- factor(ycoarse, levels=c("1", "2"), ordered=FALSE)

ylogit1 <- y
ylogit1[ylogit1==3] <- 2
ylogit1 <- factor(ylogit1, levels=c("1", "2"), ordered=FALSE)

rm(ret)

convertToSlopeIntcpt <- function(intcpt, beta){
	stopifnot(length(intcpt) == 1)
	stopifnot(length(beta) == 2)

	intcpt_ret <- -intcpt/beta[2]
	slope_ret <- -beta[1]/beta[2]

	return(list(intcpt=intcpt_ret, slope=slope_ret))
}

# Bayes Decision boundaries

intcp_1 <- convertToSlopeIntcpt(intercepts[1], beta_mat[, 1])$intcpt
slope_1 <- convertToSlopeIntcpt(intercepts[1], beta_mat[, 1])$slope

intcp_2 <- convertToSlopeIntcpt(intercepts[2], beta_mat[, 2])$intcpt
slope_2 <- convertToSlopeIntcpt(intercepts[2], beta_mat[, 2])$slope

# Estimate decision boundaries via logistic regression
df1 <- data.frame(X, ycoarse)
colnames(df1) <- c("x1", "x2", "y")

df2 <- data.frame(X, y)
colnames(df2) <- c("x1", "x2", "y")

dflogit2 <- data.frame(X, ylogit1)
colnames(dflogit2) <- c("x1", "x2", "y") 

logit1 <- glm(y ~., family=binomial, data=df1)

logit2 <- glm(y ~., family=binomial, data=dflogit2)

est_intcp_1 <- convertToSlopeIntcpt(coef(logit1)[1], coef(logit1)[2:3])$intcpt
est_slope_1 <- convertToSlopeIntcpt(coef(logit1)[1], coef(logit1)[2:3])$slope

est_intcp_2 <- convertToSlopeIntcpt(coef(logit2)[1], coef(logit2)[2:3])$intcpt
est_slope_2 <- convertToSlopeIntcpt(coef(logit2)[1], coef(logit2)[2:3])$slope

# Proportional odds model
df_po <- data.frame(X, y)
colnames(df_po) <- c("x1", "x2", "y")

prop_odds <- polr(y ~., df_po)

po_est_intcp_1 <- -convertToSlopeIntcpt(prop_odds$zeta[1],
	coef(prop_odds))$intcpt
po_est_slope_1 <- convertToSlopeIntcpt(prop_odds$zeta[1], coef(prop_odds))$slope

po_est_intcp_2 <- -convertToSlopeIntcpt(prop_odds$zeta[2],
	coef(prop_odds))$intcpt
po_est_slope_2 <- convertToSlopeIntcpt(prop_odds$zeta[2], coef(prop_odds))$slope

# Plot data and decision boundaries

plot1 <- ggplot(df1, aes(x=x1, y=x2, color=y)) + geom_point() +
	geom_abline(slope=slope_2, intercept=intcp_2, color=cat_colors[2],
		linetype="dashed", size=1) + xlim(-1.5, 1) + ylim(-1.2, 1) +
	scale_color_manual(values = c("blue", cat_colors[2]))
print(plot1)

plot1_est <- ggplot(df1, aes(x=x1, y=x2, color=y)) + geom_point() +
	geom_abline(slope=slope_2, intercept=intcp_2, color=cat_colors[2],
		linetype="dashed", size=1) + xlim(-1.5, 1) + ylim(-1.2, 1) +
	scale_color_manual(values = c("blue", cat_colors[2])) +
	geom_abline(slope=est_slope_1, intercept=est_intcp_1, color="orange",
		linetype="dashed", size=1)
print(plot1_est)

plot2 <- ggplot(df2, aes(x=x1, y=x2, color=y)) + geom_point() +
	geom_abline(slope=slope_1, intercept=intcp_1, color="black",
		linetype="dashed", size=1) +
	geom_abline(slope=slope_2, intercept=intcp_2, color=cat_colors[2],
		linetype="dashed", size=1) + xlim(-1.5, 1) + ylim(-1.2, 1) +
	scale_color_manual(values = c("blue", cat_colors))
print(plot2)

plot2_est <- ggplot(df2, aes(x=x1, y=x2, color=y)) + geom_point() +
	geom_abline(slope=slope_1, intercept=intcp_1, color="black",
		linetype="dashed", size=1) +
	geom_abline(slope=slope_2, intercept=intcp_2, color=cat_colors[2],
		linetype="dashed", size=1) + xlim(-1.5, 1) + ylim(-1.2, 1) +
	scale_color_manual(values = c("blue", cat_colors)) +
	geom_abline(slope=est_slope_2, intercept=est_intcp_2, color="orange",
		linetype="dashed", size=1)
print(plot2_est)

plot2_po <- ggplot(df2, aes(x=x1, y=x2, color=y)) + geom_point() +
	geom_abline(slope=slope_1, intercept=intcp_1, color="black",
		linetype="dashed", size=1) +
	geom_abline(slope=slope_2, intercept=intcp_2, color=cat_colors[2],
		linetype="dashed", size=1) + xlim(-2, 1) + ylim(-1.6, 1) +
	scale_color_manual(values = c("blue", cat_colors)) +
	# geom_abline(slope=est_slope_2, intercept=est_intcp_2, color="orange",
	# 	linetype="dashed", size=1) +
	geom_abline(slope=est_slope_1, intercept=est_intcp_1, color="orange",
		linetype="dashed", size=1) +
	geom_abline(slope=po_est_slope_1, intercept=po_est_intcp_1, color="green",
		linetype="dashed", size=1) +
	geom_abline(slope=po_est_slope_2, intercept=po_est_intcp_2, color="green",
		linetype="dashed", size=1)
print(plot2_po)

# TODO(gfaletto): create another figure estimating the proportional odds model
# on all three classes and show that it leads to a better estimate of the 
# decision boundary between classes 2 and 3.







