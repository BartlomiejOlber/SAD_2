#### set params of experiment
set.seed(2137)
alpha <- 9
beta <- 6
start <- 10
step <- 10
no_of_steps <- 100
no_of_simulations_in_single_step <- 1000


#### simualte alpha and beta parameters
simulate_bias_and_variance <- function(alpha, beta, sample_size){
  
  # 1. calculate maximum likelihood estimator
  # according to wikipedia estimator of alpha = sample_mean of x * beta
  alpha_estimators = c()
  beta_estimators = c()
  for (i in 1:no_of_simulations_in_single_step){
    random_sample <- rgamma(sample_size, shape=alpha, rate=beta)

    alpha_estimator <- mean(random_sample) * beta
    alpha_estimators <- c(alpha_estimators, alpha_estimator)
    
    beta_estimator <- alpha / mean(random_sample) 
    beta_estimators <- c(beta_estimators, beta_estimator)
  }
  final_alpha_est <- mean(alpha_estimators)
  final_beta_est <- mean(beta_estimators)
  
  # 2. calculate biases
  alpha_bias <- final_alpha_est - alpha
  beta_bias <- final_beta_est - beta
  
  # 3. calculate cr bound
  fisher_inf_alpha <- trigamma(final_alpha_est)
  cr_bound_alpha <- 1 / (sample_size * fisher_inf_alpha)
  
  fisher_inf_beta <- alpha / final_beta_est^2
  cr_bound_beta <- 1 / (sample_size * fisher_inf_beta)
  
  # 4. calculate variance
  variance_alpha <- mean((alpha_estimators - alpha)^2)
  variance_beta <- mean((beta_estimators - beta)^2)
  
  return (c(sample_size, final_alpha_est, final_beta_est, alpha_bias, beta_bias, variance_alpha, variance_beta, cr_bound_alpha, cr_bound_beta))
}

column_names <- c("sample_size", 'final_alpha_est', 'final_beta_est', 'alpha_bias', 'beta_bias', 'variance_alpha', 'variance_beta', 'cr_bound_alpha', 'cr_bound_beta')
results <- data.frame(matrix(ncol = length(column_names), nrow = 0))
colnames(results) <- column_names

for (i in 1:no_of_steps){
  sample_size <- start + step * (i-1)
  result <- simulate_bias_and_variance(alpha, beta, sample_size)

  result <- as.data.frame(t(result))
  colnames(result) <- column_names
  results <- rbind(results, result)
}

# plot bias
plot(results$alpha_bias ~ results$sample_size,
    main = 'Alpha over sample size', xlab = 'sample size', ylab = 'bias')
abline(h = 0, col = 'red', lty=2, lwd=2)

plot(results$beta_bias ~ results$sample_size,
     main = 'Beta Bias over sample size', xlab = 'sample size', ylab = 'bias')
abline(h = 0, col = 'red', lty=2, lwd=2)

# plot cr-bound and variance
plot(
  results$cr_bound_alpha ~ results$sample_size,
  main = 'alpha Variance and CR bound over sample size',
  xlab = 'sample size', ylab = 'value',
  type = 'l', col = 'red', lwd = 2, ylim=c(0, 0.5), log = 'x'
)
lines(results$sample_size, results$variance_alpha, col='blue', lwd = 1)
legend(
  x = 'topright', c('variance', 'cr-bound'), col = c('blue', 'red'),
  lty=c(1, 1), lwd = 2
)

plot(
  results$cr_bound_beta ~ results$sample_size,
  main = 'Beta Variance and CR bound over sample size',
  xlab = 'sample size', ylab = 'value',
  type = 'l', col = 'red', lwd = 2, ylim=c(0, 0.5), log = 'x'
)
lines(results$sample_size, results$variance_beta, col='blue', lwd = 1)
legend(
  x = 'topright', c('variance', 'cr-bound'), col = c('blue', 'red'),
  lty=c(1, 1), lwd = 2
)


#### simualte estimators normality
simulate_normality <- function(alpha, beta, number_of_estimators, sample_size){
  
  alpha_estimators = c()
  beta_estimators = c()
  for (i in 1:number_of_estimators){
    random_sample <- rgamma(sample_size, shape=alpha, rate=beta)
    
    alpha_estimator <- mean(random_sample) * beta
    alpha_estimators <- c(alpha_estimators, alpha_estimator)
    
    beta_estimator <- alpha / mean(random_sample) 
    beta_estimators <- c(beta_estimators, beta_estimator)
  }
  
  
  norm_alpha_test = ks.test(alpha_estimators, "pnorm", mean=mean(alpha_estimators), sd=sd(alpha_estimators))  
  
  norm_beta_test = ks.test(beta_estimators, "pnorm", mean=mean(beta_estimators), sd=sd(beta_estimators))  
  
  
  normalized_alpha <- sqrt(sample_size) * (alpha_estimators - alpha)
  hist(normalized_alpha, prob=TRUE, xlab="Value", main=paste('Normalized alpha for sample size = ', as.character(sample_size)), xlim=c(-10,10), ylim=c(0, 0.5), breaks = 100)
  curve(dnorm(x, mean=0, sd=sd(normalized_alpha)), add=TRUE, col = 'red', lwd = 2)

  normalized_beta <- sqrt(sample_size) * (beta_estimators - beta)
  hist(normalized_beta, prob=TRUE, xlab="Value", main=paste('Normalized beta for sample size = ', as.character(sample_size)), xlim=c(-6,6), ylim=c(0, 0.5), breaks = 100)
  curve(dnorm(x, mean=0, sd=sd(normalized_beta)), add=TRUE, col = 'red', lwd = 2)
  
  return (c(number_of_estimators, sample_size, norm_alpha_test$p.value, norm_beta_test$p.value, norm_alpha_test$statistic, norm_beta_test$statistic))
}

column_names <- c("number_of_estimators", "sample_size", "normality_pvalue_alpha", "normality_pvalue_beta")
results <- data.frame(matrix(ncol = length(column_names), nrow = 0))
colnames(results) <- column_names

start <- 10
step <- 1000
no_of_steps <- 3

for (i in 1:no_of_steps){
  sample_size <- start + step * (i-1)
  result <- simulate_normality(alpha, beta, number_of_estimators=2000, sample_size=sample_size)
  
  result <- as.data.frame(t(result))
  colnames(result) <- column_names
  results <- rbind(results, result)
}

plot(
  results$normality_pvalue_alpha ~ results$sample_size, 
  main = 'alpha normality (p > 0.05 signifies normality)', 
  xlab = 'sample size', ylab = 'Kolmogorov-Smirnov test p-value', 
  type = 'l', col = 'red', lwd = 2, ylim=c(0, 1),
)

abline(h = 0.05, col="blue")

plot(
  results$normality_pvalue_beta ~ results$sample_size, 
  main = 'beta normality (p > 0.05 signifies normality)', 
  xlab = 'sample size', ylab = 'Kolmogorov-Smirnov test p-value', 
  type = 'l', col = 'red', lwd = 2, ylim=c(0, 1),
)
abline(h = 0.05, col="blue")