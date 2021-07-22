A_prior <- 9/1.5
sigma <- 2.2
B <- 100
n <- 100000
set.seed(1)
sim_three_comp_normal <- three_component_normal_sim(n, B, sigma, A_prior)
sim_three_comp_unif   <- three_component_uniform_sim(n, B, sigma, A_prior)
sim_three_comp_laplace   <- three_component_laplace_sim(n, B, sigma, A_prior)

for (sim_tmp in list(sim_three_comp_normal, sim_three_comp_unif, sim_three_comp_laplace)){
  expect_equal(sim_tmp$A, A_prior)
  expect_equal(sim_tmp$sigma, sigma)
  expect_equal(sim_tmp$B, B)
  expect_equal(length(sim_tmp$true_mus), n)
  expect_equal(sim_tmp$m, n)
  expect_equal(dim(sim_tmp$Zs), c(n,B))
  expect_equal(max(sim_tmp$true_mus), 3)
  expect_equal(min(sim_tmp$true_mus), -3)
  expect_true( abs(A_prior - var(sim_tmp$true_mus)) < 0.01*A_prior)
  expect_true( abs(sigma - mean(apply(sim_three_comp_normal$Zs, 1, sd))/sqrt(B)) < 0.01*sigma)
}

set.seed(2)
B_2 = 30
sigma_2 = 1
sim_three_comp_unif_2   <- three_component_uniform_sim(n*10, B_2, sigma_2, 1)
expect_true( abs(range(sim_three_comp_unif_2$Zs - sim_three_comp_unif_2$true_mus)[2] - sqrt(3*sigma_2^2*B_2)) < 0.001)

mu_midrange <- AuroraR::location_midrange(sim_three_comp_unif_2$Zs)
mu_mean <- AuroraR::location_mean(sim_three_comp_unif_2$Zs)

mse_midrange <- mean( (mu_midrange - sim_three_comp_unif_2$true_mus)^2)
mse_mean <- mean( (mu_mean - sim_three_comp_unif_2$true_mus)^2)

expect_true(abs( mse_midrange/mse_mean - 6*B_2/( (B_2+1)*(B_2+2))) < 0.001)

unif_max_squared <- 3*sigma_2^2*B_2

expect_true(abs(mse_midrange - 2*unif_max_squared/( (B_2+1)*(B_2+2))) < 0.001)
expect_true(abs(mse_midrange - 6*sigma_2^2*B_2/( (B_2+1)*(B_2+2))) < 0.001)





# Test Pareto

set.seed(1)
tmp_alpha <- 6
tmp_B <- 5000
tmp_m <- 10000
pareto_sim <- uniform_pareto_sim(m=tmp_m, B=tmp_B, pareto_alpha=tmp_alpha, unif_lb=0.3, unif_ub=5)

expect_true(range(pareto_sim$true_mus)[1] >= 0.3)
expect_true(range(pareto_sim$true_mus)[2] <= 5)

pareto_mles <- pareto_mle(pareto_sim$Zs)
pareto_means <- location_mean(pareto_sim$Zs)

mean(pareto_mles - pareto_sim$true_mus)
mean(pareto_means - pareto_sim$true_mus)

pareto_mle_error <- mean( (pareto_mles - pareto_sim$true_mus)^2)

pareto_mean_error <- mean( (pareto_means - pareto_sim$true_mus)^2)
expect_true(pareto_mean_error > pareto_mle_error)
expect_true( (pareto_mean_error - mean(pareto_sim$sigmas^2 / 100 ))/pareto_mean_error <= 0.01)

true_xs <- pareto_sim$true_mus*(tmp_alpha - 1)/tmp_alpha
true_medians <- true_xs * 2^(1/tmp_alpha)

true_pdfs <- tmp_alpha*true_xs^(tmp_alpha)/( true_medians^(tmp_alpha+1))
medians_vars <- 1/(tmp_B*4*true_pdfs^2)
  
  
pareto_medians <- location_median(pareto_sim$Zs)
expect_true(length(pareto_medians) == tmp_m)

expect_true( abs(mean(pareto_medians - true_medians)) < sqrt(mean(medians_vars)/tmp_m)*4) 
abs(mean(pareto_medians - pareto_sim$true_mus))

pareto_median_error <- mean( (pareto_medians - pareto_sim$true_mus)^2)

pareto_median_median_error <- mean( (pareto_medians -true_medians)^2)
expect_true( pareto_median_median_error <= 1.1*mean(medians_vars))
expect_true( pareto_median_median_error >= 0.9*mean(medians_vars))

expect_true(pareto_median_error > pareto_median_median_error)
expect_true(pareto_mean_error < pareto_median_error)
