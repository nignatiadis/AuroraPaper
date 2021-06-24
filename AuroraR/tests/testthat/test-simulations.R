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



