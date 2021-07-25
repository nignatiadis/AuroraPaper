A_prior <- 9/1.5
sigma <- 2.2
B <- 100
n <- 100000
set.seed(1)

# Test homoskedastic
sim_three_comp_normal  <- three_component_normal_sim(n, B, sigma, A_prior)
sim_three_comp_unif    <- three_component_uniform_sim(n, B, sigma, A_prior)
sim_three_comp_laplace <- three_component_laplace_sim(n, B, sigma, A_prior)
sim_normal_normal  <- normal_normal_sim(n, B, sigma, A_prior)
sim_normal_unif    <- normal_uniform_sim(n, B, sigma, A_prior)
sim_normal_laplace <- normal_laplace_sim(n, B, sigma, A_prior)


for (sim_tmp in list(sim_three_comp_normal, sim_three_comp_unif, sim_three_comp_laplace,  
                     sim_normal_normal, sim_normal_unif, sim_normal_laplace)){
  expect_equal(sim_tmp$A, A_prior)
  expect_equal(sim_tmp$sigma, sigma)
  expect_equal(sim_tmp$B, B)
  expect_equal(length(sim_tmp$true_mus), n)
  expect_equal(sim_tmp$m, n)
  expect_equal(dim(sim_tmp$Zs), c(n,B))
  expect_true( abs(A_prior - var(sim_tmp$true_mus)) < 0.01*A_prior)
  expect_true( abs(sigma - mean(apply(sim_tmp$Zs, 1, sd))/sqrt(B)) < 0.01*sigma)
  expect_true( abs( sigma^2 - mean((location_mean(sim_tmp$Zs) - sim_tmp$true_mus)^2)) < 20*sigma/sqrt(n))
}

for (sim_tmp in list(sim_three_comp_normal, sim_three_comp_unif, sim_three_comp_laplace)){
  expect_equal(max(sim_tmp$true_mus), 3)
  expect_equal(min(sim_tmp$true_mus), -3)
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



sim_normal_normal2  <- normal_normal_sim(n, B, sigma, A_prior, prior_mu=10)
expect_true(abs(mean(sim_normal_normal2$true_mus) - 10) < 4*sqrt(A_prior)/sqrt(n) )
sim_normal_unif2    <- normal_uniform_sim(n, B, sigma, A_prior, prior_mu=10)
expect_true(abs(mean(sim_normal_unif2$true_mus) - 10) < 4*sqrt(A_prior)/sqrt(n) )
sim_normal_laplace2 <- normal_laplace_sim(n, B, sigma, A_prior, prior_mu=10)
expect_true(abs(mean(sim_normal_laplace2$true_mus) - 10) < 4*sqrt(A_prior)/sqrt(n) )
expect_true(sim_normal_normal2$sigma == sigma)

#sanity check for normal normal sim

js_res <- james_stein(sim_normal_normal2$Zs, sigma=sim_normal_normal2$sigma)
mean_res <- location_mean(sim_normal_normal2$Zs)
sure_res <- sure_grandmean(sim_normal_normal2$Zs, sim_normal_normal2$sigma)  
spher_shrink_res <- grouplinear(sim_normal_normal2$Zs, sim_normal_normal2$sigma)  

julia_dir = "/Applications/Julia-1.5.app/Contents/Resources/julia/bin"
aurora_setup(JULIA_HOME=julia_dir)

cc_res <- ebcc(sim_normal_normal2$Zs)
auroral_res <- auroral(sim_normal_normal2$Zs)


bayes_risk <- A_prior*sigma^2/(A_prior+sigma^2)
expect_true(abs( mean( (js_res - sim_normal_normal2$true_mus)^2) - bayes_risk) < 0.01)
expect_true(abs( mean( (sure_res - sim_normal_normal2$true_mus)^2) - bayes_risk) < 0.01)
expect_true(abs( mean( (spher_shrink_res - sim_normal_normal2$true_mus)^2) - bayes_risk) < 0.01)
expect_false(abs( mean( (mean_res - sim_normal_normal2$true_mus)^2) - bayes_risk) < 0.01)
expect_true(abs( mean( (mean_res - sim_normal_normal2$true_mus)^2) - sigma^2) < 0.03)

expect_true(abs( mean( (cc_res - sim_normal_normal2$true_mus)^2) - bayes_risk) < 0.01)
expect_true(abs( mean( (auroral_res - sim_normal_normal2$true_mus)^2) - bayes_risk) < 0.1)

# code below requires Mosek
#gmleb_res <- gmleb(sim_normal_normal2$Zs, sim_normal_normal2$sigma)  
#expect_true(abs( mean( (gmleb_res - sim_normal_normal2$true_mus)^2) - bayes_risk) < 0.03)





# Test heteroskedastic

set.seed(1)
l <- 1
u <- 2.1
B <- 50

sim_wmbz_a<- wmbz_example_a_sim(n, B, l, u, A_prior)
sim_wmbz_c<- wmbz_example_c_sim(n, B, l, u)
sim_wmbz_f<- wmbz_example_f_sim(n, B, l, u)





for (sim_tmp in list(sim_wmbz_a, sim_wmbz_c, sim_wmbz_f)){
  expect_equal(sim_tmp$sigma_squared_lower, l)
  expect_equal(sim_tmp$sigma_squared_upper, u)
  expect_true(all(   (sim_tmp$sigmas <= sqrt(u)) & (sim_tmp$sigmas > sqrt(l))    ))
  expect_true(  abs( mean(sim_tmp$sigmas) - (u+l)/2 ) < (l-u)^2 )
  
  expect_equal(sim_tmp$B, B)
  expect_equal(length(sim_tmp$true_mus), n)
  expect_equal(sim_tmp$m, n)
  expect_equal(dim(sim_tmp$Zs), c(n,B))

  expect_true( mean(abs(sim_tmp$sigmas - mean(apply(sim_tmp$Zs, 1, sd))/sqrt(B))) < 4*sigma/sqrt(B))
  expect_true( abs(mean(sim_tmp$sigmas^2 - mean(apply(sim_tmp$Zs, 1, var))/B)) < 4*sigma/sqrt(B)/sqrt(n))
}


expect_true( abs(A_prior - var(sim_wmbz_a$true_mus)) < 4*A_prior/sqrt(n))

# WMBZ A

wmbz_A_bayes_risk_fun <- function(sigma_squared_upper, sigma_squared_lower=0.1, A=0.5){
  l = sigma_squared_lower
  u <- sigma_squared_upper
  A - A^2*(log(A+u)-log(A+l))/(u-l)
}

wmbz_A_risk <- wmbz_A_bayes_risk_fun(u, sigma_squared_lower=l, A=A_prior)
sure_res <- sure_grandmean(sim_wmbz_a$Zs, sim_wmbz_a$sigmas)  
spher_shrink_res <- grouplinear(sim_wmbz_a$Zs, sim_wmbz_a$sigmas)  
mean_res <- location_mean(sim_wmbz_a$Zs)


expect_true(abs( mean( (sure_res - sim_wmbz_a$true_mus)^2) - wmbz_A_risk) < 0.01)
expect_true(abs( mean( (spher_shrink_res - sim_wmbz_a$true_mus)^2) - wmbz_A_risk) < 0.01)
expect_true(abs( mean( (mean_res - sim_wmbz_a$true_mus)^2) - (l+u)/2) < 0.01)


#gmleb_res <- gmleb(sim_wmbz_a$Zs, sim_wmbz_a$sigmas)  
#gmleb_res_dd <- gmleb(sim_wmbz_a$Zs, sigma=apply(sim_wmbz_a$Zs, 1 , sd)/sqrt(ncol(sim_wmbz_a$Zs)))  
#expect_true(abs( mean( (gmleb_res - sim_wmbz_a$true_mus)^2) - wmbz_A_risk) < 0.03)
#mean( (gmleb_res_dd - sim_wmbz_a$true_mus)^2)


# WMBZ C & F

expect_equal( sim_wmbz_c$true_mus, sim_wmbz_c$sigmas^2)
expect_equal( sim_wmbz_f$true_mus, sim_wmbz_f$sigmas^2)


expect_true(abs( mean( (location_mean(sim_wmbz_c$Zs) - sim_wmbz_c$true_mus)^2) - (l+u)/2) < 0.01)
expect_true(abs( mean( (location_mean(sim_wmbz_f$Zs) - sim_wmbz_f$true_mus)^2) - (l+u)/2) < 0.01)


cc_risk_wmbc_pred <- function(sigma_squared_upper, sigma_squared_lower=0.1){
  l <- sigma_squared_lower
  u <- sigma_squared_upper
  lambda <- (u+l)/( u +l + (u-l)^2/6)
  lambda^2*(u-l)^2/12 + (1-lambda)^2*(u+l)/2
}


sure_res_c <- sure_grandmean(sim_wmbz_c$Zs, sim_wmbz_c$sigmas)  
cc_res_c <- ebcc(sim_wmbz_c$Zs)
auroral_res_c <- auroral(sim_wmbz_c$Zs)

mse_sure_c <- mean( (sure_res_c - sim_wmbz_c$true_mus)^2)
mse_cc_c <- mean( (cc_res_c - sim_wmbz_c$true_mus)^2)
mse_auroral_c <- mean( (auroral_res_c - sim_wmbz_c$true_mus)^2)

expect_true(mse_sure_c < mse_cc_c)
expect_true(mse_auroral_c < mse_cc_c)

cc_risk_theory <- cc_risk_wmbc_pred(u, sigma_squared_lower=l)
expect_true(abs(cc_risk_theory - mse_cc_c) < 4*cc_risk_theory/sqrt(n))



sure_res_f <- sure_grandmean(sim_wmbz_f$Zs, sim_wmbz_f$sigmas)  
cc_res_f <- ebcc(sim_wmbz_f$Zs)
auroral_res_f <- auroral(sim_wmbz_f$Zs)

mse_sure_f <- mean( (sure_res_f - sim_wmbz_c$true_mus)^2)
mse_cc_f <- mean( (cc_res_f - sim_wmbz_c$true_mus)^2)
mse_auroral_f <- mean( (auroral_res_f - sim_wmbz_c$true_mus)^2)

expect_true(mse_sure_f < mse_cc_f)
expect_true(mse_auroral_f < mse_cc_f)

expect_true(abs(cc_risk_theory - mse_cc_f) < 4*cc_risk_theory/sqrt(n))


expect_true( mse_auroral_f < mse_auroral_c)

expect_true( abs(mse_sure_f - mse_sure_c) < 4*mse_sure_c/sqrt(n))
expect_true( abs(mse_cc_f - mse_cc_c) < 4*mse_cc_c/sqrt(n))
expect_true( abs(mse_cc_f - mse_cc_c) < 4*mse_cc_c/sqrt(n))

expect_true( max( sim_wmbz_f$Zs - sim_wmbz_f$true_mus ) < sqrt(u)*sqrt(3*B))
expect_true( min( sim_wmbz_f$Zs - sim_wmbz_f$true_mus ) >- sqrt(u)*sqrt(3*B))

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
