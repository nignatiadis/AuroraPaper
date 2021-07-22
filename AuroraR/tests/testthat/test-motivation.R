# test code for motivation.Rmd

# three component prior
set.seed(1)
n <- 500000
sigma_3c <- 1
effect_sizes_3c <- 3*c(-1:1)
three_component_two_reps <- three_component_normal_sim(n, 2,
                                                       sigma=sigma_3c/sqrt(2), 
                                                       A = 9/1.5)


expect_true( abs(mean((three_component_two_reps$Zs[,1] - three_component_two_reps$true_mus)^2) - 1) < 0.01)
expect_true( abs(mean((three_component_two_reps$Zs[,2] - three_component_two_reps$true_mus)^2) - 1) < 0.01)


bayes_three_component <- function(zs,  sigma=1, effect_sizes = 3*(-1:1)){
  num <- sapply(zs, function(z) sum( effect_sizes*dnorm(z, effect_sizes, sigma)))
  denom <- sapply(zs, function(z) sum(dnorm(z, effect_sizes, sigma)))
  num/denom
}

three_component_df_two_reps <- with(three_component_two_reps,
                                    tibble(xs=Zs[,1],
                                           ys=Zs[,2], 
                                           mus=true_mus,
                                           bayes_rule = bayes_three_component(
                                             Zs[,1], 
                                             sigma=sigma_3c, 
                                             effect_sizes = effect_sizes_3c), 
                                           prior= "Three-point prior"))



expect_true( abs(mean((three_component_df_two_reps$xs - three_component_df_two_reps$mus)^2) - 1) < 0.01)
expect_true( abs(mean((three_component_df_two_reps$ys - three_component_df_two_reps$mus)^2) - 1) < 0.01)


fisher_info_integrand <- function(z, effect_sizes=effect_sizes_3c){
  num_sqrt <- sum( (z-effect_sizes)*dnorm(z, mean=effect_sizes))/3
  denom <- sum(dnorm(z, mean=effect_sizes))/3
  num_sqrt^2/denom
}

fisher_info <-stats::integrate(Vectorize(fisher_info_integrand), -30, +30)

expect_true( abs(with(three_component_df_two_reps, mean( (bayes_rule - mus)^2)) - 1 + fisher_info$value) < 0.005)



# test normal prior
set.seed(2)
A_nn <- 4
sigma_nn <- 1
prior_mu <- 0.5
normal_normal_two_reps <- normal_normal_sim(m=n, B=2,
                                            sigma=sigma_nn/sqrt(2),
                                            A=A_nn,
                                            prior_mu=0.5)

bayes_normal_normal <- function(zs, A=1, sigma=1, prior_mu=prior_mu){
  lambda <- A/(A+sigma^2)
  zs*lambda + prior_mu*(1-lambda)
}

normal_normal_two_reps_df <- with(normal_normal_two_reps,
                                  tibble(xs=Zs[,1],
                                         ys=Zs[,2], 
                                         mus=true_mus,
                                         bayes_rule=bayes_normal_normal(Zs[,1],
                                                                        A=A_nn,
                                                                        sigma=sigma_nn,
                                                                        prior_mu=prior_mu),
                                         prior="Normal prior"))





expect_true( abs(mean((normal_normal_two_reps_df$xs - normal_normal_two_reps_df$mus)^2) - 1) < 0.01)
expect_true( abs(mean((normal_normal_two_reps_df$ys - normal_normal_two_reps_df$mus)^2) - 1) < 0.01)

expect_true( abs(mean(normal_normal_two_reps_df$mus) - prior_mu) < 5*sqrt(A_nn)/sqrt(n))

true_bayes_risk <- A_nn/(1+A_nn)

expect_true( abs(with(normal_normal_two_reps_df, mean( (bayes_rule - mus)^2)) - true_bayes_risk) < 0.005)
