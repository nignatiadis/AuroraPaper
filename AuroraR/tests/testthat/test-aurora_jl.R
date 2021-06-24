test_that("multiplication works", {
  expect_equal(2 * 2, 4)
})


aurora_setup(JULIA_HOME="/Applications/Julia-1.5.app/Contents/Resources/julia/bin")

tmp_sim <- normal_normal_sim(10000, 10, sqrt(0.4), 0.5, prior_mu=0.5)

auroraknn_fit <- aurora_knn(tmp_sim$Zs, predictions_only=FALSE)

a <- apply(tmp_sim$Zs, 1, sd)/sqrt(1000)
mean(a) - sqrt(0.4)

