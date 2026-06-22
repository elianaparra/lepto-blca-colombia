// ============================================================
// blca_stan.stan  —  BLCA in Stan (robustness check)
// ============================================================
// 
install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))

file <- file.path(cmdstan_path(), "examples", "bernoulli", "bernoulli.stan")
mod  <- cmdstan_model(file)

data_list <- list(N = 10, y = c(1,0,1,1,0,1,0,1,1,0))
fit <- mod$sample(data = data_list, chains = 4, parallel_chains = 4)

install.packages("cmdstanr", repos = c("https://mc-stan.org/r-packages/", getOption("repos")))
library(cmdstanr)
check_cmdstan_toolchain(fix = TRUE)
print(fit$summary())
install_cmdstan()
cmdstanr::cmdstan_version()
