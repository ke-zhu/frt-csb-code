library(tidyverse)
library(parallel)
RNGkind("L'Ecuyer-CMRG")
set.seed(2024)

n_rep <- 500
n_cores <- 8


# 1 setup -----------------------------------------------

mag_bias <- 0:8

prop_unbias <- 0.5

model_spec <- tribble(
  ~sm, ~om,
  T, T,
  #T, F,
  #F, T,
  #F, F
)

ep <- "continuous"

bias_stru <- c("const")

sample_size <- tribble(
  ~n1, ~n0,  ~m,
  50, 25, 50,
  50, 25, 300
)

setup <- expand_grid(
  ep,
  model_spec,
  sample_size,
  bias_stru,
  mag_bias,
  prop_unbias
) %>% 
  mutate(
    prop_unbias = ifelse(mag_bias == 0, 0, prop_unbias)
  ) %>% 
  distinct() %>% 
  filter(!(bias_stru == "heter" & (mag_bias > 3 | mag_bias == 0 | m != 300))) %>% 
  mutate(case = row_number(), .before = everything()) %>% 
  group_by(ep, n1, n0, m) %>% 
  mutate(main_case = case[1]) %>% 
  ungroup

if (!dir.exists("data")) {dir.create("data")}

save(setup, file = "data/setup.RData")



# 2 generate simulation data -------------------------------------------------

for (case in 1:nrow(setup)) {
  # 2.1 load setup parameters
  list2env(setup[case,], envir = .GlobalEnv)
  N <- n1 + n0 + m
  sd_ec <- 0.5
  
  # 2.2 generate (X, S, A, Y0, Y1)
  # if mag_bias == 0, generate data; otherwise, skip
  if (mag_bias == 0) {
    simdata_main <- mclapply(1:n_rep, function(rep) {
      # X
      p <- 2
      X <- matrix(runif(N * p, -2, 2), N, p)
      Xt <- scale(X^3, center = F)
      
      # S
      if (sm) {
        X_pi <- X
      } else {
        X_pi <- Xt
      }
      eta <- rep(0.1, p)
      X_pi_eta <- X_pi %*% eta
      eta0 <- uniroot(function(.x) {
        mean(plogis(.x + X_pi_eta)) - (n1 + n0) / N
      }, interval = c(-100, 100)
      )$root
      pi_S <- plogis(eta0 + X_pi_eta)
      S <- rbinom(N, size = 1, prob = pi_S)
      #S <- rbinom(N, 1, (n1 + n0) / (n1 + n0 + m))
      
      # A
      A <- S
      A[S == 1] <- rbinom(sum(S == 1), 1, n1 / (n1 + n0))
      
      # Y0, Y1
      if (om) {
        X_mu <- X
      } else {
        X_mu <- Xt
      }
      beta0 <- rep(1, p)
      beta1 <- rep(2, p)
      if (ep == "continuous") {
        mu0 <- as.vector(X_mu %*% beta0)
        mu1 <- as.vector(X_mu %*% beta1) + 0.4
        Y0 <- mu0 + rnorm(N)
        Y1 <- mu1 + rnorm(N)
        R2 <- mean(c(var(mu0) / var(Y0), var(mu1) / var(Y1)))
      }
      tau <- mean((Y1 - Y0)[S == 1])
      
      # output
      lst(X, Xt, S, A, beta0, beta1, mu0, mu1, Y0, Y1, R2, tau)
    }, mc.cores = n_cores)
    # save data
    save(simdata_main, file = str_glue("data/simdata_main_{case}.RData"))
  } else {
    load(str_glue("data/simdata_main_{main_case}.RData"))
  }
  
  # 2.3 generate (Y00, Y, Ynull)
  simdata <- mclapply(1:n_rep, function(rep) {
    # load required quantities
    mu0 <- simdata_main[[rep]]$mu0
    S <- simdata_main[[rep]]$S
    A <- simdata_main[[rep]]$A
    Y0 <- simdata_main[[rep]]$Y0
    Y1 <- simdata_main[[rep]]$Y1
    X <- simdata_main[[rep]]$X
    beta0 <- simdata_main[[rep]]$beta0
    # Y00
    if (ep == "continuous") {
      if (mag_bias == 0) {
        id_unbias <- which(S == 0)
        mu00 <- mu0
      } else {
        id_unbias <- sample(
          which(S == 0), round(sum(S == 0) * prop_unbias)
        )
        if (bias_stru == "const") {
          mu00 <- mu0 - mag_bias
          mu00[id_unbias] <- mu0[id_unbias]
        } else if (bias_stru == "heter") {
          mu00 <- as.vector(X %*% (beta0 + mag_bias))
          mu00[id_unbias] <- mu0[id_unbias]
        }
      }
      Y00 <- mu00 + rnorm(N, sd = sd_ec)
    }
    # (Y, Ynull)
    Y <- S * A * Y1 + S * (1 - A) * Y0 + (1 - S) * Y00
    Ynull <- S * Y0 + (1 - S) * Y00
    # output
    lst(mu00, Y00, Y, Ynull, id_unbias)
  }, mc.cores = n_cores)
  save(simdata, file = str_glue("data/simdata_{case}.RData"))
}




