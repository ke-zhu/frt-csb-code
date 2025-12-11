library(tidyverse)
library(tictoc)
library(parallel)
library(SelectiveIntegrative)
library(caret)
library(quantreg)
library(intFRT)
RNGkind("L'Ecuyer-CMRG")
set.seed(2025)

n_c <- 8
n_f <- 10000
n_rep_gamma <- 100

# data --------------------------------------------------------------------

load("data/C9633+NCDB.RData")

Y <- dat$y
A <- dat$treat
S <- dat$sample
X <- dat %>% 
  select(sex, age, race, hist, tsize) %>% 
  as.matrix()
fam <- "gaussian"

# analysis -------------------------------------------------------------

tic()
ada_gamma <- compute_ada_gamma(
  Y, A, S, X, fam,
  n_rep_gamma = n_rep_gamma,
  parallel = T,
  cf_score = "CQR",
  cf = "jackknife+"
)
toc()
# 82.482 sec elapsed
cat(ada_gamma)


tic()
result <- list(
  ec_borrow(Y, A, S, X, "No Borrow DiM", fam, n_f, parallel = T),
  ec_borrow(Y, A, S, X, "No Borrow AIPW", fam, n_f, parallel = T),
  ec_borrow(Y, A, S, X, "Borrow AIPW", fam, n_f, parallel = T),
  ec_borrow(Y, A, S, X, "AdaLasso Selective Borrow ACW", fam, n_f, parallel = T),
  ec_borrow(Y, A, S, X, "Conformal Selective Borrow AIPW", fam, n_f, 
            gamma_sel = 0.6, cf_score = "CQR", cf = "jackknife+",
            parallel = T) %>% 
    intFRT:::add_name("(g=0.6)"),
  ec_borrow(Y, A, S, X, "Conformal Selective Borrow AIPW", fam, n_f, 
            gamma_sel = ada_gamma, cf_score = "CQR", cf = "jackknife+",
            parallel = T) %>% 
    intFRT:::add_name("(g=ada)")
)
total_time0 <- toc()
total_time <- unname(total_time0$toc - total_time0$tic)

#1601.595 sec elapsed

#map_dbl(result, "gamma_sel")

tab <- map_dfr(result, 1) %>%
  mutate(
    across(c(est, se, ci_l, ci_u, p_value), ~ round(.,3)),
    across(c(n_sel, ess_sel), round)
  ) %>% select(-ess_sel)


# save results
if (!dir.exists("output")) {dir.create("output")}
save(result, tab, total_time, file = "output/result.RData")


