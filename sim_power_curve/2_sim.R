library(tidyverse)
library(parallel)
library(tictoc)
# library(SelectiveIntegrative) # for method = "AdaLasso Selective Borrow ACW"
# library(caret) # for method = "AdaLasso Selective Borrow ACW"
# library(nleqslv) # for method = "AdaLasso Selective Borrow ACW"
library(intFRT)

RNGkind("L'Ecuyer-CMRG")
set.seed(2024)

if (is.na(as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID")))) {
  test <- TRUE
  if (!exists("case")) {
    case <- 3
  }
  n_c <- 8
  n_rep <- 40
  n_f <- 2
  n_rep_g <- 2
} else {
  test <- FALSE
  case <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID"))
  n_c <- as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK"))
  n_rep <- 500
  n_f <- 5000
  n_rep_g <- 100
}

print(str_glue("test = {test}"))
print(str_glue("case = {case}"))
print(str_glue("n_cores = {n_c}"))
print(str_glue("n_rep = {n_rep}"))


# 1 setup -----------------------------------------------

load("data/setup.RData")
list2env(setup[case,], envir = .GlobalEnv)
if (ep == "continuous") {
  fam <- "gaussian"
} else if (ep == "binary") {
  fam <- "binomial"
}
load(str_glue("data/simdata_main_{main_case}.RData"))
load(str_glue("data/simdata_{case}.RData"))

# 2 simulation -------------------------------------------------

if (!dir.exists("output")) {dir.create("output")}
if (!dir.exists("output/cache")) {dir.create("output/cache")}
writeLines("", str_glue("output/cache/prog_{case}.txt"))

# i_rep <- 1


tic()
raw_result <- mclapply(1:n_rep, function(i_rep) {
  # save oracle information
  oracle_info <- tibble(
    tau = simdata_main[[i_rep]]$tau, 
    R2 = simdata_main[[i_rep]]$R2, 
    id_unbias = list(simdata[[i_rep]]$id_unbias)
  )
  # load data
  A <- simdata_main[[i_rep]]$A
  S <- simdata_main[[i_rep]]$S
  X <- simdata_main[[i_rep]]$X
  Y <- simdata[[i_rep]]$Y
  Ynull <- simdata[[i_rep]]$Ynull
  
  # save one sample data & output_frt
  if (i_rep == 1) {
    save(
      Y, A, S, X, Ynull, oracle_info,
      file = str_glue("output/cache/dt_{case}_{i_rep}.RData")
    )
    o_f <- T
  } else {
    o_f <- F
  }
  
  # adaptive selection threshold
  ada_gamma <- compute_ada_gamma(
    Y, A, S, X, fam,
    n_rep_gamma = n_rep_g
  )
  
  # inference
  res_alt <- res_null <- NULL
  
  if (ep == "continuous") {
    tryCatch({
      # result under alternative
      res_alt <- list(
        ec_borrow(Y, A, S, X, "No Borrow AIPW", fam, n_f, output_frt = o_f),
        ec_borrow(Y, A, S, X, "Borrow AIPW", fam, n_f, output_frt = o_f),
        ec_borrow(Y, A, S, X, "Conformal Selective Borrow AIPW", fam, n_f, 
                  gamma_sel = ada_gamma, output_frt = o_f)
      )
      # result under null
      res_null <- res_alt # these results are not needed; just for convenience
    }, error = function(e) {
      # output error
      error_msg <- conditionMessage(e)
      save(
        error_msg, Y, A, S, X, Ynull,
        file = str_glue("output/cache/err_{case}_{i_rep}.RData")
      )
    })
  }

  # output
  cat(i_rep, "\n", file = str_glue("output/cache/prog_{case}.txt"), 
      append = TRUE)
  lst(res_alt, res_null, oracle_info)
}, mc.cores = n_c)
total_time0 <- toc()
total_time <- unname(total_time0$toc - total_time0$tic)

save(raw_result, total_time, file = str_glue("output/cache/raw_{case}.RData"))

# organize raw result
inf_result <- imap_dfr(raw_result, function(raw_i, i) {
  if (!is.null(raw_i$res_alt)) {
    # each replication
    map2_dfr(raw_i$res_alt, raw_i$res_null, function(res_alt_j, res_null_j) {
      # each method
      id_ec_i <- res_alt_j$dat_info$id_ec[[1]]
      id_unbias_i <- raw_i$oracle_info$id_unbias[[1]]
      id_sel_j <- res_alt_j$out$id_sel[[1]]
      n_ec <- length(id_ec_i)
      n_unbias <- length(id_unbias_i)
      n_sel <- length(id_sel_j)
      n_TU <- intersect(id_sel_j, id_unbias_i) %>% length()
      n_FU <- n_sel - n_TU
      n_FB <- n_unbias - n_TU
      n_TB <- n_ec - n_TU - n_FU - n_FB
      #c(n_TU, n_FU, n_FB, n_TB) %>% matrix(2, 2, byrow = T)
      res_alt_j$res %>% 
        mutate(
          p_value_null = res_null_j$res$p_value,
          `FU/Sel` = ifelse(n_sel == 0, 0, n_FU / n_sel),
          `TU/Unbias` = n_TU / n_unbias,
          gamma_sel = res_alt_j$gamma_sel
        )
    }) %>%
      mutate(
        rep = i, 
        tau = raw_i$oracle_info$tau,
        R2 = raw_i$oracle_info$R2,
        .before = everything()
      )
  } else {
    NULL
  }
}) %>%
  mutate(case = case, .before = everything()) %>% 
  mutate(method = as_factor(method))

save(inf_result, file = str_glue("output/cache/inf_{case}.RData"))


# 3 summary ----------------------------------------------------------------

metrics <- inf_result %>% 
  mutate(
    R2 = mean(R2),
    tau = mean(tau)
  ) %>%
  group_by(
    case, 
    R2, 
    tau, 
    method
  ) %>%
  summarise(
    Bias = mean(est - tau),
    SD = sd(est),
    Var = var(est),
    MSE = mean((est - tau)^2),
    CP = mean(ci_l <= tau & tau <= ci_u),
    Width = mean(ci_u - ci_l),
    `Type I` = mean(p_value_null <= 0.05), 
    Power = mean(p_value <= 0.05),
    n_sel = mean(n_sel),
    ess_sel = round(mean(ess_sel)),
    `FU/Sel` = ifelse(is.na(n_sel), NA, mean(`FU/Sel`)),
    `TU/Unbias` = ifelse(is.na(n_sel), NA, mean(`TU/Unbias`)),
    runtime = mean(runtime),
    gamma_sel = mean(gamma_sel),
    .groups = "drop"
  ) %>% 
  mutate(`Bias/SD%` = round(abs(Bias / SD) * 100), .after = `Bias`) %>%
  select(-SD) %>% 
  mutate(`Var%` = round((Var / Var[1]) * 100), .after = Var) %>%
  mutate(`MSE%` = round((MSE / MSE[1]) * 100), .after = MSE) %>%
  mutate(`Width%` = round((Width / Width[1]) * 100), .after = Width) %>%
  mutate(`Pow%` = round((Power / Power[1]) * 100), .after = Power)

save(metrics, file = str_glue("output/metrics_{case}.RData"))



# estimated time

# total_time * n_c / n_rep / (8 * (n_f + 1) + (n_rep_g + 1) * 10) /
#   (80 / 200 / (8 * (5000 + 1) + (100 + 1) * 10)) /
#   3600 / 24






