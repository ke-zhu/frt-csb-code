library(tidyverse)
library(scales)
library(cowplot)
library(latex2exp)
theme_set(theme_bw())
set.seed(2024)


case_seq <- 1:12


method_all <- c("NB", "FB", "ALSB", "CSB0.8", "CSB0.6", "CSB0.4", "CSBada")
prop_met_id <- 7

lab_all <- c(
  TeX("No Borrow ($\\gamma=1$)"),
  TeX("Full Borrow ($\\gamma=0$)"),
  TeX("Adaptive Lasso Selective Borrow"),
  TeX("Conformal Selective Borrow ($\\gamma=0.8$)"),
  TeX("Conformal Selective Borrow ($\\gamma=0.6$)"),
  TeX("Conformal Selective Borrow ($\\gamma=0.4$)"),
  TeX("Conformal Selective Borrow ($\\hat{\\gamma}$)")
)

# show_col(hue_pal()(18))
color_all <- hue_pal()(18)[c(7, 1, 2, 3, 14, 18, 12)]
#show_col(color_all)

load("data/setup.RData")

# extract all results -----------------------------------------------------

res <- map_dfr(case_seq, function(case) {
  if (case %in% c(1, 7)) {
    load(str_glue("../sim_main/output/metrics_{case}.RData"))
    load(str_glue("../sim_main/output/cache/inf_{case}.RData"))
  } else {
    load(str_glue("output/metrics_{case}.RData"))
    load(str_glue("output/cache/inf_{case}.RData"))
  }
  
  n_rep <- inf_result$rep %>% max
  
  met_frt <- metrics %>% 
    filter(str_detect(method, "FRT")) %>% 
    mutate(method = str_remove(method, "\\+FRT")) %>% 
    select(method, `Type I`, Power)
  
  out <- metrics %>%
    filter(
      !str_detect(method, "FRT"),
      method %in% c(
        "No Borrow AIPW",
        "Borrow AIPW",
        "Conformal Selective Borrow AIPW",
        "Conformal Selective Borrow AIPW (g=ada)"
      )
      #method != "AdaLasso Selective Borrow ACW"
    ) %>% 
    left_join(met_frt, by = "method") %>% 
    mutate(
      Method = case_when(
        method == "No Borrow AIPW" ~ "NB",
        method == "Borrow AIPW" ~ "FB",
        # method == "AdaLasso Selective Borrow ACW" ~ "ALSB",
        # method == "Conformal Selective Borrow AIPW (g=0.8)" ~ "CSB0.8",
        # method == "Conformal Selective Borrow AIPW (g=0.6)" ~ "CSB0.6",
        # method == "Conformal Selective Borrow AIPW (g=0.4)" ~ "CSB0.4",
        method == "Conformal Selective Borrow AIPW" |
          method == "Conformal Selective Borrow AIPW (g=ada)" ~ "CSBada"
      )
    ) %>% 
    mutate(
      b = setup$mag_bias[case], 
      tau = setup$tau[case], 
      .before = everything()
    ) %>% 
    mutate(
      `Type I.y` = ifelse(Method == "ALSB", `Type I.x`, `Type I.y`),
      `Power.y` = ifelse(Method == "ALSB", `Power.x`, `Power.y`),
      `Type I.y.low` = map_dbl(`Type I.y`, ~{
        binom.test(. * n_rep, n_rep, conf.level = 0.95)$conf.int[1]
      }),
      `Type I.y.up` = map_dbl(`Type I.y`, ~{
        binom.test(. * n_rep, n_rep, conf.level = 0.95)$conf.int[2]
      })
    )
  
  if (case %in% c(1, 7)) {
    out$Power.y = out$`Type I.y`
    out
  } else {
    out
  }
}) %>% 
  mutate(
    Method = Method %>% 
      as_factor %>% 
      fct_relevel(
        "NB", 
        "FB", 
        #"ALSB", "CSB0.8", "CSB0.6", "CSB0.4", 
        "CSBada"
      ),
    proposed = ifelse(str_detect(Method, "CSB"), "y", "n") %>% 
      as_factor %>% 
      fct_relevel("y", "n")
  )

met_id <- c(1, 2, 7)

if (length(met_id) <= 3) {
  lrow = 1
} else {
  lrow = 2
}

legend_lab <- scale_color_manual(
  values = color_all[met_id],
  labels = lab_all[met_id]
) 

my_theme <- theme(legend.position = "bottom", 
                  #axis.title.y = element_blank(),
                  plot.title = element_text(hjust = 0.5))

my_line <- geom_line(aes(group = Method, linetype = proposed),show.legend = FALSE)

my_x <- "Magnitude of Hidden Bias" # Magnitude of Hidden Bias

tImax <- max(
  res %>% 
    filter(Method %in% method_all[met_id]) %>% 
    pull(`Type I.y.up`),
  0.5
)

yp <- range(
  res %>% 
    filter(Method %in% method_all[met_id]) %>% 
    pull(Power.y), na.rm=T
)+c(-0.1,0.1)
scale_power <- ylim(yp[1], yp[2])

res %>% 
  filter(Method %in% method_all[met_id]) %>%
  mutate(b = ifelse(b == 0, "No Hidden Bias", "Half of ECs Exhibit Hidden Bias") %>% as_factor) %>% 
  ggplot(aes(tau, Power.y, color = Method)) +
  my_line +
  geom_point() +
  my_theme +
  legend_lab +
  scale_power +
  facet_grid(.~b) +
  xlab(TeX("$\\tau$")) +
  ylab(TeX("Power"))

ggsave(str_glue("chart/sim_power.pdf"), width = 6, height = 4)


lab_all <- c(
  TeX("No Borrow AIPW ($\\gamma=1$)"),
  TeX("Borrow AIPW ($\\gamma=0$)"),
  TeX("Adaptive Lasso Selective Borrow"),
  TeX("Conformal Selective Borrow ($\\gamma=0.8$)"),
  TeX("Conformal Selective Borrow ($\\gamma=0.6$)"),
  TeX("Conformal Selective Borrow ($\\gamma=0.4$)"),
  TeX("CSB AIPW ($\\hat{\\gamma}$)")
)

legend_lab <- scale_color_manual(
  values = color_all[met_id],
  labels = lab_all[met_id]
) 

res %>% 
  filter(Method %in% method_all[met_id]) %>%
  mutate(b = ifelse(b == 0, "No Hidden Bias (b = 0)", "With Hidden Bias (b = 8)") %>% as_factor) %>% 
  ggplot(aes(tau, Power.y, color = Method)) +
  my_line +
  geom_point() +
  my_theme +
  legend_lab +
  scale_power +
  facet_grid(.~b) +
  xlab(TeX("$\\tau$")) +
  ylab(TeX("Power"))+
  labs(color="") +
  guides(color = guide_legend(nrow = 1))

ggsave(str_glue("chart/sim_power2.pdf"), width = 5, height = 5)

ggsave(str_glue("chart/sim_power3.pdf"), width = 6, height = 4)




lab_all <- c(
  TeX("No Borrow AIPW"),
  TeX("Borrow AIPW"),
  TeX("Adaptive Lasso Selective Borrow"),
  TeX("Conformal Selective Borrow ($\\gamma=0.8$)"),
  TeX("Conformal Selective Borrow ($\\gamma=0.6$)"),
  TeX("Conformal Selective Borrow ($\\gamma=0.4$)"),
  TeX("Selective Borrow AIPW")
)

legend_lab <- scale_color_manual(
  values = color_all[met_id],
  labels = lab_all[met_id]
) 

res %>% 
  filter(Method %in% method_all[met_id]) %>%
  mutate(b = ifelse(b == 0, "No Outcome Incomparability", "Half With Outcome Incomparability") %>% as_factor) %>% 
  ggplot(aes(tau, Power.y, color = Method)) +
  my_line +
  geom_point() +
  my_theme +
  legend_lab +
  scale_power +
  facet_grid(.~b) +
  xlab(TeX("$\\tau$")) +
  ylab(TeX("Power of FRT"))+
  labs(color="") +
  guides(color = guide_legend(nrow = 1))

ggsave(str_glue("chart/sim_power_grant.pdf"), width = 5, height = 3)




