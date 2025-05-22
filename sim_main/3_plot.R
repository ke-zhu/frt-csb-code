library(tidyverse)
library(scales)
library(cowplot)
library(latex2exp)
theme_set(theme_bw())
set.seed(2024)


case_pre <- "const_ne50"
case_seq <- 1:9

# case_pre <- "const_ne300"
# case_seq <- 10:18

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
  load(str_glue("output/metrics_{case}.RData"))
  load(str_glue("output/cache/inf_{case}.RData"))
  n_rep <- inf_result$rep %>% max
  
  met_frt <- metrics %>% 
    filter(str_detect(method, "FRT")) %>% 
    mutate(method = str_remove(method, "\\+FRT")) %>% 
    select(method, `Type I`, Power)
  
  metrics %>%
    filter(
      !str_detect(method, "FRT"),
      #method != "AdaLasso Selective Borrow ACW"
    ) %>% 
    left_join(met_frt, by = "method") %>% 
    mutate(
      Method = case_when(
        method == "No Borrow AIPW" ~ "NB",
        method == "Borrow AIPW" ~ "FB",
        method == "AdaLasso Selective Borrow ACW" ~ "ALSB",
        method == "Conformal Selective Borrow AIPW (g=0.8)" ~ "CSB0.8",
        method == "Conformal Selective Borrow AIPW (g=0.6)" ~ "CSB0.6",
        method == "Conformal Selective Borrow AIPW (g=0.4)" ~ "CSB0.4",
        method == "Conformal Selective Borrow AIPW (g=ada)" ~ "CSBada"
      )
    ) %>% 
    mutate(b = setup$mag_bias[case], .before = everything()) %>% 
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
}) %>% 
  mutate(
    Method = Method %>% 
      as_factor %>% 
      fct_relevel(
        "NB", 
        "FB", "ALSB", "CSB0.8", "CSB0.6", "CSB0.4", 
        "CSBada"
      ),
    proposed = ifelse(str_detect(Method, "CSB"), "y", "n") %>% 
      as_factor %>% 
      fct_relevel("y", "n")
  )


# fig:sim_pdist -----------------------------------------------------------

met_id <- c(1, 2, prop_met_id)

case_show <- range(case_seq)

p_all0 <- map(case_show, function(case) {
  load(str_glue("output/cache/inf_{case}.RData"))
  load(str_glue("output/metrics_{case}.RData"))
  
  inf_tab <- inf_result %>% 
    mutate(Method = case_when(
      method == "No Borrow AIPW+FRT" ~ "NB",
      method == "Borrow AIPW+FRT" ~ "FB",
      method == "AdaLasso Selective Borrow ACW+FRT" ~ "ALSB",
      method == "Conformal Selective Borrow AIPW+FRT (g=0.8)" ~ "CSB0.8",
      method == "Conformal Selective Borrow AIPW+FRT (g=0.6)" ~ "CSB0.6",
      method == "Conformal Selective Borrow AIPW+FRT (g=0.4)" ~ "CSB0.4",
      method == "Conformal Selective Borrow AIPW+FRT (g=ada)" ~ "CSBada"
    ) %>% as_factor %>% fct_relevel("NB", "FB")) %>% 
    filter(Method %in% method_all[met_id]) %>% 
    select(Method, p_value, p_value_null)
  
  sum_tab <- inf_tab %>% 
    group_by(Method) %>% 
    summarise(
      `Type I` = mean(p_value_null <= 0.05),
      Power = mean(p_value <= 0.05)
    )
  
  if (case == case_show[1]) {
    t1 <- ggtitle(TeX("(A) Distribution under $H_0$ (No Hidden Bias in ECs)"))
    t2 <- ggtitle(TeX("(B) Distribution under $H_1$ (No Hidden Bias in ECs)"))
  } else {
    t1 <- ggtitle(TeX("(C) Distribution under $H_0$ (Hidden Bias in ECs)"))
    t2 <- ggtitle(TeX("(D) Distribution under $H_1$ (Hidden Bias in ECs)"))
  }
  
  p1 <- inf_tab %>% 
    ggplot() +
    geom_histogram(aes(p_value_null, fill = Method, color = Method), 
                   alpha = 0.5) +
    geom_vline(xintercept = 0.05, linetype = 2) +
    geom_label(data = sum_tab, ##
               aes(x = Inf, y = Inf, 
                   label = sprintf("Type I error = %.3f", `Type I`), ##
                   color = Method,
                   vjust = 1.5
               ),
               hjust = 1.1, size = 3.5, show.legend = F) +
    scale_fill_manual(values = color_all[met_id], labels = lab_all[met_id]) +
    scale_color_manual(values = color_all[met_id], labels = lab_all[met_id]) +
    facet_grid(Method~.)+
    ylab(NULL) +
    xlab("Exact p-value") +
    theme(strip.text = element_blank()) +
    theme(legend.position = "bottom") +
    t1
  
  p2 <- inf_tab %>% 
    ggplot() +
    geom_histogram(aes(p_value, fill = Method, color = Method), alpha = 0.5) + ##
    geom_vline(xintercept = 0.05, linetype = 2) +
    geom_label(data = sum_tab, ##
               aes(x = Inf, y = Inf, 
                   label = sprintf("Power = %.3f", Power), ##
                   color = Method,
                   vjust = 1.5
               ),
               hjust = 1.1, size = 3.5, show.legend = F) +
    scale_fill_manual(values = color_all[met_id], labels = lab_all[met_id]) +
    scale_color_manual(values = color_all[met_id], labels = lab_all[met_id]) +
    facet_grid(Method~.)+
    ylab(NULL) +
    xlab("Exact p-value") +
    theme(strip.text = element_blank()) +
    theme(legend.position = "bottom") +
    t2
  lst(p1, p2)
})

p_all <- list(p_all0[[1]][[1]], p_all0[[1]][[2]], 
              p_all0[[2]][[1]], p_all0[[2]][[2]])

legend <- get_plot_component(
  p_all[[1]], "guide-box", return_all = T
)[[3]]

combined_plot <- do.call(
  plot_grid, 
  map(p_all, ~{. + theme(legend.position="none")})
)

plot_grid(combined_plot, legend, ncol = 1, rel_heights = c(1, 0.1))
ggsave(str_glue("chart/sim_{case_pre}_pdist.pdf"), width = 9, height = 6)


# fig:sim_pdist2 -----------------------------------------------------------

met_id <- c(1, 2)

case_show <- range(case_seq)

p_all0 <- map(case_show, function(case) {
  load(str_glue("output/cache/inf_{case}.RData"))
  load(str_glue("output/metrics_{case}.RData"))
  
  inf_tab <- inf_result %>% 
    mutate(Method = case_when(
      method == "No Borrow AIPW+FRT" ~ "NB",
      method == "Borrow AIPW+FRT" ~ "FB",
      method == "AdaLasso Selective Borrow ACW+FRT" ~ "ALSB",
      method == "Conformal Selective Borrow AIPW+FRT (g=0.8)" ~ "CSB0.8",
      method == "Conformal Selective Borrow AIPW+FRT (g=0.6)" ~ "CSB0.6",
      method == "Conformal Selective Borrow AIPW+FRT (g=0.4)" ~ "CSB0.4",
      method == "Conformal Selective Borrow AIPW+FRT (g=ada)" ~ "CSBada"
    ) %>% as_factor %>% fct_relevel("NB", "FB")) %>% 
    filter(Method %in% method_all[met_id]) %>% 
    select(Method, p_value, p_value_null)
  
  sum_tab <- inf_tab %>% 
    group_by(Method) %>% 
    summarise(
      `Type I` = mean(p_value_null <= 0.05),
      Power = mean(p_value <= 0.05)
    )
  
  if (case == case_show[1]) {
    t1 <- ggtitle(TeX("(A) $H_0$, unbiased ECs"))
    t2 <- ggtitle(TeX("(B) $H_1$, unbiased ECs"))
  } else {
    t1 <- ggtitle(TeX("(C) $H_0$, biased ECs"))
    t2 <- ggtitle(TeX("(D) $H_1$, biased ECs"))
  }
  
  p1 <- inf_tab %>% 
    ggplot() +
    geom_histogram(aes(p_value_null, fill = Method, color = Method), 
                   alpha = 0.5) +
    geom_vline(xintercept = 0.05, linetype = 2) +
    geom_label(data = sum_tab, ##
               aes(x = Inf, y = Inf, 
                   label = sprintf("Type I error = %.3f", `Type I`), ##
                   color = Method,
                   vjust = 1.5
               ),
               hjust = 1.1, size = 3.5, show.legend = F) +
    scale_fill_manual(values = color_all[met_id], labels = lab_all[met_id]) +
    scale_color_manual(values = color_all[met_id], labels = lab_all[met_id]) +
    facet_grid(Method~.)+
    ylab(NULL) +
    xlab("Exact p-value") +
    theme(strip.text = element_blank()) +
    theme(legend.position = "bottom") +
    t1
  
  p2 <- inf_tab %>% 
    ggplot() +
    geom_histogram(aes(p_value, fill = Method, color = Method), alpha = 0.5) + ##
    geom_vline(xintercept = 0.05, linetype = 2) +
    geom_label(data = sum_tab, ##
               aes(x = Inf, y = Inf, 
                   label = sprintf("Power = %.3f", Power), ##
                   color = Method,
                   vjust = 1.5
               ),
               hjust = 1.1, size = 3.5, show.legend = F) +
    scale_fill_manual(values = color_all[met_id], labels = lab_all[met_id]) +
    scale_color_manual(values = color_all[met_id], labels = lab_all[met_id]) +
    facet_grid(Method~.)+
    ylab(NULL) +
    xlab("Exact p-value") +
    theme(strip.text = element_blank()) +
    theme(legend.position = "bottom") +
    t2
  lst(p1, p2)
})

p_all <- list(p_all0[[1]][[1]], p_all0[[1]][[2]], 
              p_all0[[2]][[1]], p_all0[[2]][[2]])

legend <- get_plot_component(
  p_all[[1]], "guide-box", return_all = T
)[[3]]

combined_plot <- do.call(
  plot_grid, 
  c(
    map(p_all, ~{. + theme(legend.position="none")}),
    ncol = 4
  )
)

plot_grid(combined_plot, legend, ncol = 1, rel_heights = c(1, 0.1))
ggsave(str_glue("chart/sim_{case_pre}_pdist2.pdf"), width = 10, height = 3)



# fig:sim_main ------------------------------------------------------------

#met_id <- c(1,2,3,5)

# custom_trans <- function(x) {
#   ifelse(x > 0.3 & x < 0.7, 0.3 + 4 * (x - 0.3), x)
# }
# 
# inverse_custom_trans <- function(x) {
#   ifelse(x > 0.3 & x < 2.1, 0.3 + (x - 0.3) / 4, x)
# }
# scale_power = scale_y_continuous(trans = scales::trans_new(
#   name = "custom",
#   transform = custom_trans,
#   inverse = inverse_custom_trans
# ))

plot_sim <- function(met_id) {
  
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
                    axis.title.y = element_blank(),
                    plot.title = element_text(hjust = 0.5))
  
  my_line <- geom_line(aes(group = Method, linetype = proposed),show.legend = FALSE)
  
  my_x <- "Magnitude of Hidden Bias" # Magnitude of Hidden Bias
  
  tImax <- max(
    res %>% 
      filter(Method %in% method_all[met_id]) %>% 
      pull(`Type I.y.up`),
    0.6
  )
  
  yp <- range(
    res %>% 
      filter(Method %in% method_all[met_id]) %>% 
      pull(Power.y), na.rm=T
  )+c(-0.1,0.1)
  scale_power <- ylim(yp[1], yp[2])
  
  # my_line <- geom_smooth(
  #   aes(group = Method), method = "lm", 
  #   formula = y ~ poly(x, 6), se = F, linewidth = 0.5
  # )
  
  p1 <- res %>% 
    filter(Method %in% method_all[met_id]) %>%
    mutate(`Absolute Bias` = abs(Bias)) %>% 
    ggplot(aes(b, `Absolute Bias`, color = Method)) +
    my_line +
    geom_point() +
    my_theme +
    labs(x = my_x,
         title = "(A) Absolute Bias") +
    legend_lab
  
  p2 <- res %>% 
    filter(Method %in% method_all[met_id]) %>%
    ggplot(aes(b, Var, color = Method)) +
    my_line +
    geom_point() +
    my_theme +
    labs(x = my_x,
         title = "(B) Variance") +
    legend_lab
  
  p3 <- res %>% 
    filter(Method %in% method_all[met_id]) %>%
    ggplot(aes(b, MSE, color = Method)) +
    my_line +
    geom_point() +
    my_theme +
    labs(x = my_x,
         title = "(C) MSE") +
    legend_lab +
    scale_y_log10()
  
  p4 <- res %>% 
    filter(Method %in% method_all[met_id]) %>%
    ggplot(aes(b, `Type I.y`, color = Method)) +
    geom_hline(yintercept = 0.05) +
    #my_line +
    geom_point(position = position_dodge(width = 0.5)) +
    geom_errorbar(aes(ymin = `Type I.y.low`, ymax = `Type I.y.up`), 
                  width = 0.4, position = position_dodge(width = 0.5)) +
    my_theme +
    ylim(0, tImax) +
    labs(x = my_x, title = "(D) Type I Error Rate") +
    legend_lab
  
  p5 <- res %>% 
    filter(Method %in% method_all[met_id]) %>%
    ggplot(aes(b, Power.y, color = Method)) +
    my_line +
    geom_point() +
    my_theme +
    labs(x = my_x, title = "(E) Power") +
    legend_lab +
    scale_power
  
  p6 <- res %>% 
    filter(Method %in% method_all[met_id]) %>%
    ggplot(aes(b, `FU/Sel`, color = Method)) +
    my_line +
    geom_point() +
    my_theme +
    labs(x = my_x, title = "(F) #Biased / #Selected") +
    legend_lab
  
  p_all <- map(list(p1, p2, p3, p4, p5, p6), ~ {
    .x + theme(legend.position="none")
  })
  p_all$ncol <- 3
  legend <- get_plot_component(
    p1 + guides(color = guide_legend(nrow = lrow)), "guide-box", return_all = T
  )[[3]]
  combined_plot <- do.call(plot_grid, p_all)
  
  plot_grid(combined_plot, legend, ncol = 1, rel_heights = c(1, 0.1))
}

plot_sim(c(1, 2, prop_met_id))
ggsave(str_glue("chart/sim_{case_pre}_main.pdf"), width = 9, height = 6)

# fig:sim_alasso ------------------------------------------------------------

plot_sim(c(1, 3, prop_met_id))
ggsave(str_glue("chart/sim_{case_pre}_alasso.pdf"), width = 9, height = 6)


# fig:sim_difg -----------------------------------------------------------

plot_sim(
  c(1, 4, 5, 6, 7)
)
ggsave(str_glue("chart/sim_{case_pre}_difg.pdf"), width = 9, height = 6)



# fig:sim_main2 ------------------------------------------------------------

#met_id <- c(1,2,3,5)

# custom_trans <- function(x) {
#   ifelse(x > 0.3 & x < 0.7, 0.3 + 4 * (x - 0.3), x)
# }
# 
# inverse_custom_trans <- function(x) {
#   ifelse(x > 0.3 & x < 2.1, 0.3 + (x - 0.3) / 4, x)
# }
# scale_power = scale_y_continuous(trans = scales::trans_new(
#   name = "custom",
#   transform = custom_trans,
#   inverse = inverse_custom_trans
# ))

plot_sim <- function(met_id, leg_h = 0.1) {
  
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
                    axis.title.y = element_blank(),
                    plot.title = element_text(hjust = 0.5))
  
  my_line <- geom_line(aes(group = Method, linetype = proposed),show.legend = FALSE)
  
  my_x <- "Magnitude of Hidden Bias" # Magnitude of Hidden Bias
  
  tImax <- max(
    res %>% 
      filter(Method %in% method_all[met_id]) %>% 
      pull(`Type I.y.up`),
    0.2
  )
  
  yp <- range(
    res %>% 
      filter(Method %in% method_all[met_id]) %>% 
      pull(Power.y), na.rm=T
  )+c(-0.1,0.1)
  scale_power <- ylim(yp[1], yp[2])
  
  # my_line <- geom_smooth(
  #   aes(group = Method), method = "lm", 
  #   formula = y ~ poly(x, 6), se = F, linewidth = 0.5
  # )
  
  p1 <- res %>% 
    filter(Method %in% method_all[met_id]) %>%
    mutate(`Absolute Bias` = abs(Bias)) %>% 
    ggplot(aes(b, `Absolute Bias`, color = Method)) +
    my_line +
    geom_point() +
    my_theme +
    labs(x = my_x,
         title = "(A) Absolute Bias") +
    legend_lab
  
  p2 <- res %>% 
    filter(Method %in% method_all[met_id]) %>%
    ggplot(aes(b, Var, color = Method)) +
    my_line +
    geom_point() +
    my_theme +
    labs(x = my_x,
         title = "(B) Variance") +
    legend_lab
  
  p3 <- res %>% 
    filter(Method %in% method_all[met_id]) %>%
    ggplot(aes(b, MSE, color = Method)) +
    my_line +
    geom_point() +
    my_theme +
    labs(x = my_x,
         title = "(C) MSE") +
    legend_lab +
    scale_y_log10()
  
  p4 <- res %>% 
    filter(Method %in% method_all[met_id]) %>%
    ggplot(aes(b, `Type I.y`, color = Method)) +
    geom_hline(yintercept = 0.05) +
    #my_line +
    geom_point(position = position_dodge(width = 0.5)) +
    geom_errorbar(aes(ymin = `Type I.y.low`, ymax = `Type I.y.up`), 
                  width = 0.4, position = position_dodge(width = 0.5)) +
    my_theme +
    ylim(0, tImax) +
    labs(x = my_x, title = "(D) Type I Error Rate") +
    legend_lab
  
  p5 <- res %>% 
    filter(Method %in% method_all[met_id]) %>%
    ggplot(aes(b, Power.y, color = Method)) +
    my_line +
    geom_point() +
    my_theme +
    labs(x = my_x, title = "(E) Power") +
    legend_lab +
    scale_power
  
  # p6 <- res %>% 
  #   filter(Method %in% method_all[met_id]) %>%
  #   ggplot(aes(b, `FU/Sel`, color = Method)) +
  #   my_line +
  #   geom_point() +
  #   my_theme +
  #   labs(x = my_x, title = "(F) #Biased / #Selected") +
  #   legend_lab
  
  p_all <- map(list(p1, p2, p3, p4, p5), ~ {
    .x + theme(legend.position="none")
  })
  p_all$ncol <- 5
  legend <- get_plot_component(
    p1 + guides(color = guide_legend(nrow = lrow)), "guide-box", return_all = T
  )[[3]]
  combined_plot <- do.call(plot_grid, p_all)
  
  plot_grid(combined_plot, legend, ncol = 1, rel_heights = c(1, leg_h))
}

plot_sim(c(1, 2, prop_met_id))
ggsave(str_glue("chart/sim_{case_pre}_main2.pdf"), width = 12, height = 3)






# fig:sim_alasso2 ------------------------------------------------------------

plot_sim(c(1, 3, prop_met_id))
ggsave(str_glue("chart/sim_{case_pre}_alasso2.pdf"), width = 12, height = 3)


# fig:sim_difg2 -----------------------------------------------------------

plot_sim(
  c(1, 4, 5, 6, 7), leg_h = 0.3
)
ggsave(str_glue("chart/sim_{case_pre}_difg2.pdf"), width = 12, height = 3)


# fig:sim_bias ----------------------------------------------------------------

dat_fig <- map_dfr(case_seq, function(case) {
  load("data/setup.RData")
  load(str_glue("output/cache/dt_{case}_1.RData"))
  load(str_glue("output/cache/raw_{case}.RData"))
  dat <- tibble(Y, A, S, X)
  id_sel <- raw_result[[1]]$res_alt[[prop_met_id]]$out$id_sel[[1]]
  fit <- glm(S ~ X, data = dat)
  dat %>%
    mutate(
      #y = ifelse(y == 3, y + runif(n(), 0, 0.1), y),
      `Sampling Score` = predict(fit, type = "response"),
      sel = ifelse(row_number() %in% id_sel, 1, 0),
      Sample = case_when(
        A == 1 & S == 1 ~ "RCT Treated",
        A == 0 & S == 1 ~ "RCT Controlled",
        sel   == 0 & S == 0 ~ "EC (Unselected)",
        sel   == 1 & S == 0 ~ "EC (Selected)"
      ) %>% 
        as_factor() %>% 
        fct_relevel(
          "RCT Treated",
          "RCT Controlled",
          "EC (Selected)",
          "EC (Unselected)"
        )
    ) %>% 
    filter(Sample != "RCT Treated") %>% 
    arrange(desc(Sample)) %>% 
    mutate(b = setup$mag_bias[case] %>% as_factor, .before = everything())
})

dat_fig %>% 
  ggplot(aes(
    `Sampling Score`, Y, 
    color = Sample,
    #shape = Sample
  )) +
  geom_point(alpha = 0.9, size = 2) +
  # hue_pal()(18)[c(7, 12, 1)]
  scale_color_manual(values = c("#5A5A5A",hue_pal()(18)[c(12, 1)]), name = "") +
  #scale_shape_manual(values = c(19,19,19))+
  theme(legend.position = "bottom") +
  # theme(plot.title = element_text(hjust = 0.5))+
  # ggtitle(str_glue("Magnitude of Hidden Bias = {setup$mag_bias[case]}"), )
  labs(
    x = "Sampling Score",
    y = "Outcome"
  ) +
  facet_wrap(vars(b), labeller = "label_both") +
  guides(
    #shape = guide_legend(""), 
    color = guide_legend("")
  )

ggsave(str_glue("chart/sim_{case_pre}_bias.pdf"), width = 9, height = 6)


# fig:sim_bias2 ----------------------------------------------------------------

dat_fig <- map_dfr(case_seq, function(case) {
  load("data/setup.RData")
  load(str_glue("output/cache/dt_{case}_1.RData"))
  load(str_glue("output/cache/raw_{case}.RData"))
  dat <- tibble(Y, A, S, X)
  id_sel <- raw_result[[1]]$res_alt[[prop_met_id]]$out$id_sel[[1]]
  fit <- glm(S ~ X, data = dat)
  dat %>%
    mutate(
      #y = ifelse(y == 3, y + runif(n(), 0, 0.1), y),
      `Sampling Score` = predict(fit, type = "response"),
      sel = ifelse(row_number() %in% id_sel, 1, 0),
      Sample = case_when(
        A == 1 & S == 1 ~ "RCT Treated",
        A == 0 & S == 1 ~ "RCT Controlled",
        sel   == 0 & S == 0 ~ "EC (Unselected)",
        sel   == 1 & S == 0 ~ "EC (Selected)"
      ) %>% 
        as_factor() %>% 
        fct_relevel(
          "RCT Treated",
          "RCT Controlled",
          "EC (Selected)",
          "EC (Unselected)"
        )
    ) %>% 
    filter(Sample != "RCT Treated") %>% 
    arrange(desc(Sample)) %>% 
    mutate(b = setup$mag_bias[case] %>% as_factor, .before = everything())
})

dat_fig %>% 
  ggplot(aes(
    `Sampling Score`, Y, 
    color = Sample,
    #shape = Sample
  )) +
  geom_point(alpha = 0.9, size = 2) +
  # hue_pal()(18)[c(7, 12, 1)]
  scale_color_manual(values = c("#5A5A5A",hue_pal()(18)[c(12, 1)]), name = "") +
  #scale_shape_manual(values = c(19,19,19))+
  theme(legend.position = "bottom") +
  # theme(plot.title = element_text(hjust = 0.5))+
  # ggtitle(str_glue("Magnitude of Hidden Bias = {setup$mag_bias[case]}"), )
  labs(
    x = "Sampling Score",
    y = "Outcome"
  ) +
  facet_grid(~b, labeller = "label_both") +
  guides(
    #shape = guide_legend(""), 
    color = guide_legend("")
  )

ggsave(str_glue("chart/sim_{case_pre}_bias2.pdf"), width = 12, height = 3)




# fig:sim_adag ---------------------------------------------------------------

res_inf <- map_dfr(case_seq, function(case) {
  load(str_glue("output/cache/inf_{case}.RData"))
  inf_result %>% 
    filter(method == "Conformal Selective Borrow AIPW (g=ada)") %>% 
    mutate(b = setup$mag_bias[case] %>% as_factor, .before = everything())
})

# p1 <- res_inf %>% 
#   ggplot(aes(b, gamma_sel)) + 
#   geom_boxplot() +
#   ylab(TeX("$\\gamma$"))

res_inf %>% 
  ggplot() + 
  geom_bar(aes(y = gamma_sel)) +
  facet_grid(. ~ b, labeller = label_both) +
  ylab(TeX("$\\hat{\\gamma}$"))

#plot_grid(p1, p2, nrow = 2)

ggsave(str_glue("chart/sim_{case_pre}_adag.pdf"), width = 9, height = 3)




# fig:sim_pdist_pre1 -------------------------------------------------------

met_id <- c(1, 2)

case_show <- 1

p_all0 <- map(case_show, function(case) {
  load(str_glue("output/cache/inf_{case}.RData"))
  load(str_glue("output/metrics_{case}.RData"))
  
  inf_tab <- inf_result %>% 
    mutate(Method = case_when(
      method == "No Borrow AIPW+FRT" ~ "NB",
      method == "Borrow AIPW+FRT" ~ "FB",
      method == "AdaLasso Selective Borrow ACW+FRT" ~ "ALSB",
      method == "Conformal Selective Borrow AIPW+FRT (g=0.8)" ~ "CSB0.8",
      method == "Conformal Selective Borrow AIPW+FRT (g=0.6)" ~ "CSB0.6",
      method == "Conformal Selective Borrow AIPW+FRT (g=0.4)" ~ "CSB0.4",
      method == "Conformal Selective Borrow AIPW+FRT (g=ada)" ~ "CSBada"
    ) %>% as_factor %>% fct_relevel("NB", "FB")) %>% 
    filter(Method %in% method_all[met_id]) %>% 
    select(Method, p_value, p_value_null)
  
  sum_tab <- inf_tab %>% 
    group_by(Method) %>% 
    summarise(
      `Type I` = mean(p_value_null <= 0.05),
      Power = mean(p_value <= 0.05)
    )
  
  t1 <- ggtitle(TeX("Distribution under $H_0$ (No Hidden Bias in ECs)"))
  t2 <- ggtitle(TeX("Distribution under $H_1$ (No Hidden Bias in ECs)"))
  
  p1 <- inf_tab %>% 
    ggplot() +
    geom_histogram(aes(p_value_null, fill = Method, color = Method), 
                   alpha = 0.5) +
    geom_vline(xintercept = 0.05, linetype = 2) +
    geom_label(data = sum_tab, ##
               aes(x = Inf, y = Inf, 
                   label = sprintf("Type I error = %.3f", `Type I`), ##
                   color = Method,
                   vjust = 1.5
               ),
               hjust = 1.1, size = 3.5, show.legend = F) +
    scale_fill_manual(values = color_all[met_id], labels = lab_all[met_id]) +
    scale_color_manual(values = color_all[met_id], labels = lab_all[met_id]) +
    facet_grid(Method~.)+
    ylab(NULL) +
    xlab("Exact p-value") +
    theme(strip.text = element_blank()) +
    theme(legend.position = "bottom") +
    t1
  
  p2 <- inf_tab %>% 
    ggplot() +
    geom_histogram(aes(p_value, fill = Method, color = Method), alpha = 0.5) + ##
    geom_vline(xintercept = 0.05, linetype = 2) +
    geom_label(data = sum_tab, ##
               aes(x = Inf, y = Inf, 
                   label = sprintf("Power = %.3f", Power), ##
                   color = Method,
                   vjust = 1.5
               ),
               hjust = 1.1, size = 3.5, show.legend = F) +
    scale_fill_manual(values = color_all[met_id], labels = lab_all[met_id]) +
    scale_color_manual(values = color_all[met_id], labels = lab_all[met_id]) +
    facet_grid(Method~.)+
    ylab(NULL) +
    xlab("Exact p-value") +
    theme(strip.text = element_blank()) +
    theme(legend.position = "bottom") +
    t2
  lst(p1, p2)
})

p_all <- list(p_all0[[1]][[1]], p_all0[[1]][[2]])

legend <- get_plot_component(
  p_all[[1]], "guide-box", return_all = T
)[[3]]


combined_plot <- do.call(
  plot_grid, 
  map(p_all, ~{. + theme(legend.position="none")})
)

plot_grid(combined_plot, legend, ncol = 1, rel_heights = c(1, 0.1))
ggsave(str_glue("chart/sim_{case_pre}_pdist_pre1.pdf"), width = 9, height = 3)



# fig:sim_pdist_pre2 -------------------------------------------------------

met_id <- c(1, 2)

case_show <- 9

p_all0 <- map(case_show, function(case) {
  load(str_glue("output/cache/inf_{case}.RData"))
  load(str_glue("output/metrics_{case}.RData"))
  
  inf_tab <- inf_result %>% 
    mutate(Method = case_when(
      method == "No Borrow AIPW+FRT" ~ "NB",
      method == "Borrow AIPW+FRT" ~ "FB",
      method == "AdaLasso Selective Borrow ACW+FRT" ~ "ALSB",
      method == "Conformal Selective Borrow AIPW+FRT (g=0.8)" ~ "CSB0.8",
      method == "Conformal Selective Borrow AIPW+FRT (g=0.6)" ~ "CSB0.6",
      method == "Conformal Selective Borrow AIPW+FRT (g=0.4)" ~ "CSB0.4",
      method == "Conformal Selective Borrow AIPW+FRT (g=ada)" ~ "CSBada"
    ) %>% as_factor %>% fct_relevel("NB", "FB")) %>% 
    filter(Method %in% method_all[met_id]) %>% 
    select(Method, p_value, p_value_null)
  
  sum_tab <- inf_tab %>% 
    group_by(Method) %>% 
    summarise(
      `Type I` = mean(p_value_null <= 0.05),
      Power = mean(p_value <= 0.05)
    )
  
  t1 <- ggtitle(TeX("Distribution under $H_0$ (Hidden Bias in ECs)"))
  t2 <- ggtitle(TeX("Distribution under $H_1$ (Hidden Bias in ECs)"))
  
  p1 <- inf_tab %>% 
    ggplot() +
    geom_histogram(aes(p_value_null, fill = Method, color = Method), 
                   alpha = 0.5) +
    geom_vline(xintercept = 0.05, linetype = 2) +
    geom_label(data = sum_tab, ##
               aes(x = Inf, y = Inf, 
                   label = sprintf("Type I error = %.3f", `Type I`), ##
                   color = Method,
                   vjust = 1.5
               ),
               hjust = 1.1, size = 3.5, show.legend = F) +
    scale_fill_manual(values = color_all[met_id], labels = lab_all[met_id]) +
    scale_color_manual(values = color_all[met_id], labels = lab_all[met_id]) +
    facet_grid(Method~.)+
    ylab(NULL) +
    xlab("Exact p-value") +
    theme(strip.text = element_blank()) +
    theme(legend.position = "bottom") +
    t1
  
  p2 <- inf_tab %>% 
    ggplot() +
    geom_histogram(aes(p_value, fill = Method, color = Method), alpha = 0.5) + ##
    geom_vline(xintercept = 0.05, linetype = 2) +
    geom_label(data = sum_tab, ##
               aes(x = Inf, y = Inf, 
                   label = sprintf("Power = %.3f", Power), ##
                   color = Method,
                   vjust = 1.5
               ),
               hjust = 1.1, size = 3.5, show.legend = F) +
    scale_fill_manual(values = color_all[met_id], labels = lab_all[met_id]) +
    scale_color_manual(values = color_all[met_id], labels = lab_all[met_id]) +
    facet_grid(Method~.)+
    ylab(NULL) +
    xlab("Exact p-value") +
    theme(strip.text = element_blank()) +
    theme(legend.position = "bottom") +
    t2
  lst(p1, p2)
})

p_all <- list(p_all0[[1]][[1]], p_all0[[1]][[2]])

legend <- get_plot_component(
  p_all[[1]], "guide-box", return_all = T
)[[3]]


combined_plot <- do.call(
  plot_grid, 
  map(p_all, ~{. + theme(legend.position="none")})
)

plot_grid(combined_plot, legend, ncol = 1, rel_heights = c(1, 0.1))
ggsave(str_glue("chart/sim_{case_pre}_pdist_pre2.pdf"), width = 9, height = 3)



# fig:sim_bias_pre ----------------------------------------------------------------

dat_fig <- map_dfr(case_seq, function(case) {
  load("data/setup.RData")
  load(str_glue("output/cache/dt_{case}_1.RData"))
  load(str_glue("output/cache/raw_{case}.RData"))
  dat <- tibble(Y, A, S, X)
  id_sel <- raw_result[[1]]$res_alt[[prop_met_id]]$out$id_sel[[1]]
  fit <- glm(S ~ X, data = dat)
  dat %>%
    mutate(
      #y = ifelse(y == 3, y + runif(n(), 0, 0.1), y),
      `Sampling Score` = predict(fit, type = "response"),
      sel = ifelse(row_number() %in% id_sel, 1, 0),
      Sample = case_when(
        A == 1 & S == 1 ~ "RCT Treated",
        A == 0 & S == 1 ~ "RCT Controlled",
        S == 0 ~ "EC"
      ) %>% 
        as_factor() %>% 
        fct_relevel(
          "RCT Treated",
          "RCT Controlled",
          "EC"
        )
    ) %>% 
    filter(Sample != "RCT Treated") %>% 
    arrange(desc(Sample)) %>% 
    mutate(b = setup$mag_bias[case] %>% as_factor, .before = everything())
})

dat_fig %>% 
  ggplot(aes(
    `Sampling Score`, Y, 
    color = Sample,
    shape = Sample
  )) +
  geom_point(alpha = 0.9, size = 2, color = "#5A5A5A") +
  # hue_pal()(18)[c(7, 12, 1)]
  #scale_color_manual(values = c("#5A5A5A",hue_pal()(18)[c(12, 1)]), name = "") +
  scale_shape_manual(values = c(19, 1), name = "") +
  theme(legend.position = "bottom") +
  # theme(plot.title = element_text(hjust = 0.5))+
  # ggtitle(str_glue("Magnitude of Hidden Bias = {setup$mag_bias[case]}"), )
  labs(
    x = "Sampling Score",
    y = "Outcome"
  ) +
  facet_wrap(vars(b), labeller = "label_both") +
  guides(
    #shape = guide_legend(""), 
    color = guide_legend("")
  )

ggsave(str_glue("chart/sim_{case_pre}_bias_pre.pdf"), width = 9, height = 6)

# fig:sim_main_pre ------------------------------------------------------------

#met_id <- c(1,2,3,5)

# custom_trans <- function(x) {
#   ifelse(x > 0.3 & x < 0.7, 0.3 + 4 * (x - 0.3), x)
# }
# 
# inverse_custom_trans <- function(x) {
#   ifelse(x > 0.3 & x < 2.1, 0.3 + (x - 0.3) / 4, x)
# }
# scale_power = scale_y_continuous(trans = scales::trans_new(
#   name = "custom",
#   transform = custom_trans,
#   inverse = inverse_custom_trans
# ))

plot_sim <- function(met_id) {
  
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
                    axis.title.y = element_blank(),
                    plot.title = element_text(hjust = 0.5))
  
  my_line <- geom_line(aes(group = Method, linetype = proposed),show.legend = FALSE)
  
  my_x <- "Magnitude of Hidden Bias" # Magnitude of Hidden Bias
  
  tImax <- max(
    res %>% 
      filter(Method %in% method_all[met_id]) %>% 
      pull(`Type I.y.up`),
    0.6
  )
  
  yp <- range(
    res %>% 
      filter(Method %in% method_all[met_id]) %>% 
      pull(Power.y), na.rm=T
  )+c(-0.1,0.1)
  scale_power <- ylim(yp[1], yp[2])
  
  # my_line <- geom_smooth(
  #   aes(group = Method), method = "lm", 
  #   formula = y ~ poly(x, 6), se = F, linewidth = 0.5
  # )
  
  # p1 <- res %>% 
  #   filter(Method %in% method_all[met_id]) %>%
  #   mutate(`Absolute Bias` = abs(Bias)) %>% 
  #   ggplot(aes(b, `Absolute Bias`, color = Method)) +
  #   my_line +
  #   geom_point() +
  #   my_theme +
  #   labs(x = my_x,
  #        title = "(A) Absolute Bias") +
  #   legend_lab
  
  # p2 <- res %>% 
  #   filter(Method %in% method_all[met_id]) %>%
  #   ggplot(aes(b, Var, color = Method)) +
  #   my_line +
  #   geom_point() +
  #   my_theme +
  #   labs(x = my_x,
  #        title = "(B) Variance") +
  #   legend_lab
  
  p3 <- res %>% 
    filter(Method %in% method_all[met_id]) %>%
    ggplot(aes(b, MSE, color = Method)) +
    my_line +
    geom_point() +
    my_theme +
    labs(x = my_x,
         title = "MSE") +
    legend_lab +
    scale_y_log10()
  
  p4 <- res %>% 
    filter(Method %in% method_all[met_id]) %>%
    ggplot(aes(b, `Type I.y`, color = Method)) +
    geom_hline(yintercept = 0.05) +
    #my_line +
    geom_point(position = position_dodge(width = 0.5)) +
    geom_errorbar(aes(ymin = `Type I.y.low`, ymax = `Type I.y.up`), 
                  width = 0.4, position = position_dodge(width = 0.5)) +
    my_theme +
    ylim(0, tImax) +
    labs(x = my_x, title = "Type I Error Rate") +
    legend_lab
  
  p5 <- res %>% 
    filter(Method %in% method_all[met_id]) %>%
    ggplot(aes(b, Power.y, color = Method)) +
    my_line +
    geom_point() +
    my_theme +
    labs(x = my_x, title = "Power") +
    legend_lab +
    scale_power
  
  # p6 <- res %>% 
  #   filter(Method %in% method_all[met_id]) %>%
  #   ggplot(aes(b, `FU/Sel`, color = Method)) +
  #   my_line +
  #   geom_point() +
  #   my_theme +
  #   labs(x = my_x, title = "(F) #Biased / #Selected") +
  #   legend_lab
  
  p_all <- map(list(p3, p4, p5), ~ {
    .x + theme(legend.position="none")
  })
  p_all$ncol <- 3
  legend <- get_plot_component(
    p3 + guides(color = guide_legend(nrow = lrow)), "guide-box", return_all = T
  )[[3]]
  combined_plot <- do.call(plot_grid, p_all)
  
  plot_grid(combined_plot, legend, ncol = 1, rel_heights = c(1, 0.1))
}

plot_sim(c(1, 2, prop_met_id))
ggsave(str_glue("chart/sim_{case_pre}_main_pre.pdf"), width = 8, height = 3)






# fig:sim_main_poster ------------------------------------------------------------

#met_id <- c(1,2,3,5)

# custom_trans <- function(x) {
#   ifelse(x > 0.3 & x < 0.7, 0.3 + 4 * (x - 0.3), x)
# }
# 
# inverse_custom_trans <- function(x) {
#   ifelse(x > 0.3 & x < 2.1, 0.3 + (x - 0.3) / 4, x)
# }
# scale_power = scale_y_continuous(trans = scales::trans_new(
#   name = "custom",
#   transform = custom_trans,
#   inverse = inverse_custom_trans
# ))

plot_sim <- function(met_id, leg_h = 0.1) {
  
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
                    axis.title.y = element_blank(),
                    plot.title = element_text(hjust = 0.5))
  
  my_line <- geom_line(aes(group = Method, linetype = proposed),show.legend = FALSE)
  
  my_x <- "b" # Magnitude of Hidden Bias
  
  tImax <- max(
    res %>% 
      filter(Method %in% method_all[met_id]) %>% 
      pull(`Type I.y.up`),
    0.6
  )
  
  yp <- range(
    res %>% 
      filter(Method %in% method_all[met_id]) %>% 
      pull(Power.y), na.rm=T
  )+c(-0.1,0.1)
  scale_power <- ylim(yp[1], yp[2])
  
  # my_line <- geom_smooth(
  #   aes(group = Method), method = "lm", 
  #   formula = y ~ poly(x, 6), se = F, linewidth = 0.5
  # )
  
  p1 <- res %>% 
    filter(Method %in% method_all[met_id]) %>%
    mutate(`Absolute Bias` = abs(Bias)) %>% 
    ggplot(aes(b, `Absolute Bias`, color = Method)) +
    my_line +
    geom_point() +
    my_theme +
    labs(x = my_x,
         title = "(A) Absolute Bias") +
    legend_lab
  
  p2 <- res %>% 
    filter(Method %in% method_all[met_id]) %>%
    ggplot(aes(b, Var, color = Method)) +
    my_line +
    geom_point() +
    my_theme +
    labs(x = my_x,
         title = "(B) Variance") +
    legend_lab
  
  p3 <- res %>% 
    filter(Method %in% method_all[met_id]) %>%
    ggplot(aes(b, MSE, color = Method)) +
    my_line +
    geom_point() +
    my_theme +
    labs(x = my_x,
         title = "(C) MSE") +
    legend_lab +
    scale_y_log10()
  
  p4 <- res %>% 
    filter(Method %in% method_all[met_id]) %>%
    ggplot(aes(b, `Type I.y`, color = Method)) +
    geom_hline(yintercept = 0.05) +
    #my_line +
    geom_point(position = position_dodge(width = 0.5)) +
    geom_errorbar(aes(ymin = `Type I.y.low`, ymax = `Type I.y.up`), 
                  width = 0.4, position = position_dodge(width = 0.5)) +
    my_theme +
    ylim(0, tImax) +
    labs(x = my_x, title = "(D) Type I Error Rate") +
    legend_lab
  
  p5 <- res %>% 
    filter(Method %in% method_all[met_id]) %>%
    ggplot(aes(b, Power.y, color = Method)) +
    my_line +
    geom_point() +
    my_theme +
    labs(x = my_x, title = "(E) Power") +
    legend_lab +
    scale_power
  
  # p6 <- res %>% 
  #   filter(Method %in% method_all[met_id]) %>%
  #   ggplot(aes(b, `FU/Sel`, color = Method)) +
  #   my_line +
  #   geom_point() +
  #   my_theme +
  #   labs(x = my_x, title = "(F) #Biased / #Selected") +
  #   legend_lab
  
  p_all <- map(list(p1, p2, p3, p4, p5), ~ {
    .x + theme(legend.position="none")
  })
  p_all$ncol <- 5
  legend <- get_plot_component(
    p1 + guides(color = guide_legend(nrow = lrow)) + labs(color = ""),
    "guide-box", return_all = T
  )[[3]]
  combined_plot <- do.call(plot_grid, p_all)
  
  plot_grid(combined_plot, legend, ncol = 1, rel_heights = c(1, leg_h))
}

lab_all <- c(
  TeX("No Borrow AIPW ($\\gamma=1$)"),
  TeX("Borrow AIPW ($\\gamma=0$)"),
  TeX("Adaptive Lasso Selective Borrow"),
  TeX("Conformal Selective Borrow ($\\gamma=0.8$)"),
  TeX("Conformal Selective Borrow ($\\gamma=0.6$)"),
  TeX("Conformal Selective Borrow ($\\gamma=0.4$)"),
  TeX("Conformal Selective Borrow AIPW ($\\hat{\\gamma}$)")
)

plot_sim(c(1, 2, prop_met_id))
ggsave(str_glue("chart/sim_{case_pre}_poster.pdf"), width = 12, height = 3)


# fig:sim_bias_poster ----------------------------------------------------------------

dat_fig <- map_dfr(case_seq, function(case) {
  load("data/setup.RData")
  load(str_glue("output/cache/dt_{case}_1.RData"))
  load(str_glue("output/cache/raw_{case}.RData"))
  dat <- tibble(Y, A, S, X)
  id_sel <- raw_result[[1]]$res_alt[[prop_met_id]]$out$id_sel[[1]]
  fit <- glm(S ~ X, data = dat)
  dat %>%
    mutate(
      #y = ifelse(y == 3, y + runif(n(), 0, 0.1), y),
      `Sampling Score` = predict(fit, type = "response"),
      sel = ifelse(row_number() %in% id_sel, 1, 0),
      Sample = case_when(
        A == 1 & S == 1 ~ "RCT-Treated",
        A == 0 & S == 1 ~ "RCT-Control",
        sel   == 0 & S == 0 ~ "EC (Unselected)",
        sel   == 1 & S == 0 ~ "EC (Selected)"
      ) %>% 
        as_factor() %>% 
        fct_relevel(
          "RCT-Treated",
          "RCT-Control",
          "EC (Selected)",
          "EC (Unselected)"
        )
    ) %>% 
    filter(Sample != "RCT-Treated") %>% 
    arrange(desc(Sample)) %>% 
    mutate(b = setup$mag_bias[case] %>% as_factor, .before = everything())
})

dat_fig %>% 
  ggplot(aes(
    `Sampling Score`, Y, 
    color = Sample,
    #shape = Sample
  )) +
  geom_point(alpha = 0.9, size = 1) +
  # hue_pal()(18)[c(7, 12, 1)]
  scale_color_manual(values = c("#5A5A5A",hue_pal()(18)[c(12, 1)]), name = "") +
  #scale_shape_manual(values = c(19,19,19))+
  theme(legend.position = "bottom") +
  # theme(plot.title = element_text(hjust = 0.5))+
  # ggtitle(str_glue("Magnitude of Hidden Bias = {setup$mag_bias[case]}"), )
  labs(
    x = "Sampling Score",
    y = "Outcome"
  ) +
  #facet_wrap(vars(b), labeller = "label_both") +
  facet_grid(~b, labeller = "label_both") +
  guides(
    #shape = guide_legend(""), 
    color = guide_legend("")
  )+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave(str_glue("chart/sim_{case_pre}_bias_poster.pdf"), width = 8, height = 2.5)


# fig:sim_csb_asym --------------------------------------------------------

lab_all <- c(
  TeX("No Borrow ($\\gamma=1$) + Asym Inf"),
  TeX("Full Borrow ($\\gamma=0$)"),
  TeX("Adaptive Lasso Selective Borrow + Asym Inf"),
  TeX("Conformal Selective Borrow ($\\gamma=0.8$)"),
  TeX("Conformal Selective Borrow ($\\gamma=0.6$)"),
  TeX("Conformal Selective Borrow ($\\gamma=0.4$)"),
  TeX("Conformal Selective Borrow ($\\hat{\\gamma}$) + Asym Inf")
)


res <- map_dfr(case_seq, function(case) {
  load(str_glue("output/metrics_{case}.RData"))
  load(str_glue("output/cache/inf_{case}.RData"))
  n_rep <- inf_result$rep %>% max
  
  met_frt <- metrics %>% 
    filter(str_detect(method, "FRT")) %>% 
    mutate(method = str_remove(method, "\\+FRT")) %>% 
    select(method, `Type I`, Power)
  
  metrics %>%
    filter(
      !str_detect(method, "FRT"),
      #method != "AdaLasso Selective Borrow ACW"
    ) %>% 
    left_join(met_frt, by = "method") %>% 
    mutate(
      Method = case_when(
        method == "No Borrow AIPW" ~ "NB",
        method == "Borrow AIPW" ~ "FB",
        method == "AdaLasso Selective Borrow ACW" ~ "ALSB",
        method == "Conformal Selective Borrow AIPW (g=0.8)" ~ "CSB0.8",
        method == "Conformal Selective Borrow AIPW (g=0.6)" ~ "CSB0.6",
        method == "Conformal Selective Borrow AIPW (g=0.4)" ~ "CSB0.4",
        method == "Conformal Selective Borrow AIPW (g=ada)" ~ "CSBada"
      )
    ) %>% 
    mutate(b = setup$mag_bias[case], .before = everything()) %>% 
    mutate(
      `Type I.y` = ifelse(Method %in% c("NB", "ALSB", "CSBada"), 
                          `Type I.x`, `Type I.y`), ###
      `Power.y` = ifelse(Method %in% c("NB", "ALSB", "CSBada"),
                         `Power.x`, `Power.y`), ###
      `Type I.y.low` = map_dbl(`Type I.y`, ~{
        binom.test(. * n_rep, n_rep, conf.level = 0.95)$conf.int[1]
      }),
      `Type I.y.up` = map_dbl(`Type I.y`, ~{
        binom.test(. * n_rep, n_rep, conf.level = 0.95)$conf.int[2]
      })
    )
}) %>% 
  mutate(
    Method = Method %>% 
      as_factor %>% 
      fct_relevel(
        "NB", 
        "FB", "ALSB", "CSB0.8", "CSB0.6", "CSB0.4", 
        "CSBada"
      ),
    proposed = ifelse(str_detect(Method, "CSB"), "y", "n") %>% 
      as_factor %>% 
      fct_relevel("y", "n")
  )

plot_sim_inf <- function(met_id, leg_h = 0.1) {
  
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
                    axis.title.y = element_blank(),
                    plot.title = element_text(hjust = 0.5))
  
  my_line <- geom_line(aes(group = Method, linetype = proposed),show.legend = FALSE)
  
  my_x <- "Magnitude of Hidden Bias" # Magnitude of Hidden Bias
  
  tImax <- max(
    res %>% 
      filter(Method %in% method_all[met_id]) %>% 
      pull(`Type I.y.up`),
    0.2
  )
  
  yp <- range(
    res %>% 
      filter(Method %in% method_all[met_id]) %>% 
      pull(Power.y), na.rm=T
  )+c(-0.1,0.1)
  scale_power <- ylim(yp[1], yp[2])
  
  # my_line <- geom_smooth(
  #   aes(group = Method), method = "lm", 
  #   formula = y ~ poly(x, 6), se = F, linewidth = 0.5
  # )
  
  p1 <- res %>% 
    filter(Method %in% method_all[met_id]) %>%
    mutate(`Absolute Bias` = abs(Bias)) %>% 
    ggplot(aes(b, `Absolute Bias`, color = Method)) +
    my_line +
    geom_point() +
    my_theme +
    labs(x = my_x,
         title = "(A) Absolute Bias") +
    legend_lab
  
  p2 <- res %>% 
    filter(Method %in% method_all[met_id]) %>%
    ggplot(aes(b, Var, color = Method)) +
    my_line +
    geom_point() +
    my_theme +
    labs(x = my_x,
         title = "(B) Variance") +
    legend_lab
  
  p3 <- res %>% 
    filter(Method %in% method_all[met_id]) %>%
    ggplot(aes(b, MSE, color = Method)) +
    my_line +
    geom_point() +
    my_theme +
    labs(x = my_x,
         title = "(C) MSE") +
    legend_lab +
    scale_y_log10()
  
  p4 <- res %>% 
    filter(Method %in% method_all[met_id]) %>%
    ggplot(aes(b, `Type I.y`, color = Method)) +
    geom_hline(yintercept = 0.05) +
    #my_line +
    geom_point(position = position_dodge(width = 0.5)) +
    geom_errorbar(aes(ymin = `Type I.y.low`, ymax = `Type I.y.up`), 
                  width = 0.4, position = position_dodge(width = 0.5)) +
    my_theme +
    ylim(0, tImax) +
    labs(x = my_x, title = "Type I Error Rate") +
    legend_lab
  
  p5 <- res %>% 
    filter(Method %in% method_all[met_id]) %>%
    ggplot(aes(b, Power.y, color = Method)) +
    my_line +
    geom_point() +
    my_theme +
    labs(x = my_x, title = "Power") +
    legend_lab +
    scale_power
  
  # p6 <- res %>% 
  #   filter(Method %in% method_all[met_id]) %>%
  #   ggplot(aes(b, `FU/Sel`, color = Method)) +
  #   my_line +
  #   geom_point() +
  #   my_theme +
  #   labs(x = my_x, title = "(F) #Biased / #Selected") +
  #   legend_lab
  
  p_all <- map(list(p4, p5), ~ {
    .x + theme(legend.position="none")
  })
  p_all$ncol <- 2
  legend <- get_plot_component(
    p1 + guides(color = guide_legend(nrow = lrow)), "guide-box", return_all = T
  )[[3]]
  combined_plot <- do.call(plot_grid, p_all)
  
  plot_grid(combined_plot, legend, ncol = 1, rel_heights = c(1, leg_h))
}


plot_sim_inf(c(1, 3, prop_met_id))
ggsave(str_glue("chart/sim_csb_asym.pdf"), width = 10, height = 5)

