library(tidyverse)
library(scales)
library(latex2exp)
theme_set(theme_bw())
load(file = "output/result.RData")

tab_frt <- tab %>% 
  filter(str_detect(method, "FRT")) %>% 
  mutate(method = str_remove(method, "\\+FRT")) %>% 
  select(method, p_value, runtime)

res <- tab %>%
  filter(
    !str_detect(method, "FRT"),
    method != "Conformal Selective Borrow AIPW (g=0.6)",
    method != "AdaLasso Selective Borrow ACW"
  ) %>% 
  left_join(tab_frt, by = "method") %>% 
  mutate(
    Method = case_when(
      method == "No Borrow DiM" ~ "No Borrow DifMean", 
      method == "No Borrow AIPW" ~ "No Borrow AIPW",
      method == "Borrow AIPW" ~ "Borrow AIPW",
      method == "Conformal Selective Borrow AIPW (g=ada)" ~ "CSB AIPW"
        #"Conformal Selective Borrow AIPW"
    ) %>% as_factor,
    Est = est, 
    SE = se, 
    CI_l = ci_l,
    CI_u = ci_u,
    asyp = ifelse(p_value.x < 0.001, "<0.001",sprintf("%.3f", p_value.x)), 
    frtp = ifelse(p_value.y < 0.001, "<0.001",sprintf("%.3f", p_value.y)), 
    label1 = "Asymptotic p-value",
    label2 = "FRT p-value",
    .keep = "none"
  )

res %>% 
  ggplot(aes(x = Method, color = Method)) +
  geom_label(
    aes(x = Method, y = max(CI_u) + 0.1, label = asyp, fill = label1),
    inherit.aes = FALSE,
    color = "black",  
    vjust = 0.8,
    size = 3
  ) +
  geom_label(
    aes(x = Method, y = max(CI_u) + 0.05, label = frtp, fill = label2),
    inherit.aes = FALSE,
    color = "black",  
    vjust = 0.8,
    size = 3
  ) +
  geom_point(aes(y = Est)) +
  geom_errorbar(aes(ymin = CI_l, ymax = CI_u), width = 0.2) +
  geom_hline(yintercept = 0, linetype = 2) +
  scale_color_manual(values = c("#5A5A5A", hue_pal()(18)[c(7, 1, 12)]), name = "") +
  scale_fill_manual(values = c(
    "Asymptotic p-value" = "gray",
    "FRT p-value" = "white"
  ), name = "") +
  theme(axis.text.x = element_blank(), axis.title.x = element_blank()) +
  ylab(TeX("$\\hat{\\tau}$"))
  #theme(legend.position = "bottom", legend.box = "vertical")+
  # guides(
  #   fill = guide_legend(order = 1, nrow = 2),
  #   color = guide_legend(order = 2, nrow = 1)
  # )

ggsave("chart/res.pdf", width = 4.5, height = 2.5)  


res <- tab %>%
  filter(
    !str_detect(method, "FRT"),
    method != "Conformal Selective Borrow AIPW (g=0.6)",
    method != "AdaLasso Selective Borrow ACW"
  ) %>% 
  left_join(tab_frt, by = "method") %>% 
  mutate(
    Method = case_when(
      method == "No Borrow DiM" ~ "No Borrow DifMean", 
      method == "No Borrow AIPW" ~ "No Borrow AIPW",
      method == "Borrow AIPW" ~ "Full Borrow AIPW",
      method == "Conformal Selective Borrow AIPW (g=ada)" ~ "CSB AIPW"
      #"Conformal Selective Borrow AIPW"
    ) %>% as_factor,
    Est = est, 
    SE = se, 
    CI_l = ci_l,
    CI_u = ci_u,
    asyp = ifelse(p_value.x < 0.001, "<0.001",sprintf("%.3f", p_value.x)), 
    frtp = ifelse(p_value.y < 0.001, "<0.001",sprintf("%.3f", p_value.y)), 
    label1 = "Asymptotic p-value",
    label2 = "FRT p-value",
    .keep = "none"
  )

res %>% 
  ggplot(aes(x = Method, color = Method)) +
  geom_label(
    aes(x = Method, y = max(CI_u) + 0.1, label = asyp, fill = label1),
    inherit.aes = FALSE,
    color = "black",  
    vjust = 0.8,
    size = 3
  ) +
  geom_label(
    aes(x = Method, y = max(CI_u) + 0.05, label = frtp, fill = label2),
    inherit.aes = FALSE,
    color = "black",  
    vjust = 0.8,
    size = 3
  ) +
  geom_point(aes(y = Est)) +
  geom_errorbar(aes(ymin = CI_l, ymax = CI_u), width = 0.2) +
  geom_hline(yintercept = 0, linetype = 2) +
  scale_color_manual(values = c("#5A5A5A", hue_pal()(18)[c(7, 1, 12)]), name = "") +
  scale_fill_manual(values = c(
    "Asymptotic p-value" = "gray",
    "FRT p-value" = "white"
  ), name = "") +
  theme(axis.text.x = element_blank(), axis.title.x = element_blank()) +
  ylab(TeX("$\\hat{\\tau}$"))
#theme(legend.position = "bottom", legend.box = "vertical")+
# guides(
#   fill = guide_legend(order = 1, nrow = 2),
#   color = guide_legend(order = 2, nrow = 1)
# )

ggsave("chart/res2.pdf", width = 4.5, height = 2.5)  



