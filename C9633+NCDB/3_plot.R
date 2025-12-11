library(tidyverse)
library(RColorBrewer)
library(scales)
library(latex2exp)
library(quantreg)
theme_set(theme_bw())
load("data/C9633+NCDB.RData")
load(file = "output/result.RData")


id_sel <- result[[6]]$out$id_sel[[1]]
#id_sel <- result[[4]]$out$id_sel[[1]]

fit <- glm(sample ~ sex + age + race + hist + tsize, data = dat)

dat_fig <- dat %>%
  mutate(
    #y = ifelse(y == 3, y + runif(n(), 0, 0.1), y),
    `Sampling Score` = predict(fit, type = "response"),
    sel = ifelse(row_number() %in% id_sel, 1, 0),
    Sample = case_when(
      treat == 1 & sample == 1 ~ "CALGB 9633 Treated",
      treat == 0 & sample == 1 ~ "CALGB 9633 Controlled",
      sel   == 0 & sample == 0 ~ "NCDB Controlled (Unselected)",
      sel   == 1 & sample == 0 ~ "NCDB Controlled (Selected)"
    ) %>% 
      as_factor() %>% 
      fct_relevel(
        "CALGB 9633 Treated", 
        "CALGB 9633 Controlled",
        "NCDB Controlled (Selected)",
        "NCDB Controlled (Unselected)",
      )
  ) %>% 
  filter(Sample != "CALGB 9633 Treated") %>% 
  arrange(desc(Sample))


# rd_ss_init ------------------------------------------------------------

dat_fig0 <- dat_fig %>% 
  mutate(
    Sample = case_when(
      treat == 0 & sample == 1 ~ "CALGB 9633 Control (S = 1, A = 0)",
      sample == 0 ~ "Matched NCDB Control (S = 0, A = 0)"
    ) %>% 
      as_factor() %>% 
      fct_relevel(
        "CALGB 9633 Control (S = 1, A = 0)"
      )
  ) 

qu <- 0.975
ql <- 0.025


fit975 <- rq(
  y ~ `Sampling Score`, tau = qu, 
  data = dat_fig0,
  subset = dat_fig0$Sample == "CALGB 9633 Control (S = 1, A = 0)"
)
dat_fig0$pred975 <- predict(fit975, newdata = dat_fig0)

fit025 <- rq(
  y ~ `Sampling Score`, tau = ql, 
  data = dat_fig0,
  subset = dat_fig0$Sample == "CALGB 9633 Control (S = 1, A = 0)"
)
dat_fig0$pred025 <- predict(fit025, newdata = dat_fig0)


dat_fig0 %>% 
  ggplot() +
  geom_ribbon(
    aes(`Sampling Score`, ymin = pred025, ymax = pred975), 
    fill = "grey80", alpha = 0.5
  ) +
  geom_point(
    aes(
      `Sampling Score`, y, 
      #color = Sample,
      shape = Sample
    ), 
    alpha = 0.8, size = 2, color = "#5A5A5A"
  ) +
  scale_shape_manual(values = c(19, 1), name = "") +
  #scale_color_manual(values = c("#5A5A5A",hue_pal()(18)[c(12)]), name = "") +
  theme(legend.position = "bottom") +
  labs(x = "Sampling Score", y = "Outcome (3-year RMST)")

ggsave("chart/rd_ss_init.pdf", width = 8, height = 5.5)

# rd_ss_sel ----------------------------------------------------------

dat_fig$pred975 <- predict(fit975, newdata = dat_fig)
dat_fig$pred025 <- predict(fit025, newdata = dat_fig)

# dat_fig %>% 
#   filter(y == 3) %>% 
#   ggplot(aes(`Sampling Score`, Sample, color = Sample)) +
#   geom_point(alpha = 0.8, size = 2, position = "jitter") +
#   scale_color_manual(values = c("#5A5A5A",hue_pal()(18)[c(12, 1)]), name = "") +
#   theme(legend.position = "bottom") +
#   labs(x = "Sampling Score", y = "Outcome (3-year RMST)")

dat_fig %>% 
  ggplot() +
  geom_ribbon(
    aes(`Sampling Score`, ymin = pred025, ymax = pred975), 
    fill = "grey80", alpha = 0.5
  ) +
  geom_point(
    aes(`Sampling Score`, y, color = Sample), 
    alpha = 0.8, size = 2
  ) +
  scale_color_manual(values = c("#5A5A5A",hue_pal()(18)[c(12, 1)]), name = "") +
  theme(legend.position = "bottom") +
  labs(x = "Sampling Score", y = "Outcome (3-year RMST)")

ggsave("chart/rd_ss_sel.pdf", width = 8, height = 5.5)

# rd_ss_sel2 ----------------------------------------------------------

dat_fig$pred975 <- predict(fit975, newdata = dat_fig)
dat_fig$pred025 <- predict(fit025, newdata = dat_fig)

# dat_fig %>% 
#   filter(y == 3) %>% 
#   ggplot(aes(`Sampling Score`, Sample, color = Sample)) +
#   geom_point(alpha = 0.8, size = 2, position = "jitter") +
#   scale_color_manual(values = c("#5A5A5A",hue_pal()(18)[c(12, 1)]), name = "") +
#   theme(legend.position = "bottom") +
#   labs(x = "Sampling Score", y = "Outcome (3-year RMST)")

dat_fig %>% 
  ggplot() +
  geom_ribbon(
    aes(`Sampling Score`, ymin = pred025, ymax = pred975), 
    fill = "grey80", alpha = 0.5
  ) +
  geom_point(
    aes(`Sampling Score`, y, color = Sample), 
    alpha = 0.8, size = 2
  ) +
  scale_color_manual(values = c("#5A5A5A",hue_pal()(18)[c(12, 1)]), name = "") +
  theme(legend.position = "bottom") +
  labs(x = "Sampling Score", y = "Outcome (3-year RMST)")

ggsave("chart/rd_ss_sel2.pdf", width = 10, height = 4)

# rd_ts_init (tumor size) --------------------------------------------------------------

dat_fig0ts <- dat_fig %>% 
  mutate(
    Sample = case_when(
      treat == 0 & sample == 1 ~ "CALGB 9633 Controlled",
      sample == 0 ~ "NCDB Controlled"
    ) %>% 
      as_factor() %>% 
      fct_relevel(
        "CALGB 9633 Controlled",
        "NCDB Controlled",
      )
  ) 

# qu <- 0.95
# ql <- 0.05

qu <- 0.975
ql <- 0.025


fit975 <- rq(
  y ~ tsize, tau = qu, 
  data = dat_fig0ts,
  subset = dat_fig0ts$Sample == "CALGB 9633 Controlled"
)
dat_fig0ts$pred975 <- predict(fit975, newdata = dat_fig0ts)

fit025 <- rq(
  y ~ tsize, tau = ql, 
  data = dat_fig0ts,
  subset = dat_fig0ts$Sample == "CALGB 9633 Controlled"
)
dat_fig0ts$pred025 <- predict(fit025, newdata = dat_fig0ts)

dat_fig0ts %>% 
  ggplot() +
  geom_ribbon(
    aes(tsize, ymin = pred025, ymax = pred975), 
    fill = "grey80", alpha = 0.5
  ) +
  geom_point(
    aes(
      tsize, y, 
      #color = Sample,
      shape = Sample
    ), 
    alpha = 0.8, size = 2, color = "#5A5A5A"
  ) +
  scale_shape_manual(values = c(19, 1), name = "") +
  #scale_color_manual(values = c("#5A5A5A",hue_pal()(18)[c(12)]), name = "") +
  theme(legend.position = "bottom") +
  labs(x = "Tumor Size", y = "Outcome (3-year RMST)")

ggsave("chart/rd_ts_init.pdf", width = 8, height = 5.5)

# rd_ts_sel (tumor size) ----------------------------------------------------------

dat_fig$pred975 <- predict(fit975, newdata = dat_fig)
dat_fig$pred025 <- predict(fit025, newdata = dat_fig)

# dat_fig %>% 
#   filter(y == 3) %>% 
#   ggplot(aes(`Sampling Score`, Sample, color = Sample)) +
#   geom_point(alpha = 0.8, size = 2, position = "jitter") +
#   scale_color_manual(values = c("#5A5A5A",hue_pal()(18)[c(12, 1)]), name = "") +
#   theme(legend.position = "bottom") +
#   labs(x = "Sampling Score", y = "Outcome (3-year RMST)")

dat_fig %>% 
  ggplot() +
  geom_ribbon(
    aes(tsize, ymin = pred025, ymax = pred975), 
    fill = "grey80", alpha = 0.5
  ) +
  geom_point(
    aes(tsize, y, color = Sample), 
    alpha = 0.8, size = 2
  ) +
  scale_color_manual(values = c("#5A5A5A",hue_pal()(18)[c(12, 1)]), name = "") +
  theme(legend.position = "bottom") +
  labs(x = "Tumor Size", y = "Outcome (3-year RMST)")


ggsave("chart/rd_ts_sel.pdf", width = 8, height = 5.5)



# rd_ss_init_poster ------------------------------------------------------------

dat_fig0 <- dat_fig %>% 
  mutate(
    Sample = case_when(
      treat == 0 & sample == 1 ~ "CALGB 9633 Control (S = 1, A = 0)",
      sample == 0 ~ "Matched NCDB Control (S = 0, A = 0)"
    ) %>% 
      as_factor() %>% 
      fct_relevel(
        "CALGB 9633 Control (S = 1, A = 0)"
      )
  ) 

qu <- 0.975
ql <- 0.025


fit975 <- rq(
  y ~ `Sampling Score`, tau = qu, 
  data = dat_fig0,
  subset = dat_fig0$Sample == "CALGB 9633 Control (S = 1, A = 0)"
)
dat_fig0$pred975 <- predict(fit975, newdata = dat_fig0)

fit025 <- rq(
  y ~ `Sampling Score`, tau = ql, 
  data = dat_fig0,
  subset = dat_fig0$Sample == "CALGB 9633 Control (S = 1, A = 0)"
)
dat_fig0$pred025 <- predict(fit025, newdata = dat_fig0)


dat_fig0 %>% 
  ggplot() +
  geom_ribbon(
    aes(`Sampling Score`, ymin = pred025, ymax = pred975), 
    fill = "grey80", alpha = 0.5
  ) +
  geom_point(
    aes(
      `Sampling Score`, y, 
      #color = Sample,
      shape = Sample
    ), 
    alpha = 0.8, size = 2, color = "#5A5A5A"
  ) +
  scale_shape_manual(values = c(19, 1), name = "") +
  #scale_color_manual(values = c("#5A5A5A",hue_pal()(18)[c(12)]), name = "") +
  theme(legend.position = "bottom") +
  labs(x = "Sampling Score", y = "Outcome (3-year RMST)")

ggsave("chart/rd_ss_init_poster.pdf", width = 8, height = 3)


# rd_ss_sel_poster ----------------------------------------------------------


dat_fig <- dat %>%
  mutate(
    #y = ifelse(y == 3, y + runif(n(), 0, 0.1), y),
    `Sampling Score` = predict(fit, type = "response"),
    sel = ifelse(row_number() %in% id_sel, 1, 0),
    Sample = case_when(
      treat == 1 & sample == 1 ~ "CALGB 9633 Treated",
      treat == 0 & sample == 1 ~ "CALGB 9633 Control",
      sel   == 0 & sample == 0 ~ "NCDB Control (Unselected)",
      sel   == 1 & sample == 0 ~ "NCDB Control (Selected)"
    ) %>% 
      as_factor() %>% 
      fct_relevel(
        "CALGB 9633 Treated", 
        "CALGB 9633 Control",
        "NCDB Control (Selected)",
        "NCDB Control (Unselected)",
      )
  ) %>% 
  filter(Sample != "CALGB 9633 Treated") %>% 
  arrange(desc(Sample))

dat_fig$pred975 <- predict(fit975, newdata = dat_fig)
dat_fig$pred025 <- predict(fit025, newdata = dat_fig)

# dat_fig %>% 
#   filter(y == 3) %>% 
#   ggplot(aes(`Sampling Score`, Sample, color = Sample)) +
#   geom_point(alpha = 0.8, size = 2, position = "jitter") +
#   scale_color_manual(values = c("#5A5A5A",hue_pal()(18)[c(12, 1)]), name = "") +
#   theme(legend.position = "bottom") +
#   labs(x = "Sampling Score", y = "Outcome (3-year RMST)")

dat_fig %>% 
  ggplot() +
  geom_ribbon(
    aes(`Sampling Score`, ymin = pred025, ymax = pred975), 
    fill = "grey80", alpha = 0.5
  ) +
  geom_point(
    aes(`Sampling Score`, y, color = Sample), 
    alpha = 0.8, size = 2
  ) +
  scale_color_manual(values = c("#5A5A5A",hue_pal()(18)[c(12, 1)]), name = "") +
  theme(legend.position = "bottom") +
  labs(x = "Sampling Score", y = "Outcome (3-year RMST)")

ggsave("chart/rd_ss_sel_poster.pdf", width = 8, height = 3)
