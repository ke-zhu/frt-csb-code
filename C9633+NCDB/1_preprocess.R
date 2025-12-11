library(tidyverse)
library(tictoc)
library(haven)
#library(summarytools)
library(survival)
library(eventglm)
library(MatchIt)
library(cobalt)
library(cowplot)
library(scales)
set.seed(2024)
theme_set(theme_bw())

dat_raw <- read_sas("data/all.sas7bdat", NULL)

# dat_raw %>%
#   filter(cohort == "C9633") %>%
#   count(treat, hist)

# dat_raw %>% 
#   filter(! ((treat == 1) & (cohort == "NCDB"))) %>% 
#   count(cohort, survcens)
# 
# 141 / (141 + 194)
# 5584 / (5584+6116)

# # EDA
# dat_raw %>%
#   dfSummary %>%
#   view
# 
# dat_raw %>% 
#   group_by(patid) %>%
#   filter(n() > 1) %>% 
#   View
# 
# dat_raw %>% 
#   filter(cohort == "C9633") %>% 
#   dfSummary %>% 
#   view


# 1 impute missing X ------------------------------------------

dat_rct <- dat_raw %>% 
  filter(cohort == "C9633") %>% 
  select(
    treat,
    sex, age, race, hist, tsize,
    survtime, survcens
  ) %>% 
  mutate(
    # impute 8 missing samples with median
    tsize = ifelse(is.na(tsize), median(tsize, na.rm = T), tsize),
    sample = 1
  )

# dat_rct %>% 
#   dfSummary %>% 
#   view

dat_ec <- dat_raw %>% 
  filter(cohort == "NCDB", treat == 0) %>% 
  select(
    treat,
    sex, age, race, hist, tsize,
    survtime, survcens
  ) %>% 
  drop_na() %>% 
  mutate(sample = 0)

# dat_ec %>% 
#   dfSummary %>% 
#   view



# 2 transform to pseudo-observation ---------------------------------------

res_time <- 3

pse_obs_rct <- rmeanglm(
  Surv(survtime, survcens) ~ 1,
  time = res_time,
  model.censoring = "stratified",
  formula.censoring = ~ sex + race + hist, 
  data = dat_rct
)$y

pse_obs_ec <- rmeanglm(
  Surv(survtime, survcens) ~ 1,
  time = res_time,
  model.censoring = "stratified",
  formula.censoring = ~ sex + race + hist, 
  data = dat_ec
)$y

dat_rct_trans <- dat_rct %>% 
  mutate(y = pse_obs_rct) %>% 
  mutate(y = pmin(y, res_time))
dat_ec_trans <- dat_ec %>% 
  mutate(y = pse_obs_ec) %>% 
  mutate(y = pmin(y, res_time))


# figure 1: pseudo-observation vs. censored time ----

dat_fig1 <- rbind(dat_rct_trans, dat_ec_trans) %>%
  mutate(
    `Censored Time` = survtime,
    `Pseudo-observation` = y,
    Sample = ifelse(sample == 1, "CALGB 9633", "NCDB") %>% as_factor(),
    Censored = ifelse(survcens == 1, "No", "Yes") %>% as_factor()
  ) 

dat_fig1 %>% 
  ggplot(aes(`Censored Time`, `Pseudo-observation`, color = Censored)) +
  geom_point(alpha=0.5) +
  geom_abline(slope = 1, intercept = 0, linetype=2) +
  facet_grid(Sample~.) +
  theme(legend.position = "bottom")+
  coord_fixed(ratio = 1, ylim = c(-0.2, res_time+0.3))
  #ggtitle("Pseudo-observation vs. Censored Time"))
ggsave(filename = "chart/rd_pseudo.pdf", width = 7, height = 5)

# dat_fig1 %>%
#   mutate(Sample = case_when(
#     treat == 1 & sample == 1 ~ "CALGB 9633 (treat)",
#     treat == 0 & sample == 1 ~ "CALGB 9633 (control)",
#     treat == 0 & sample == 0 ~ "NCDB (external control)"
#   ) %>% as_factor() %>% fct_relevel("CALGB 9633 (treat)")) %>% 
#   ggplot(aes(`Pseudo-observation`)) +
#   geom_histogram() +
#   scale_y_sqrt() +
#   facet_grid(Sample~., scales = "free_y")
#   #ggtitle("Distribution of Pseudo-observation"))
# 
# ggsave(filename = "chart/pseudo_dist.pdf", width = 10, height = 6)





# 3 matching ----------------------------------------------------------------


# remove ECs that are out of range
dat_trans <- rbind(
  dat_rct_trans,
  dat_ec_trans %>% 
    filter(
      age >= min(dat_rct_trans$age) & age <= max(dat_rct_trans$age),
      tsize >= min(dat_rct_trans$tsize) & tsize <= max(dat_rct_trans$tsize)
    )
)




# match
m.out <- matchit(
  sample ~ sex + age + race + hist + tsize,
  data = dat_trans,
  method = "nearest", distance = "glm", replace = FALSE,
  exact = c("hist")
)
summary(m.out)

# matched data
dat <- match.data(m.out) %>% select(-weights, -subclass)

# summary statistics before/after matching

dat %>% 
  group_by(treat, sample) %>% 
  summarise(m = mean(y))

# save data -----------------------------------------------------

save(dat_rct_trans, dat_ec_trans, dat, 
     file = "data/C9633+NCDB.RData")



# figure 2: X balance, RCT vs. EC ----
dat_fig2 <- dat_trans %>% 
  rename(
    Sex = sex, Age = age, Race = race,
    Histology = hist, `Tumor Size` = tsize
  )
  
m.out <- matchit(
  sample ~ Sex + Age + Race + Histology + `Tumor Size`,
  data = dat_fig2,
  method = "nearest", distance = "glm", replace = FALSE,
  exact = c("Histology")
)

legend_plot <- bal.plot(
  m.out, "Age", which = "both", 
  sample.names = c("Unmatched", "Matched")
) + 
  scale_fill_discrete(name = "Sample",labels = c( "NCDB (S = 0)", "CALGB 9633 (S = 1)")) +
  theme_void() + 
  theme(legend.position = "bottom")
legend <- get_plot_component(legend_plot, "guide-box", return_all = T)[[3]]

p_all <- map(c("Sex", "Age", "Race", "Histology", "Tumor Size"), ~ {
  bal.plot(
    m.out, .x, which = "both", 
    sample.names = c("Unmatched", "Matched")
  ) + 
    scale_fill_discrete(name = "Sample") +
    theme(legend.position = "none") +
    ggtitle(NULL)
})
p_all[[6]] <- bal.plot(
  m.out, "distance", which = "both", 
  sample.names = c("Unmatched", "Matched")
) + 
  scale_fill_discrete(name = "Sample") +
  theme(legend.position = "none") +
  xlab("Sampling Score") +
  ggtitle(NULL)
  #ggtitle("Distributional Balance for \"Sampling Score\"")

combined_plot <- do.call(plot_grid, p_all)

plot_grid(combined_plot, legend, ncol = 1, rel_heights = c(1, 0.1))


ggsave(filename = "chart/rd_cov.pdf", width = 8, height = 5)





# plot for poster ----
dat_ori <- dat_trans %>% 
  mutate(
    Sex = factor(sex) %>% 
      recode_factor(`1` = "Male", `0` = "Female"), 
    Age = age,
    Race = factor(race) %>% 
      recode_factor(`1` = "White", `0` = "Other"), 
    Histology = factor(hist) %>% 
      recode_factor(`1` = "Squamous", `0` = "Other"), 
    `Tumor Size` = tsize
  )
m.out <- matchit(
  sample ~ Sex + Age + Race + Histology + `Tumor Size`,
  data = dat_ori,
  method = "nearest", replace = FALSE
)
p_all <- map(c("Sex", "Age", "Race", "Histology", "Tumor Size"), ~ {
  bal.plot(
    m.out, ., which = "unadjusted", sample.names = c("Unmatched"),
    colors = hue_pal()(18)[c(1, 12)]
  ) + 
    theme(legend.position = "none", strip.text = element_blank()) +
    ggtitle(NULL)
})
legend_plot <- bal.plot(
  m.out, "Age", which = "unadjusted", sample.names = c("Unmatched")
) + 
  scale_fill_manual(
    name = "Sample",
    labels = c( "NCDB (S = 0)", "CALGB 9633 (S = 1)"),
    values = hue_pal()(18)[c(1, 12)]
  ) +
  theme_void() + 
  theme(legend.position = "bottom")
legend <- get_plot_component(legend_plot, "guide-box", return_all = T)[[3]]

combined_plot <- do.call(plot_grid, c(p_all, ncol = 5))

plot_grid(combined_plot, legend, ncol = 1, rel_heights = c(1, 0.1))
ggsave(filename = "chart/rd_cov_poster.pdf", width = 10, height = 3)




# for slides --------------------------------------------------------------
dat_fig2 <- dat_trans %>% 
  rename(
    Sex = sex, Age = age, Race = race,
    Histology = hist, `Tumor Size` = tsize
  )

m.out <- matchit(
  sample ~ Sex + Age + Race + Histology + `Tumor Size`,
  data = dat_fig2,
  method = "nearest", distance = "glm", replace = FALSE,
  exact = c("Histology")
)

legend_plot <- bal.plot(
  m.out, "Age", which = "both", 
  sample.names = c("Unmatched", "Matched")
) + 
  scale_fill_manual(
    name = "Sample",labels = c( "NCDB (S = 0)", "CALGB 9633 (S = 1)"),
    values = hue_pal()(18)[c(1, 12)]
  ) +
  theme_void() + 
  theme(legend.position = "bottom")
legend <- get_plot_component(legend_plot, "guide-box", return_all = T)[[3]]

p_all <- map(c("Sex", "Age", "Race", "Histology", "Tumor Size"), ~ {
  bal.plot(
    m.out, .x, which = "both", 
    sample.names = c("Unmatched", "Matched")
  ) + 
    scale_fill_manual(name = "Sample",  values = hue_pal()(18)[c(1, 12)]) +
    theme(legend.position = "none") +
    ggtitle(NULL)
})

combined_plot <- do.call(plot_grid, p_all)

plot_grid(combined_plot, legend, ncol = 1, rel_heights = c(1, 0.1))


ggsave(filename = "chart/rd_cov_slides.pdf", width = 8, height = 5)


