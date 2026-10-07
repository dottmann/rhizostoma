#####################################################################
## Project: Early growth and morphological transition of
##          Rhizostoma pulmo ephyrae in culture
##
## Authors: Daniel Ottmann
## Email:   daniel.ottmann.riera@gmail.com
##
## Date created: October 2026
## Last update:  October 2026
##
## Design notes
## - 19 ephyrae (one per polyp), followed individually from release
##   (day 0, < 24 h old) to day 30. Released in three batches on
##   different calendar days and measured Mon/Wed/Fri. Neither feeding
##   nor measuring took place at weekends, so on any given day-age only
##   a subset of individuals was measured (n = 7-18).
## - Because of this unbalanced, repeated-measures design, growth is
##   estimated from each individual's fitted trajectory (mixed-effects
##   models with individual as random effect), NOT from consecutive
##   day-to-day differences or from regressions on daily means.
## - Consecutive-interval growth is kept only as a diagnostic, to
##   quantify the effect of unfed weekends (supplementary figure).
#####################################################################

# ---- Packages -------------------------------------------------------
library(tidyverse)
library(lme4)
library(patchwork)

# ---- Working directory (project root = parent of /scripts) ---------
if (Sys.getenv("RSTUDIO") == "1") {
  script_dir <- dirname(rstudioapi::getActiveDocumentContext()$path)
  setwd(path.expand(script_dir))
  setwd("../")
}
dir.create("plots",  showWarnings = FALSE)
dir.create("output", showWarnings = FALSE)

# Metric labels. TBL in the raw files is renamed TBD (total body
# diameter) to follow Straehler-Pohl & Jarms (2010).
metric_levels <- c("TBD", "CDD", "TMLL", "LW")

# Colour-blind-safe palette (Okabe-Ito) whose hues match the arrows in
# Fig. 1 (blue = TBD, vermillion = CDD, green = TMLL, orange = LW).
metric_cols <- c(TBD = "#0072B2", CDD = "#D55E00", TMLL = "#009E73", LW = "#E69F00")

lmm_ctl <- lmerControl(optimizer = "bobyqa")


# ==============================================================================
# 1. DATA LOADING AND CLEANING
# ==============================================================================

read_metric <- function(file_path, metric_name) {
  read.delim(file_path, header = TRUE, check.names = FALSE) %>%
    pivot_longer(-id, names_to = "day", values_to = "size") %>%
    mutate(day = as.numeric(str_remove(day, "Day ")),
           metric = metric_name) %>%
    filter(!is.na(size))
}

df_long <- bind_rows(
  read_metric("data/tbl.txt",  "TBD"),
  read_metric("data/cdd.txt",  "CDD"),
  read_metric("data/tmll.txt", "TMLL"),
  read_metric("data/lw.txt",   "LW")
)

# Data-entry correction: for ephyra 8, TBD, CDD and TMLL were entered
# at day 12 but LW at day 13. In that batch's calendar day 12 falls on
# a Sunday (no measurements); all other ephyrae of the batch were
# measured on day 13. The record is therefore re-assigned to day 13.
df_long <- df_long %>%
  mutate(day = if_else(id == 8 & day == 12, 13, day))

df_long <- df_long %>%
  mutate(metric = factor(metric, levels = metric_levels)) %>%
  arrange(metric, id, day)

# ---- Feeding calendar ------------------------------------------------
# Measurements were taken Mon/Wed/Fri, so a 3-day gap between
# consecutive measurements always spans a weekend (Fri -> Mon). The
# first such gap fixes each individual's weekday calendar; Saturdays
# and Sundays (no food) then follow every 7 days.
calendar <- df_long %>%
  filter(metric == "TBD") %>%
  group_by(id) %>%
  arrange(day, .by_group = TRUE) %>%
  summarise(friday = day[which(diff(day) == 3)[1]], .groups = "drop") %>%
  mutate(friday_mod = friday %% 7,
         batch = LETTERS[as.integer(factor(friday_mod, levels = c(3, 0, 2)))])

is_unfed_day <- function(day, friday_mod) ((day - friday_mod) %% 7) %in% c(1, 2)

df_long <- df_long %>% left_join(calendar %>% select(id, batch, friday_mod), by = "id")

cat("\n--- Release batches ---\n")
print(calendar %>% count(batch, name = "n_ephyrae"))


# ==============================================================================
# 2. DESCRIPTIVE SUMMARY BY DAY-AGE (n, mean, SD, CV)
# ==============================================================================
# Daily CV = 100 * SD / mean across the individuals measured that day.
# It describes between-individual relative variability at a given age
# (biological variability + measurement error), not measurement
# precision per se. n varies by day and is always reported with it.

summary_df <- df_long %>%
  group_by(metric, day) %>%
  summarise(n_samples = n(),
            avg_size  = mean(size),
            sd_size   = sd(size),
            .groups   = "drop") %>%
  mutate(daily_cv = 100 * sd_size / avg_size)

cat("\n--- n measured per day-age (TBD) ---\n")
print(summary_df %>% filter(metric == "TBD") %>% select(day, n_samples), n = Inf)

write_csv(summary_df, "output/daily_summary.csv")


# ==============================================================================
# 3. INTERNAL CONSISTENCY OF THE FOUR MEASUREMENTS
# ==============================================================================
# Geometrically, TBD (lappet tip to opposite lappet tip) should equal
# CDD + 2 * TMLL. Because all three come from the same image, the
# scatter of this ratio indicates how consistently landmarks were
# placed, and any drift with age shows changes in body geometry.

wide <- df_long %>%
  select(id, batch, day, metric, size) %>%
  pivot_wider(names_from = metric, values_from = size) %>%
  mutate(closure   = (CDD + 2 * TMLL) / TBD,
         cdd_ratio = CDD / TBD)

closure_fit <- lmer(closure ~ day + (1 | id), data = wide, control = lmm_ctl)

cat("\n--- Closure (CDD + 2*TMLL) / TBD ---\n")
print(wide %>% summarise(mean = mean(closure, na.rm = TRUE),
                         sd   = sd(closure, na.rm = TRUE),
                         cv   = 100 * sd / mean,
                         min  = min(closure, na.rm = TRUE),
                         max  = max(closure, na.rm = TRUE)))
cat("Slope of closure with age (per day):\n")
print(summary(closure_fit)$coefficients)


# ==============================================================================
# 4. GROWTH FROM INDIVIDUAL FITTED TRAJECTORIES (mixed-effects models)
# ==============================================================================
# Candidate trajectory forms for each metric, all with random intercept
# and slope per individual, fitted by maximum likelihood:
#   linear       size      = a + b * day
#   quadratic    size      = a + b * day + c * day^2
#   exponential  log(size) = a + g * day   (g = specific growth rate)
# AIC of the log-scale model is put on the original scale by adding the
# Jacobian term 2 * sum(log(size)), so the three are comparable.

fit_candidates <- function(dat) {
  list(
    linear      = lmer(size ~ day + (1 + day | id), data = dat, REML = FALSE, control = lmm_ctl),
    quadratic   = lmer(size ~ day + I(day^2) + (1 + day | id), data = dat, REML = FALSE, control = lmm_ctl),
    exponential = lmer(log(size) ~ day + (1 + day | id), data = dat, REML = FALSE, control = lmm_ctl)
  )
}

model_table <- list()
best_fits   <- list()

for (m in metric_levels) {
  dat  <- filter(df_long, metric == m)
  fits <- suppressMessages(fit_candidates(dat))
  aic  <- c(linear      = AIC(fits$linear),
            quadratic   = AIC(fits$quadratic),
            exponential = AIC(fits$exponential) + 2 * sum(log(dat$size)))
  model_table[[m]] <- tibble(metric = m, model = names(aic), AIC = aic,
                             dAIC = aic - min(aic))
  best <- names(which.min(aic))
  # Prefer the simpler linear form when it is within 2 AIC units of the best
  if (best != "linear" && aic["linear"] - min(aic) < 2) best <- "linear"
  best_fits[[m]] <- list(form = best, fit = fits[[best]], data = dat)
}

model_table <- bind_rows(model_table)
cat("\n--- Trajectory model comparison (AIC, ML) ---\n")
print(model_table %>% mutate(across(c(AIC, dAIC), ~ round(.x, 1))), n = Inf)
write_csv(model_table, "output/trajectory_model_aic.csv")

# ---- Fixed-effect (population) equations ---------------------------
equation_table <- imap_dfr(best_fits, function(b, m) {
  cf <- fixef(b$fit); se <- sqrt(diag(vcov(b$fit)))
  tibble(metric = m, form = b$form, term = names(cf),
         estimate = unname(cf), se = unname(se))
})
cat("\n--- Population trajectory equations (fixed effects) ---\n")
print(equation_table %>% mutate(across(c(estimate, se), ~ signif(.x, 3))), n = Inf)
write_csv(equation_table, "output/trajectory_equations.csv")

# ---- Individual fitted trajectories ----------------------------------
# Predicted size for each individual over its own observed age range
# (no extrapolation), back-transformed for the exponential form.
predict_individual <- function(b) {
  b$data %>%
    group_by(id) %>%
    summarise(day_start = min(day), day_end = max(day), .groups = "drop") %>%
    mutate(day = map2(day_start, day_end, ~ seq(.x, .y, by = 0.25))) %>%
    unnest(day) %>%
    mutate(pred = predict(b$fit, newdata = ., re.form = NULL),
           pred = if (b$form == "exponential") exp(pred) else pred)
}

ind_traj <- imap_dfr(best_fits, ~ predict_individual(.x) %>% mutate(metric = .y))

# Growth estimates per individual, from its fitted trajectory between
# the first and last day it was measured:
#   AGR = absolute growth rate (mm d-1)
#   SGR = specific growth rate (% d-1) = 100 * ln(L_end / L_start) / dt
growth_individual <- ind_traj %>%
  group_by(metric, id) %>%
  summarise(day_start  = first(day), day_end = last(day),
            size_start = first(pred), size_end = last(pred),
            .groups = "drop") %>%
  mutate(AGR = (size_end - size_start) / (day_end - day_start),
         SGR = 100 * log(size_end / size_start) / (day_end - day_start))

growth_summary <- growth_individual %>%
  group_by(metric) %>%
  summarise(n = n(),
            AGR_mean = mean(AGR), AGR_sd = sd(AGR),
            SGR_mean = mean(SGR), SGR_sd = sd(SGR),
            .groups = "drop") %>%
  left_join(map_dfr(best_fits, ~ tibble(form = .x$form), .id = "metric"), by = "metric")

cat("\n--- Growth from individual fitted trajectories (mean +- SD across ephyrae) ---\n")
print(growth_summary %>% mutate(across(where(is.double), ~ round(.x, 3))))
write_csv(growth_individual, "output/growth_individual.csv")
write_csv(growth_summary,    "output/growth_summary.csv")

# TMLL follows a dome: age at which each individual's lappets reach
# their maximum projected length (only meaningful for the quadratic form)
if (best_fits$TMLL$form == "quadratic") {
  tmll_fit  <- best_fits$TMLL$fit
  tmll_coef <- coef(tmll_fit)$id
  tmll_peak <- tibble(id   = as.integer(rownames(tmll_coef)),
                      peak = -tmll_coef[["day"]] / (2 * tmll_coef[["I(day^2)"]]))
  pop_peak  <- -fixef(tmll_fit)[["day"]] / (2 * fixef(tmll_fit)[["I(day^2)"]])
  cat("\n--- Age of maximum TMLL ---\n")
  cat(sprintf("Population: day %.1f; individuals: %.1f +- %.1f (range %.1f-%.1f)\n",
              pop_peak, mean(tmll_peak$peak), sd(tmll_peak$peak),
              min(tmll_peak$peak), max(tmll_peak$peak)))
}


# ==============================================================================
# 5. MORPHOLOGICAL TRANSITION: BREAKPOINT IN CDD/TBD
# ==============================================================================
# Piecewise-linear ("broken-stick") mixed model:
#   ratio = a + b1 * day + b2 * max(day - bp, 0) + random(id)
# The breakpoint bp is estimated by profile likelihood over a grid; its
# approximate 95% CI contains all bp within 3.84 deviance units of the
# minimum. Breakpoints are also estimated separately for each
# individual to describe developmental asynchrony.

ratio_dat <- wide %>% filter(!is.na(cdd_ratio))
bp_grid   <- seq(8, 24, by = 0.1)

fit_bp <- function(bp, dat) {
  dat$hinge <- pmax(dat$day - bp, 0)
  suppressMessages(lmer(cdd_ratio ~ day + hinge + (1 + hinge | id),
                        data = dat, REML = FALSE, control = lmm_ctl))
}

bp_profile <- tibble(bp = bp_grid) %>%
  mutate(deviance = map_dbl(bp, ~ deviance(fit_bp(.x, ratio_dat))))

bp_hat  <- bp_profile$bp[which.min(bp_profile$deviance)]
bp_ci   <- range(bp_profile$bp[bp_profile$deviance - min(bp_profile$deviance) <= qchisq(0.95, 1)])
bp_fit  <- fit_bp(bp_hat, ratio_dat)

# Compare against a single straight line (no transition)
lin_ratio_fit <- suppressMessages(lmer(cdd_ratio ~ day + (1 + day | id), data = ratio_dat,
                                       REML = FALSE, control = lmm_ctl))

cat("\n--- CDD/TBD breakpoint (population) ---\n")
cat(sprintf("Breakpoint: day %.1f (95%% profile CI %.1f-%.1f)\n", bp_hat, bp_ci[1], bp_ci[2]))
cat(sprintf("Slope before: %.2f %% d-1; slope after: %.2f %% d-1\n",
            100 * fixef(bp_fit)[["day"]],
            100 * (fixef(bp_fit)[["day"]] + fixef(bp_fit)[["hinge"]])))
cat(sprintf("Ratio at breakpoint: %.1f %%\n",
            100 * (fixef(bp_fit)[["(Intercept)"]] + fixef(bp_fit)[["day"]] * bp_hat)))
cat(sprintf("AIC broken-stick: %.1f; AIC straight line: %.1f\n", AIC(bp_fit), AIC(lin_ratio_fit)))

# Is the change abrupt or a smooth acceleration? Compare with smooth
# polynomial curves. The breakpoint was estimated, so it is counted as
# one extra parameter (+2 AIC) in the broken-stick model.
quad_ratio_fit  <- suppressMessages(lmer(cdd_ratio ~ day + I(day^2) + (1 + day | id), data = ratio_dat,
                                         REML = FALSE, control = lmm_ctl))
cubic_ratio_fit <- suppressMessages(lmer(cdd_ratio ~ day + I(day^2) + I(day^3) + (1 + day | id),
                                         data = ratio_dat, REML = FALSE, control = lmm_ctl))
ratio_model_table <- tibble(
  model = c("broken-stick", "straight line", "quadratic", "cubic"),
  AIC   = c(AIC(bp_fit) + 2, AIC(lin_ratio_fit), AIC(quad_ratio_fit), AIC(cubic_ratio_fit))
) %>% mutate(dAIC = AIC - min(AIC))
cat("\n--- CDD/TBD: abrupt vs smooth change (AIC; breakpoint counted as a parameter) ---\n")
print(ratio_model_table %>% mutate(across(where(is.double), ~ round(.x, 1))))
write_csv(ratio_model_table, "output/ratio_model_aic.csv")

# Individual breakpoints (least squares, at least 3 points each side)
individual_bp <- ratio_dat %>%
  group_by(id) %>%
  group_modify(function(d, key) {
    grid <- bp_grid[map_lgl(bp_grid, ~ sum(d$day < .x) >= 3 & sum(d$day > .x) >= 3)]
    rss  <- map_dbl(grid, function(bp) {
      deviance(lm(cdd_ratio ~ day + pmax(day - bp, 0), data = d))
    })
    tibble(bp = grid[which.min(rss)])
  }) %>%
  ungroup()

cat(sprintf("Individual breakpoints: %.1f +- %.1f d (range %.1f-%.1f)\n",
            mean(individual_bp$bp), sd(individual_bp$bp),
            min(individual_bp$bp), max(individual_bp$bp)))

if (exists("tmll_peak")) {
  bp_vs_peak <- individual_bp %>% left_join(tmll_peak, by = "id")
  cat(sprintf("Correlation individual breakpoint vs TMLL peak: r = %.2f\n",
              cor(bp_vs_peak$bp, bp_vs_peak$peak)))
  write_csv(bp_vs_peak, "output/individual_transition.csv")
}

# Daily CV before and after the transition
cv_by_phase <- summary_df %>%
  mutate(phase = if_else(day < bp_hat, "before", "after")) %>%
  group_by(metric, phase) %>%
  summarise(mean_daily_cv = weighted.mean(daily_cv, n_samples), .groups = "drop") %>%
  pivot_wider(names_from = phase, values_from = mean_daily_cv)


# ==============================================================================
# 6. VARIABILITY OF EACH MEASUREMENT
# ==============================================================================
# (a) Between-individual variability at a given age: daily CV, averaged
#     across days weighted by the number of ephyrae measured.
# (b) Within-individual scatter: residual SD around each ephyra's own
#     fitted trajectory, relative to its fitted size. This combines
#     measurement error and short-term fluctuations (e.g. contraction
#     state, unfed weekends) and is the closer proxy for measurement
#     repeatability available without repeated imaging.

residual_cv <- imap_dfr(best_fits, function(b, m) {
  fitted_size <- fitted(b$fit)
  obs_size    <- b$data$size
  if (b$form == "exponential") fitted_size <- exp(fitted_size)
  tibble(metric = m, residual_cv = 100 * sd((obs_size - fitted_size) / fitted_size))
})

variability_table <- summary_df %>%
  group_by(metric) %>%
  summarise(mean_daily_cv   = weighted.mean(daily_cv, n_samples),
            min_daily_cv    = min(daily_cv),
            max_daily_cv    = max(daily_cv),
            .groups = "drop") %>%
  left_join(cv_by_phase %>% rename(daily_cv_before = before, daily_cv_after = after),
            by = "metric") %>%
  left_join(residual_cv, by = "metric") %>%
  arrange(mean_daily_cv)

cat("\n--- Variability of each measurement (%) ---\n")
print(variability_table %>% mutate(across(where(is.double), ~ round(.x, 1))))
write_csv(variability_table, "output/variability_table.csv")


# ==============================================================================
# 7. DIAGNOSTIC: EFFECT OF UNFED WEEKENDS ON INTERVAL GROWTH
# ==============================================================================
# SGR between consecutive measurements of the same individual, flagged
# by whether the interval included unfed weekend days. Used only to
# document the weekly feeding cycle (supplementary figure); headline
# growth estimates come from section 4.
# With Mon/Wed/Fri measurements, every interval that contains unfed days
# runs from a Friday to the following Monday (3 days, 2 of them unfed);
# all other intervals fall within the week (Mon-Wed, Wed-Fri).

interval_growth <- df_long %>%
  group_by(metric, id) %>%
  arrange(day, .by_group = TRUE) %>%
  mutate(day_prev  = lag(day),
         size_prev = lag(size),
         dt        = day - day_prev,
         sgr       = 100 * log(size / size_prev) / dt) %>%
  filter(!is.na(sgr)) %>%
  rowwise() %>%
  mutate(unfed_days = sum(is_unfed_day(seq(day_prev + 1, day), friday_mod)),
         interval   = if_else(unfed_days > 0, "Friday to Monday (unfed weekend)", "Within the week (fed daily)")) %>%
  ungroup()

weekend_effect <- map_dfr(metric_levels, function(m) {
  dat <- filter(interval_growth, metric == m, day_prev >= 1)  # exclude release -> first measurement
  f   <- suppressMessages(lmer(sgr ~ interval + day + (1 | id), data = dat, control = lmm_ctl))
  ct  <- summary(f)$coefficients["intervalWithin the week (fed daily)", ]
  tibble(metric = m,
         sgr_weekend_interval  = mean(dat$sgr[dat$interval == "Friday to Monday (unfed weekend)"]),
         sgr_weekday_interval  = mean(dat$sgr[dat$interval == "Within the week (fed daily)"]),
         difference_estimate   = ct[["Estimate"]],
         difference_se         = ct[["Std. Error"]],
         t_value               = ct[["t value"]])
})

cat("\n--- Interval SGR (% d-1): weekend vs weekday intervals ---\n")
print(weekend_effect %>% mutate(across(where(is.double), ~ round(.x, 2))))
write_csv(weekend_effect, "output/weekend_effect.csv")

first_interval <- interval_growth %>% filter(day_prev == 0)
cat(sprintf("\nRelease to first measurement, TBD SGR: %.1f +- %.1f %% d-1 (n = %d)\n",
            mean(first_interval$sgr[first_interval$metric == "TBD"]),
            sd(first_interval$sgr[first_interval$metric == "TBD"]),
            sum(first_interval$metric == "TBD")))


# ==============================================================================
# 8. FIGURES
# ==============================================================================

theme_paper <- theme_bw(base_size = 9) +
  theme(panel.grid = element_blank(),
        legend.position = "none",
        plot.tag = element_text(face = "bold"))

# ---- Fig. 2A: individual trajectories and population fits -----------
pop_traj <- imap_dfr(best_fits, function(b, m) {
  nd <- tibble(day = seq(0, 30, by = 0.25))
  p  <- predict(b$fit, newdata = nd, re.form = NA)
  tibble(metric = m, day = nd$day, pred = if (b$form == "exponential") exp(p) else p)
}) %>% mutate(metric = factor(metric, levels = metric_levels))

label_pos <- pop_traj %>% filter(day == 30) %>%
  mutate(pred = case_when(metric == "TMLL" ~ pred - 0.22,
                          metric == "LW"   ~ pred + 0.25,
                          TRUE             ~ pred))

fig2a <- ggplot() +
  geom_line(data = df_long, aes(day, size, group = interaction(id, metric), colour = metric),
            linewidth = 0.25, alpha = 0.35) +
  geom_point(data = df_long, aes(day, size, colour = metric), size = 0.6, alpha = 0.45) +
  geom_line(data = pop_traj, aes(day, pred, colour = metric), linewidth = 1) +
  geom_text(data = label_pos, aes(30.8, pred, label = metric, colour = metric),
            hjust = 0, size = 2.8, fontface = "bold") +
  scale_colour_manual(values = metric_cols) +
  scale_x_continuous(breaks = seq(0, 30, 5), limits = c(0, 35.5), expand = c(0.01, 0)) +
  labs(x = "Age (days after release)", y = "Size (mm)") +
  theme_paper

# ---- Fig. 2B: CDD/TBD with breakpoint -------------------------------
ratio_pop <- tibble(day = seq(0, 30, by = 0.25)) %>%
  mutate(hinge = pmax(day - bp_hat, 0))
ratio_pop$pred <- predict(bp_fit, newdata = ratio_pop, re.form = NA)

fig2b <- ggplot(ratio_dat, aes(day, 100 * cdd_ratio)) +
  annotate("rect", xmin = bp_ci[1], xmax = bp_ci[2], ymin = -Inf, ymax = Inf,
           fill = "grey85", alpha = 0.6) +
  geom_vline(xintercept = bp_hat, linetype = "dashed", colour = "grey40", linewidth = 0.4) +
  geom_line(aes(group = id), colour = metric_cols[["CDD"]], linewidth = 0.25, alpha = 0.35) +
  geom_point(colour = metric_cols[["CDD"]], size = 0.6, alpha = 0.5) +
  geom_line(data = ratio_pop, aes(day, 100 * pred), colour = "black", linewidth = 1) +
  geom_rug(data = individual_bp, aes(x = bp), inherit.aes = FALSE,
           sides = "b", colour = "grey30", length = unit(0.03, "npc")) +
  scale_x_continuous(breaks = seq(0, 30, 5), limits = c(0, 30.5), expand = c(0.01, 0)) +
  labs(x = "Age (days after release)", y = "CDD / TBD (%)") +
  theme_paper

fig2 <- fig2a + fig2b + plot_annotation(tag_levels = "A")
ggsave("plots/figure_2.png", fig2, width = 175, height = 75, units = "mm", dpi = 600)
ggsave("plots/figure_2.tiff", fig2, width = 175, height = 75, units = "mm", dpi = 600,
       compression = "lzw")

# ---- Fig. S1: interval growth, weekend vs weekday intervals ---------
figS1 <- interval_growth %>%
  filter(day_prev >= 1) %>%
  ggplot(aes(day, sgr, colour = interval)) +
  geom_hline(yintercept = 0, colour = "grey60", linewidth = 0.3) +
  geom_point(size = 0.8, alpha = 0.6) +
  geom_smooth(method = "loess", se = FALSE, linewidth = 0.7, formula = y ~ x) +
  facet_wrap(~ metric, nrow = 1) +
  scale_colour_manual(values = c("Friday to Monday (unfed weekend)" = "#7F7F7F", "Within the week (fed daily)" = "#0072B2"),
                      name = NULL) +
  labs(x = "Age at end of interval (days)", y = expression("SGR (% d"^-1*")")) +
  theme_bw(base_size = 9) +
  theme(panel.grid = element_blank(), legend.position = "bottom")
ggsave("plots/figure_S1.png", figS1, width = 175, height = 70, units = "mm", dpi = 300)

cat("\nDone. Figures in plots/, tables in output/.\n")

#        END OF SCRIPT
