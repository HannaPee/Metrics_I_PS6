###########################
########## PS6 ############
###########################

#Loading packages # 

library(dplyr)
library(tidyr)
library(tidyverse)
library(ggplot2)
library(readr)
library(xtable)
library(fixest)
library(broom)
library(sandwich)
library(purrr)
library(knitr)
library(patchwork)
library(readxl)
library(stargazer)



####### Question 1 ########¨

# b # 


# Generate Monte Carlo data # 

n_sim = 10000
n_obs_1 = 100 
n_obs_2 = 2000


gen_data_1 <- function(n_sim, n) {

  data <- map(1:n_sim, ~ {
    tibble(
    x = rchisq(n, df = 2) - 2,
    e = rchisq(n, df = 2) - 2,
    y = 0 + 0 * x + e   
   )
  })
  
  return(data)

}

#Gen data for b # 

data_100 <- gen_data_1(n_sim = n_sim, n = n_obs_1)
data_2000 <- gen_data_1(n_sim = n_sim, n = n_obs_2)

#Create functin for estimating model# 

run_estim <- function(n, data) {
  
  
  # Estimate models# 
 
  model <- map(data, ~ lm(y ~ x, data = .x))
  
  
  # Estimate standard errors # 

  vcov_classical <- map(model, vcov)
  vcov_HC1 <-  map(model, ~ vcovHC(.x, type = "HC1"))
  vcov_HC2 <-  map(model, ~ vcovHC(.x, type = "HC2"))
  vcov_HC3 <- map(model, ~ vcovHC(.x, type = "HC3"))
  
  # Helper to extract se for coefficient "x" from a vcov matrix safely
  extract_se <- function(vcov_mat, coef_name = "x") {
    # return NA if something odd happens (robust to missing names)
    if (is.null(rownames(vcov_mat)) || !(coef_name %in% rownames(vcov_mat))) return(NA_real_)
    sqrt(vcov_mat[coef_name, coef_name])
  }
  
  # Extract beta1_hat # 
  beta1_hat <- map_dbl(model, ~ coef(.x)["x"])
  
  # 5. Extract  SE(beta1)
  se_classical <- map_dbl(vcov_classical, extract_se)
  se_HC1       <- map_dbl(vcov_HC1, extract_se)
  se_HC2       <- map_dbl(vcov_HC2, extract_se)
  se_HC3       <- map_dbl(vcov_HC3, extract_se)
  
  
  # Compute t-statistic # 
  t_classical <- beta1_hat / se_classical
  t_HC1       <- beta1_hat / se_HC1
  t_HC2       <- beta1_hat / se_HC2
  t_HC3       <- beta1_hat / se_HC3
  
  # bundle results in a tibble for easy summaries / plotting #
  results <- tibble(
    beta1_hat = beta1_hat,
    se_classical = se_classical,
    se_HC1 = se_HC1,
    se_HC2 = se_HC2,
    se_HC3 = se_HC3,
    t_classical = t_classical,
    t_HC1 = t_HC1,
    t_HC2 = t_HC2,
    t_HC3 = t_HC3
  )
  
  out <- list(
    n = n,
    models = models,
    vcov = list(classical = vcov_classical, HC1 = vcov_HC1, HC2 = vcov_HC2, HC3 = vcov_HC3),
    results = results
  )
  
  return(out)
}

# Function for summarizing the results# 

summarize_for <- function(mc) {
  df <- mc$results
  tibble(
    n = mc$n,
    mean_beta1 = mean(df$beta1_hat),
    sd_beta1   = sd(df$beta1_hat),
    se_classical = mean(df$se_classical),
    se_HCI = mean(df$se_HC1), 
    rej_5_classical = mean(abs(df$t_classical) > 1.96, na.rm = TRUE),
    rej_5_HC1       = mean(abs(df$t_HC1) > 1.96, na.rm = TRUE)
  )

}


#Results# 

result_b_100 <- run_estim(n = n_obs_1, data = data_100)
result_b_2000 <- run_estim(n = n_obs_2, data = data_2000)

sum_b_100 <- summarize_for(result_b_100)
sum_b_2000 <- summarize_for(result_b_2000)

#Present the result# 


print(kable(sum_b_100,
            format = "latex",
            booktabs = TRUE,
            digits = 4,
            caption = paste0("Results")))

print(kable(sum_b_2000,
            format = "latex",
            booktabs = TRUE,
            digits = 4,
            caption = paste0("Results")))

#Plot the distribution of beta_1_hat using classical errors# 

plot_1 <- ggplot(result_b_100$results, aes(x = beta1_hat)) +
  geom_histogram(aes(y = after_stat(density)),
                 bins = 60,
                 fill = "steelblue",
                 alpha = 0.6) +
  geom_density(color = "darkblue", linewidth = 1) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  labs(x = expression(hat(beta)[1]),
       y = "Density") +
  theme_minimal()

plot_2 <- ggplot(result_b_2000$results, aes(x = beta1_hat)) +
  geom_histogram(aes(y = after_stat(density)),
                 bins = 60,
                 fill = "steelblue",
                 alpha = 0.6) +
  geom_density(color = "darkblue", linewidth = 1) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  labs(x = expression(hat(beta)[1]),
       y = "Density") +
  theme_minimal()

print(plot_2)

comined_beta_plot <- plot_1 + plot_2

ggsave(
  filename = "beta_1_hat_dist.pdf",
  plot     = comined_beta_plot,
  width    = 8,
  height   = 6
)


# Plot t-statistcs whith classical SE # 

plot_3 <- ggplot(result_b_100$results, aes(x = t_classical)) +
  geom_histogram(aes(y = after_stat(density)),
                 bins = 60,
                 fill = "seagreen",
                 alpha = 0.6) +
  geom_density(color = "darkgreen", linewidth = 1) +
  stat_function(fun = dnorm,
                args = list(mean = 0, sd = 1),
                color = "red",
                linewidth = 1,
                linetype = "dashed") +
  labs(x = "t-statistic",
       y = "Density") +
  theme_minimal()

print(plot_3)


plot_4 <- ggplot(result_b_2000$results, aes(x = t_classical)) +
  geom_histogram(aes(y = after_stat(density)),
                 bins = 60,
                 fill = "seagreen",
                 alpha = 0.6) +
  geom_density(color = "darkgreen", linewidth = 1) +
  stat_function(fun = dnorm,
                args = list(mean = 0, sd = 1),
                color = "red",
                linewidth = 1,
                linetype = "dashed") +
  labs(x = "t-statistic",
       y = "Density") +
  theme_minimal()

print(plot_4)

combined_plot_t_stats <- plot_3 + plot_4

ggsave(
  filename = "t_stat_dist.pdf",
  plot     = combined_plot_t_stats,
  width    = 8,
  height   = 6
)



### d #### 

# Generate new data # 

n_sim = 10000
n_obs = 100 

# i # 

gen_data_i <- function(n_sim, n) {
  
  data <- map(1:n_sim, ~ {
    tibble(
      x = rnorm(n, mean = 0, sd=2),
      e = rchisq(n, df = 2) - 2,
      y = 0 + 0 * x + e   
    )
  })
  
  return(data)
  
}

#ii#

gen_data_ii <- function(n_sim, n) {
  
  data <- map(1:n_sim, ~ {
    tibble(
      x = rchisq(n, df = 2) - 2,
      e = rnorm(n, mean = 0, sd = 2),
      y = 0 + 0 * x + e   
    )
  })
  
  return(data)
  
}

 #iii#
gen_data_iii <- function(n_sim, n) {
  
  data <- map(1:n_sim, ~ {
    tibble(
      x = rnorm(n, mean = 0, sd = 2),
      e = rnorm(n, mean = 0, sd = 2),
      y = 0 + 0 * x + e   
    )
  })
  
  return(data)
  
}

# Generate the data # 

data_i <- gen_data_i(n_sim, n_obs)Mean
data_ii <- gen_data_ii(n_sim, n_obs)
data_iii <- gen_data_iii(n_sim, n_obs)

# Run estimations # 

results_c_i <- run_estim(n_obs, data_i)
results_c_ii <- run_estim(n_obs, data_ii)
results_c_iii <- run_estim(n_obs, data_iii)

# Summary stats# 

sum_c_i <- summarize_for(results_c_i)
sum_c_ii <- summarize_for(results_c_ii)
sum_c_iii <- summarize_for(results_c_iii)

# Print in one table # 

sum_c_combined <- rbind(sum_c_i,sum_c_ii, sum_c_iii)
sum_c_combined <- sum_c_combined%>%
  mutate(Case = c("i", "ii", "iii"))

print(kable(sum_c_combined,
            format = "latex",
            booktabs = TRUE,
            digits = 4,
            caption = paste0("Results")))


# Plot the distribution of beta_hat# 

datalist <- list(results_c_i, results_c_ii, results_c_iii)

data_plotting <- imap_dfr(datalist, ~ {
  # .x is a single mc object, .y is its name
  tibble(
    design    = .y,
    n         = .x$n,
    beta1_hat = .x$results$beta1_hat,
    t_classical = .x$results$t_classical, 
    t_HC1 = .x$results$t_HC1
  )
}) %>%
  mutate(design = as.factor(design),
         n = as.factor(n))

# Plot betas#

p_beta <- ggplot(data_plotting, aes(x = beta1_hat, color = design, fill = design)) +
  geom_density(alpha = 0.25, linewidth = 0.8) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "red") +
  labs(x = expression(hat(beta)[1]),
       y = "Density",
       color = "Case", fill = "Case") +
  theme_minimal()

print(p_beta)

ggsave(
  filename = "beta_1_dist_c.pdf",
  plot     = p_beta,
  width    = 8,
  height   = 6
)


# Plot t_statistics with classical errors#


p_t <- ggplot(data_plotting, aes(x = t_classical, color = design, fill = design)) +
  geom_density(alpha = 0.25, linewidth = 0.8) +
  stat_function(fun = dnorm, args = list(mean = 0, sd = 1),
                linetype = "dashed", color = "black", inherit.aes = FALSE) +
  labs(x = "t",
       y = "Density",
       color = "Case", fill = "Case") +
  theme_minimal()

print(p_t)

ggsave(
  filename = "t_stats_dist_c.pdf",
  plot     = p_t,
  width    = 8,
  height   = 6
)


# Plot t_statistics with HC1 errors#
p_t_HC1 <- ggplot(data_plotting, aes(x = t_HC1, color = design, fill = design)) +
  geom_density(alpha = 0.25, linewidth = 0.8) +
  stat_function(fun = dnorm, args = list(mean = 0, sd = 1),
                linetype = "dashed", color = "black", inherit.aes = FALSE) +
  labs(x = "t",
       y = "Density",
       color = "Case", fill = "Case") +
  theme_minimal()

print(p_t_HC1)

ggsave(
  filename = "t_statsHC1__dist_c.pdf",
  plot     = p_t_HC1,
  width    = 8,
  height   = 6
)

##### Question 4 ######### 

# Import data # 
cps09mar <- read_excel("cps09mar.xlsx")

summary(cps09mar)

# Cleaning data # 

df_4 <- cps09mar %>%
  mutate(wage = earnings/(hours*week), 
         ln_wage = log(wage)) %>%
  filter(wage > 1)

summary(df_4)


# Plot the data # 

plot_y <- ggplot(df_4, aes(x = ln_wage)) +
  geom_density(fill = "steelblue", alpha = 0.5) +
  labs(
    x = "Log Wage",
    y = "Density"
  ) +
  theme_minimal()

plot(plot_y)

ggsave(
  filename = "dist_ln_wage.pdf",
  plot     = plot_y,
  width    = 8,
  height   = 6
)

# Estimation #   

model_4 <- lm(ln_wage ~ age + hours + female, data = df_4)
summary(model_4)

#Print results# 
stargazer(model_4, type = "latex")

## Plot confidence intervals with different SE:s ##

# Extract standard errors # 
coefs <- coef(model_4)
vcov_classical <- vcov(model_4)
vcov_hc1       <- vcovHC(model_4, type = "HC1")

# build a table with estimate + both SEs and CIs
coef_tbl <- tibble(
  term = names(coefs),
  estimate = as.numeric(coefs),
  se_classical = sqrt(diag(vcov_classical)),
  se_HC1       = sqrt(diag(vcov_hc1))
) %>%
  mutate(
    ci_lo_classical = estimate - 1.96 * se_classical,
    ci_hi_classical = estimate + 1.96 * se_classical,
    ci_lo_HC1 = estimate - 1.96 * se_HC1,
    ci_hi_HC1 = estimate + 1.96 * se_HC1
  ) %>%
  # optional: remove intercept for plotting
  filter(term != "(Intercept)")

# Print coef table # 

print(kable(coef_tbl,
            format = "latex",
            booktabs = TRUE,
            digits = 3,
            caption = paste0("Coefficients")))

# Plot one graph# 
plot_data <- coef_tbl %>%
  pivot_longer(cols = c(ci_lo_classical, ci_hi_classical, ci_lo_HC1, ci_hi_HC1), 
                   names_to = c(".value", "method"), names_pattern = "(lo|hi)_(.*)")

# 4) Plot: point + error bars, dodged so the two methods sit side-by-side
plot_4 <-ggplot(plot_data, aes(x = method, y = estimate)) + geom_point(size = 3) + 
  geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  facet_wrap(~ term, scales = "free_y") +
  labs(
         x = "",
         y = "Coefficient estimate"
       ) +   theme_minimal()


# show plot
print(plot_4)

ggsave(
  filename = "ci_combined.pdf",
  plot     = plot_4,
  width    = 8,
  height   = 6
)

#Plot 3 separate graphs# 

dir.create("coef_plots", showWarnings = FALSE)  # folder for plots

for (i in seq_len(nrow(coef_tbl))) {
  row <- coef_tbl[i, ]
  term_i <- row$term
  
  plot_df <- tibble(
    method = c("Classical", "HC1"),
    estimate = c(row$estimate, row$estimate),
    lo = c(row$ci_lo_classical, row$ci_lo_HC1),
    hi = c(row$ci_hi_classical, row$ci_hi_HC1)
  )
  
  p <- ggplot(plot_df, aes(x = method, y = estimate, color = method)) +
    geom_point(size = 3) +
    geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.15, size = 1) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
    labs(x = "", y = "Estimate (ln wage)") +
    theme_minimal() +
    theme(legend.position = "none")
  
    print(p)
  
  fname <- file.path("coef_plots", paste0("coef_plot_", term_i, ".pdf"))
  ggsave(filename = fname, plot = p, width = 5, height = 4)
}

# Plot the residuals over age, female and hours #

df_4$residuals <- resid(model_4)


plot_resid_data <- df_4 %>%
  select(residuals, age, hours, female) %>%
  pivot_longer(
    cols = c(age, hours, female),
    names_to = "variable",
    values_to = "value"
  )

plot_resid<- ggplot(plot_resid_data, aes(x = value, y = residuals)) +
  geom_point(alpha = 0.4) +
  geom_smooth(method = "loess", se = FALSE, color = "red") +
  facet_wrap(~ variable, scales = "free_x") +
  labs(
    x = "",
    y = "Residuals"
  ) +
  theme_minimal()

print(plot_resid)

ggsave(
  filename = "residuals.pdf",
  plot     = plot_resid,
  width    = 8,
  height   = 6
)

# c - T-test for female # 
t_test_beta_3 <- coef_tbl %>%
  filter(term == "female") %>%
  mutate(
    t_beta_3_classical = (estimate - (-0.25)) / se_classical,
    t_beta_3_HC1 = (estimate - (-0.25)) / se_HC1,
    p_value_classical = 2 * pt(abs(t_beta_3_classical),
                               df = df.residual(model_4),
                               lower.tail = FALSE),
    p_value_HC1 = 2 * pt(abs(t_beta_3_HC1),
                         df = df.residual(model_4),
                         lower.tail = FALSE)
  )%>%
  select(term, t_beta_3_classical, t_beta_3_HC1, p_value_classical, p_value_HC1)

print(kable(t_test_beta_3,
            format = "latex",
            booktabs = TRUE,
            digits = 3,
            caption = "t-test result"))

# d - Print the variance -covariance matrix #

print(vcov_classical)

print(kable(vcov_classical,
            format = "latex",
            booktabs = TRUE,
            digits = 10,
            caption = "Variance -covariance matrix"))

# e - t test for beta_1 = beta_2#








      