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
library(car)
library(modelsummary)
library(data.table)



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

t_stat_test_2 <- (coefs["age"] + coefs["hours"]) /
  sqrt(vcov_classical["age","age"] + vcov_classical["hours","hours"] + 2*vcov_classical["age","hours"])

print(t_stat_test_2) # Because t stat is above 

### f ###

df_4 <- df_4 %>%
  mutate(part_time = if_else(hours < 40, 1, 0), 
         greedy = if_else(hours > 40, 1, 0), 
         int_female_part_time = female * part_time, 
         int_female_greedy = female*greedy)

# Estimate the new model # 

model_4_b <- lm(ln_wage ~ age + part_time + greedy + female + int_female_part_time + int_female_greedy, data = df_4)
# Print result# 

stargazer(model_4_b, type = "latex")

# Extract standard errors #
test_f <- linearHypothesis(model_4_b,
                 c("int_female_part_time = 0",
                   "int_female_greedy = 0"))

F_stat <- test_f$F[2]
p_val  <- test_f$`Pr(>F)`[2]
df1    <- test_f$Df[2]
df2    <- test_f$Res.Df[2]

latex_f <- sprintf("\\[ F(%d,%d) = %.2f, \\; p = %.3f \\]", df1, df2, F_stat, p_val)

cat(latex_f)


########## Question 6 ############

#Simulate data # 

set.seed(12345)   

# Parameters from the problem
n_industry   <- 4       # number of industries (g)
firms_per_ind <- 25     # firms per industry (f)
emp_per_firm  <- 100    # employees per firm (e)
N_firms       <- n_industry * firms_per_ind
N_employees   <- N_firms * emp_per_firm

# DGP parameters (change these if you want other cases)
beta  <- 0      # β (for Y = beta * X + phi_f + eps)
lambda <- 0     # λ (weight on U_f when building X)
omega  <- 0     # ω (weight on W_g when building X)
# note: the weight on eta is (1 - lambda - omega)

# check weights sensible
if ( (1 - lambda - omega) < 0 ) stop("lambda + omega must be <= 1")

# sds implied by distributions in the screenshot:
sd_U  <- sqrt(5)    # U_f ~ N(0,5)  (variance = 5)
sd_W  <- sqrt(5)    # W_g ~ N(0,5)
sd_eta <- sqrt(5)   # eta_e ~ N(0,5)
sd_phi <- sqrt(25)  # phi_f ~ N(0,25)
sd_eps <- sqrt(25)  # eps_e ~ N(0,25)

# Build mapping: industry -> firms -> employees
industry_id <- rep(1:n_industry, each = firms_per_ind * emp_per_firm)
# create firm ids from 1..N_firms
firm_seq <- rep(1:N_firms, each = emp_per_firm)
# alternatively create firm-to-industry mapping
firm_to_industry <- rep(1:n_industry, each = firms_per_ind)

# Generate firm-level and industry-level random variables
# U_f for each firm:
U_f <- rnorm(N_firms, mean = 0, sd = sqrt(5))

# W_g for each industry:
W_g <- rnorm(n_industry, mean = 0, sd = sqrt(5) )

# phi_f (firm fixed effect) for each firm:
phi_f <- rnorm(N_firms, mean = 0, sd = sqrt(25))

# Employee-level shocks eta_e and eps_e
eta_e <- rnorm(N_employees, mean = 0, sd = sqrt(5))
eps_e <- rnorm(N_employees, mean = 0, sd = sqrt(25))

# Now assemble the full dataset
df <- data.frame(
  emp_id = 1:N_employees,
  firm_id = firm_seq,
  industry_id = industry_id
)

# attach firm-level and industry-level variables to employees
df$U_f <- U_f[df$firm_id]                    # firm-level U
df$W_g <- W_g[df$industry_id]                # industry-level W
df$phi_f <- phi_f[df$firm_id]                # firm FE
df$eta_e <- eta_e
df$eps_e <- eps_e

# Construct X_e according to X_e = lambda * U_f + omega * W_g + (1-lambda-omega)*eta_e
df$X <- lambda * df$U_f + omega * df$W_g + (1 - lambda - omega) * df$eta_e

# Construct outcome Y_e = beta * X + phi_f + eps
df$Y <- beta * df$X + df$phi_f + df$eps_e

# quick sanity checks
cat("Total employees:", nrow(df), "\n")
cat("Total firms:", length(unique(df$firm_id)), "\n")
cat("Firms per industry (table):\n"); print(table(firm_to_industry))
cat("\nSummary of key vars:\n")
print(summary(df[, c("X","Y","phi_f","U_f","W_g","eta_e","eps_e")]))

# e.g., verify that each firm has emp_per_firm employees
print(table(df$firm_id)[1:6])   # show counts for first 6 firms

## Run OLS regression ## 

model_6 <- feols(Y ~ X, data = df, cluster = "firm_id")

# Print to overleaf # 
etable(model_6, file = "cluster_table.tex", title = "OLS result - Clustered errors")

##### b #####

#MC simulations and estimations of data # 


run_mc_compare_se <- function(lambda = 0, omega = 0, beta = 0,
                              n_industry = 4, firms_per_ind = 25, emp_per_firm = 100,
                              reps = 1000, seed = 123, show_progress = TRUE) {
  set.seed(seed)
  if (lambda + omega > 1) stop("lambda + omega must be <= 1")
  
  N_firms <- n_industry * firms_per_ind
  N_employees <- N_firms * emp_per_firm
  
  results <- tibble::tibble(
    rep = seq_len(reps),
    est = as.numeric(NA),          # coefficient estimate (same for both)
    se_classic = as.numeric(NA),   # classical (lm) se
    p_classic  = as.numeric(NA),   # p-value using classical se
    se_cluster = as.numeric(NA),   # cluster-robust se (firm level)
    p_cluster  = as.numeric(NA)    # p-value using clustered se
  )
  
  if (show_progress) pb <- txtProgressBar(min = 0, max = reps, style = 3)
  
  for (r in seq_len(reps)) {
    # ---------- DGP ----------
    U_f   <- rnorm(N_firms, 0, sqrt(5))
    W_g   <- rnorm(n_industry, 0, sqrt(5))
    phi_f <- rnorm(N_firms, 0, sqrt(25))
    
    eta_e <- rnorm(N_employees, 0, sqrt(5))
    eps_e <- rnorm(N_employees, 0, sqrt(25))
    
    firm_id     <- rep(1:N_firms, each = emp_per_firm)
    industry_id <- rep(rep(1:n_industry, each = firms_per_ind), each = emp_per_firm)
    
    df <- tibble::tibble(
      emp_id = seq_len(N_employees),
      firm_id = firm_id,
      industry_id = industry_id
    ) %>%
      dplyr::mutate(
        U_f = U_f[firm_id],
        W_g = W_g[industry_id],
        phi_f = phi_f[firm_id],
        eta_e = eta_e,
        eps_e = eps_e,
        X = lambda * U_f + omega * W_g + (1 - lambda - omega) * eta_e,
        Y = beta * X + phi_f + eps_e
      )
    
    # ---------- Estimation ----------
    # classical OLS (lm) for i.i.d. SE
    fit_lm <- tryCatch(lm(Y ~ X, data = df), error = function(e) e)
    if (inherits(fit_lm, "error")) {
      # leave NA and continue; print a message for the first/last iteration
      if (r == 1 || r == reps) message("lm error at rep ", r, ": ", fit_lm$message)
      if (show_progress) setTxtProgressBar(pb, r)
      next
    }
    
    # clustered estimator using fixest (cluster at firm level)
    fit_fe <- tryCatch(feols(Y ~ X, data = df, cluster = "firm_id"),
                       error = function(e) e)
    if (inherits(fit_fe, "error")) {
      if (r == 1 || r == reps) message("feols error at rep ", r, ": ", fit_fe$message)
      if (show_progress) setTxtProgressBar(pb, r)
      next
    }
    
    # extract coefficient (should be identical for lm and feols without FE)
    coef_lm <- coef(fit_lm)
    if (!"X" %in% names(coef_lm)) {
      if (r == 1 || r == reps) message("No coef X in lm at rep ", r, ". Names: ", paste(names(coef_lm), collapse = ", "))
      if (show_progress) setTxtProgressBar(pb, r)
      next
    }
    est <- as.numeric(coef_lm["X"])
    
    # classical se from lm summary
    lm_sum <- summary(fit_lm)
    se_classic <- as.numeric(lm_sum$coefficients["X", "Std. Error"])
    t_classic  <- est / se_classic
    p_classic  <- 2 * pnorm(-abs(t_classic))   # normal approx
    
    # clustered se from fixest vcov
    vc_cluster <- tryCatch(vcov(fit_fe), error = function(e) e)
    if (inherits(vc_cluster, "error")) {
      if (r == 1 || r == reps) message("vcov error at rep ", r, ": ", vc_cluster$message)
      if (show_progress) setTxtProgressBar(pb, r)
      next
    }
    if (!("X" %in% rownames(vc_cluster))) {
      if (r == 1 || r == reps) message("vcov missing 'X' row at rep ", r, ". Row names: ", paste(rownames(vc_cluster), collapse = ", "))
      if (show_progress) setTxtProgressBar(pb, r)
      next
    }
    se_cluster <- sqrt(as.numeric(vc_cluster["X","X"]))
    t_cluster  <- est / se_cluster
    p_cluster  <- 2 * pnorm(-abs(t_cluster))
    
    # store
    results[r, "est"] <- est
    results[r, "se_classic"] <- se_classic
    results[r, "p_classic"]  <- p_classic
    results[r, "se_cluster"] <- se_cluster
    results[r, "p_cluster"]  <- p_cluster
    
    if (show_progress) setTxtProgressBar(pb, r)
  } # end reps
  
  if (show_progress) close(pb)
  
  # ---------- Summaries ----------
  # valid counts
  n_valid_classic <- sum(is.finite(results$se_classic))
  n_valid_cluster <- sum(is.finite(results$se_cluster))
  n_valid_est     <- sum(is.finite(results$est))
  
  # compute rejection rates at alpha = 0.05 (two-sided)
  rej_classic <- mean(results$p_classic < 0.05, na.rm = TRUE)
  rej_cluster <- mean(results$p_cluster < 0.05, na.rm = TRUE)
  
  # other MC stats for est
  mean_est <- mean(results$est, na.rm = TRUE)
  bias     <- mean_est - beta
  emp_sd   <- sd(results$est, na.rm = TRUE)
  mean_se_classic <- mean(results$se_classic, na.rm = TRUE)
  mean_se_cluster <- mean(results$se_cluster, na.rm = TRUE)
  
  summary_tbl <- tibble::tibble(
    lambda = lambda,
    omega  = omega,
    beta   = beta,
    n_industry = n_industry,
    firms_per_ind = firms_per_ind,
    emp_per_firm = emp_per_firm,
    reps = reps,
    n_valid_est = n_valid_est,
    n_valid_classic = n_valid_classic,
    n_valid_cluster = n_valid_cluster,
    mean_est = mean_est,
    bias = bias,
    emp_sd = emp_sd,
    mean_se_classic = mean_se_classic,
    mean_se_cluster = mean_se_cluster,
    rej_classic = rej_classic,
    rej_cluster = rej_cluster
  )
  
  list(results = results, summary = summary_tbl)
}


result_6_b <- run_mc_compare_se(lambda = 0, omega = 0, beta = 0,
                                n_industry = 4, firms_per_ind = 25, emp_per_firm = 100,
                                reps = 1000, seed = 2026, show_progress = TRUE)

# View summary:
print(result_6_b$summary)

# Print in overleaf # 
summary_6_b <- result_6_b$summary
summary_6_b <- summary_6_b %>% 
  select(lambda, omega, mean_est, mean_se_classic, mean_se_cluster, rej_classic, rej_cluster)

print(kable(summary_6_b,
            format = "latex",
            booktabs = TRUE,
            digits = 3,
            caption = "Monte Carlo simulation - Result"))


#### c #### 

#Change lambda and omega # 


result_6_c <- run_mc_compare_se(lambda = 0.2, omega = 0.7, beta = 0,
                                n_industry = 4, firms_per_ind = 25, emp_per_firm = 100,
                                reps = 1000, seed = 2026, show_progress = TRUE)

#Summarise and present the results# 

summary_6_c <- result_6_c$summary
summary_6_c <- summary_6_c %>% 
  select(lambda, omega, mean_est, mean_se_classic, mean_se_cluster, rej_classic, rej_cluster)

print(kable(summary_6_c,
            format = "latex",
            booktabs = TRUE,
            digits = 3,
            caption = "Monte Carlo simulation - Result 2"))
      