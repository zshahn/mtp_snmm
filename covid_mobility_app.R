# Packages
library(data.table)
library(dplyr)
library(lubridate)
library(zoo)

# ---- 1) Download raw data ----
# Google mobility: global CSV (contains US rows with census_fips_code)
mob_url <- "https://www.gstatic.com/covid19/mobility/Global_Mobility_Report.csv"
mob <- fread(mob_url)

# NYT county cases
nyt_url <- "https://raw.githubusercontent.com/nytimes/covid-19-data/master/us-counties.csv"
nyt <- fread(nyt_url)

# ---- 2) Keep US counties and select needed columns ----
mob_us <- mob %>%
  filter(country_region_code == "US",
         !is.na(census_fips_code),                 # county rows
         !is.na(date)) %>%
  transmute(
    fips = sprintf("%05d", as.integer(census_fips_code)),
    date = as.Date(date),
    workplaces = workplaces_percent_change_from_baseline,
    retail = retail_and_recreation_percent_change_from_baseline,
    grocery = grocery_and_pharmacy_percent_change_from_baseline,
    parks = parks_percent_change_from_baseline,
    transit = transit_stations_percent_change_from_baseline,
    residential = residential_percent_change_from_baseline
  )

# Clean NYT (cumulative counts)
nyt_clean <- nyt %>%
  mutate(fips = sprintf("%05s", fips), date = as.Date(date)) %>%
  filter(!is.na(fips)) %>%
  arrange(fips, date) %>%
  group_by(fips) %>%
  mutate(
    # daily new cases; guard against small negative corrections
    new_cases = pmax(cases - lag(cases, default = first(cases)), 0),
    pop = NA_real_  # fill later if you add a population table
  ) %>% ungroup()

# If you want per-100k, merge a county population file here (e.g., Census ACS).
# For a quick demo, approximate per-100k by within-county standardization or skip scaling.

# ---- 3) Join mobility to cases by fips & date ----
dt <- mob_us %>%
  inner_join(nyt_clean, by = c("fips", "date"))

# ---- 4) Construct features for a point-exposure episode at t0 ----
make_episodes <- function(df,
                          t0_start = as.Date("2020-06-01"),
                          t0_end   = as.Date("2021-03-31"),
                          lag_hist_days     = 14,
                          outcome_start_lag = 14,
                          outcome_window    = 7) {
  df %>%
    arrange(fips, date) %>%
    group_by(fips) %>%
    mutate(
      # A0: 7-day avg workplaces ending at t0 (right-aligned)
      A0 = zoo::rollapply(workplaces, 7, mean, align = "right", fill = NA_real_),
      
      # Optional history features (7-day avgs of other categories)
      R7_retail = zoo::rollapply(retail, 7, mean, align = "right", fill = NA_real_),
      R7_grocery = zoo::rollapply(grocery, 7, mean, align = "right", fill = NA_real_),
      R7_transit = zoo::rollapply(transit, 7, mean, align = "right", fill = NA_real_),
      R7_resid   = zoo::rollapply(residential, 7, mean, align = "right", fill = NA_real_),
      
      # History: past 14 days cases (sum/mean)
      lag_cases_14      = zoo::rollapply(new_cases, lag_hist_days, sum,  align = "right", fill = NA_real_),
      lag_cases_14_mean = zoo::rollapply(new_cases, lag_hist_days, mean, align = "right", fill = NA_real_),
      
      # Outcome Y: sum of new cases in [t0+outcome_start_lag, ..., +outcome_start_lag+outcome_window-1]
      # Do it vectorized: first a forward (left-aligned) rolling 7-day sum, then shift forward by the lag.
      fwd_sum_outcome = zoo::rollapply(new_cases, outcome_window, sum, align = "left", fill = NA_real_),
      Y = dplyr::lead(fwd_sum_outcome, outcome_start_lag)
    ) %>%
    ungroup() %>%
    # keep the decision-time window and rows with all needed variables
    filter(date >= t0_start, date <= t0_end,
           !is.na(A0), !is.na(Y), !is.na(lag_cases_14)) %>%
    mutate(t0 = date)
}
episodes <- make_episodes(dt)

#add some covariates and convert outcome to per 100k
library(haven)
ses <- read_dta("~/snmm/medicaid/master_area_dataset_4MEPS.dta")
sesdic <- labelled::generate_dictionary(ses)
names(ses)[1] = 'fips'
ses = ses[ses$year==2017,]
ses = ses[!is.na(ses$co_tot_pop),]
episodes = merge(episodes,ses[,c('fips','co_tot_pop','co_ruc_code13','co_per_black','co_per_hisp','co_rep_voteshr_pres','co_per_poor_saipe','co_unemploymentrate','co_sahie1864_inspct_all','co_land_area','co_pop_den','co_pcpi','st_mw_fred','st_mw_super','st_h_rep_prct')])
episodes$log_pop = log(episodes$co_tot_pop)
episodes$cases = episodes$Y
episodes$Y = episodes$cases/(episodes$co_tot_pop/100000)
# episodes now has: fips, t0 (decision day), exposure A0 (workplaces 7-day avg),
# history H0 (lag_cases_14, other mobility avgs), and outcome Y (future 7-day cases after 14d lag).

# ---- deps ----
suppressPackageStartupMessages({
  library(data.table)
  library(Matrix)
  library(xgboost)
})

# ---- helpers ----
dense_model_matrix <- function(df) {
  # quick sparse->dense builder (xgboost wants numeric matrix)
  as.matrix(data.frame(lapply(df, function(v) as.numeric(v))))
}

fit_xgb_reg <- function(X, y, nrounds = 300, max_depth = 5, eta = 0.05, subsample = 0.8, colsample_bytree = 0.8) {
  dtr <- xgb.DMatrix(X, label = y)
  params <- list(
    objective = "reg:squarederror",
    max_depth = max_depth,
    eta = eta,
    subsample = subsample,
    colsample_bytree = colsample_bytree,
    nthread = 0,
    verbosity = 0
  )
  xgb.train(params = params, data = dtr, nrounds = nrounds, watchlist = list(), verbose = 0)
}

predict_xgb <- function(model, Xnew) {
  predict(model, xgb.DMatrix(Xnew))
}

# Normal density ratio p(a - delta | H)/p(a | H) using mean m(H), sd sigma
dens_ratio_normal_shift <- function(a, m, sigma, delta) {
  dnorm(a - delta, mean = m, sd = sigma) / dnorm(a, mean = m, sd = sigma)
}

# ---- main estimator ----
# data: data.frame with columns Y, A, and H_cols
# H_cols: character vector of history feature names (numeric)
# s_fun: function(H_df_row) -> numeric vector s(H), low-dim (e.g., c(1, H$lag_cases_14))
# delta: shift (e.g., -5 for -5pp)
# K: number of cross-fitting folds
mtp_point_shift_cf <- function(data, H_cols, s_fun,
                               delta = -5, K = 2,
                               ratio_cap = 50,
                               xgb_mu = list(nrounds=400, max_depth=6, eta=0.05),
                               xgb_m  = list(nrounds=300, max_depth=5, eta=0.05),
                               seed = 1) {
  set.seed(seed)
  stopifnot(all(c("Y","A") %in% names(data)))
  df <- data[complete.cases(data[, c("Y","A", H_cols), drop=FALSE]), c("Y","A", H_cols), drop=FALSE]
  n  <- nrow(df)
  folds <- sample(rep_len(1:K, n))
  
  # storage
  s_mat <- NULL            # n x p_s (we'll fill per fold)
  qmt   <- NULL            # n x p_s  (q - tilde q)
  mu_AL <- mu_g <- e <- rep(NA_real_, n)
  
  # linear system M ψ = rhs  (p_s x p_s)
  first_i <- s_fun(df[1, H_cols, drop=FALSE])     # get dimension of s(H)
  p_s <- length(first_i)
  M   <- matrix(0, p_s, p_s)
  rhs <- rep(0, p_s)
  
  # Pre-allocate
  s_mat <- matrix(NA_real_, n, p_s)
  qmt   <- matrix(NA_real_, n, p_s)
  
  for (k in 1:K) {
    idx_tr <- which(folds != k); idx_te <- which(folds == k)
    dtr <- df[idx_tr, , drop=FALSE]
    dte <- df[idx_te, , drop=FALSE]
    
    # ----- Fit nuisances on train -----
    Xtr_mu <- dense_model_matrix(dtr[, c("A", H_cols), drop=FALSE])
    mu_fit <- fit_xgb_reg(Xtr_mu, dtr$Y,
                          nrounds = xgb_mu$nrounds, max_depth = xgb_mu$max_depth, eta = xgb_mu$eta)
    
    # mean model for A|H
    Xtr_m <- dense_model_matrix(dtr[, H_cols, drop=FALSE])
    m_fit <- fit_xgb_reg(Xtr_m, dtr$A,
                         nrounds = xgb_m$nrounds, max_depth = xgb_m$max_depth, eta = xgb_m$eta)
    
    # estimate residual sd on train for A|H
    m_tr  <- predict_xgb(m_fit, Xtr_m)
    sigA  <- sd(dtr$A - m_tr)
    
    # ----- Predict on test -----
    Xte_mu_now <- dense_model_matrix(dte[, c("A", H_cols), drop=FALSE])
    mu_now <- predict_xgb(mu_fit, Xte_mu_now)
    # shift A by delta for mu_g
    A_shift <- dte$A + delta
    Xte_mu_shift <- dense_model_matrix(data.frame(A = A_shift, dte[, H_cols, drop=FALSE]))
    mu_shift <- predict_xgb(mu_fit, Xte_mu_shift)
    
    # m(H) and density ratios
    Xte_m <- dense_model_matrix(dte[, H_cols, drop=FALSE])
    m_te  <- predict_xgb(m_fit, Xte_m)
    ratio <- dens_ratio_normal_shift(dte$A, m_te, sigA, delta)
    ratio <- pmin(pmax(ratio, 1/ratio_cap), ratio_cap)
    
    # build s(H) and pullback tilde{q} = s(H) * ratio  (since q(H,A)=s(H))
    s_list <- lapply(seq_len(nrow(dte)), function(i) s_fun(dte[i, H_cols, drop=FALSE]))
    S      <- do.call(rbind, s_list)              # n_te x p_s
    tQ     <- S * ratio                           # elementwise multiply by ratio
    
    # residuals e = Y - mu(H,A)
    e_te <- dte$Y - mu_now
    
    # accumulate system and stash per-obs pieces
    for (j in seq_len(nrow(dte))) {
      i <- idx_te[j]
      s_i <- S[j, ]
      q_i <- s_i
      tq_i <- tQ[j, ]
      
      s_mat[i, ] <- s_i
      qmt[i, ]   <- (q_i - tq_i)
      mu_AL[i]   <- mu_now[j]
      mu_g[i]    <- mu_shift[j]
      e[i]       <- e_te[j]
      
      M   <- M   + as.numeric(delta) * (s_i %o% s_i)
      rhs <- rhs - ( s_i * (mu_now[j] - mu_shift[j]) + (q_i - tq_i) * e_te[j] )
    }
  }
  
  # ----- Solve for psi -----
  psi_hat <- as.numeric(solve(M, rhs))
  names(psi_hat) <- paste0("psi_", seq_len(p_s))
  
  # ----- IF-based sandwich (recommended, stable for intercept) -----
  # Jacobian G = (delta / n) * sum s s^T
  G_hat <- (as.numeric(delta) / n) * crossprod(s_mat, s_mat)
  Vinv  <- solve(G_hat)
  
  # per-observation scores φ_i(ψ̂) = (q-~q)e + s{ μ-μ_g + δ s^T ψ̂ }
  phi <- qmt * e + s_mat * ((mu_AL - mu_g) + as.numeric(delta) * as.numeric(s_mat %*% psi_hat))
  
  # Influence functions IF_i = G^{-1} φ_i
  IF  <- phi %*% t(Vinv)
  IFc <- IF - matrix(colMeans(IF), n, ncol(IF), byrow = TRUE)
  Sigma_IF <- crossprod(IFc) / n
  V_hat <- Sigma_IF / n
  se <- sqrt(diag(V_hat))
  
  list(
    psi = psi_hat,
    se  = se,
    ci  = cbind(lower = psi_hat - 1.96*se, upper = psi_hat + 1.96*se),
    delta = delta,
    dim_s = p_s,
    G = G_hat, V = V_hat
  )
}

## Suppose your episodes df has:
##   Y   : future 7-day cases (or per-100k) after a 14-day lag
##   A   : 7-day avg workplaces mobility ending at t0
##   H_cols: c("lag_cases_14", "R7_retail", "R7_transit", "R7_resid", "dow", "time_spline1", ...)

episodes$lag_cases_14_100k = episodes$lag_cases_14/(episodes$co_tot_pop/100000)
episodes$log_lag_cases = log(episodes$lag_cases_14_100k+1)
H_cols <- c("log_lag_cases",'co_ruc_code13','log_pop','co_per_black','co_per_hisp','co_rep_voteshr_pres','co_per_poor_saipe','co_unemploymentrate','co_sahie1864_inspct_all','co_land_area','co_pop_den','co_pcpi')#,'st_mw_fred','st_mw_super','st_h_rep_prct')

library(dplyr)

target <- as.Date("2020-06-01")
win    <- 14L  # days

episodes_nearest <- episodes %>%
  mutate(dist_days = abs(as.integer(t0 - target))) %>%   # distance to target in days
  filter(dist_days <= win) %>%                           # keep only within ±14 days
  group_by(fips) %>%
  slice_min(order_by = dist_days, with_ties = FALSE) %>% # pick nearest; ties broken arbitrarily
  ungroup() %>%
  select(-dist_days)

# quick audit:
kept    <- n_distinct(episodes_nearest$fips)
total   <- n_distinct(episodes$fips)
dropped <- total - kept
cat(sprintf("Kept %d counties; dropped %d with no t0 within ±%d days of %s.\n",
            kept, dropped, win, target))

episodes = episodes_nearest

## Choose a small blip/test basis s(H). Keep it low-dim!
s_fun <- function(Hrow) {
  c(1, as.numeric(Hrow$co_ruc_code13))  # intercept + disease activity
  # You can expand modestly, e.g., c(1, Hrow$lag_cases_14, Hrow$R7_transit)
}

## Fit with a -5 pp shift, 5-fold cross-fitting
fit <- mtp_point_shift_cf(
  data = episodes %>% dplyr::rename(A = A0),   # rename exposure to A
  H_cols = H_cols,
  s_fun = s_fun,
  delta = -5,
  K = 5
)

fit$psi
fit$se
fit$ci

(-5)*sum(fit$psi * c(1,quantile(episodes$log_lag_cases,.75)))
pred = function(x){
  -5*sum(fit$psi*c(1,x))
}
plot(1:9,sapply(1:9,pred))

delta_effect <- function(psi_hat, V_hat, sH, delta, z = qnorm(0.975), scale = 1) {
  est <- scale * delta * as.numeric(crossprod(sH, psi_hat))
  var <- (scale^2) * (delta^2) * as.numeric(t(sH) %*% V_hat %*% sH)
  se  <- sqrt(var)
  c(estimate = est, se = se, lo = est - z*se, hi = est + z*se)
}

H_ex <- data.frame(co_ruc_code13 = 9)
sH   <- s_fun(H_ex)
out  <- delta_effect(psi_hat = fit$psi, V_hat = fit$V, sH = sH, delta = -5)
out
