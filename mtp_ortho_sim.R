expit <- function(x) 1/(1+exp(-x))
clip  <- function(x, lo=1e-8, hi=1-1e-8) pmax(lo, pmin(hi, x))
solve_ridge <- function(M, b, ridge = 1e-8) {
  p <- ncol(M); solve(M + diag(ridge, p), b, tol = 1e-12)
}

## Normal conditional density and ratio for a shift by delta:
## A|L ~ N(m(L), s^2). Then:
##   r = p(a - delta | L) / p(a | L) = dnorm(a - delta; m, s) / dnorm(a; m, s)
dens_ratio_normal_shift <- function(a, m, s, delta) {
  dnorm(a - delta, mean = m, sd = s) / dnorm(a, mean = m, sd = s)
}

## ---------- DGP: point-exposure, continuous A, shift delta ----------
simulate_point_shift <- function(n, delta = 0.5,
                                 theta = c(0.2, 0.8, -0.4),  # m(L) = th0 + th1*L + th2*L^2
                                 sigA  = 1.0,
                                 xi    = c(0.3, 0.5, 0.2),   # baseline in Y
                                 beta  = c(0.8, -0.6),       # TRUE psi = beta
                                 seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  L <- rnorm(n, 0, 1)
  mL <- theta[1] + theta[2]*L + theta[3]*L^2
  A  <- rnorm(n, mean = mL, sd = sigA)
  Y  <- xi[1] + xi[2]*L + xi[3]*L^2 + (beta[1] + beta[2]*L)*A + rnorm(n, 0, 1)
  data.frame(Y=Y, A=A, L=L, mL_true=mL)
}

## ---------- Estimator: cross-fitted orthogonal, continuous shift ----------
estimate_point_shift_fixed <- function(dat, delta = 0.5, K = 2,
                                       dens_floor = 1e-6, ratio_cap = 50) {
  n <- nrow(dat)
  folds <- sample(rep_len(1:K, n))
  s_of_L <- function(L) c(1, L)  # basis s(L)
  p <- 2
  
  M   <- matrix(0, p, p)
  rhs <- rep(0, p)
  
  for (k in 1:K) {
    idx_tr <- which(folds != k); idx_te <- which(folds == k)
    dtr <- dat[idx_tr, , drop=FALSE]
    dte <- dat[idx_te, , drop=FALSE]
    
    ## mu(a,l) ≈ E[Y|A,L]
    fit_mu <- lm(Y ~ A + L + I(L^2) + A:L, data = dtr)
    mu_hat <- function(a, l) predict(fit_mu, newdata = data.frame(Y=NA, A=a, L=l))
    
    ## p(A|L): Normal m(L), s; then use ratio p(a-δ|l)/p(a|l)
    fit_mA <- lm(A ~ L + I(L^2), data = dtr)
    m_te   <- predict(fit_mA, newdata = dte)
    s_te   <- rep(sd(residuals(fit_mA)), length(idx_te))
    
    A <- dte$A; L <- dte$L; Y <- dte$Y
    sL_list <- lapply(L, s_of_L)
    q_list  <- sL_list
    
    ## pullback: tilde q = q(L, A-δ) * ratio; here q doesn't depend on A, so same s(L)
    ratio <- dnorm(A - delta, mean = m_te, sd = s_te) /
      pmax(dnorm(A,         mean = m_te, sd = s_te), dens_floor)
    ## clip extreme ratios for stability
    ratio <- pmin(pmax(ratio, 1/ratio_cap), ratio_cap)
    tq_list <- Map(function(sv, r) sv * r, sL_list, as.list(ratio))
    
    mu_AL <- mu_hat(A, L)
    p_AL  <- pmax(dnorm(A, mean = m_te, sd = s_te), dens_floor)
    r_i   <- (Y - mu_AL) / p_AL
    
    mu_g  <- mu_hat(A + delta, L)
    
    for (i in seq_along(A)) {
      s_i  <- sL_list[[i]]
      q_i  <- q_list[[i]]
      tq_i <- tq_list[[i]]
      
      ## LEFT: δ * s s^T
      M   <- M + delta * (s_i %*% t(s_i))
      
      ## RIGHT:  -[ s(μ - μ_g) + (q - \tilde q) * r ]
      rhs <- rhs - ( s_i * (mu_AL[i] - mu_g[i]) + (q_i - tq_i) * r_i[i] )
    }
  }
  
  as.numeric(solve(M + diag(1e-8, p), rhs))  # tiny ridge
}

estimate_point_shift_oracle <- function(dat, delta=0.5,
                                        xi=c(0.3,0.5,0.2), beta=c(0.8,-0.6),
                                        theta=c(0.2,0.8,-0.4), sigA=1) {
  s_of_L <- function(L) c(1, L)
  mu_true <- function(a, l) xi[1] + xi[2]*l + xi[3]*l^2 + (beta[1] + beta[2]*l)*a
  mL <- theta[1] + theta[2]*dat$L + theta[3]*dat$L^2
  
  M <- matrix(0, 2, 2); rhs <- c(0,0)
  for (i in seq_len(nrow(dat))) {
    L <- dat$L[i]; A <- dat$A[i]; Y <- dat$Y[i]
    s  <- s_of_L(L)
    mu <- mu_true(A, L)
    mu_g <- mu_true(A + delta, L)
    pA  <- dnorm(A, mean = mL[i], sd = sigA)
    r   <- (Y - mu) / pA
    ratio <- dnorm(A - delta, mean = mL[i], sd = sigA) / pA
    tq <- s * ratio; q <- s
    M   <- M + delta * (s %*% t(s))
    rhs <- rhs - ( s * (mu - mu_g) + (q - tq) * r )
  }
  as.numeric(solve(M + diag(1e-12,2), rhs))
}

estimate_point_shift_oracle_stable <- function(dat, delta=0.5,
                                               xi=c(0.3,0.5,0.2),
                                               beta=c(0.8,-0.6),
                                               theta=c(0.2,0.8,-0.4),
                                               sigA=1) {
  s_of_L <- function(L) c(1, L)
  mu_true <- function(a, l) xi[1] + xi[2]*l + xi[3]*l^2 + (beta[1] + beta[2]*l)*a
  mL <- theta[1] + theta[2]*dat$L + theta[3]*dat$L^2
  
  M <- matrix(0, 2, 2); rhs <- c(0,0)
  for (i in seq_len(nrow(dat))) {
    L <- dat$L[i]; A <- dat$A[i]; Y <- dat$Y[i]
    s  <- s_of_L(L)
    mu <- mu_true(A, L)
    mu_g <- mu_true(A + delta, L)
    e  <- Y - mu
    ratio <- dnorm(A - delta, mean = mL[i], sd = sigA) /
      dnorm(A,         mean = mL[i], sd = sigA)
    tq <- s * ratio; q <- s
    M   <- M + delta * (s %*% t(s))
    rhs <- rhs - ( s * (mu - mu_g) + (q - tq) * e )
  }
  as.numeric(solve(M + diag(1e-12,2), rhs))
}

estimate_point_shift_stable <- function(dat, delta = 0.5, K = 2,
                                        ratio_cap = 50, ridge = 1e-8) {
  n <- nrow(dat)
  folds <- sample(rep_len(1:K, n))
  s_of_L <- function(L) c(1, L)  # basis in L only
  p <- 2
  M <- matrix(0, p, p)
  rhs <- rep(0, p)
  
  for (k in 1:K) {
    idx_tr <- which(folds != k); idx_te <- which(folds == k)
    dtr <- dat[idx_tr, , drop=FALSE]
    dte <- dat[idx_te, , drop=FALSE]
    
    ## Outcome regression mu(a,l)
    fit_mu <- lm(Y ~ A + L + I(L^2) + A:L, data = dtr)
    mu_hat <- function(a, l) predict(fit_mu, newdata = data.frame(Y=NA, A=a, L=l))
    
    ## Conditional density ratio for the shift: p(a-δ|l)/p(a|l)
    fit_mA <- lm(A ~ L + I(L^2), data = dtr)
    m_te   <- predict(fit_mA, newdata = dte)
    s_te   <- rep(sd(residuals(fit_mA)), length(idx_te))
    
    A <- dte$A; L <- dte$L; Y <- dte$Y
    sL <- lapply(L, s_of_L)   # q(L,A) = s(L)
    
    ## Pullback:  tilde q = q(L, A-δ) * p(A-δ|L)/p(A|L)  (here q doesn't depend on A)
    ratio <- dnorm(A - delta, mean = m_te, sd = s_te) /
      dnorm(A,         mean = m_te, sd = s_te)
    ratio <- pmin(pmax(ratio, 1/ratio_cap), ratio_cap)  # tame extremes
    tq    <- Map(function(sv, r) sv * r, sL, as.list(ratio))
    
    ## Residual WITHOUT dividing by a density
    e <- Y - mu_hat(A, L)
    
    ## Counterfactual regression at shifted A
    mu_g <- mu_hat(A + delta, L)
    
    for (i in seq_along(A)) {
      s_i  <- sL[[i]]
      tq_i <- tq[[i]]
      q_i  <- s_i
      
      ## δ * s s^T  ψ  =  - [ s (μ - μ_g) + (q - \tilde q) * e ]
      M   <- M + delta * (s_i %*% t(s_i))
      rhs <- rhs - ( s_i * ( (mu_hat(A[i], L[i]) - mu_g[i]) ) + (q_i - tq_i) * e[i] )
    }
  }
  as.numeric(solve(M + diag(ridge, p), rhs))
}



mc_point_shift <- function(ns = c(400, 1000, 3000), R = 200, delta = 0.5, seed = 123) {
  set.seed(seed)
  TRUE_psi <- c(0.8, -0.6)
  out <- list()
  for (n in ns) {
    ests <- matrix(NA_real_, R, 2)
    for (r in 1:R) {
      dat <- simulate_point_shift(n, delta = delta)
      ests[r, ] <- estimate_point_shift_stable(dat, delta = delta)#, K = 2)
    }
    colnames(ests) <- c("psi0","psi1")
    bias <- colMeans(ests) - TRUE_psi
    rmse <- sqrt(colMeans((t(t(ests) - TRUE_psi))^2))
    out[[as.character(n)]] <- list(ests=ests, bias=bias, rmse=rmse)
    cat(sprintf("[Point/shift] n=%d  bias=(%.3f, %.3f)  rmse=(%.3f, %.3f)\n",
                n, bias[1], bias[2], rmse[1], rmse[2]))
  }
  invisible(out)
}

estimate_point_shift_with_se <- function(
    dat, delta = 0.5, K = 2, ratio_cap = 50, ridge = 1e-8
) {
  n <- nrow(dat); p <- 2L
  folds <- sample(rep_len(1:K, n))
  s_of_L <- function(L) c(1, L)
  
  ## store per-observation ingredients for SEs
  S         <- matrix(NA_real_, n, p)     # s(L_i)
  mu        <- numeric(n)                 # mu(A_i, L_i)
  mu_g      <- numeric(n)                 # mu(A_i+delta, L_i)
  e         <- numeric(n)                 # residual Y - mu
  qminustq  <- matrix(NA_real_, n, p)     # q_i - tilde q_i
  
  ## linear system for psi
  M   <- matrix(0, p, p)
  rhs <- numeric(p)
  
  for (k in 1:K) {
    idx_tr <- which(folds != k); idx_te <- which(folds == k)
    dtr <- dat[idx_tr, , drop = FALSE]
    dte <- dat[idx_te, , drop = FALSE]
    
    ## mu(a,l) ≈ E[Y|A,L]
    fit_mu <- lm(Y ~ A + L + I(L^2) + A:L, data = dtr)
    mu_hat <- function(a, l) predict(fit_mu, newdata = data.frame(Y=NA, A=a, L=l))
    
    ## density ratio for shift via normal A|L
    fit_mA <- lm(A ~ L + I(L^2), data = dtr)
    m_te   <- predict(fit_mA, newdata = dte)
    s_te   <- rep(sd(residuals(fit_mA)), length(idx_te))
    
    A <- dte$A; L <- dte$L; Y <- dte$Y
    
    ratio <- dens_ratio_normal_shift(A, m_te, s_te, delta)
    ratio <- pmin(pmax(ratio, 1/ratio_cap), ratio_cap)  # tame extremes
    
    mu_AL <- mu_hat(A, L)
    mu_gL <- mu_hat(A + delta, L)
    e_te  <- Y - mu_AL
    
    for (j in seq_along(A)) {
      i   <- idx_te[j]
      s_i <- s_of_L(L[j])
      q_i <- s_i
      tq_i <- s_i * ratio[j]
      qmt_i <- q_i - tq_i
      
      ## store for SE
      S[i, ]        <- s_i
      mu[i]         <- mu_AL[j]
      mu_g[i]       <- mu_gL[j]
      e[i]          <- e_te[j]
      qminustq[i, ] <- qmt_i
      
      ## accumulate linear system:  δ s s^T ψ = -[ s(μ-μ_g) + (q-tilde q)e ]
      M   <- M + delta * (s_i %*% t(s_i))
      rhs <- rhs - ( s_i * (mu_AL[j] - mu_gL[j]) + qmt_i * e_te[j] )
    }
  }
  
  psi_hat <- as.numeric(solve_ridge(M, rhs, ridge))
  
  ## Sandwich SEs
  ## G_hat = (1/n) sum δ s s^T  ;  φ_i(ψ̂) = (q-tilde q)e + s{ μ-μ_g + δ s^T ψ̂ }
  G_hat <- (delta / n) * crossprod(S, S)
  phi   <- matrix(NA_real_, n, p)
  for (i in 1:n) {
    phi[i, ] <- qminustq[i, ] * e[i] + S[i, ] * ( (mu[i] - mu_g[i]) + delta * as.numeric(S[i, ] %*% psi_hat) )
  }
  ## center φ for numerical stability
  phi_bar <- colMeans(phi)
  Sigma_hat <- crossprod(phi - matrix(phi_bar, n, p, byrow = TRUE)) / n
  
  Vinv <- solve(G_hat)
  V_hat <- Vinv %*% Sigma_hat %*% t(Vinv) / n  # divide by n for estimator variance
  se <- sqrt(diag(V_hat))
  
  list(psi = psi_hat, se = se,
       ci = cbind(psi_hat - 1.96*se, psi_hat + 1.96*se),
       V = V_hat, G = G_hat, Sigma = Sigma_hat)
}

mc_point_shift_with_se <- function(ns = c(400, 1000, 3000), R = 500,
                                   delta = 0.5, seed = 123) {
  set.seed(seed)
  TRUE_psi <- c(0.8, -0.6)
  results <- list()
  for (n in ns) {
    est  <- matrix(NA_real_, R, 2)
    se   <- matrix(NA_real_, R, 2)
    covg <- matrix(NA_integer_, R, 2)
    for (r in 1:R) {
      dat <- simulate_point_shift(n, delta = delta)
      fit <- estimate_point_shift_with_se(dat, delta = delta, K = 2)
      est[r, ] <- fit$psi
      se[r, ]  <- fit$se
      CI       <- fit$ci
      covg[r, ] <- as.integer((TRUE_psi >= CI[,1]) & (TRUE_psi <= CI[,2]))
    }
    colnames(est) <- colnames(se) <- colnames(covg) <- c("psi0","psi1")
    
    bias <- colMeans(est) - TRUE_psi
    rmse <- sqrt(colMeans((t(t(est) - TRUE_psi))^2))
    emp_sd <- apply(est, 2, sd)
    mean_se <- colMeans(se)
    cover <- colMeans(covg)  # 95% coverage
    
    cat(sprintf(
      "[Point/shift] n=%d  bias=(%.3f, %.3f)  rmse=(%.3f, %.3f)  empSD=(%.3f, %.3f)  meanSE=(%.3f, %.3f)  cov95=(%.2f, %.2f)\n",
      n, bias[1], bias[2], rmse[1], rmse[2], emp_sd[1], emp_sd[2], mean_se[1], mean_se[2], cover[1], cover[2]
    ))
    results[[as.character(n)]] <- list(
      est=est, se=se, bias=bias, rmse=rmse,
      emp_sd=emp_sd, mean_se=mean_se, cover=cover
    )
  }
  invisible(results)
}

mc_point_shift_with_se()
