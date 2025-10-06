## ---------- utilities ----------
solve_ridge <- function(M, b, ridge = 1e-8) {
  p <- ncol(M); solve(M + diag(ridge, p), b, tol = 1e-12)
}
dens_ratio_normal_shift <- function(a, m, s, delta) {
  dnorm(a - delta, mean = m, sd = s) / dnorm(a, mean = m, sd = s)
}

## ---------- DGP parameters you control ----------
## --- parameters you can tweak coherently ---
## --- configuration with time-1 heterogeneity (psi1_L1 != 0) ---
dgp_cfg_het <- list(
  delta0 = 0.4, delta1 = 0.5,
  # time-1 blip params (true):
  psi1   = c(0.5, 0.3),                 # (psi1_0, psi1_L1)  <-- L1 heterogeneity on A1 effect
  # forward laws:
  rho    = c(0.1, 0.6, 0.8),            # L1 <- rho0 + rho1*L0 + rho2*A0 + noise
  sdL1   = 0.5,
  kappa  = c(0.2, 0.35),                # A1 <- k0 + k2*L0 + noise   (vector (k0, k2))
  sdA1   = 1.0,
  # outcome baselines:
  b1     = function(L0,L1) 0.25 + 0.5*L1 + 0.2*L0,  # b1_L1 = 0.5
  b0     = function(L0) 0.1 + 0.25*L0,              # free baseline under g
  # A0|L0:
  A0_mean = function(L0) 0.4*L0, sdA0 = 1.0,
  sdY   = 1.0
)

## implied psi0 from calibration
true_psi_from_cfg_het <- function(cfg = dgp_cfg_het){
  rho2 <- cfg$rho[3]; k0 <- cfg$kappa[1]; k2 <- cfg$kappa[2]
  psi1_1 <- cfg$psi1[2]; b1_L1 <- 0.5   # must match coefficient of L1 in cfg$b1
  psi0_0 <- rho2*b1_L1 + rho2*psi1_1*(k0 + cfg$delta1)
  psi0_1 <- rho2*psi1_1*k2
  list(psi0 = c(psi0_0, psi0_1), psi1 = cfg$psi1,
       delta0 = cfg$delta0, delta1 = cfg$delta1)
}

## simulator (coherent with the calibration)
simulate_long_mtpSNMM_het <- function(n, cfg = dgp_cfg_het){
  L0 <- rnorm(n)
  A0 <- rnorm(n, mean = cfg$A0_mean(L0), sd = cfg$sdA0)
  L1 <- rnorm(n, mean = cfg$rho[1] + cfg$rho[2]*L0 + cfg$rho[3]*A0, sd = cfg$sdL1)
  A1 <- rnorm(n, mean = cfg$kappa[1] + cfg$kappa[2]*L0, sd = cfg$sdA1)
  
  psi1 <- cfg$psi1
  mu1  <- function(a1, L0, L1) cfg$b1(L0,L1) + (psi1[1] + psi1[2]*L1) * a1
  
  Y <- mu1(A1, L0, L1) + rnorm(n, sd = cfg$sdY)
  data.frame(L0=L0, A0=A0, L1=L1, A1=A1, Y=Y)
}

## ---------- the analytic truth you’ll compare to ----------
true_psi <- function(cfg = dgp_defaults) {
  list(delta0 = cfg$delta0, delta1 = cfg$delta1, psi0 = cfg$psi0, psi1 = cfg$psi1)
}


## ---------- stable longitudinal estimator for continuous shifts ----------
## Uses: q_t(H_t, A_t) = s_t(H_t) (no A_t term), so pullback is s_t(H_t) * ratio
## helpers
solve_ridge <- function(M, b, ridge=1e-8){ p <- ncol(M); solve(M + diag(ridge,p), b, tol=1e-12) }
dens_ratio_normal_shift <- function(a, m, s, d) dnorm(a-d, m, s) / dnorm(a, m, s)

estimate_long_shift_with_se_basic <- function(dat, delta0, delta1, K=2, ratio_cap=50, ridge=1e-8){
  n <- nrow(dat); s0 <- function(L0) c(1, L0); s1 <- function(L1) c(1, L1)
  folds <- sample(rep_len(1:K, n))
  
  S0 <- matrix(NA_real_, n, 2); S1 <- matrix(NA_real_, n, 2)
  mu0 <- mu0_g <- e0 <- numeric(n); mu1 <- mu1_g <- e1 <- numeric(n)
  qmt0 <- matrix(NA_real_, n, 2); qmt1 <- matrix(NA_real_, n, 2)
  M <- matrix(0,4,4); rhs <- numeric(4)
  
  for (k in 1:K){
    idx_tr <- which(folds != k); idx_te <- which(folds == k)
    dtr <- dat[idx_tr, , drop=FALSE]; dte <- dat[idx_te, , drop=FALSE]
    
    ## t=1: E[Y|H1,A1]
    fit_mu1 <- lm(Y ~ A1 + L1 + L0 + A0 + A1:L1, data=dtr)     # includes the true A1×L1 effect-modifier
    mu1_hat <- function(a1, L1, L0, A0){
      predict(fit_mu1, newdata=data.frame(Y=NA, A1=a1, L1=L1, L0=L0, A0=A0))
    }
    ## A1|H1 density model (mean in L0 only per DGP; we allow L1 too harmlessly)
    fit_m1 <- lm(A1 ~ L0 + L1, data=dtr)
    m1_te  <- predict(fit_m1, newdata=dte)
    s1_te  <- rep(sd(residuals(fit_m1)), length(idx_te))
    
    L1 <- dte$L1; L0 <- dte$L0; A0 <- dte$A0; A1 <- dte$A1; Y <- dte$Y
    mu1_AL <- mu1_hat(A1, L1, L0, A0); e1_te <- Y - mu1_AL
    mu1_gL <- mu1_hat(A1 + delta1, L1, L0, A0)
    
    s1_list <- lapply(L1, s1)
    ratio1  <- dens_ratio_normal_shift(A1, m1_te, s1_te, delta1)
    ratio1  <- pmin(pmax(ratio1, 1/ratio_cap), ratio_cap)
    tq1     <- Map(function(sv, r) sv * r, s1_list, as.list(ratio1))
    
    for (j in seq_along(A1)){
      i <- idx_te[j]; s1i <- s1_list[[j]]; q1i <- s1i; tq1i <- tq1[[j]]
      S1[i,] <- s1i; mu1[i] <- mu1_AL[j]; mu1_g[i] <- mu1_gL[j]; e1[i] <- e1_te[j]; qmt1[i,] <- (q1i - tq1i)
      M[3:4,3:4] <- M[3:4,3:4] + delta1 * (s1i %*% t(s1i))
      rhs[3:4]   <- rhs[3:4] - ( s1i * (mu1_AL[j] - mu1_gL[j]) + (q1i - tq1i) * e1_te[j] )
    }
    
    ## V1 for t=0
    V1_tr <- mu1_hat(dtr$A1 + delta1, dtr$L1, dtr$L0, dtr$A0)
    V1_te <- mu1_gL
    
    ## t=0: E[V1|H0,A0]
    fit_mu0 <- lm(V1_tr ~ A0 + L0 + A0:L0, data=cbind(dtr, V1_tr=V1_tr))
    mu0_hat <- function(a0, L0){ predict(fit_mu0, newdata=data.frame(V1_tr=NA, A0=a0, L0=L0)) }
    fit_m0  <- lm(A0 ~ L0, data=dtr)
    m0_te   <- predict(fit_m0, newdata=dte)
    s0_te   <- rep(sd(residuals(fit_m0)), length(idx_te))
    
    A0_te <- dte$A0; L0_te <- dte$L0
    mu0_AL <- mu0_hat(A0_te, L0_te); e0_te <- V1_te - mu0_AL
    mu0_gL <- mu0_hat(A0_te + delta0, L0_te)
    
    s0_list <- lapply(L0_te, s0)
    ratio0  <- dens_ratio_normal_shift(A0_te, m0_te, s0_te, delta0)
    ratio0  <- pmin(pmax(ratio0, 1/ratio_cap), ratio_cap)
    tq0     <- Map(function(sv, r) sv * r, s0_list, as.list(ratio0))
    
    for (j in seq_along(A0_te)){
      i <- idx_te[j]; s0i <- s0_list[[j]]; q0i <- s0i; tq0i <- tq0[[j]]
      S0[i,] <- s0i; mu0[i] <- mu0_AL[j]; mu0_g[i] <- mu0_gL[j]; e0[i] <- e0_te[j]; qmt0[i,] <- (q0i - tq0i)
      M[1:2,1:2] <- M[1:2,1:2] + delta0 * (s0i %*% t(s0i))
      rhs[1:2]   <- rhs[1:2] - ( s0i * (mu0_AL[j] - mu0_gL[j]) + (q0i - tq0i) * e0_te[j] )
    }
  }
  
  psi_hat <- as.numeric(solve_ridge(M, rhs, ridge))
  names(psi_hat) <- c("psi0_0","psi0_L0","psi1_0","psi1_L1")
  
  ## sandwich
  G_hat <- matrix(0,4,4)
  G_hat[1:2,1:2] <- (delta0/ n) * crossprod(S0,S0)
  G_hat[3:4,3:4] <- (delta1/ n) * crossprod(S1,S1)
  phi <- matrix(NA_real_, n, 4)
  psi0 <- psi_hat[1:2]; psi1 <- psi_hat[3:4]
  for (i in 1:n){
    phi0 <- qmt0[i,]*e0[i] + S0[i,]*((mu0[i]-mu0_g[i]) + delta0*as.numeric(S0[i,] %*% psi0))
    phi1 <- qmt1[i,]*e1[i] + S1[i,]*((mu1[i]-mu1_g[i]) + delta1*as.numeric(S1[i,] %*% psi1))
    phi[i,] <- c(phi0, phi1)
  }
  phi_bar <- colMeans(phi); Sigma_hat <- crossprod(phi - matrix(phi_bar, n, 4, byrow=TRUE)) / n
  Vinv <- solve(G_hat); V_hat <- Vinv %*% Sigma_hat %*% t(Vinv)/n 
  list(psi=psi_hat, se=sqrt(diag(V_hat)), ci=cbind(psi_hat-1.96*sqrt(diag(V_hat)), psi_hat+1.96*sqrt(diag(V_hat))))
}



mc_long_mtpSNMM_het <- function(ns=c(800,2000,6000), R=300, cfg=dgp_cfg_het, seed=7){
  set.seed(seed)
  tr <- true_psi_from_cfg_het(cfg); truth <- c(tr$psi0, tr$psi1)
  for (n in ns){
    est <- matrix(NA_real_, R, 4); se <- matrix(NA_real_, R, 4); covg <- matrix(NA_integer_, R, 4); boot_se = matrix(NA_real_, R, 4); boot_covg <- matrix(NA_integer_, R, 4)
    for (r in 1:R){
      dat <- simulate_long_mtpSNMM_het(n)
      fit <- estimate_long_shift_with_se_basic(dat, delta0=cfg$delta0, delta1=cfg$delta1, K=2)
      est[r,] <- fit$psi; se[r,] <- fit$se
      est_boots <- matrix(NA_real_, 200, 4)
      for(b in 1:200){
        dat_boot = dat[sample(1:n,n,replace=T),]
        fit_boot = estimate_long_shift_with_se_basic(dat_boot, delta0=cfg$delta0, delta1=cfg$delta1, K=2)
        est_boots[b,]=fit_boot$psi
      }
      boot_se[r,]=apply(est_boots,2,sd,na.rm=T)
      CI <- fit$ci; covg[r,] <- as.integer((truth >= CI[,1]) & (truth <= CI[,2]))
      CI_boot <- cbind(fit$psi-1.96*boot_se[r,],fit$psi+1.96*boot_se[r,]); boot_covg[r,] <- as.integer((truth >= CI_boot[,1]) & (truth <= CI_boot[,2]))
    }
    bias <- colMeans(est) - truth
    rmse <- sqrt(colMeans((t(t(est)-truth))^2))
    empSD <- apply(est, 2, sd); meanSE <- colMeans(se); cover <- colMeans(covg); boot_cover = colMeans(boot_covg)
    cat(sprintf("[Long/MTP het] n=%d  bias=(%.3f, %.3f, %.3f, %.3f)  rmse=(%.3f, %.3f, %.3f, %.3f)  empSD=(%.3f, %.3f, %.3f, %.3f)  meanSE=(%.3f, %.3f, %.3f, %.3f)  cov95=(%.2f, %.2f, %.2f, %.2f) bootcov95=(%.2f, %.2f, %.2f, %.2f) \n",
                n, bias[1],bias[2],bias[3],bias[4], rmse[1],rmse[2],rmse[3],rmse[4],
                empSD[1],empSD[2],empSD[3],empSD[4], meanSE[1],meanSE[2],meanSE[3],meanSE[4],
                cover[1],cover[2],cover[3],cover[4],boot_cover[1],boot_cover[2],boot_cover[3],boot_cover[4]))
  }
  invisible(NULL)
}



Sys.time()
mc_long_mtpSNMM_het(ns=1000,R=500)
Sys.time()
