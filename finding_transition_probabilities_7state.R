# Authors: Created by Marcus Moen, MS and edited by James Dickersonm MD, MS.
# Oversight: Jeremy Goldhaber-Feibert, PhD and Fernando Alarid-Escobedo, PhD

# ---- Libraries and Options ----
options(scipen = 999)  # Avoid scientific notation
library(expm)          # Matrix exponentials
library(ggplot2)       # Plotting

#The target variables in line 161 are the values we are calibrating to. Change this if you want to change the calibration
#The sd corresponding to the targets is found at line 296

#Function for calculating the standard error of a survival curve given 
#number at risk, n_risk, number of deaths at each timestep, n_deaths, and the survival probabilities, surv_prob



# ---- Transition Matrix Generator ----

#' Generate Transition Matrices for 7-State Model
#' @param input_var Numeric vector of 9 parameters used in the model
#' @param calibrate Boolean, if TRUE uses fewer time steps for faster calibration
#' @return A list of two lists of transition matrices: one for chemo and one for T-DXd
transition_matrices_7state <- function(input_var, calibrate=T){
  # Extract rates and parameters from input
  r_pf2prog <- input_var[1]
  r_prog2death <- input_var[2]
  hr_pf2prog_chemo2tdxd <- input_var[3]
  gamma_chemo <- input_var[4]; alpha_chemo <- input_var[5]
  gamma_tdxd  <- input_var[6]; alpha_tdxd  <- input_var[7]
  gamma_ILD   <- input_var[8]; alpha_ILD   <- input_var[9]
  
  # Background mortality rates
  r_OC_death <- rep(c(0.000455117, 0.000667712, 0.000939403, 0.001487732, 0.002481213,
                      0.004368355, 0.008008381, 0.014505686, 0.024617632), each = 60)
  
  
  n <- if (calibrate) 40 else 121
  
  # Initialize vectors for AE and ILD rates
  Ft_chemo <- Ft_tdxd <- Ft_tdxd_ILD <- Ht_chemo <- Ht_tdxd <- Ht_tdxd_ILD <- c()
  r_AE <- r_AE_tdxd <- r_ILD <- r_ILD_tdxd <- c()
  
  for(t in 1:(n+1)){
    Ft_chemo <- append(Ft_chemo, alpha_chemo*(1-exp(-gamma_chemo*t)))
    Ft_tdxd <- append(Ft_tdxd, alpha_tdxd*(1-exp(-gamma_tdxd*t)))
    Ft_tdxd_ILD <- append(Ft_tdxd_ILD, alpha_ILD*(1-exp(-gamma_ILD*t)))
    
    Ht_chemo <- append(Ht_chemo, -log(1-Ft_chemo[t]))
    Ht_tdxd <- append(Ht_tdxd, -log(1-Ft_tdxd[t]))
    Ht_tdxd_ILD <- append(Ht_tdxd_ILD, -log(1-Ft_tdxd_ILD[t]))
    
    if(t > 1){
      r_AE <- append(r_AE, Ht_chemo[t]-Ht_chemo[t-1])
      r_AE_tdxd <- append(r_AE_tdxd, Ht_tdxd[t]-Ht_tdxd[t-1])
      r_ILD <- append(r_ILD, 0)
      r_ILD_tdxd <- append(r_ILD_tdxd, Ht_tdxd_ILD[t]-Ht_tdxd_ILD[t-1])
    }
  }
  
  
  # Hazard ratios for various transitions (may be modified later)
  hr_defaults <- rep(1, 7)
  names(hr_defaults) <- c("hr_pfAE2progAE", "hr_progAE2death", "hr_pfae2prog_chemo2tdxd",
                          "hr_pf2prog2death_chemo2tdxd", "hr_pfAE2progAE2death_chemo2tdxd",
                          "hr_pf2progILD_chemo2tdxd", "hr_pf2progILD2death_chemo2tdxd")
  
  
  with(as.list(hr_defaults), {
    # Adjusted rates for each treatment arm using hazard ratios
    r_pfAE2progAE <- r_pf2prog * hr_pfAE2progAE
    r_pf2progILD <- r_pf2prog * 0
    r_pf2progILD_death <- r_prog2death * 0
    r_progAE2death <- r_prog2death * hr_progAE2death
    
    r_pf2prog_tdxd <- r_pf2prog * hr_pf2prog_chemo2tdxd
    r_pfAE2progAE_tdxd <- r_pf2prog * hr_pfae2prog_chemo2tdxd
    r_prog2death_tdxd <- r_prog2death * hr_pf2prog2death_chemo2tdxd
    r_progAE2death_tdxd <- r_prog2death * hr_pfAE2progAE2death_chemo2tdxd
    r_pf2progILD_tdxd <- r_pf2prog * hr_pf2progILD_chemo2tdxd
    r_pf2progILD_death_tdxd <- r_prog2death * hr_pf2progILD2death_chemo2tdxd
    
    P_chemo <- list()
    P_tdxd <- list()
    
    for (t in 1:n) {
      # Generator matrix for chemo
      G_chemo <- matrix(c(
        -(r_pf2prog + r_OC_death[t] + r_AE[t] + r_ILD[t]), r_AE[t], r_ILD[t], r_pf2prog, 0, 0, r_OC_death[t],
        0, -(r_pfAE2progAE + r_OC_death[t]), 0, 0, r_pfAE2progAE, 0, r_OC_death[t],
        0, 0, -(r_pf2progILD + r_OC_death[t]), 0, 0, r_pf2progILD, r_OC_death[t],
        0, 0, 0, -(r_OC_death[t] + r_prog2death), 0, 0, r_OC_death[t] + r_prog2death,
        0, 0, 0, 0, -(r_OC_death[t] + r_progAE2death), 0, r_OC_death[t] + r_progAE2death,
        0, 0, 0, 0, 0, -(r_OC_death[t] + r_pf2progILD_death), r_OC_death[t] + r_pf2progILD_death,
        0, 0, 0, 0, 0, 0, 0),
        dimnames = list(c("PF", "PF_AE", "PF_ILD", "Prog", "Prog_AE", "Prog_ILD", "Death"),
                        c("PF", "PF_AE", "PF_ILD", "Prog", "Prog_AE", "Prog_ILD", "Death")),
        ncol = 7, byrow = TRUE)
      
      # Generator matrix for TDxd
      G_tdxd <- matrix(G_chemo, 
                       dimnames = list(c("PF", "PF_AE", "PF_ILD", "Prog", "Prog_AE", "Prog_ILD", "Death"),
                                       c("PF", "PF_AE", "PF_ILD", "Prog", "Prog_AE", "Prog_ILD", "Death")),
                       ncol = 7)
      G_tdxd[1, ] <- c(-(r_pf2prog_tdxd + r_OC_death[t] + r_AE_tdxd[t] + r_ILD_tdxd[t]), r_AE_tdxd[t], r_ILD_tdxd[t], r_pf2prog_tdxd, 0, 0, r_OC_death[t])
      G_tdxd[2, ] <- c(0, -(r_pfAE2progAE_tdxd + r_OC_death[t]), 0, 0, r_pfAE2progAE_tdxd, 0, r_OC_death[t])
      G_tdxd[3, ] <- c(0, 0, -(r_pf2progILD_tdxd + r_OC_death[t]), 0, 0, r_pf2progILD_tdxd, r_OC_death[t])
      G_tdxd[4, 4:7] <- c(-(r_OC_death[t] + r_prog2death_tdxd), 0, 0, r_OC_death[t] + r_prog2death_tdxd)
      G_tdxd[5, 5:7] <- c(-(r_OC_death[t] + r_progAE2death_tdxd), 0, r_OC_death[t] + r_progAE2death_tdxd)
      G_tdxd[6, 6:7] <- c(-(r_OC_death[t] + r_pf2progILD_death_tdxd), r_OC_death[t] + r_pf2progILD_death_tdxd)
      
      P_chemo[[t]] <- expm(G_chemo)
      P_tdxd[[t]] <- expm(G_tdxd)
    }
    
    return(list(P_chemo = P_chemo, P_tdxd = P_tdxd))
  })
}




# ---- Evaluation Function ----

#' Evaluate Parameter Fit to Target Values
#'
#' @param input_var Vector of model parameters (same 9 values as above)
#' @param target_var Named vector of calibration targets (median OS, PFS, AE rates)
#' @return Negative log-likelihood value for optimization
#' 
eval_variables_7state <- function(input_var, target_var = c(median_os_chemo = 16.8, 
                                                            median_pf_chemo = 5.1, 
                                                            median_os_tdxd = 23.4, 
                                                            median_pf_tdxd = 9.9, 
                                                            target_ae_tdxd_chemo = 0.081, 
                                                            target_ae_tdxd = 0.0757, 
                                                            target_ild_tdxd = 0.0863)){
  
  
  matrix_list <- transition_matrices_7state(input_var)
  A_chemo_list <- matrix_list$P_chemo
  A_tdxd_list <- matrix_list$P_tdxd
  
  
  
  # Initialize survival vectors
  s0_chemo <- s0_tdxd <- c(1, 0, 0, 0, 0, 0, 0)
  km_os_chemo <- km_pf_chemo <- km_os_tdxd <- km_pf_tdxd <- c()
  ae_chemo <- ae_tdxd <- ild_tdxd <- c(0)
  
  # Initialize at time 0
  km_os_chemo[1] <- sum(s0_chemo[c(1, 4)])
  km_pf_chemo[1] <- s0_chemo[1]
  km_os_tdxd[1] <- sum(s0_tdxd[c(1, 4)])
  km_pf_tdxd[1] <- s0_tdxd[1]
  
  n_cycles <- 30
  
  for (t in 1:n_cycles) {
    A_chemo <- A_chemo_list[[t]]
    A_tdxd <- A_tdxd_list[[t]]
    s1_chemo <- s0_chemo %*% A_chemo
    s1_tdxd <- s0_tdxd %*% A_tdxd
    
    km_os_chemo[t + 1] <- km_os_chemo[t] * sum(s1_chemo[c(1, 4)]) /
      (sum(s1_chemo[c(1, 4)]) + s0_chemo[1]*A_chemo[1,7] + s0_chemo[4]*A_chemo[4,7])
    km_pf_chemo[t + 1] <- km_pf_chemo[t] * s1_chemo[1] /
      (s1_chemo[1] + s0_chemo[1]*(A_chemo[1,4] + A_chemo[1,7]))
    
    km_os_tdxd[t + 1] <- km_os_tdxd[t] * sum(s1_tdxd[c(1, 4)]) /
      (sum(s1_tdxd[c(1, 4)]) + s0_tdxd[1]*A_tdxd[1,7] + s0_tdxd[4]*A_tdxd[4,7])
    km_pf_tdxd[t + 1] <- km_pf_tdxd[t] * s1_tdxd[1] /
      (s1_tdxd[1] + s0_tdxd[1]*(A_tdxd[1,4] + A_tdxd[1,7]))
    
    ae_chemo[t + 1] <- ae_chemo[t] + s0_chemo[1] * A_chemo[1,2]
    ae_tdxd[t + 1] <- ae_tdxd[t] + s0_tdxd[1] * A_tdxd[1,2]
    ild_tdxd[t + 1] <- ild_tdxd[t] + s0_tdxd[1] * A_tdxd[1,3]
    
    s0_chemo <- s1_chemo
    s0_tdxd <- s1_tdxd
  }
  
  # Extract model outputs
  model_targets <- c(
    which(km_os_chemo < 0.5)[1] - 1,
    which(km_pf_chemo < 0.5)[1] - 1,
    which(km_os_tdxd < 0.5)[1] - 1,
    which(km_pf_tdxd < 0.5)[1] - 1,
    ae_chemo[12], ae_tdxd[12], ild_tdxd[12]
  )
  
  # Target standard deviations
  v_target_sd <- c(1.4, 0.66, 1.22, 0.59, 0.0203, 0.0189, 0.0216)
  
  # Negative log-likelihood for optimizer
  -sum(dnorm(x = target_var, mean = model_targets, sd = v_target_sd, log = TRUE))
}




#Find the optimal values (Here we could do a wider search for starting values since we might have missed the global optimum.)
fit_out <- optim(c(0.15, 0.08, 0.5, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2), 
                 eval_variables_7state,
                 hessian = T)

m_cov <- solve(fit_out$hessian)
m_cor <- cov2cor(m_cov)
v_se <- sqrt(diag(m_cov))

#Extract the optimal values for our three variables
opt_var <- fit_out$par

opt_var


# ---- Optimization Wrapper ----

#' Fit Optimal Parameters to Match Target Outcomes
#' @param target Named vector of target summary statistics (optional)
#' @return Vector of optimal parameters
find_opt_var_7state <- function(target = c(median_os_chemo = 16.8, 
                                           median_pf_chemo = 5.1, 
                                           median_os_tdxd = 23.4, 
                                           median_pf_tdxd = 9.9, 
                                           target_ae_tdxd_chemo = 0.081, 
                                           target_ae_tdxd = 0.0757, 
                                           target_ild_tdxd = 0.0863)){
  fit_out <- optim(c(0.15, 0.08, 0.5, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2),
                   eval_variables_7state,
                   target_var = target)
  opt_var <- fit_out$par
  
  return(opt_var)
}



#Find the optimal transition matrices
tm <- transition_matrices_7state(opt_var, calibrate = F)







# ---- Validation Plot Function ----

#' Plot Modeled vs. Observed KM Curves for OS and PFS
#'
#' Compares survival curves simulated from transition matrices with observed Kaplan-Meier curves.
#' Useful for visually assessing model calibration.
#'
#' @param input_var A numeric vector of 9 model parameters
plot_function <- function(input_var) {
  matrix_list <- transition_matrices_7state(input_var)
  A_chemo_list <- matrix_list$P_chemo
  A_tdxd_list <- matrix_list$P_tdxd
  
  # Observed KM curves (trimmed)
  km_os_chemo <- c(1, 0.986, 0.981, 0.957, 0.94, 0.927, 0.884, 0.841, 0.793, 0.738, 0.712, 0.683, 0.669, 0.634)
  km_pf_chemo <- c(1, 0.983, 0.756, 0.62, 0.603, 0.501, 0.439, 0.398, 0.34, 0.277, 0.268, 0.226, 0.217, 0.185)
  km_os_tdxd  <- c(1, 0.992, 0.985, 0.967, 0.958, 0.939, 0.928, 0.899, 0.87, 0.857, 0.829, 0.815, 0.793, 0.774)
  km_pf_tdxd  <- c(1, 0.988, 0.896, 0.82, 0.809, 0.765, 0.681, 0.633, 0.592, 0.554, 0.491, 0.457, 0.422, 0.386)
  
  # Initialize vectors for modeled survival
  s0_chemo <- s0_tdxd <- c(1, 0, 0, 0, 0, 0, 0)
  km_os_chemo_model <- km_pf_chemo_model <- km_os_tdxd_model <- km_pf_tdxd_model <- c()
  ae_test_chemo <- ae_test_tdxd <- ae_test_ild <- c(0)
  
  km_os_chemo_model[1] <- sum(s0_chemo[c(1, 4)])
  km_pf_chemo_model[1] <- s0_chemo[1]
  km_os_tdxd_model[1] <- sum(s0_tdxd[c(1, 4)])
  km_pf_tdxd_model[1] <- s0_tdxd[1]
  
  n <- length(km_os_chemo)
  
  for (t in 1:n) {
    A_chemo <- A_chemo_list[[t]]
    A_tdxd  <- A_tdxd_list[[t]]
    s1_chemo <- s0_chemo %*% A_chemo
    s1_tdxd  <- s0_tdxd %*% A_tdxd
    
    km_os_chemo_model[t + 1] <- km_os_chemo_model[t] * sum(s1_chemo[c(1, 4)]) /
      (sum(s1_chemo[c(1, 4)]) + s0_chemo[1]*A_chemo[1, 7] + s0_chemo[4]*A_chemo[4, 7])
    km_pf_chemo_model[t + 1] <- km_pf_chemo_model[t] * s1_chemo[1] /
      (s1_chemo[1] + s0_chemo[1]*(A_chemo[1, 4] + A_chemo[1, 7]))
    
    km_os_tdxd_model[t + 1] <- km_os_tdxd_model[t] * sum(s1_tdxd[c(1, 4)]) /
      (sum(s1_tdxd[c(1, 4)]) + s0_tdxd[1]*A_tdxd[1, 7] + s0_tdxd[4]*A_tdxd[4, 7])
    km_pf_tdxd_model[t + 1] <- km_pf_tdxd_model[t] * s1_tdxd[1] /
      (s1_tdxd[1] + s0_tdxd[1]*(A_tdxd[1, 4] + A_tdxd[1, 7]))
    
    ae_test_chemo[t + 1] <- ae_test_chemo[t] + s0_chemo[1]*A_chemo[1, 2]
    ae_test_tdxd[t + 1]  <- ae_test_tdxd[t] + s0_tdxd[1]*A_tdxd[1, 2]
    ae_test_ild[t + 1]   <- ae_test_ild[t] + s0_tdxd[1]*A_tdxd[1, 3]
    
    s0_chemo <- s1_chemo
    s0_tdxd <- s1_tdxd
  }
  
  # Plot AE accumulation
  time <- 0:(length(ae_test_chemo) - 1)
  par(mfrow = c(1, 3))
  plot(time, ae_test_chemo, type = "l", col = "blue", main = "AE: Chemo", xlab = "Month", ylab = "Prob")
  plot(time, ae_test_tdxd, type = "l", col = "red", main = "AE: T-DXd", xlab = "Month", ylab = "Prob")
  plot(time, ae_test_ild, type = "l", col = "purple", main = "ILD: T-DXd", xlab = "Month", ylab = "Prob")
}



plot_function(opt_var)




