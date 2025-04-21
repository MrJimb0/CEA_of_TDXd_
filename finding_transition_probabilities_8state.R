# Authors: Created by Marcus Moen, MS and edited by James Dickerson, MD, MS
# Oversight: Jeremy Goldhaber-Feibert, PhD and Fernando Alarid-Escobedo, PhD

# ---- Libraries and Global Options ----
options(scipen = 999)  # Avoid scientific notation
library(expm)          # Matrix exponentials
library(ggplot2)       # Plotting

# ---- Notation Guide ----
# PF = Progression-Free first-line
# Prog_drug = Progressed, receiving second-line drug
# Prog_nodrug = Progressed, no second-line drug (palliative care)
# AE = Adverse Event (non-ILD)
# ILD = Interstitial Lung Disease

# ---- Configuration Notes ----
# - Change proportion receiving second-line therapy at: line with 'proportion_r_pf2prog_drug'
# - PFS for second-line therapy configured near line 439-447

# ---- Import and Initialization ----
source("finding_transition_probabilities_7state.R")  # Uses opt_var from prior script

# Global model parameters from the 7-state optimization
seven_state_var <- c(
  r_pf2prog              = opt_var[1],
  r_prog2death           = opt_var[2],
  hr_pf2prog_chemo2tdxd  = opt_var[3],
  gamma_chemo            = opt_var[4],
  alpha_chemo            = opt_var[5],
  gamma_tdxd             = opt_var[6],
  alpha_tdxd             = opt_var[7],
  gamma_ILD              = opt_var[8],
  alpha_ILD              = opt_var[9]
)

# Proportion of patients progressing to second-line therapy vs palliative care
# Important sensitivity parameter
proportion_r_pf2prog_drug <- 0.6

# ---- Helper: Chemo-Chemo Transition Variables ----

#' Define Variables for Chemo-Chemo Treatment Path
#'
#' Splits PF progression rates into those receiving chemo again vs. palliative care.
#' Also returns AE and ILD parameters.
#'
#' @param seven_state_var Vector of core model rates/parameters
#' @param proportion_r_pf2prog_drug Proportion progressing to drug (not palliation)
#' @return Named vector of model parameters used in chemo-chemo transition generation
chemo_chemo_var <- function(seven_state_var, proportion_r_pf2prog_drug) {
  
  # Progression rates split between second-line drug and no-drug
  r_pf2prog_drug   <- seven_state_var[1] * proportion_r_pf2prog_drug
  r_pf2prog_nodrug <- seven_state_var[1] * (1 - proportion_r_pf2prog_drug)
  
  # Assume progression under second-line chemo mirrors first-line
  r_prog_drug2prog_nodrug <- seven_state_var[1]
  
  # AE Weibull parameters from 7-state input
  gamma_AE <- seven_state_var[4]
  alpha_AE <- seven_state_var[5]
  
  # No ILD assumed for chemo (set to 0)
  gamma_ILD <- 0
  alpha_ILD <- 0
  
  return(c(
    r_pf2prog_drug, r_pf2prog_nodrug, r_prog_drug2prog_nodrug,
    gamma_AE, alpha_AE, gamma_ILD, alpha_ILD
  ))
}

# ---- Transition Matrix Generator: Chemo-Chemo ----

#' Generate Transition Matrices for Chemo-Chemo Sequence (8-state model)
#'
#' Constructs monthly transition matrices accounting for AE, ILD, progression,
#' background mortality, and second-line treatment outcomes.
#'
#' @param r_prog_nodrug2death Rate of death for patients not receiving second-line drug
#' @param proportion_r_pf2prog_drug Proportion of PF patients receiving second-line drug
#' @param seven_state_var Vector of global parameters from 7-state optimization
#' @param calibrate Boolean: Use shortened horizon for calibration (TRUE = 40 steps, FALSE = 121)
#' @param firstline_2line Character: treatment sequence descriptor (default "Chemo-Chemo")
#' @return List of monthly transition probability matrices (length n)
transition_matrices_chemo2chemo <- function(r_prog_nodrug2death,
                                            proportion_r_pf2prog_drug,
                                            seven_state_var,
                                            calibrate = TRUE,
                                            firstline_2line = "Chemo-Chemo") {
  
  # Extract model parameters based on treatment path
  if (firstline_2line == "Chemo-Chemo") {
    variables <- chemo_chemo_var(seven_state_var, proportion_r_pf2prog_drug)
  }
  
  r_pf2prog_drug   <- variables[1]
  r_pf2prog_nodrug <- variables[2]
  r_prog_drug2prog_nodrug <- variables[3]
  gamma_AE <- variables[4]; alpha_AE <- variables[5]
  gamma_ILD <- variables[6]; alpha_ILD <- variables[7]
  
  # Background mortality by age group (repeated per month)
  r_OC_death <- rep(c(0.000455117, 0.000667712, 0.000939403, 0.001487732,
                      0.002481213, 0.004368355, 0.008008381, 0.014505686, 0.024617632),
                    each = 60)
  
  n <- if (calibrate) 40 else 121
  
  # Time-varying transition rates for AE and ILD (Weibull-based)
  Ft_AE <- alpha_AE * (1 - exp(-gamma_AE * (1:(n + 1))))
  Ft_ILD <- alpha_ILD * (1 - exp(-gamma_ILD * (1:(n + 1))))
  Ht_AE <- -log(1 - Ft_AE)
  Ht_ILD <- -log(1 - Ft_ILD)
  
  r_pf2AE <- r_pf2ILD <- c()
  for (t in 2:(n+1)) {
    r_pf2AE  <- append(r_pf2AE, Ht_AE[t] - Ht_AE[t - 1])
    r_pf2ILD <- append(r_pf2ILD, Ht_ILD[t] - Ht_ILD[t - 1])
  }
  
  # Transition rates sourced from 7-state model (fixed)
  r_pfAE2progAE <- 0.15362337
  r_progAE2death <- 0.07171597
  r_pfILD2progILD <- r_pfAE2progAE
  r_pf2progILD_death <- r_progAE2death
  
  # Initialize list of transition matrices
  P <- vector("list", n)
  
  for (t in 1:n) {
    min_pf <- r_OC_death[t] + r_pf2prog_nodrug + r_pf2prog_drug + r_pf2ILD[t] + r_pf2AE[t]
    
    G <- matrix(c(
      -min_pf, r_pf2AE[t], r_pf2ILD[t], r_pf2prog_drug, r_pf2prog_nodrug, 0, 0, r_OC_death[t],
      0, -(r_OC_death[t] + r_pfAE2progAE), 0, 0, 0, r_pfAE2progAE, 0, r_OC_death[t],
      0, 0, -(r_OC_death[t] + r_pfILD2progILD), 0, 0, 0, r_pfILD2progILD, r_OC_death[t],
      0, 0, 0, -(r_OC_death[t] + r_prog_drug2prog_nodrug), r_prog_drug2prog_nodrug, 0, 0, r_OC_death[t],
      0, 0, 0, 0, -(r_OC_death[t] + r_prog_nodrug2death), 0, 0, r_OC_death[t] + r_prog_nodrug2death,
      0, 0, 0, 0, 0, -(r_OC_death[t] + r_progAE2death), 0, r_OC_death[t] + r_progAE2death,
      0, 0, 0, 0, 0, 0, -(r_OC_death[t] + r_pf2progILD_death), r_OC_death[t] + r_pf2progILD_death,
      0, 0, 0, 0, 0, 0, 0, 0
    ),
    ncol = 8, byrow = TRUE,
    dimnames = list(
      c("PF", "PF_AE", "PF_ILD", "Prog_drug", "Prog_nodrug", "Prog_AE", "Prog_ILD", "Death"),
      c("PF", "PF_AE", "PF_ILD", "Prog_drug", "Prog_nodrug", "Prog_AE", "Prog_ILD", "Death")
    ))
    
    P[[t]] <- expm(G)
  }
  
  return(P)
}




# ---- Evaluation Function ----

#' Evaluate Model Fit for Chemo-Chemo Scenario
#'
#' Simulates survival under the Chemo-Chemo treatment path and compares predicted
#' median overall survival (OS) to a target value using a likelihood function.
#'
#' @param r_prog_nodrug2death Death rate for patients receiving no second-line drug
#' @param proportion_r_pf2prog_drug Proportion of PF patients receiving second-line drug
#' @param seven_state_var Vector of core model parameters
#' @param drug_sequence String identifying treatment pathway (default = "Chemo-Chemo")
#' @return Negative log-likelihood of model fit to OS target

eval_variables <- function(r_prog_nodrug2death,
                           proportion_r_pf2prog_drug,
                           seven_state_var,
                           drug_sequence = "Chemo-Chemo") {
  
  if (drug_sequence == "Chemo-Chemo") {
    median_os <- 16.8
    A_list <- transition_matrices_chemo2chemo(
      r_prog_nodrug2death,
      proportion_r_pf2prog_drug,
      seven_state_var
    )
  }
  
  n_cycles <- 30
  s0 <- c(1, 0, 0, 0, 0, 0, 0, 0)  # Initial state vector
  km_os_model <- numeric(n_cycles + 1)
  km_os_model[1] <- sum(s0[c(1, 4, 5)])  # PF, Prog_drug, Prog_nodrug
  
  model_median_os <- 0
  
  for (t in 1:n_cycles) {
    A <- A_list[[t]]
    s1 <- s0 %*% A
    
    survival <- km_os_model[t] * sum(s1[c(1, 4, 5)]) /
      (sum(s1[c(1, 4, 5)]) + s0[1]*A[1,8] + s0[4]*A[4,8] + s0[5]*A[5,8])
    
    km_os_model[t + 1] <- survival
    
    if ((km_os_model[t + 1] < 0.5) && (model_median_os == 0)) {
      model_median_os <- t
      break
    }
    
    s0 <- s1
  }
  
  # Calculate likelihood from target
  v_target <- c(median_os)
  v_target_sd <- c(1)
  v_output <- c(model_median_os)
  
  nllik <- -sum(dnorm(x = v_target, mean = v_output, sd = v_target_sd, log = TRUE))
  return(nllik)
}




# ---- Wrapper for Evaluation: Optimization-Ready ----

#' Wrapper for Optimization to Evaluate OS Fit
#'
#' Simplifies `eval_variables()` interface to support use in optimization routines.
#'
#' @param x Single parameter: r_prog_nodrug2death (numeric)
#' @return Negative log-likelihood from model fit

eval_variables2 <- function(x) {
  eval_variables(
    r_prog_nodrug2death = x,
    proportion_r_pf2prog_drug = proportion_r_pf2prog_drug,
    seven_state_var = seven_state_var
  )
}


# ---- Optimization to Calibrate OS ----

#' Run Optimization to Fit Death Rate for No-Drug Progression
#'
#' Uses `optimize()` to calibrate `r_prog_nodrug2death` such that model-predicted
#' OS matches target (16.8 months). Result is stored in `opt_var`.

fit_out <- optimize(eval_variables2, c(0.01, 0.7))

# Optimal parameter value
opt_var <- fit_out$minimum
r_prog_nodrug2death <- opt_var



# ---- Helper Functions: Transition Variables by Sequence ----

#' TDxD-Chemo Variable Configuration
#'
#' Constructs rate parameters assuming TDxD followed by Chemo.
#' Includes progression, AE, and ILD rates.
#'
#' @param r_prog_nodrug2death Death rate for palliative patients
#' @param seven_state_var Global model parameters
#' @param proportion_r_pf2prog_drug Proportion receiving second-line drug
#' @return Vector of model input parameters

tdxd_chemo_var <- function(r_prog_nodrug2death, seven_state_var, proportion_r_pf2prog_drug) {
  r_pf2prog_drug   <- seven_state_var[1] * seven_state_var[3] * proportion_r_pf2prog_drug
  r_pf2prog_nodrug <- seven_state_var[1] * seven_state_var[3] * (1 - proportion_r_pf2prog_drug)
  
  gamma_AE <- seven_state_var[6]; alpha_AE <- seven_state_var[7]
  gamma_ILD <- seven_state_var[8]; alpha_ILD <- seven_state_var[9]
  
  return(c(r_pf2prog_drug, r_pf2prog_nodrug, r_prog_nodrug2death,
           gamma_AE, alpha_AE, gamma_ILD, alpha_ILD))
}


#' Chemo-TDxD Variable Configuration
#'
#' Constructs parameters for Chemo followed by TDxD.
#'
#' @return Vector of model input parameters

chemo_tdxd_var <- function(r_prog_nodrug2death, seven_state_var, proportion_r_pf2prog_drug) {
  r_pf2prog_drug   <- seven_state_var[1] * proportion_r_pf2prog_drug
  r_pf2prog_nodrug <- seven_state_var[1] * (1 - proportion_r_pf2prog_drug)
  gamma_AE <- seven_state_var[4]; alpha_AE <- seven_state_var[5]
  gamma_ILD <- 0; alpha_ILD <- 0
  
  return(c(r_pf2prog_drug, r_pf2prog_nodrug, r_prog_nodrug2death,
           gamma_AE, alpha_AE, gamma_ILD, alpha_ILD))
}


#' TDxD-SG Variable Configuration
#'
#' Constructs parameters for TDxD followed by SG treatment.
#'
#' @return Vector of model input parameters

tdxd_sg_var <- function(r_prog_nodrug2death, seven_state_var, proportion_r_pf2prog_drug) {
  r_pf2prog_drug   <- seven_state_var[1] * seven_state_var[3] * proportion_r_pf2prog_drug
  r_pf2prog_nodrug <- seven_state_var[1] * seven_state_var[3] * (1 - proportion_r_pf2prog_drug)
  gamma_AE <- seven_state_var[6]; alpha_AE <- seven_state_var[7]
  gamma_ILD <- seven_state_var[8]; alpha_ILD <- seven_state_var[9]
  
  return(c(r_pf2prog_drug, r_pf2prog_nodrug, r_prog_nodrug2death,
           gamma_AE, alpha_AE, gamma_ILD, alpha_ILD))
}






# ---- Transition Matrix Generator (Generalized for 8-State Model) ----

#' Generate Transition Matrices for Various 8-State Sequences
#'
#' Builds monthly transition probability matrices depending on specified first- and second-line treatment.
#'
#' @param r_prog_drug2prog_nodrug Rate of progression from second-line drug to palliative care
#' @param proportion_r_pf2prog_drug Proportion of first-line PF patients who receive second-line drug
#' @param seven_state_var Vector of model parameters from 7-state optimization
#' @param calibrate Boolean: Use shorter time horizon for calibration
#' @param firstline_2line Character string (e.g., "TDxD-Chemo", "Chemo-TDxD", "TDxD-SG")
#' @return List of transition matrices over time steps

transition_matrices <- function(r_prog_drug2prog_nodrug,
                                proportion_r_pf2prog_drug,
                                seven_state_var,
                                calibrate = TRUE,
                                firstline_2line = "TDxD-Chemo") {
  
  # Select appropriate transition rates based on treatment sequence
  variables <- switch(firstline_2line,
                      "TDxD-Chemo" = tdxd_chemo_var(r_prog_nodrug2death, seven_state_var, proportion_r_pf2prog_drug),
                      "Chemo-TDxD" = chemo_tdxd_var(r_prog_nodrug2death, seven_state_var, proportion_r_pf2prog_drug),
                      "TDxD-SG"    = tdxd_sg_var(r_prog_nodrug2death, seven_state_var, proportion_r_pf2prog_drug)
  )
  
  r_pf2prog_drug   <- variables[1]
  r_pf2prog_nodrug <- variables[2]
  r_prog_nodrug2death <- variables[3]
  gamma_AE <- variables[4]; alpha_AE <- variables[5]
  gamma_ILD <- variables[6]; alpha_ILD <- variables[7]
  
  # Background mortality (monthly)
  r_OC_death <- rep(
    c(0.000455117, 0.000667712, 0.000939403, 0.001487732,
      0.002481213, 0.004368355, 0.008008381, 0.014505686, 0.024617632),
    each = 60
  )
  
  n <- if (calibrate) 40 else 121
  
  # Time-varying rates from Weibull hazards
  Ft_AE <- alpha_AE * (1 - exp(-gamma_AE * (1:(n + 1))))
  Ft_ILD <- alpha_ILD * (1 - exp(-gamma_ILD * (1:(n + 1))))
  Ht_AE <- -log(1 - Ft_AE)
  Ht_ILD <- -log(1 - Ft_ILD)
  
  r_pf2AE <- r_pf2ILD <- c()
  for (t in 2:(n+1)) {
    r_pf2AE  <- append(r_pf2AE, Ht_AE[t] - Ht_AE[t - 1])
    r_pf2ILD <- append(r_pf2ILD, Ht_ILD[t] - Ht_ILD[t - 1])
  }
  
  # Constants from 7-state model
  r_pfAE2progAE <- 0.15362337
  r_progAE2death <- 0.07171597
  r_pfILD2progILD <- r_pfAE2progAE
  r_pf2progILD_death <- r_progAE2death
  
  # Store transition matrices
  P <- vector("list", n)
  
  for (t in 1:n) {
    min_pf <- r_OC_death[t] + r_pf2prog_nodrug + r_pf2prog_drug + r_pf2ILD[t] + r_pf2AE[t]
    
    G <- matrix(c(
      -min_pf, r_pf2AE[t], r_pf2ILD[t], r_pf2prog_drug, r_pf2prog_nodrug, 0, 0, r_OC_death[t],
      0, -(r_OC_death[t] + r_pfAE2progAE), 0, 0, 0, r_pfAE2progAE, 0, r_OC_death[t],
      0, 0, -(r_OC_death[t] + r_pfILD2progILD), 0, 0, 0, r_pfILD2progILD, r_OC_death[t],
      0, 0, 0, -(r_OC_death[t] + r_prog_drug2prog_nodrug), r_prog_drug2prog_nodrug, 0, 0, r_OC_death[t],
      0, 0, 0, 0, -(r_OC_death[t] + r_prog_nodrug2death), 0, 0, r_OC_death[t] + r_prog_nodrug2death,
      0, 0, 0, 0, 0, -(r_OC_death[t] + r_progAE2death), 0, r_OC_death[t] + r_progAE2death,
      0, 0, 0, 0, 0, 0, -(r_OC_death[t] + r_pf2progILD_death), r_OC_death[t] + r_pf2progILD_death,
      0, 0, 0, 0, 0, 0, 0, 0
    ),
    ncol = 8, byrow = TRUE,
    dimnames = list(
      c("PF", "PF_AE", "PF_ILD", "Prog_drug", "Prog_nodrug", "Prog_AE", "Prog_ILD", "Death"),
      c("PF", "PF_AE", "PF_ILD", "Prog_drug", "Prog_nodrug", "Prog_AE", "Prog_ILD", "Death")
    ))
    
    P[[t]] <- expm(G)
  }
  
  return(P)
}




# ---- Progression Rate Calibration (Second-Line to Palliative) ----

#' Calibrate Progression Rate from Drug to Palliative Care
#'
#' Simulates a simplified 3-state model (PF -> Palliative -> Death) to match
#' median progression-free survival (PFS) target.
#'
#' @param r_prog_drug2prog_nodrug Rate of progression to palliation
#' @param median_pf Target PFS duration (in months)
#' @return Negative log-likelihood of fit to median PFS
calc_r_prog_drug2prog_nodrug <- function(r_prog_drug2prog_nodrug, median_pf = 3.7) {
  
  # Monthly background mortality rates
  r_OC_death <- rep(
    c(0.000455117, 0.000667712, 0.000939403, 0.001487732, 0.002481213,
      0.004368355, 0.008008381, 0.014505686, 0.024617632),
    each = 60
  )
  
  n <- 40  # Time horizon (months)
  P <- vector("list", n)
  
  # Build generator matrix for 3-state model (PF, Palliative, Death)
  for (t in 1:n) {
    G <- matrix(c(
      -(r_prog_drug2prog_nodrug + r_OC_death[t]), r_prog_drug2prog_nodrug, r_OC_death[t],
      0, -0.04853368, 0.04853368,
      0, 0, 0
    ),
    ncol = 3, byrow = TRUE,
    dimnames = list(c("PF", "P", "D"), c("PF", "P", "D")))
    
    P[[t]] <- expm(G)
  }
  
  s0 <- c(1, 0, 0)
  km_pf_model <- numeric(n + 1)
  km_pf_model[1] <- s0[1]
  model_pf <- 0
  
  for (t in 1:n) {
    A <- P[[t]]
    s1 <- s0 %*% A
    km_pf_model[t + 1] <- km_pf_model[t] * s1[1] / (s1[1] + s0[1] * A[1, 2] + s0[1] * A[1, 3])
    
    if ((km_pf_model[t + 1] < 0.5) && (model_pf == 0)) {
      model_pf <- t
    }
    
    s0 <- s1
  }
  
  # Likelihood comparison
  v_target <- c(median_pf)
  v_target_sd <- c(1)
  v_output <- c(model_pf)
  nllik <- -sum(dnorm(x = v_target, mean = v_output, sd = v_target_sd, log = TRUE))
  return(nllik)
}


# ---- Evaluation Wrappers by Sequence ----

#' Evaluate Fit for Chemo -> TDxD Sequence
#' @param r_prog_drug2prog_nodrug Progression rate to palliation
#' @param median_pf Target median PFS (default: 9.9 months)
#' @return Negative log-likelihood from model fit
eval_chemo_tdxd <- function(r_prog_drug2prog_nodrug, median_pf = 9.9) {
  calc_r_prog_drug2prog_nodrug(r_prog_drug2prog_nodrug, median_pf)
}

#' Evaluate Fit for TDxD -> Chemo Sequence
#' @param r_prog_drug2prog_nodrug Progression rate to palliation
#' @param median_pf Target median PFS (default: 4 months)
#' @return Negative log-likelihood from model fit
eval_tdxd_chemo <- function(r_prog_drug2prog_nodrug, median_pf = 4) {
  calc_r_prog_drug2prog_nodrug(r_prog_drug2prog_nodrug, median_pf)
}

#' Evaluate Fit for TDxD -> SG Sequence
#' @param r_prog_drug2prog_nodrug Progression rate to palliation
#' @param median_pf Target median PFS (default: 7 months)
#' @return Negative log-likelihood from model fit
eval_tdxd_sg <- function(r_prog_drug2prog_nodrug, median_pf = 7) {
  calc_r_prog_drug2prog_nodrug(r_prog_drug2prog_nodrug, median_pf)
}


# ---- Optimization for Each Sequence ----

# Optimize to fit TDxD -> Chemo scenario
fit_out_tdxd_chemo <- optimize(eval_tdxd_chemo, c(0.01, 0.5))
opt_var_tdxd_chemo <- fit_out_tdxd_chemo$minimum

# Optimize to fit Chemo -> TDxD scenario
fit_out_chemo_tdxd <- optimize(eval_chemo_tdxd, c(0.01, 0.5))
opt_var_chemo_tdxd <- fit_out_chemo_tdxd$minimum

# Optimize to fit TDxD -> SG scenario
fit_out_tdxd_sg <- optimize(eval_tdxd_sg, c(0.01, 0.5))
opt_var_tdxd_sg <- fit_out_tdxd_sg$minimum




# ---- Transition Matrix Outputs for CEA File ----

#' Final Transition Matrices by Treatment Sequence
#'
#' These matrices are generated using calibrated parameters and exported for use in the
#' cost-effectiveness analysis (CEA) file.

# Chemo -> Chemo sequence
tm_chemo_chemo <- transition_matrices_chemo2chemo(
  r_prog_nodrug2death,
  proportion_r_pf2prog_drug,
  seven_state_var,
  calibrate = FALSE
)

# TDxD -> Chemo sequence
tm_tdxd_chemo <- transition_matrices(
  opt_var_tdxd_chemo,
  proportion_r_pf2prog_drug,
  seven_state_var,
  calibrate = FALSE
)

# Chemo -> TDxD sequence
tm_chemo_tdxd <- transition_matrices(
  opt_var_chemo_tdxd,
  proportion_r_pf2prog_drug,
  seven_state_var,
  calibrate = FALSE,
  firstline_2line = "Chemo-TDxD"
)

# TDxD -> SG sequence
tm_tdxd_sg <- transition_matrices(
  opt_var_tdxd_sg,
  proportion_r_pf2prog_drug,
  seven_state_var,
  calibrate = FALSE,
  firstline_2line = "TDxD-SG"
)


