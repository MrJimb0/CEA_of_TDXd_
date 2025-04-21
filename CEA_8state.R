# Authors: Created by Marcus Moen, MS and edited by James Dickersonm MD, MS.
# Oversight: Jeremy Goldhaber-Feibert, PhD and Fernando Alarid-Escobedo, PhD

# ---- Libraries and Options ----
options(scipen = 999)  # Avoid scientific notation

# ---- Cost and Utility Inputs ----

#' Cost and QALY Inputs for Each Treatment Strategy
#'
#' These vectors define per-cycle costs and utilities by health state for each treatment path.

# TDxD → Chemo
cost_tdxd_chemo <- c(
  cost_pf = 14113.690,
  cost_p_drug = 7203.56,
  cost_p_nodrug = 10349.86159,
  cost_pfAE = 5093.37,
  cost_pAE = 10349.86159,
  cost_pfILD = 5093.37,
  cost_pILD = 10349.86159,
  additional_cost = 7461.730261
)

qaly_tdxd_chemo <- c(
  qaly_pf = 0.65,
  qaly_p_drug = 0.54,
  qaly_p_nodrug = 0.40,
  qaly_pfAE = 0.547,
  qaly_pAE = 0.54,
  qaly_pfILD = 0.547,
  qaly_pILD = 0.54,
  decrement_qaly_ae = 0.05604
)

# Chemo → TDxD
cost_chemo_tdxd <- c(
  cost_pf = 7203.56,
  cost_p_drug = 14113.69,
  cost_p_nodrug = 10349.86159,
  cost_pfAE = 5093.37,
  cost_pAE = 10349.86159,
  cost_pfILD = 5093.37,
  cost_pILD = 10349.86159,
  additional_cost = 10435.36588
)

qaly_chemo_tdxd <- c(
  qaly_pf = 0.65,
  qaly_p_drug = 0.54,
  qaly_p_nodrug = 0.40,
  qaly_pfAE = 0.547,
  qaly_pAE = 0.54,
  qaly_pfILD = 0.547,
  qaly_pILD = 0.54,
  decrement_qaly_ae = 0.0714327
)

# TDxD → SG
cost_tdxd_sg <- c(
  cost_pf = 14113.69,
  cost_p_drug = 23970.46,
  cost_p_nodrug = 10349.86159,
  cost_pfAE = 5093.37,
  cost_pAE = 10349.86159,
  cost_pfILD = 5093.37,
  cost_pILD = 10349.86159,
  additional_cost = 7461.730261
)

qaly_tdxd_sg <- c(
  qaly_pf = 0.65,
  qaly_p_drug = 0.54,
  qaly_p_nodrug = 0.40,
  qaly_pfAE = 0.547,
  qaly_pAE = 0.54,
  qaly_pfILD = 0.547,
  qaly_pILD = 0.54,
  decrement_qaly_ae = 0.05604
)

# Chemo → Chemo
cost_chemo_chemo <- c(
  cost_pf = 7203.56,
  cost_p_drug = 5093.37,
  cost_p_nodrug = 10349.86159,
  cost_pfAE = 5093.37,
  cost_pAE = 10349.86159,
  cost_pfILD = 5093.37,
  cost_pILD = 10349.86159,
  additional_cost = 10435.36588
)

qaly_chemo_chemo <- c(
  qaly_pf = 0.65,
  qaly_p_drug = 0.54,
  qaly_p_nodrug = 0.40,
  qaly_pfAE = 0.547,
  qaly_pAE = 0.54,
  qaly_pfILD = 0.547,
  qaly_pILD = 0.54,
  decrement_qaly_ae = 0.0714327
)



# ---- 8-State Markov Model Function ----

#' Run 8-State Markov Simulation
#'
#' Simulates patient transitions across 8 health states using a time-dependent transition matrix list.
#' Returns both wide-format state probabilities and a long-format version for plotting.
#'
#' @param A_list List of transition matrices (one per month)
#' @param n_cycles Integer, number of cycles (months) to simulate
#' @return List: [1] state probabilities, [2] plotting data, [3] discount rates
cea8state <- function(A_list, n_cycles) {
  # Initialize storage
  df <- data.frame(
    cycle = 0:n_cycles,
    ProgressionFree = 0:n_cycles,
    ProgressionFreeAE = 0:n_cycles,
    ProgressionFreeILD = 0:n_cycles,
    Progressed_drug = 0:n_cycles,
    Progressed_nodrug = 0:n_cycles,
    ProgressedAE = 0:n_cycles,
    ProgressedILD = 0:n_cycles,
    Dead = 0:n_cycles
  )
  
  df_plot <- data.frame(
    cycle = rep(0:n_cycles, each = 8),
    state = rep(c("ProgressionFree", "ProgressionFreeAE", "ProgressionFreeILD",
                  "Progressed_drug", "Progressed_nodrug", "ProgressedAE",
                  "ProgressedILD", "Dead"), times = n_cycles + 1),
    value = numeric((n_cycles + 1) * 8)
  )
  
  # Discounting (3% annually)
  cycle_length <- 1 / 12
  d_c <- d_e <- 0.03
  v_dwc <- 1 / ((1 + d_c * cycle_length) ^ (0:n_cycles))  # cost
  dr_v <- 1 / ((1 + d_e * cycle_length) ^ (0:n_cycles))   # QALYs
  
  # Initial state
  s0 <- matrix(c(1, 0, 0, 0, 0, 0, 0, 0), ncol = 8)
  df[1, 2:9] <- s0
  df_plot$value[1:8] <- s0
  
  # Run model
  for (t in 1:n_cycles) {
    s1 <- s0 %*% A_list[[t]]
    df[t + 1, 2:9] <- s1
    df_plot$value[((t * 8) + 1):((t + 1) * 8)] <- s1
    s0 <- s1
  }
  
  return(list(df, df_plot, dr_v))
}


# ---- Import and Initialization ----
source("finding_transition_probabilities_8state.R") # Find transition matrices from previous code


# ---- Simulation Setup ----

#' Import Transition Matrices and Run Simulations
#'
#' This section loads the model matrices from the previous script and runs each of the
#' four treatment strategies for a 10-year (120-cycle) horizon.

n_cycles <- 120

# Store imported transition matrices
A_tdxd_chemo <- tm_tdxd_chemo
A_chemo_chemo <- tm_chemo_chemo
A_chemo_tdxd <- tm_chemo_tdxd
A_tdxd_sg <- tm_tdxd_sg

# Run model for each treatment arm
df_list_tdxd_chemo <- cea8state(A_tdxd_chemo, n_cycles)
df_tdxd_chemo <- df_list_tdxd_chemo[[1]]
df_plot_tdxd_chemo <- df_list_tdxd_chemo[[2]]

df_list_chemo_chemo <- cea8state(A_chemo_chemo, n_cycles)
df_chemo_chemo <- df_list_chemo_chemo[[1]]
df_plot_chemo_chemo <- df_list_chemo_chemo[[2]]

df_list_chemo_tdxd <- cea8state(A_chemo_tdxd, n_cycles)
df_chemo_tdxd <- df_list_chemo_tdxd[[1]]
df_plot_chemo_tdxd <- df_list_chemo_tdxd[[2]]

df_list_tdxd_sg <- cea8state(A_tdxd_sg, n_cycles)
df_tdxd_sg <- df_list_tdxd_sg[[1]]
df_plot_tdxd_sg <- df_list_tdxd_sg[[2]]



# ---- Summary Function for Cost and QALY Calculation ----

#' Calculate Total Costs, QALYs, and Life Years
#'
#' Applies state-specific cost and QALY weights to Markov state probabilities
#' and incorporates discounting. Also includes one-time additional costs and
#' subtracts AE-related QALY decrements.
#'
#' @param df DataFrame of state probabilities per cycle (from `cea8state`)
#' @param cost_data Named vector of per-state and additional costs
#' @param qaly_data Named vector of per-state utilities and AE decrement
#' @param dr_v Discount factor vector for QALYs
#' @return Named vector with total and discounted costs, QALYs, and life-years
calc_summary_data <- function(df, cost_data, qaly_data, dr_v) {
  # Extract costs
  cost_pf <- cost_data[1]
  cost_p_drug <- cost_data[2]
  cost_p_nodrug <- cost_data[3]
  cost_pfAE <- cost_data[4]
  cost_pAE <- cost_data[5]
  cost_pfILD <- cost_data[6]
  cost_pILD <- cost_data[7]
  additional_cost <- cost_data[8]
  
  # Extract QALYs
  qaly_pf <- qaly_data[1]
  qaly_p_drug <- qaly_data[2]
  qaly_p_nodrug <- qaly_data[3]
  qaly_pfAE <- qaly_data[4]
  qaly_pAE <- qaly_data[5]
  qaly_pfILD <- qaly_data[6]
  qaly_pILD <- qaly_data[7]
  decrement_qaly_ae <- qaly_data[8]
  
  # Cost calculations (undiscounted + discounted)
  cost <- (sum(df$ProgressionFree) * cost_pf +
             sum(df$Progressed_drug) * cost_p_drug +
             sum(df$Progressed_nodrug) * cost_p_nodrug +
             sum(df$ProgressionFreeAE) * cost_pfAE +
             sum(df$ProgressedAE) * cost_pAE +
             sum(df$ProgressionFreeILD) * cost_pfILD +
             sum(df$ProgressedAE) * cost_pILD) + additional_cost
  
  cost_d <- (sum(df$ProgressionFree * dr_v) * cost_pf +
               sum(df$Progressed_drug * dr_v) * cost_p_drug +
               sum(df$Progressed_nodrug * dr_v) * cost_p_nodrug +
               sum(df$ProgressionFreeAE * dr_v) * cost_pfAE +
               sum(df$ProgressedAE * dr_v) * cost_pAE +
               sum(df$ProgressionFreeILD * dr_v) * cost_pfILD +
               sum(df$ProgressedAE * dr_v) * cost_pILD) + additional_cost
  
  # QALY calculations (undiscounted + discounted)
  qaly <- (sum(df$ProgressionFree) * qaly_pf +
             sum(df$Progressed_drug) * qaly_p_drug +
             sum(df$Progressed_nodrug) * qaly_p_nodrug +
             sum(df$ProgressionFreeAE) * qaly_pfAE +
             sum(df$ProgressedAE) * qaly_pAE +
             sum(df$ProgressionFreeILD) * qaly_pfILD +
             sum(df$ProgressedAE) * qaly_pILD) / 12 - decrement_qaly_ae
  
  qaly_d <- (sum(df$ProgressionFree * dr_v) * qaly_pf +
               sum(df$Progressed_drug * dr_v) * qaly_p_drug +
               sum(df$Progressed_nodrug * dr_v) * qaly_p_nodrug +
               sum(df$ProgressionFreeAE * dr_v) * qaly_pfAE +
               sum(df$ProgressedAE * dr_v) * qaly_pAE +
               sum(df$ProgressionFreeILD * dr_v) * qaly_pfILD +
               sum(df$ProgressedAE * dr_v) * qaly_pILD) / 12 - decrement_qaly_ae
  
  # Life-years (undiscounted)
  ly <- (sum(df$ProgressionFree) +
           sum(df$Progressed_drug) +
           sum(df$Progressed_nodrug) +
           sum(df$ProgressionFreeAE) +
           sum(df$ProgressedAE) +
           sum(df$ProgressionFreeILD) +
           sum(df$ProgressedILD)) / 12
  
  return(c(cost = cost, cost_d = cost_d, qaly = qaly, qaly_d = qaly_d, ly = ly))
}


# ---- Run Summary for Each Strategy ----

res_vec_tdxd_chemo <- calc_summary_data(df_tdxd_chemo, cost_tdxd_chemo, qaly_tdxd_chemo, df_list_tdxd_chemo[[3]])
res_vec_chemo_tdxd <- calc_summary_data(df_chemo_tdxd, cost_chemo_tdxd, qaly_chemo_tdxd, df_list_chemo_tdxd[[3]])
res_vec_chemo_chemo <- calc_summary_data(df_chemo_chemo, cost_chemo_chemo, qaly_chemo_chemo, df_list_chemo_chemo[[3]])
res_vec_tdxd_sg <- calc_summary_data(df_tdxd_sg, cost_tdxd_sg, qaly_tdxd_sg, df_list_tdxd_sg[[3]])


# ---- Compile Results Table ----

df_res <- data.frame(
  Strategy = c("Chemo-Chemo:", "Chemo-TDxd:", "TDxD-Chemo:", "TDxd-SG:"),
  Costs = c(res_vec_chemo_chemo[1], res_vec_chemo_tdxd[1], res_vec_tdxd_chemo[1], res_vec_tdxd_sg[1]),
  DiscountedCost = c(res_vec_chemo_chemo[2], res_vec_chemo_tdxd[2], res_vec_tdxd_chemo[2], res_vec_tdxd_sg[2]),
  LifeYears = c(res_vec_chemo_chemo[5], res_vec_chemo_tdxd[5], res_vec_tdxd_chemo[5], res_vec_tdxd_sg[5]),
  QALYs = c(res_vec_chemo_chemo[3], res_vec_chemo_tdxd[3], res_vec_tdxd_chemo[3], res_vec_tdxd_sg[3]),
  DiscountedQALY = c(res_vec_chemo_chemo[4], res_vec_chemo_tdxd[4], res_vec_tdxd_chemo[4], res_vec_tdxd_sg[4]),
  Incremental_dCosts = c(0, (res_vec_chemo_tdxd[2]-res_vec_chemo_chemo[2]), 
                         (res_vec_tdxd_chemo[2]-res_vec_chemo_tdxd[2]), 
                         (res_vec_tdxd_sg[2]-res_vec_tdxd_chemo[2])),
  Incremental_dQALYs = c(0, (res_vec_chemo_tdxd[4]-res_vec_chemo_chemo[4]), 
                         (res_vec_tdxd_chemo[4]-res_vec_chemo_tdxd[4]), 
                         (res_vec_tdxd_sg[4]-res_vec_tdxd_chemo[4])),
  ICER = c(0, (res_vec_chemo_tdxd[2] - res_vec_chemo_chemo[2])/(res_vec_chemo_tdxd[4] - res_vec_chemo_chemo[4]), 
           (res_vec_tdxd_chemo[2]-res_vec_chemo_tdxd[2])/(res_vec_tdxd_chemo[4] - res_vec_chemo_tdxd[4]), 
           (res_vec_tdxd_sg[2] - res_vec_tdxd_chemo[2])/(res_vec_tdxd_sg[4] - res_vec_tdxd_chemo[4]))
)
