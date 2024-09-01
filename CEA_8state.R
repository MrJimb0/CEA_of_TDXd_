#This code was created by Marcus Moen, MS and edited by James Dickerson, MD, MS.
#Oversight from Jeremy Goldhaber-Feibert, PhD and Fernando Alarid-Escudero, Ph.D.

options(scipen=999)

#INPUT DATA

tdxd_chemo_cost <- c(cost_pf = 14113.690, 
                     cost_p_drug = 7203.56, 
                     cost_p_nodrug = 10349.86159, 
                     cost_pfAE = 5093.37, 
                     cost_pAE = 10349.86159, 
                     cost_pfILD = 5093.37, 
                     cost_pILD = 10349.86159,
                     additional_cost = 7461.730261)
tdxd_chemo_qaly <- c(qaly_pf = 0.65, 
                     qaly_p_drug = 0.54, 
                     qaly_p_nodrug = 0.40, 
                     qaly_pfAE = 0.547, 
                     qaly_pAE = 0.54, 
                     qaly_pfILD = 0.547, 
                     qaly_pILD = 0.54, 
                     decrement_qaly_ae = 0.05604)

chemo_tdxd_cost <- c(cost_pf = 7203.56, 
                     cost_p_drug = 14113.69, 
                     cost_p_nodrug = 10349.86159, 
                     cost_pfAE = 5093.37, 
                     cost_pAE = 10349.86159, 
                     cost_pfILD = 5093.37, 
                     cost_pILD = 10349.86159,
                     additional_cost = 10435.36588)
chemo_tdxd_qaly <- c(qaly_pf = 0.65, 
                     qaly_p_drug = 0.54, 
                     qaly_p_nodrug = 0.40, 
                     qaly_pfAE = 0.547, 
                     qaly_pAE = 0.54, 
                     qaly_pfILD = 0.547, 
                     qaly_pILD = 0.54, 
                     decrement_qaly_ae = 0.0714327)

tdxd_sg_cost <- c(cost_pf = 14113.69, 
                  cost_p_drug = 23970.46, 
                  cost_p_nodrug = 10349.86159, 
                  cost_pfAE = 5093.37, 
                  cost_pAE = 10349.86159, 
                  cost_pfILD = 5093.37, 
                  cost_pILD = 10349.86159,
                  additional_cost = 7461.730261)
tdxd_sg_qaly <- c(qaly_pf = 0.65, 
                  qaly_p_drug = 0.54, 
                  qaly_p_nodrug = 0.40, 
                  qaly_pfAE = 0.547, 
                  qaly_pAE = 0.54, 
                  qaly_pfILD = 0.547, 
                  qaly_pILD = 0.54, 
                  decrement_qaly_ae = 0.05604)

chemo_chemo_cost <- c(cost_pf = 7203.56, 
                      cost_p_drug = 5093.37, 
                      cost_p_nodrug = 10349.86159, 
                      cost_pfAE = 5093.37, 
                      cost_pAE = 10349.86159, 
                      cost_pfILD = 5093.37, 
                      cost_pILD = 10349.86159,
                      additional_cost = 10435.36588)
chemo_chemo_qaly <- c(qaly_pf = 0.65, 
                      qaly_p_drug = 0.54, 
                      qaly_p_nodrug = 0.40, 
                      qaly_pfAE = 0.547, 
                      qaly_pAE = 0.54, 
                      qaly_pfILD = 0.547, 
                      qaly_pILD = 0.54, 
                      decrement_qaly_ae = 0.0714327)





#Input:
#A = a 8x8 transition matrix
#n_cycles = number of cycles
#Function for doing n_cycles given a 8x8 transition matrix
#Output:
#A list of two DataFrames. One with the probability values for each of the three states at each cycle. And one which is used for plotting.
cea8state <- function(A_list, n_cycles){
  
  #Make dataframes to store values
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
    cycle = rep(0:n_cycles, 8),
    state = 0:(n_cycles*8+7),
    value = 0:(n_cycles*8+7)
  )
  
  #The discount rate per cycle (assuming an annual discount rate of 3%)
  #A discount rate vector (For storing the discount at each timestep)
  # FAE comments:
  cycle_length <- 1/12
  d_c <- 0.03
  d_e <- 0.03
  v_dwc  <- 1 / ((1 + (d_c * cycle_length)) ^ (0:n_cycles))
  dr_v  <- 1 / ((1 + (d_e * cycle_length)) ^ (0:n_cycles))
  
  
  ly <- 1:n_cycles
  
  #The initial state
  s0 <- matrix(c(1,0,0,0,0,0,0,0), ncol = 8)

  
  #The Markov simulation
  for (t in 1:(n_cycles)) {
    
    #Update the discount rate vector
    #dr_v[t+1] <- (1/(1+0.03))^((t)*cycle_length)
    
    #Do one cycle
    s1 <- s0 %*% A_list[[t]]
    #Update the dataframes
    df["ProgressionFree"][df["ProgressionFree"] == t] <- s1[1,1]
    df["ProgressionFreeAE"][df["ProgressionFreeAE"] == t] <- s1[1,2]
    df["ProgressionFreeILD"][df["ProgressionFreeILD"] == t] <- s1[1,3]
    df["Progressed_drug"][df["Progressed_drug"] == t] <- s1[1,4]
    df["Progressed_nodrug"][df["Progressed_nodrug"] == t] <- s1[1,5]
    df["ProgressedAE"][df["ProgressedAE"] == t] <- s1[1,6]
    df["ProgressedILD"][df["ProgressedILD"] == t] <- s1[1,7]
    df["Dead"][df["Dead"] == t] <- s1[1,8] 
    
    df_plot["state"][df_plot["value"] == t] <- "ProgressionFree"
    df_plot["state"][df_plot["value"] == t + n_cycles+1] <- "ProgressionFreeAE"
    df_plot["state"][df_plot["value"] == t + n_cycles*2+2] <- "ProgressionFreeILD"
    df_plot["state"][df_plot["value"] == t + n_cycles*3+3] <- "Progressed_drug"
    df_plot["state"][df_plot["value"] == t + n_cycles*4+4] <- "Progressed_nodrug"
    df_plot["state"][df_plot["value"] == t + n_cycles*5+5] <- "ProgressedAE"
    df_plot["state"][df_plot["value"] == t + n_cycles*6+6] <- "ProgressedILD"
    df_plot["state"][df_plot["value"] == t + n_cycles*7+7] <- "Dead"
    df_plot["value"][df_plot["value"] == t] <- s1[1,1]
    df_plot["value"][df_plot["value"] == t + n_cycles+1] <- s1[1,2]
    df_plot["value"][df_plot["value"] == t + n_cycles*2+2] <- s1[1,3]
    df_plot["value"][df_plot["value"] == t + n_cycles*3+3] <- s1[1,4]
    df_plot["value"][df_plot["value"] == t + n_cycles*4+4] <- s1[1,5]
    df_plot["value"][df_plot["value"] == t + n_cycles*5+5] <- s1[1,6]
    df_plot["value"][df_plot["value"] == t + n_cycles*6+6] <- s1[1,7]
    df_plot["value"][df_plot["value"] == t + n_cycles*7+7] <- s1[1,8]
    
    s0 <- s1
    
  }
  
  df["ProgressionFree"][df["ProgressionFree"] == 0] <- 1
  df_plot["value"][df_plot["state"] == 0] <- 1
  df_plot["value"][df_plot["state"] == n_cycles+1] <- 0
  df_plot["value"][df_plot["state"] == n_cycles*2+2] <- 0
  df_plot["value"][df_plot["state"] == n_cycles*3+3] <- 0
  df_plot["value"][df_plot["state"] == n_cycles*4+4] <- 0
  df_plot["value"][df_plot["state"] == n_cycles*5+5] <- 0
  df_plot["value"][df_plot["state"] == n_cycles*6+6] <- 0
  df_plot["value"][df_plot["state"] == n_cycles*7+7] <- 0
  
  
  df_plot["state"][df_plot["state"] == 0] <- "ProgressionFree"
  df_plot["state"][df_plot["state"] == n_cycles+1] <- "ProgressionFreeAE"
  df_plot["state"][df_plot["state"] == n_cycles*2+2] <- "ProgressionFreeILD"
  df_plot["state"][df_plot["state"] == n_cycles*3+3] <- "Progressed_drug"
  df_plot["state"][df_plot["state"] == n_cycles*4+4] <- "Progressed_nodrug"
  df_plot["state"][df_plot["state"] == n_cycles*5+5] <- "ProgressedAE"
  df_plot["state"][df_plot["state"] == n_cycles*6+6] <- "ProgressedILD"
  df_plot["state"][df_plot["state"] == n_cycles*7+7] <- "Dead"
  
  
  return(list(df, df_plot, dr_v))
}




#Import the transition matrices from finding_transition_probabilities.R
source("finding_transition_probabilities_8state.R")

#Number of cycles. This is 120 months or 10 years. 
n_cycles <- 120

#Store the transition matrices from the "finding_transition_probabilities_8state.R" file
A_tdxd_chemo <- tm_tdxd_chemo
A_chemo_chemo <- tm_chemo_chemo
A_chemo_tdxd <- tm_chemo_tdxd
A_tdxd_sg <- tm_tdxd_sg

#Run the CEA on all the four options
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




calc_summary_data <- function(df, cost_data, qaly_data, dr_v){
  #The cost values
  cost_pf <- cost_data[1]
  cost_p_drug <- cost_data[2]
  cost_p_nodrug <- cost_data[3]
  cost_pfAE <- cost_data[4]
  cost_pAE <- cost_data[5]
  cost_pfILD <- cost_data[6]
  cost_pILD <- cost_data[7]
  additional_cost <- cost_data[8]
  
  #The QALY values
  qaly_pf <- qaly_data[1]
  qaly_p_drug <- qaly_data[2]
  qaly_p_nodrug <- qaly_data[3]
  qaly_pfAE <- qaly_data[4]
  qaly_pAE <- qaly_data[5]
  qaly_pfILD <- qaly_data[6]
  qaly_pILD <- qaly_data[7]
  decrement_qaly_ae <- qaly_data[8]
  
  
  #Calculating the cost
  cost <- (sum(df$ProgressionFree)*cost_pf + sum(df$Progressed_drug)*cost_p_drug +
             sum(df$Progressed_nodrug)*cost_p_nodrug +
             sum(df$ProgressionFreeAE)*cost_pfAE + sum(df$ProgressedAE)*cost_pAE +
             sum(df$ProgressionFreeILD)*cost_pfILD + sum(df$ProgressedAE)*cost_pILD) + additional_cost

  
  cost_d <- (sum(df$ProgressionFree*dr_v)*cost_pf + sum(df$Progressed_drug*dr_v)*cost_p_drug +
               sum(df$Progressed_nodrug*dr_v)*cost_p_nodrug +
               sum(df$ProgressionFreeAE*dr_v)*cost_pfAE + sum(df$ProgressedAE*dr_v)*cost_pAE +
               sum(df$ProgressionFreeILD*dr_v)*cost_pfILD + sum(df$ProgressedAE*dr_v)*cost_pILD) + additional_cost
  
  
  #Calculate the qalys
  qaly <- (sum(df$ProgressionFree)*qaly_pf + sum(df$Progressed_drug)*qaly_p_drug +
             sum(df$Progressed_nodrug)*qaly_p_nodrug +
             sum(df$ProgressionFreeAE)*qaly_pfAE + sum(df$ProgressedAE)*qaly_pAE +
             sum(df$ProgressionFreeILD)*qaly_pfILD + sum(df$ProgressedAE)*qaly_pILD)/12 - decrement_qaly_ae
  
  qaly_d <- (sum(df$ProgressionFree*dr_v)*qaly_pf + sum(df$Progressed_drug*dr_v)*qaly_p_drug +
               sum(df$Progressed_nodrug*dr_v)*qaly_p_nodrug +
               sum(df$ProgressionFreeAE*dr_v)*qaly_pfAE + sum(df$ProgressedAE*dr_v)*qaly_pAE +
               sum(df$ProgressionFreeILD*dr_v)*qaly_pfILD + sum(df$ProgressedAE*dr_v)*qaly_pILD)/12 - decrement_qaly_ae
  
  #Calculate total life years
  ly <- (sum(df$ProgressionFree) + sum(df$Progressed_drug) + sum(df$Progressed_nodrug) +
           sum(df$ProgressionFreeAE) + sum(df$ProgressedAE) +
           sum(df$ProgressionFreeILD) + sum(df$ProgressedILD))/12
  
  v_res <- c(cost = cost, 
             cost_d = cost_d, 
             qaly = qaly, 
             qaly_d = qaly_d, 
             ly = ly)
  
  return(v_res)
}




res_vec_tdxd_chemo <- calc_summary_data(df = df_tdxd_chemo, 
                                        cost_data = tdxd_chemo_cost, 
                                        qaly_data = tdxd_chemo_qaly, 
                                        dr_v = df_list_tdxd_chemo[[3]])




res_vec_chemo_tdxd <- calc_summary_data(df_chemo_tdxd, 
                                        chemo_tdxd_cost, 
                                        chemo_tdxd_qaly, 
                                        df_list_chemo_tdxd[[3]])




res_vec_chemo_chemo <- calc_summary_data(df_chemo_chemo, 
                                         chemo_chemo_cost, 
                                         chemo_chemo_qaly, 
                                         df_list_chemo_chemo[[3]])



res_vec_tdxd_sg <- calc_summary_data(df_tdxd_sg, 
                                     tdxd_sg_cost, 
                                     tdxd_sg_qaly, 
                                     df_list_tdxd_sg[[3]])

#Summary of the results
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



