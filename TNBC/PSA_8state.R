#This code was created by Marcus Moen, MS and edited by James Dickerson, MD, MS.
#Oversight from Jeremy Goldhaber-Feibert, PhD and Fernando Alarid-Escudero, Ph.D.

options(scipen=999)

library(expm)
library(dplyr)

#For generating the QALY parameters go to line 570
#For generating the cost parameters go to line 654

#Notation used in this file:
#pf = Progression Free first-line
#prog_drug / p_drug = Progressed from first line and given new drug in second line
#prog_nodrug / p_nodrug = Progressed from first line and put in palliative care
#AE = Adverse Event (not ILD)
#ILD = Interstitial lung disease

#Import the results and functions from other files.
source("finding_transition_probabilities_7state.R") #We get opt_var from this
source("CEA_8state.R")

#The following code is from Jeremy's paper on preference order (referenced in the manuscript, table in supplement)
# Please reference the manuscript for mathematical notations and variable descriptions.
CorrelateUtils <- function(U, Q, epsilon, delta){
  n <- nrow(U) #number of PSA samples
  s <- ncol(U) #number of states
  R <- matrix(rnorm(n*s,0,1), n, s) #the reference matrix.
  
  C <- matrix(0,s,s) #a place holder for the correlation matrix
  Viol <- matrix(0,s,s) #violations matrix.
  for (j in 2:s){ # {j,k} is the selected pair of state comparisons
    for (k in 1:(j-1)){
      rho <- 1 # #bivariate correlations
      X = U[,c(j,k)] #selected columns of U
      Y = R[,c(j,k)] #selected columns of R
      viol <- 0
      while(viol<epsilon & rho>=0){ #if these conditions are met, continue
        rho <- rho - delta #reduce correltion
        Xstar = induceRankCorrelation(X, Y, rho) #correlated utilities.
        viol = mean((Q[j,k] * Xstar[,1]) < (Q[j,k] * Xstar[,2])) #compute %violations between the col. vectors.
      }
      #Viol[j,k] <- viol
      C[j,k] <- rho + delta #record the desired correlation.
    }
    #print(j) #just to show the column indices.
    #print(k)
  }
  #Fill in the other elements of C.
  C = C + t(C)
  for (j in 1:s){
    C[j,j] <- 1 # % the diagonal of ones.
  }
  ## Eigenvectors and Eigenvalues correction of C
  eigenResults <- eigen(C)
  B <- eigenResults$values
  V <- eigenResults$vectors
  B[B<=0] <- 0.0001 #to make sure C is positive definite, set eigenvalues<=0 to a very small positive number
  Cstar <- V %*% diag(B) %*% solve(V) #reconstruct C
  Ustar <- induceRankCorrelation(U, R, Cstar) #similar to above, induce the correlation.
  return(Ustar)
}
## To induce Rank correlation: inputs X: QoL vectors, Y is the reference vectors, and Sigma is the correlation matrix.
induceRankCorrelation <- function(X, Y, Sigma){
  if (length(Sigma)==1){ #if Sigma is a single value, convert it to a 2x2 matrix.
    Sigma <- matrix(c(1, Sigma,
                      Sigma, 1), 2, 2)
  }
  n <- nrow(X)
  s <- ncol(X)
  #Initialize matrices.
  Xsorted <- matrix(0, n, s)
  Yrank <- matrix(0, n, s)
  Xstar <- matrix(0, n, s)
  
  P <- chol(Sigma) #compute the upper triangular matrix
  Ystar <- Y %*% P #Sort the values in the reference vectors by multiplying by P
  cor(Ystar)
  for (j in 1:s){
    Xsorted[,j] <- sort(X[,j]) #Sort each variable
    Yrank[order(Ystar[,j]),j] <- seq(1:n) #Reverse sort
    Xstar[,j]=Xsorted[Yrank[,j],j] #sort Xsorted to have the same ranks as Ystar.
  }
  return(Xstar) #return the sorted vectors.
} 







#Generate all the QALY parameters used in the PSA. 
#Since the QALYs are in relation with each other the number of simulations we input here is not what we end up with.
generate_psa_qaly_params <- function(n_sim, seed = 0){
  #Draw values for QALYs
  qaly_pf <- rbeta(n_sim, shape1 = 2, shape2 = 2)*(0.81-0.52) + rep(0.52, n_sim)
  qaly_p_drug <- rbeta(n_sim, shape1 = 2, shape2 = 2)*(0.65-0.297) + rep(0.297, n_sim)
  qaly_p_nodrug <- rbeta(n_sim, shape1 = 2, shape2 = 2)*((0.40*1.25)-(0.40*0.75)) + rep(0.40*0.75, n_sim)
  qaly_pfAE <- rbeta(n_sim, shape1 = 2, shape2 = 2)*(0.707-0.417) + rep(0.417, n_sim)
  qaly_pAE <- rbeta(n_sim, shape1 = 2, shape2 = 2)*(0.65-0.297) + rep(0.297, n_sim)
  qaly_pfILD <- rbeta(n_sim, shape1 = 2, shape2 = 2)*(0.707-0.417) + rep(0.417, n_sim)
  qaly_pILD <- rbeta(n_sim, shape1 = 2, shape2 = 2)*(0.65-0.297) + rep(0.297, n_sim)
  
  U <- matrix(0, n_sim, 7,
              dimnames = list(seq(1,n_sim,1), c("pf", "pfAE", "pfILD", "p_drug", "p_nodrug", "pAE", "pILD")))
  
  
  U[,1] <- qaly_pf
  U[,2] <- qaly_pfAE
  U[,3] <- qaly_pfILD
  U[,4] <- qaly_p_drug
  U[,5] <- qaly_p_nodrug
  U[,6] <- qaly_pAE
  U[,7] <- qaly_pILD
  

  Q <- matrix(c( 0, 0, 0, 0, 0, 0, 0, #preference order matrix.
                 -1, 0, 0, 0, 0, 0, 0,
                 -1, 0, 0, 0, 0, 0, 0,
                 -1, -1, -1, 0, 0, 0, 0,
                 -1, -1, -1, -1, -1, 0, 0, 
                 -1, -1, -1, 0, 0, 0, 0,
                 -1, -1, -1, 0, 0, 0, 0
  ), 7, 7, byrow = TRUE) 
  
  Ustar <- CorrelateUtils(U, Q, 0.05, 0.1)
  
  
  qaly_pf <- c()
  qaly_pfAE <- c()
  qaly_pfILD <- c()
  qaly_p_drug <- c()
  qaly_p_nodrug <- c()
  qaly_pAE <- c()
  qaly_pILD <- c()
  
  for(i in 1:n_sim){
    error <- 0
    for(j in 1:7){
      a <- Ustar[i,]*Q[j,]
      if(sum(a) == 0){
        something <- 0
      }
      else if(max(a[a<0]) > Ustar[i,j]*-1){
        error <- error +1
      }
    }

    if(error == 0){
      qaly_pf <- append(qaly_pf, Ustar[i,1])
      qaly_pfAE <- append(qaly_pfAE, Ustar[i,2])
      qaly_pfILD <- append(qaly_pfILD, Ustar[i,3])
      qaly_p_drug <- append(qaly_p_drug, Ustar[i,4])
      qaly_p_nodrug <- append(qaly_p_nodrug, Ustar[i,5])
      qaly_pAE <- append(qaly_pAE, Ustar[i,6])
      qaly_pILD <- append(qaly_pILD, Ustar[i,7])
    }
  }
  
  
  df_qaly_params <- data.frame(
    
    qaly_pf = qaly_pf,
    qaly_pfAE_chemo = qaly_pfAE,
    qaly_pfILD = qaly_pfILD,
    qaly_p_drug = qaly_p_drug,
    qaly_p_nodrug = qaly_p_nodrug,
    qaly_pAE = qaly_pAE,
    qaly_pILD = qaly_pILD
  )
  
  return(df_qaly_params)
  
}



#Generate all the cost parameters used in the PSA
generate_psa_cost_params <- function(n_sim, seed = 0){
  set.seed(seed)

  #Calculate the cost of PF in chemo:
  price_cape <- ((rbeta(n_sim, shape1 = 2, shape2 = 2)*((801.2+200)-(801.2-200)) + rep((801.2-200), n_sim))*28*2325/9000)*0.201
  price_eribulin <- ((rbeta(n_sim, shape1 = 2, shape2 = 2)*((1224.17+300)-(1224.17-300)) + rep((1224.17-300), n_sim))*2*2.6/1)
  price_abraxane <- ((rbeta(n_sim, shape1 = 2, shape2 = 8)*((935.81+900)-935.81) + rep(935.81, n_sim))*483.6/100)*0.103
  price_paclitaxel <- ((rbeta(n_sim, shape1 = 8, shape2 = 2)*((155+20)-(155-95)) + rep((155-95), n_sim))*325.5/300)*0.082
  price_gemcitabine <- ((rbeta(n_sim, shape1 = 2, shape2 = 8)*((33.81+100)-(33.81-15)) + rep((33.81-15), n_sim))*2*2325/2000)*0.103
  
  admin_price <- 303.86*0.201+831.91*0.511+454.14*0.103+454.14*0.082+831.91*0.103
  cost_pf_chemo <- (price_cape + price_eribulin*0.511 + price_abraxane + price_paclitaxel + price_gemcitabine)*4/3 + rep(admin_price, n_sim)
  
  #Cost of just eribulin
  cost_pf_eribulin <- price_eribulin*4/3 + rep(831.91, n_sim)
  
  #Calculate the cost of PF in tdxd: 
  drug_price_tdxd <- (rbeta(n_sim, shape1 = 2, shape2 = 2)*((2435.71+600)-(2435.71-600)) + rep((2435.71-600), n_sim))*418.5/100
  cost_pf_tdxd <- (drug_price_tdxd)*4/3 + rep(522.43, n_sim)

  #Calculate the cost of P in chemo and tdxd:
  cost_p_chemo <- cost_pf_chemo 
  cost_p_tdxd <- cost_pf_tdxd
  
  #SG price
  drug_price_sg <- (rbeta(n_sim, shape1 = 2, shape2 = 2)*((2087.75+700)-(2087.75-400)) + rep((2087.75-500), n_sim))*775/180
  cost_p_sg <- (drug_price_sg)*4/3 + rep(831.91, n_sim)
  
  
  
  
  
  #Calculate the cost of the AE and ILD states for both chemo and tdxd:
  cost_pfAE_chemo <- rbeta(n_sim, shape1 = 2, shape2 = 2)*((4518.18*1.25)-(4518.18*0.75)) + rep(575.19, n_sim) #+ rep(518.18, n_sim)
  cost_pAE_chemo <- rgamma(n_sim, shape = (10349.861591)^2/(2114.74)^2, scale = ((2114.74)^2/10349.861591))#cost_p_chemo
  cost_pfAE_tdxd <- cost_pfAE_chemo
  cost_pAE_tdxd <- cost_pAE_chemo
  cost_pfILD_tdxd <-cost_pfAE_chemo
  cost_pILD_tdxd <-cost_pAE_chemo
  
  #ADD TERMINAL STATE COST
  cost_p_nodrug <- rbeta(n_sim, shape1 = 2, shape2 = 2)*((10349.86159*1.25)-(10349.86159*0.75)) + rep(10349.86159*0.75, n_sim)
  

  df_cost_params <- data.frame(
    cost_pf_tdxd = cost_pf_tdxd,
    cost_pf_chemo = cost_pf_chemo,
    cost_p_tdxd = cost_p_tdxd,
    cost_p_chemo = cost_p_chemo,
    cost_p_sg = cost_p_sg,
    cost_pfAE_tdxd = cost_pfAE_tdxd,
    cost_pfAE_chemo = cost_pfAE_chemo,
    cost_pAE_tdxd = cost_pAE_tdxd,
    cost_pAE_chemo = cost_pAE_chemo,
    cost_pfILD_tdxd = cost_pfILD_tdxd,
    cost_pILD_tdxd = cost_pILD_tdxd, 
    cost_p_nodrug = cost_p_nodrug
  )

  return(df_cost_params)
}






#Generate all the qaly and cost parameters used in the PSA
PSA_params <- function(n_sim, seed = 0){
  df_qaly_params <- generate_psa_qaly_params(n_sim)
  n_sim2 <- length(df_qaly_params$qaly_pf)
  df_cost_params <- generate_psa_cost_params(n_sim2)
  
  df_params <- cbind(df_qaly_params, df_cost_params)
  
  return(df_params)
}



#Draw a value from the log normal distribution
draw_log_normal <- function(mu, sigma){
  location <- log(mu^2 / sqrt(sigma^2 + mu^2))
  shape <- sqrt(log(1 + (sigma^2 / mu^2)))
  draw <- rlnorm(n=1, location, shape)
  return(draw)
}





#Main function for running the PSA
run_PSA <- function(df_params){
  
  median_os_chemo_v <- c()
  median_pf_chemo_v <- c()
  median_os_tdxd_v <- c()
  median_pf_tdxd_v <- c()
  median_os_super_tdxd_v <- c()
  median_pf_super_tdxd_v <- c()
  target_ae_chemo_v <- c()
  target_ae_tdxd_v <- c()
  target_ild_tdxd_v <- c()
  
  n_sims <- length(df_params$qaly_pf)
  
  
  #DFs to store the results:
  df_psa_chemo_chemo <- data.frame(
    DiscountedCost = 1:n_sims,
    DiscountedQALY = 1:n_sims
  )
  
  df_psa_tdxd_chemo <- data.frame(
    DiscountedCost = 1:n_sims,
    DiscountedQALY = 1:n_sims
  )
  
  df_psa_chemo_tdxd <- data.frame(
    DiscountedCost = 1:n_sims,
    DiscountedQALY = 1:n_sims
  )
  
  df_psa_tdxd_sg <- data.frame(
    DiscountedCost = 1:n_sims,
    DiscountedQALY = 1:n_sims
  )
  
  
  
  
  
  df_psa_res <- data.frame(
    sim = 1:(n_sims*5),
    DiscountedCost = 1:(n_sims*5),
    DiscountedQALY = 1:(n_sims*5),
    LY = 1:(n_sims*5),
    group = 1:(n_sims*5),
    ICER = 1:(n_sims*5)
  )
  
  df_psa_jimbo <- data.frame(
    sim = 1:n_sims,
    ChemoChemo_dQALY = 1:n_sims,
    ChemoChemo_dCOST = 1:n_sims,
    ChemoChemo_LY = 1:n_sims,
    TDxdChemo_dQALY = 1:n_sims,
    TDxdChemo_dCOST = 1:n_sims,
    TDxdChemo_LY = 1:n_sims,
    SuperTDxd_dQALY = 1:n_sims,
    SuperTDxd_dCOST = 1:n_sims,
    SuperTDxd_LY = 1:n_sims,
    ChemoTDxd_dQALY = 1:n_sims,
    ChemoTDxd_dCOST = 1:n_sims,
    ChemoTDxd_LY = 1:n_sims,
    TDxdSG_dQALY = 1:n_sims,
    TDxdSG_dCOST = 1:n_sims,
    TDxdSG_LY = 1:n_sims
  )
  
  
  
  n_cycles <- 120
  
  for(i in 1:n_sims){
    print(i)
    
    #mu is the mean and sigma is the sd found in the excel file made by Jimbo
    #The goal is to change this such that we do not have to run the whole calibration
    #every time.
    median_os_chemo <- draw_log_normal(mu = 8.3, sigma = 3.83)
    median_pf_chemo <- draw_log_normal(mu = 2.9, sigma = 0.94)
    median_os_tdxd <- draw_log_normal(mu = 18.2, sigma = 2.35)
    median_pf_tdxd <- draw_log_normal(mu = 8.5, sigma = 1.89)
    
    median_os_super_tdxd <- draw_log_normal(mu = 18.2*3, sigma = 2.35)
    median_pf_super_tdxd <- draw_log_normal(mu = 8.5*3, sigma = 1.89)
    ##Change this in your sensitivity analysis
    median_pf_sg <- draw_log_normal(mu = (5.1*1.1), sigma = (5.1*1.1)*0.15)
    
    target_ae_chemo <- draw_log_normal(mu = 0.081, sigma = 0.0203)
    target_ae_tdxd <- draw_log_normal(mu = 0.0757, sigma = 0.0189)
    target_ild_tdxd <- draw_log_normal(mu = 0.0863, sigma = 0.0216)
    
    median_os_chemo_v <- append(median_os_chemo_v, median_os_chemo)
    median_pf_chemo_v <- append(median_pf_chemo_v, median_pf_chemo)
    median_os_tdxd_v <- append(median_os_tdxd_v, median_os_tdxd)
    median_pf_tdxd_v <- append(median_pf_tdxd_v, median_pf_tdxd)
    median_os_super_tdxd_v <- append(median_os_tdxd_v, median_os_tdxd)
    median_pf_super_tdxd_v <- append(median_pf_tdxd_v, median_pf_tdxd)
    target_ae_chemo_v <- append(target_ae_chemo_v, target_ae_chemo)
    target_ae_tdxd_v <- append(target_ae_tdxd_v, target_ae_tdxd)
    target_ild_tdxd_v <- append(target_ild_tdxd_v, target_ild_tdxd)

    opt_var <- find_opt_var_7state(target = c(median_os_chemo, median_pf_chemo, median_os_tdxd, median_pf_tdxd, target_ae_chemo, target_ae_tdxd, target_ild_tdxd))
    
    opt_var_super_tdxd <- find_opt_var_7state(target = c(median_os_chemo, median_pf_chemo, median_os_super_tdxd, median_pf_super_tdxd, target_ae_chemo, target_ae_tdxd, target_ild_tdxd))

    seven_state_var <- c(r_pf2prog = opt_var[1],
                         r_prog2death = opt_var[2],
                         hr_pf2prog_chemo2tdxd = opt_var[3],
                         gamma_chemo = opt_var[4],
                         alpha_chemo = opt_var[5],
                         gamma_tdxd = opt_var[6],
                         alpha_tdxd = opt_var[7],
                         gamma_ILD = opt_var[8],
                         alpha_ILD = opt_var[9])
    
    seven_state_var_super_tdxd <- c(r_pf2prog = opt_var_super_tdxd[1],
                         r_prog2death = opt_var_super_tdxd[2],
                         hr_pf2prog_chemo2tdxd = opt_var_super_tdxd[3],
                         gamma_chemo = opt_var_super_tdxd[4],
                         alpha_chemo = opt_var_super_tdxd[5],
                         gamma_tdxd = opt_var_super_tdxd[6],
                         alpha_tdxd = opt_var_super_tdxd[7],
                         gamma_ILD = opt_var_super_tdxd[8],
                         alpha_ILD = opt_var_super_tdxd[9])

    ##Change this in your sensitivity analysis
    proportion_r_pf2prog_drug <- rbeta(1, shape1 = 57, shape2 = 38)
    
    eval_variables2 <- function(x){
      return(eval_variables(x, proportion_r_pf2prog_drug, seven_state_var))
    }
    
    eval_chemo_tdxd2 <- function(r_prog_drug2prog_nodrug){
      return(eval_chemo_tdxd(r_prog_drug2prog_nodrug, median_pf_tdxd))
    }
    eval_tdxd_chemo2 <- function(r_prog_drug2prog_nodrug){
      return(eval_tdxd_chemo(r_prog_drug2prog_nodrug, median_pf_chemo))
    }
    eval_tdxd_sg2 <- function(r_prog_drug2prog_nodrug){
      return(eval_tdxd_sg(r_prog_drug2prog_nodrug, median_pf_sg))
    }
    
    fit_out <- optimize(eval_variables2,c(0.01,0.7))
    r_prog_nodrug2death <- fit_out$minimum

    fit_out_tdxd_chemo <- optimize(eval_tdxd_chemo2, c(0.01,0.5))
    opt_var_tdxd_chemo <- fit_out_tdxd_chemo$minimum

    fit_out_chemo_tdxd <- optimize(eval_chemo_tdxd2, c(0.01,0.5))
    opt_var_chemo_tdxd <- fit_out_chemo_tdxd$minimum

    fit_out_tdxd_sg <- optimize(eval_tdxd_sg2, c(0.01,0.5))
    opt_var_tdxd_sg <- fit_out_tdxd_sg$minimum
    

    #Find the optimal transition matrices (These lists will be send to tge CEA_8state.R file)
    tm_chemo_chemo <- transition_matrices_chemo2chemo(r_prog_nodrug2death, proportion_r_pf2prog_drug, seven_state_var, calibrate=F)
    tm_tdxd_chemo <- transition_matrices(opt_var_tdxd_chemo, proportion_r_pf2prog_drug, seven_state_var, calibrate=F)
    tm_chemo_tdxd <- transition_matrices(opt_var_chemo_tdxd, proportion_r_pf2prog_drug, seven_state_var, calibrate=F, firstline_2line = "Chemo-TDxD")
    tm_tdxd_sg <- transition_matrices(opt_var_tdxd_sg, proportion_r_pf2prog_drug, seven_state_var, calibrate=F, firstline_2line = "TDxD-SG")
    
    tm_super_tdxd_chemo  <- transition_matrices(opt_var_tdxd_chemo, 
                                                proportion_r_pf2prog_drug,
                                                seven_state_var_super_tdxd,
                                                calibrate = F)

    
    #Store the transition matrices from the "finding_transition_probabilities_8state.R" file
    A_tdxd_chemo <- tm_tdxd_chemo
    A_chemo_chemo <- tm_chemo_chemo
    A_chemo_tdxd <- tm_chemo_tdxd
    A_tdxd_sg <- tm_tdxd_sg
    
    A_super_tdxd <- tm_super_tdxd_chemo
    

    #Run the CEA on all the four options
    df_list_tdxd_chemo <- cea8state(A_tdxd_chemo, n_cycles)
    df_tdxd_chemo <- df_list_tdxd_chemo[[1]]

    df_list_chemo_chemo <- cea8state(A_chemo_chemo, n_cycles)
    df_chemo_chemo <- df_list_chemo_chemo[[1]]

    df_list_chemo_tdxd <- cea8state(A_chemo_tdxd, n_cycles)
    df_chemo_tdxd <- df_list_chemo_tdxd[[1]]

    df_list_tdxd_sg <- cea8state(A_tdxd_sg, n_cycles)
    df_tdxd_sg <- df_list_tdxd_sg[[1]]
    
    df_list_super_tdxd <- cea8state(A_super_tdxd, n_cycles)
    df_super_tdxd <- df_list_super_tdxd[[1]]

    
    
    dr_v <- df_list_chemo_chemo[[3]]
    
    
    
    tdxd_chemo_cost <- c(cost_pf = df_params[i,]$cost_pf_tdxd, 
                         cost_p_drug = df_params[i,]$cost_p_chemo, 
                         cost_p_nodrug = df_params[i,]$cost_p_nodrug, 
                         cost_pfAE = df_params[i,]$cost_pfAE_tdxd, 
                         cost_pAE = df_params[i,]$cost_pAE_tdxd, 
                         cost_pfILD = df_params[i,]$cost_pfILD_tdxd, 
                         cost_pILD = df_params[i,]$cost_pILD_tdxd,
                         additional_cost = 7461.730261)
    tdxd_chemo_qaly <- c(qaly_pf = df_params[i,]$qaly_pf, 
                         qaly_p_drug = df_params[i,]$qaly_p_drug, 
                         qaly_p_nodrug = df_params[i,]$qaly_p_nodrug, 
                         qaly_pfAE = df_params[i,]$qaly_pfAE, 
                         qaly_pAE = df_params[i,]$qaly_pAE, 
                         qaly_pfILD = df_params[i,]$qaly_pfILD, 
                         qaly_pILD = df_params[i,]$qaly_pILD, 
                         decrement_qaly_ae = 0.05604)
    
    
    
    
    chemo_tdxd_cost <- c(cost_pf = df_params[i,]$cost_pf_chemo, 
                         cost_p_drug = df_params[i,]$cost_p_tdxd, 
                         cost_p_nodrug = df_params[i,]$cost_p_nodrug, 
                         cost_pfAE = df_params[i,]$cost_pfAE_chemo, 
                         cost_pAE = df_params[i,]$cost_pAE_chemo, 
                         cost_pfILD = 0, 
                         cost_pILD = 0,
                         additional_cost = 10435.36588)
    chemo_tdxd_qaly <- c(qaly_pf = df_params[i,]$qaly_pf, 
                         qaly_p_drug = df_params[i,]$qaly_p_drug, 
                         qaly_p_nodrug = df_params[i,]$qaly_p_nodrug, 
                         qaly_pfAE = df_params[i,]$qaly_pfAE, 
                         qaly_pAE = df_params[i,]$qaly_pAE, 
                         qaly_pfILD = df_params[i,]$qaly_pfILD, 
                         qaly_pILD = df_params[i,]$qaly_pILD, 
                         decrement_qaly_ae = 0.0714327)
    
    
    chemo_chemo_cost <- c(cost_pf = df_params[i,]$cost_pf_chemo, 
                          cost_p_drug = df_params[i,]$cost_p_chemo, 
                          cost_p_nodrug = df_params[i,]$cost_p_nodrug, 
                          cost_pfAE = df_params[i,]$cost_pfAE_chemo, 
                          cost_pAE = df_params[i,]$cost_pAE_chemo, 
                          cost_pfILD = 0, 
                          cost_pILD = 0,
                          additional_cost = 10435.36588)
    chemo_chemo_qaly <- c(qaly_pf = df_params[i,]$qaly_pf, 
                          qaly_p_drug = df_params[i,]$qaly_p_drug, 
                          qaly_p_nodrug = df_params[i,]$qaly_p_nodrug, 
                          qaly_pfAE = df_params[i,]$qaly_pfAE, 
                          qaly_pAE = df_params[i,]$qaly_pAE, 
                          qaly_pfILD = df_params[i,]$qaly_pfILD, 
                          qaly_pILD = df_params[i,]$qaly_pILD, 
                          decrement_qaly_ae = 0.0714327)
    
    
    tdxd_sg_cost <- c(cost_pf = df_params[i,]$cost_pf_tdxd, 
                      cost_p_drug = df_params[i,]$cost_p_sg, 
                      cost_p_nodrug = df_params[i,]$cost_p_nodrug, 
                      cost_pfAE = df_params[i,]$cost_pfAE_tdxd, 
                      cost_pAE = df_params[i,]$cost_pAE_tdxd, 
                      cost_pfILD = df_params[i,]$cost_pfILD_tdxd, 
                      cost_pILD = df_params[i,]$cost_pILD_tdxd,
                      additional_cost = 7461.730261)
    tdxd_sg_qaly <- c(qaly_pf = df_params[i,]$qaly_pf, 
                      qaly_p_drug = df_params[i,]$qaly_p_drug, 
                      qaly_p_nodrug = df_params[i,]$qaly_p_nodrug, 
                      qaly_pfAE = df_params[i,]$qaly_pfAE, 
                      qaly_pAE = df_params[i,]$qaly_pAE, 
                      qaly_pfILD = df_params[i,]$qaly_pfILD, 
                      qaly_pILD = df_params[i,]$qaly_pILD, 
                      decrement_qaly_ae = 0.05604)
    
    
    
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
    
    res_vec_super_tdxd <- calc_summary_data(df = df_super_tdxd, 
                                            cost_data = tdxd_chemo_cost, 
                                            qaly_data = tdxd_chemo_qaly, 
                                            dr_v = df_list_super_tdxd[[3]])
    
    
    
    df_psa_jimbo[i,]$sim <- i
    df_psa_jimbo[i,]$ChemoChemo_dQALY <- res_vec_chemo_chemo[[4]]
    df_psa_jimbo[i,]$ChemoChemo_dCOST <- res_vec_chemo_chemo[[2]]
    df_psa_jimbo[i,]$ChemoChemo_LY <- res_vec_chemo_chemo[[5]]
    
    df_psa_jimbo[i,]$TDxdChemo_dQALY <- res_vec_tdxd_chemo[[4]]
    df_psa_jimbo[i,]$TDxdChemo_dCOST <- res_vec_tdxd_chemo[[2]]
    df_psa_jimbo[i,]$TDxdChemo_LY <- res_vec_tdxd_chemo[[5]]
    
    df_psa_jimbo[i,]$SuperTDxd_dQALY <- res_vec_super_tdxd[[4]]
    df_psa_jimbo[i,]$SuperTDxd_dCOST <- res_vec_super_tdxd[[2]]
    df_psa_jimbo[i,]$SuperTDxd_LY <- res_vec_super_tdxd[[5]]
    
    df_psa_jimbo[i,]$ChemoTDxd_dQALY <- res_vec_chemo_tdxd[[4]]
    df_psa_jimbo[i,]$ChemoTDxd_dCOST <- res_vec_chemo_tdxd[[2]]
    df_psa_jimbo[i,]$ChemoTDxd_LY <- res_vec_chemo_tdxd[[5]]
    
    df_psa_jimbo[i,]$TDxdSG_dQALY <- res_vec_tdxd_sg[[4]]
    df_psa_jimbo[i,]$TDxdSG_dCOST <- res_vec_tdxd_sg[[2]]
    df_psa_jimbo[i,]$TDxdSG_LY <- res_vec_tdxd_sg[[5]]
    
    
    
    df_psa_res[(5*i-4),]$DiscountedCost <- res_vec_super_tdxd[[2]]
    df_psa_res[(5*i-4),]$LY <- res_vec_super_tdxd[[5]]
    df_psa_res[(5*i-4),]$DiscountedQALY <- res_vec_super_tdxd[[4]]
    df_psa_res[(5*i-4),]$group <- "Super_TDXd"
    df_psa_res[(5*i-4),]$sim <- i
    
    df_psa_res[(5*i-3),]$DiscountedCost <- res_vec_chemo_chemo[[2]]
    df_psa_res[(5*i-3),]$LY <- res_vec_chemo_chemo[[5]]
    df_psa_res[(5*i-3),]$DiscountedQALY <- res_vec_chemo_chemo[[4]]
    df_psa_res[(5*i-3),]$group <- "Chemo-Chemo"
    df_psa_res[(5*i-3),]$sim <- i
    
    df_psa_res[(5*i-2),]$DiscountedCost <- res_vec_chemo_tdxd[[2]]
    df_psa_res[(5*i-2),]$LY <- res_vec_chemo_tdxd[[5]]
    df_psa_res[(5*i-2),]$DiscountedQALY <- res_vec_chemo_tdxd[[4]]
    df_psa_res[(5*i-2),]$group <- "Chemo-TDxd"
    df_psa_res[(5*i-2),]$sim <- i
    
    df_psa_res[(5*i-1),]$DiscountedCost <- res_vec_tdxd_chemo[[2]]
    df_psa_res[(5*i-1),]$LY <- res_vec_tdxd_chemo[[5]]
    df_psa_res[(5*i-1),]$DiscountedQALY <- res_vec_tdxd_chemo[[4]]
    df_psa_res[(5*i-1),]$group <- "TDxD-Chemo"
    df_psa_res[(5*i-1),]$sim <- i
    
    df_psa_res[(5*i),]$DiscountedCost <- res_vec_tdxd_sg[[2]]
    df_psa_res[(5*i),]$LY <- res_vec_tdxd_sg[[5]]
    df_psa_res[(5*i),]$DiscountedQALY <- res_vec_tdxd_sg[[4]]
    df_psa_res[(5*i),]$group <- "TDxD-SG"
    df_psa_res[(5*i),]$sim <- i
    
    
    
  
    df_psa_res[(5*i-3),]$ICER <- "NA"
    df_psa_res[(5*i-3),]$ICER <- 0
    df_psa_res[(5*i-2),]$ICER <- (res_vec_chemo_tdxd[2]-res_vec_chemo_chemo[2])/(res_vec_chemo_tdxd[4]-res_vec_chemo_chemo[4])
    df_psa_res[(5*i-1),]$ICER <- (res_vec_tdxd_chemo[2]-res_vec_chemo_tdxd[2])/(res_vec_tdxd_chemo[4]-res_vec_chemo_tdxd[4])
    df_psa_res[(5*i),]$ICER <- (res_vec_tdxd_sg[2]-res_vec_tdxd_chemo[2])/(res_vec_tdxd_sg[4]-res_vec_tdxd_chemo[4])
    

    
    
  }
  
  df_params_full <- df_params
  df_params_full["median_os_chemo"] <- median_os_chemo_v
  df_params_full["median_pf_chemo"] <- median_pf_chemo_v
  df_params_full["median_os_tdxd"] <- median_os_tdxd_v
  df_params_full["median_pf_tdxd"] <- median_pf_tdxd_v
  df_params_full["target_ae_chemo"] <- target_ae_chemo_v
  df_params_full["target_ae_tdxd"] <- target_ae_tdxd_v
  df_params_full["target_ild_tdxd"] <- target_ild_tdxd_v
  


  return(list(df_psa_res, df_params_full, df_psa_jimbo))
}








#IMPORTANT NUMBER. Change this to change the final number of simulations
n = 750

df_param <- PSA_params(n)


print("Number of simulations:")
print(nrow(df_param))


res <- run_PSA(df_param)

View(res[[1]])
jimbo_final_df <- res[[1]]
jimbo_final_df$ICER <- NULL

jimbo_final_df <- jimbo_final_df %>%
  filter(sim <= n)

# hist(jimbo_final_df$DiscountedCost[jimbo_final_df$group == 'Chemo-Chemo'])
# hist(jimbo_final_df$DiscountedCost[jimbo_final_df$group == 'Chemo-TDxd'])
# hist(jimbo_final_df$DiscountedCost[jimbo_final_df$group == 'TDxD-Chemo'])
# hist(jimbo_final_df$DiscountedCost[jimbo_final_df$group == 'TDxD-SG'])
# 
# hist(jimbo_final_df$DiscountedQALY[jimbo_final_df$group == 'Chemo-Chemo'])
# hist(jimbo_final_df$DiscountedQALY[jimbo_final_df$group == 'Chemo-TDxd'])
# hist(jimbo_final_df$DiscountedQALY[jimbo_final_df$group == 'TDxD-Chemo'])
# hist(jimbo_final_df$DiscountedQALY[jimbo_final_df$group == 'TDxD-SG'])
# 
# mean(jimbo_final_df$DiscountedCost[jimbo_final_df$group == 'Chemo-Chemo'])
# mean(jimbo_final_df$DiscountedCost[jimbo_final_df$group == 'Chemo-TDxd'])
# mean(jimbo_final_df$DiscountedCost[jimbo_final_df$group == 'TDxD-Chemo'])
# mean(jimbo_final_df$DiscountedCost[jimbo_final_df$group == 'TDxD-SG'])
# 
# mean(jimbo_final_df$DiscountedQALY[jimbo_final_df$group == 'Chemo-Chemo'])
# mean(jimbo_final_df$DiscountedQALY[jimbo_final_df$group == 'Chemo-TDxd'])
# mean(jimbo_final_df$DiscountedQALY[jimbo_final_df$group == 'TDxD-Chemo'])
# mean(jimbo_final_df$DiscountedQALY[jimbo_final_df$group == 'TDxD-SG'])

write.csv(jimbo_final_df, file = "data/base_case_output_TNBC.csv")

#Output Table 
# Calculate mean by group
mean_by_group <- jimbo_final_df %>%
  group_by(group) %>%
  summarize(mean_DiscountedCost = mean(DiscountedCost),
            mean_DiscountedQALY = mean(DiscountedQALY),
            mean_LY = mean(LY))


