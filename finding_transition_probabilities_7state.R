#This code was created by Marcus Moen, MS and edited by James Dickerson, MD, MS.
#Oversight from Jeremy Goldhaber-Feibert, PhD and Fernando Alarid-Escudero, Ph.D.

options(scipen=999)

#The target variables in line 161 are the values we are calibrating to. Change this if you want to change the calibration
#The sd corresponding to the targets is found at line 296

#Function for calculating the standard error of a survival curve given 
#number at risk, n_risk, number of deaths at each timestep, n_deaths, and the survival probabilities, surv_prob

library(expm)
library(ggplot2)

#Function for finding the transition matrices such that they can be imported to the CEA file
#Takes as its input the optimal values for the eval_variables function and outputs a list of transition matrices
transition_matrices_7state <- function(input_var, calibrate=T){
  #print("RUNNING")
  
  r_pf2prog <- input_var[1] #Rate 1
  r_prog2death <- input_var[2] #Rate 2
  hr_pf2prog_chemo2tdxd <- input_var[3]
  
  #Variables for the AEs in chemo
  gamma_chemo <- input_var[4]
  alpha_chemo <- input_var[5]
  
  #Variables for the AEs in Tdxd
  gamma_tdxd <- input_var[6]
  alpha_tdxd <- input_var[7]
  
  #Variables for the ILDs in Tdxd
  gamma_ILD <- input_var[8]
  alpha_ILD <- input_var[9]
  
  #Rate of dying from other causes (background mortality)
  r_OC_death <- rep(c(0.000455117,0.000667712,0.000939403,0.001487732,0.002481213, 0.004368355, 0.008008381,0.014505686, 0.024617632)
                    ,each=(12*5))
  
  #Takes to long to calibrate if we use the whole length (we do not need the whole lenght to calibrate)
  if(calibrate){
    n <- 40
  }else{
    #n <- length(r_OC_death)
    n <- 121
  }
  #The rate of getting an AE and being discontinued. (Needs updating)
  Ft_chemo <- c()
  Ft_tdxd <- c()
  Ft_tdxd_ILD <- c()
  
  Ht_chemo <- c()
  Ht_tdxd <- c()
  Ht_tdxd_ILD <- c()
  
  r_AE <- c()
  r_AE_tdxd <- c()
  r_ILD <- c()
  r_ILD_tdxd <- c()
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
  
  
  #Defining some hazard ratios (Dummy variables for now, maybe changed later)
  hr_pfAE2progAE <- 1
  hr_progAE2death <- 1
  hr_pfae2prog_chemo2tdxd <- 1
  hr_pf2prog2death_chemo2tdxd <- 1
  hr_pfAE2progAE2death_chemo2tdxd <- 1
  hr_pf2progILD_chemo2tdxd <- 1
  hr_pf2progILD2death_chemo2tdxd <- 1
  
  
  #Finding the relationship between the rates
  r_pfAE2progAE <- r_pf2prog*hr_pfAE2progAE
  r_pf2progILD <- r_pf2prog*0
  r_pf2progILD_death <- r_prog2death*0
  r_progAE2death <- r_prog2death*hr_progAE2death
  
  
  r_pf2prog_tdxd <- r_pf2prog*hr_pf2prog_chemo2tdxd
  r_pfAE2progAE_tdxd <- r_pf2prog*hr_pfae2prog_chemo2tdxd
  r_prog2death_tdxd <- r_prog2death*hr_pf2prog2death_chemo2tdxd
  r_progAE2death_tdxd <- r_prog2death*hr_pfAE2progAE2death_chemo2tdxd
  r_pf2progILD_tdxd <- r_pf2prog*hr_pf2progILD_chemo2tdxd
  r_pf2progILD_death_tdxd <- r_prog2death*hr_pf2progILD2death_chemo2tdxd
  
  P_chemo <- list()
  P_tdxd <- list()
  
  for(t in 1:n){
    #Creating the transition matrices
    G_chemo <- matrix(c(-(r_pf2prog + r_OC_death[t] + r_AE[t] + r_ILD[t]), r_AE[t], r_ILD[t], r_pf2prog, 0, 0, r_OC_death[t],
                        0, -(r_pfAE2progAE + r_OC_death[t]), 0, 0, r_pfAE2progAE, 0, r_OC_death[t],
                        0, 0, -(r_pf2progILD + r_OC_death[t]), 0, 0, r_pf2progILD, r_OC_death[t],
                        0, 0, 0, -(r_OC_death[t] + r_prog2death), 0, 0, (r_OC_death[t] + r_prog2death),
                        0, 0, 0, 0, -(r_OC_death[t] + r_progAE2death), 0, (r_OC_death[t] + r_progAE2death),
                        0, 0, 0, 0, 0, -(r_OC_death[t] + r_pf2progILD_death), (r_OC_death[t] + r_pf2progILD_death),
                        0, 0, 0, 0, 0, 0, 0), 
                      ncol = 7, nrow = 7,
                      dimnames = list(c("PF", "PF_AE", "PF_ILD", "Prog", "Prog_AE", "Prog_ILD", "Death"),
                                      c("PF", "PF_AE", "PF_ILD", "Prog", "Prog_AE", "Prog_ILD", "Death")),
                      byrow = T
    )
    #print("This is G_chemo:")
    #print(G_chemo)
    
    G_tdxd <- matrix(c(-(r_pf2prog_tdxd + r_OC_death[t] + r_AE_tdxd[t] + r_ILD_tdxd[t]), r_AE_tdxd[t], r_ILD_tdxd[t], r_pf2prog_tdxd, 0, 0, r_OC_death[t],
                        0, -(r_pfAE2progAE_tdxd + r_OC_death[t]), 0, 0, r_pfAE2progAE_tdxd, 0, r_OC_death[t],
                        0, 0, -(r_pf2progILD_tdxd + r_OC_death[t]), 0, 0, r_pf2progILD_tdxd, r_OC_death[t],
                        0, 0, 0, -(r_OC_death[t] + r_prog2death_tdxd), 0, 0, (r_OC_death[t] + r_prog2death_tdxd),
                        0, 0, 0, 0, -(r_OC_death[t] + r_progAE2death_tdxd), 0, (r_OC_death[t] + r_progAE2death_tdxd),
                        0, 0, 0, 0, 0, -(r_OC_death[t] + r_pf2progILD_death_tdxd), (r_OC_death[t] + r_pf2progILD_death_tdxd),
                        0, 0, 0, 0, 0, 0, 0), 
                      ncol = 7, nrow = 7,
                      dimnames = list(c("PF", "PF_AE", "PF_ILD", "Prog", "Prog_AE", "Prog_ILD", "Death"),
                                      c("PF", "PF_AE", "PF_ILD", "Prog", "Prog_AE", "Prog_ILD", "Death")),
                      byrow = T
    )
    
    #Find paper that proves this!
    A_chemo <- expm(G_chemo)
    #print("This is A_chemo:")
    #print(A_chemo)
    A_tdxd <- expm(G_tdxd)
    
    P_chemo[[t]] <- A_chemo
    P_tdxd[[t]] <- A_tdxd
    
  }
  
  matrices <- list(P_chemo, P_tdxd)
  return(matrices)
}




#Input:
#input_var = a list of three values, the rate from PF to P, the rate from P to D, and the HR.
#Output:
#The square error when using the given inputs.
eval_variables_7state <- function(input_var, target_var = c(median_os_chemo = 16.8, 
                                                     median_pf_chemo = 5.1, 
                                                     median_os_tdxd = 23.4, 
                                                     median_pf_tdxd = 9.9, 
                                                     target_ae_tdxd_chemo = 0.081, 
                                                     target_ae_tdxd = 0.0757, 
                                                     target_ild_tdxd = 0.0863)){
  
  
  r_pf2prog <- input_var[1] #Rate 1
  r_prog2death <- input_var[2] #Rate 2
  hr_pf2prog_chemo2tdxd <- input_var[3]
  
  gamma_chemo <- input_var[4]
  alpha_chemo <- input_var[5]
  
  gamma_tdxd <- input_var[6]
  alpha_tdxd <- input_var[7]
  
  gamma_ILD <- input_var[8]
  alpha_ILD <- input_var[9]
  
  
  matrix_list <- transition_matrices_7state(input_var)
  A_chemo_list <- matrix_list[[1]]
  A_tdxd_list <- matrix_list[[2]]
  
  #target_var <- tail(input_var,7)
  

  #The extracted Kaplan-Meier values
  median_os_chemo <- target_var[1]
  median_pf_chemo <- target_var[2]
  median_os_tdxd <- target_var[3]
  median_pf_tdxd <- target_var[4]
  
  target_ae_chemo <- target_var[5]
  target_ae_tdxd <- target_var[6]
  target_ild_tdxd <- target_var[7]
  
  
  
  #Vectors for storing the estimates of the KM curves
  km_os_chemo_model <- c()
  km_pf_chemo_model <- c()
  km_os_tdxd_model <- c()
  km_pf_tdxd_model <- c()
  
  cum_AE_chemo_model <- c()
  cum_AE_tdxd_model <- c()
  cum_ILD_model <- c()
  
  
  n = 30
  
  #Start states
  s0_chemo <- c(1,0,0,0,0,0,0)
  s0_tdxd <- c(1,0,0,0,0,0,0)
  
  survival_chemo <- s0_chemo[1]+s0_chemo[4]
  km_os_chemo_model <- append(km_os_chemo_model, survival_chemo)
  km_pf_chemo_model <- append(km_pf_chemo_model, s0_chemo[1])
  
  survival_tdxd <- s0_tdxd[1]+s0_tdxd[4]
  km_os_tdxd_model <- append(km_os_tdxd_model, survival_tdxd)
  km_pf_tdxd_model <- append(km_pf_tdxd_model, s0_tdxd[1])
  
  cum_AE_chemo_model <- append(cum_AE_chemo_model, 0)
  cum_AE_tdxd_model <- append(cum_AE_tdxd_model, 0)
  cum_ILD_model <- append(cum_ILD_model, 0)
  
  
  model_os_chemo <- 0
  model_pf_chemo <- 0
  model_os_tdxd <- 0
  model_pf_tdxd <- 0
  
  #Do the simulations
  for(t in 1:n){
    
    A_chemo <- A_chemo_list[[t]]
    A_tdxd <- A_tdxd_list[[t]]
    
    
    
    s1_chemo <- s0_chemo %*% A_chemo
    s1_tdxd <- s0_tdxd %*% A_tdxd
    
    #survival_chemo <- s1_chemo[1]+s1_chemo[4]
    #km_os_chemo_model <- append(km_os_chemo_model, survival_chemo)
    #km_pf_chemo_model <- append(km_pf_chemo_model, s1_chemo[1])
    
    #survival_tdxd <- s1_tdxd[1]+s1_tdxd[4]
    #km_os_tdxd_model <- append(km_os_tdxd_model, survival_tdxd)
    #km_pf_tdxd_model <- append(km_pf_tdxd_model, s1_tdxd[1])
    
    survival_chemo <- km_os_chemo_model[t]*(s1_chemo[1]+s1_chemo[4])/(s1_chemo[1]+s1_chemo[4]+s0_chemo[1]*A_chemo[1,7]+s0_chemo[4]*A_chemo[4,7])
    km_os_chemo_model <- append(km_os_chemo_model, survival_chemo)
    km_pf_chemo_model <- append(km_pf_chemo_model, km_pf_chemo_model[t]*s1_chemo[1]/(s1_chemo[1]+s0_chemo[1]*A_chemo[1,4]+s0_chemo[1]*A_chemo[1,7]))
    
    survival_tdxd <- km_os_tdxd_model[t]*(s1_tdxd[1]+s1_tdxd[4])/(s1_tdxd[1]+s1_tdxd[4]+s0_tdxd[1]*A_tdxd[1,7]+s0_tdxd[4]*A_tdxd[4,7])
    km_os_tdxd_model <- append(km_os_tdxd_model, survival_tdxd)
    km_pf_tdxd_model <- append(km_pf_tdxd_model, km_pf_tdxd_model[t]*(s1_tdxd[1])/(s1_tdxd[1]+s0_tdxd[1]*A_tdxd[1,4]+s0_tdxd[1]*A_tdxd[1,7]))
    
    cum_AE_chemo_model <- append(cum_AE_chemo_model, cum_AE_chemo_model[t]+s0_chemo[1]*A_chemo[1,2])
    cum_AE_tdxd_model <- append(cum_AE_tdxd_model, cum_AE_tdxd_model[t]+s0_tdxd[1]*A_tdxd[1,2])
    cum_ILD_model <- append(cum_ILD_model, cum_ILD_model[t]+s0_tdxd[1]*A_tdxd[1,3])
    
    if((km_os_chemo_model[t+1]<0.5) && (model_os_chemo == 0)){
      model_os_chemo <- t
    }
    if((km_pf_chemo_model[t+1]<0.5) && (model_pf_chemo == 0)){
      model_pf_chemo <- t
    }
    if((km_os_tdxd_model[t+1]<0.5) && (model_os_tdxd == 0)){
      model_os_tdxd <- t
    }
    if((km_pf_tdxd_model[t+1]<0.5) && (model_pf_tdxd == 0)){
      model_pf_tdxd <- t
    }
    
    
    s0_chemo <- s1_chemo
    s0_tdxd <- s1_tdxd
    
  }
  
  model_ae_chemo <- cum_AE_chemo_model[12]
  model_ae_tdxd <- cum_AE_tdxd_model[12]
  model_ild_chemo <- cum_ILD_model[12]
  
  
  #Calculate the likelihood
  
  v_target <- c(median_os_chemo, median_pf_chemo, median_os_tdxd, median_pf_tdxd, target_ae_chemo, target_ae_tdxd, target_ild_tdxd)
  v_target_sd <- c(1.4,0.66,1.22,0.59,0.0203,0.0189,0.0216)
  v_output <- c(model_os_chemo, model_pf_chemo, model_os_tdxd, model_pf_tdxd, model_ae_chemo, model_ae_tdxd, model_ild_chemo)
  
  nllik <- -sum(dnorm(x = v_target, mean = v_output, sd = v_target_sd, log = TRUE))
  
  return(nllik)
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







#Plot the correct KM-curves vs the cureves from our model. This part of the file is used for validation


#The transitions from rate to probability is found in https://journals.sagepub.com/doi/pdf/10.1177/0272989X05282637 Figure 3
plot_function <- function(input_var){
  
  r_pf2p_chemo <- input_var[1] 
  r_p2d_chemo <- input_var[2] 
  hr <- input_var[3]
  
  gamma_chemo <- input_var[4]
  alpha_chemo <- input_var[5]
  
  gamma_tdxd <- input_var[6]
  alpha_tdxd <- input_var[7]
  
  gamma_ILD <- input_var[8]
  alpha_ILD <- input_var[9]
  
  matrix_list <- transition_matrices_7state(input_var)
  A_chemo_list <- matrix_list[[1]]
  A_tdxd_list <- matrix_list[[2]]
  
  
  #The extracted Kaplan-Meier values
  km_os_chemo <- c(1, 0.986, 0.981, 0.957, 0.94, 0.927, 0.884, 0.841, 0.793, 0.738, 0.712, 0.683, 0.669, 0.634)#, 0.610, 0.566, 0.520, 0.484, 0.459, 0.460, 0.457, 0.419, 0.375)
  km_pf_chemo <- c(1, 0.983, 0.756, 0.62, 0.603, 0.501, 0.439, 0.398, 0.34, 0.277, 0.268, 0.226, 0.217, 0.185) 
  km_os_tdxd <- c(1, 0.992, 0.985, 0.967, 0.958, 0.939, 0.928, 0.899, 0.87, 0.857, 0.829, 0.815, 0.793, 0.774)#, 0.742, 0.731, 0.688, 0.661, 0.623, 0.586, 0.566, 0.539, 0.507, 0.500, 0.484)
  km_pf_tdxd <- c(1, 0.988, 0.896, 0.82, 0.809, 0.765, 0.681, 0.633, 0.592, 0.554, 0.491, 0.457, 0.422, 0.386)#, 0.371, 0.363, 0.339, 0.31, 0.295, 0.274, 0.259) 
  
  #Vectors for storing the estimates of the KM curves
  km_os_chemo_model <- c()
  km_pf_chemo_model <- c()
  km_os_tdxd_model <- c()
  km_pf_tdxd_model <- c() 
  
  
  n = max(length(km_os_chemo), length(km_pf_chemo), length(km_os_tdxd), length(km_pf_tdxd))
  
  #Start states
  s0_chemo <- c(1,0,0,0,0,0,0)
  s0_tdxd <- c(1,0,0,0,0,0,0)
  
  survival_chemo <- s0_chemo[1]+s0_chemo[4]
  km_os_chemo_model <- append(km_os_chemo_model, survival_chemo)
  km_pf_chemo_model <- append(km_pf_chemo_model, s0_chemo[1])
  
  survival_tdxd <- s0_tdxd[1]+s0_tdxd[4]
  km_os_tdxd_model <- append(km_os_tdxd_model, survival_tdxd)
  km_pf_tdxd_model <- append(km_pf_tdxd_model, s0_tdxd[1])
  
  ae_test_chemo <- c(0)
  ae_test_tdxd <- c(0)
  ae_test_ild <- c(0)
  
  
  for(t in 1:n){
    A_chemo <- A_chemo_list[[t]]
    A_tdxd <- A_tdxd_list[[t]]
    
    s1_chemo <- s0_chemo %*% A_chemo
    s1_tdxd <- s0_tdxd %*% A_tdxd
    
    #survival_chemo <- s1_chemo[1]+s1_chemo[4]
    #km_os_chemo_model <- append(km_os_chemo_model, survival_chemo)
    #km_pf_chemo_model <- append(km_pf_chemo_model, s1_chemo[1])
    
    #survival_tdxd <- s1_tdxd[1]+s1_tdxd[4]
    #km_os_tdxd_model <- append(km_os_tdxd_model, survival_tdxd)
    #km_pf_tdxd_model <- append(km_pf_tdxd_model, s1_tdxd[1])
    
    survival_chemo <- km_os_chemo_model[t]*(s1_chemo[1]+s1_chemo[4])/(s1_chemo[1]+s1_chemo[4]+s0_chemo[1]*A_chemo[1,7]+s0_chemo[4]*A_chemo[4,7])
    km_os_chemo_model <- append(km_os_chemo_model, survival_chemo)
    km_pf_chemo_model <- append(km_pf_chemo_model, km_pf_chemo_model[t]*s1_chemo[1]/(s1_chemo[1]+s0_chemo[1]*A_chemo[1,4]+s0_chemo[1]*A_chemo[1,7]))
    
    survival_tdxd <- km_os_tdxd_model[t]*(s1_tdxd[1]+s1_tdxd[4])/(s1_tdxd[1]+s1_tdxd[4]+s0_tdxd[1]*A_tdxd[1,7]+s0_tdxd[4]*A_tdxd[4,7])
    km_os_tdxd_model <- append(km_os_tdxd_model, survival_tdxd)
    km_pf_tdxd_model <- append(km_pf_tdxd_model, km_pf_tdxd_model[t]*(s1_tdxd[1])/(s1_tdxd[1]+s0_tdxd[1]*A_tdxd[1,4]+s0_tdxd[1]*A_tdxd[1,7]))
    
    ae_test_chemo <- append(ae_test_chemo, ae_test_chemo[t]+s0_chemo[1]*A_chemo[1,2])
    ae_test_tdxd <- c(ae_test_tdxd, ae_test_tdxd[t]+s0_tdxd[1]*A_tdxd[1,2])
    ae_test_ild <- c(ae_test_ild, ae_test_ild[t]+s0_tdxd[1]*A_tdxd[1,3])
    
    s0_chemo <- s1_chemo
    s0_tdxd <- s1_tdxd
    
  }
  
  
  
  m = min(length(km_os_chemo), length(km_pf_chemo), length(km_os_tdxd), length(km_pf_tdxd))
  idx <- 1:m
  df <- data.frame(idx, km_os_chemo_model[1:m], km_os_tdxd_model[1:m], km_os_chemo[1:m])
  
  
  
  
  # plot1 <- ggplot(df, aes(x=idx)) + 
  #   geom_line(aes(y = km_os_chemo_model[1:m]), color = "red") + 
  #   geom_line(aes(y = km_os_chemo[1:m]), color="blue", linetype="twodash") +
  #   ylim(0, 1) +
  #   ylab("Probability") +
  #   xlab("Month") +
  #   ggtitle("OS physician's choice")+
  #   scale_x_continuous(breaks = round(seq(0, 18, by = 3),1))+theme(axis.text=element_text(size=12),
  #                                                                  axis.title=element_text(size=14))
  # 
  # plot2 <- ggplot(df, aes(x=idx)) + 
  #   geom_line(aes(y = km_os_tdxd_model[1:m]), color = "red") + 
  #   geom_line(aes(y = km_os_tdxd[1:m]), color="blue", linetype="twodash") +
  #   ylim(0, 1) +
  #   ylab("Probability")+
  #   xlab("Month")+
  #   ggtitle("OS T-DxD")+
  #   scale_x_continuous(breaks = round(seq(0, 18, by = 3),1))+theme(axis.text=element_text(size=12),
  #                                                                  axis.title=element_text(size=14))
  # 
  # plot3 <- ggplot(df, aes(x=idx)) + 
  #   geom_line(aes(y = km_pf_chemo_model[1:m]), color = "red") + 
  #   geom_line(aes(y = km_pf_chemo[1:m]), color="blue", linetype="twodash") +
  #   ylim(0, 1) +
  #   ylab("Probability")+
  #   xlab("Month")+
  #   ggtitle("PFS physician's choice")+
  #   scale_x_continuous(breaks = round(seq(0, 18, by = 3),1))+theme(axis.text=element_text(size=12),
  #                                                                  axis.title=element_text(size=14))
  # 
  # plot4 <- ggplot(df, aes(x=idx)) + 
  #   geom_line(aes(y = km_pf_tdxd_model[1:m]), color = "red") + 
  #   geom_line(aes(y = km_pf_tdxd[1:m]), color="blue", linetype="twodash") +
  #   ylim(0, 1) +
  #   ylab("Probability")+
  #   xlab("Month")+
  #   ggtitle("PFS T-DxD")+
  #   scale_x_continuous(breaks = round(seq(0, 18, by = 3),1))+theme(axis.text=element_text(size=12),
  #                                                                  axis.title=element_text(size=14))
  
  print(ae_test_chemo)
  print(ae_test_tdxd)
  print(ae_test_ild)
  plot(0:(length(ae_test_chemo)-1), ae_test_chemo, type = "l")
  plot(0:(length(ae_test_chemo)-1), ae_test_tdxd, type = "l")
  plot(0:(length(ae_test_chemo)-1), ae_test_ild, type = "l")
  
  # ggarrange(plot1, plot2, plot3, plot4,
  #           ncol = 2, nrow = 2, common.legend = TRUE,legend="bottom") 
  
}


plot_function(opt_var)






