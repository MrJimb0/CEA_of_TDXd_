#This code was created by Marcus Moen, MS and edited by James Dickerson, MD, MS.
#Oversight from Jeremy Goldhaber-Feibert, PhD and Fernando Alarid-Escudero, Ph.D.

options(scipen=999)
# setwd("/Users/jamesdickerson/Library/CloudStorage/Box-Box/Dickerson Lab/Dickerson_Lab_Github/CEA_of_TDXd/")

library(expm)
library(ggplot2)
# library(ggpubr)


#Notation used in this file:
#pf = Progression Free first-line
#prog_drug / p_drug = Progressed from first line and given new drug in second line
#prog_nodrug / p_nodrug = Progressed from first line and put in palliative care
#AE = Adverse Event (not ILD)
#ILD = Interstitial lund disease

#Go to line 40 to change the proportion that will go on to recive a second line drug instead of palliative care
#To change the PFS for the second line go to line 439-447


#In order:
#Rate PF to progressed for chemo first line
#Rate progressed to death for both
#HR for PF to prog for TDxD
#Gamma for chemo AE
#Alpha for chemo AE
#Gamma for TDxD AE
#Alpha for TDxD AE
#Gamma for TDxD ILD
#Alpha for TDxD ILD


source("finding_transition_probabilities_7state.R") #We get opt_var from this


#GLOBAL VARIABLES (START):
seven_state_var <- c(r_pf2prog = opt_var[1],
                     r_prog2death = opt_var[2],
                     hr_pf2prog_chemo2tdxd = opt_var[3],
                     gamma_chemo = opt_var[4],
                     alpha_chemo = opt_var[5],
                     gamma_tdxd = opt_var[6],
                     alpha_tdxd = opt_var[7],
                     gamma_ILD = opt_var[8],
                     alpha_ILD = opt_var[9])

#This proportion_r_pf2prog_drug is a key variable. It sets the proportion of patients who will get the next line of therapy 
#vs those that go on with palliative care alone; needs to be varied in sensitivity 

proportion_r_pf2prog_drug <- 0.6
#GLOBAL VARIABLES (END):


#Find the variables for the chemo-chemo option
chemo_chemo_var <- function(seven_state_var, proportion_r_pf2prog_drug){
  
  #Take the r_pf2prog from the seven state variables and divided it into the two groups in the second line
  r_pf2prog_drug <- seven_state_var[1]*proportion_r_pf2prog_drug
  r_pf2prog_nodrug <- seven_state_var[1]*(1 - proportion_r_pf2prog_drug)
  
  
  #Assume that the rate for progressing in the second line is the same as progressing in the first line when you get chemo both lines
  r_prog_drug2prog_nodrug <- seven_state_var[1]
  
  #Variables for the AEs in Chemo from the seven state file
  gamma_AE <- seven_state_var[4]
  alpha_AE <- seven_state_var[5]
  
  #Variables for the ILDs in Chemo
  gamma_ILD <- 0
  alpha_ILD <- 0
  
  return(c(r_pf2prog_drug, r_pf2prog_nodrug, r_prog_drug2prog_nodrug, gamma_AE, alpha_AE, gamma_ILD, alpha_ILD))
}





#Function for finding the transition matrices such that they can be imported to the CEA file
#Takes as its input the optimal values for the eval_variables function and outputs a list of transition matrices
transition_matrices_chemo2chemo <- function(r_prog_nodrug2death, 
                                            proportion_r_pf2prog_drug,
                                            seven_state_var,
                                            calibrate = T, 
                                            firstline_2line = "Chemo-Chemo"){

  
  #Extract the variables from the chemo_chemo_var function
  if (firstline_2line == "Chemo-Chemo") {
    variables <- chemo_chemo_var(seven_state_var, proportion_r_pf2prog_drug)
  }
  
  #Assign names to the different variables
  r_pf2prog_drug <- variables[1]
  r_pf2prog_nodrug <- variables[2]
  r_prog_drug2prog_nodrug <- variables[3]
  gamma_AE <- variables[4]
  alpha_AE <- variables[5]
  gamma_ILD <- variables[6]
  alpha_ILD <- variables[7]

  
  #Rate of dying from other causes (background mortality)
  r_OC_death <- rep(c(0.000455117, 0.000667712, 0.000939403, 0.001487732,
                      0.002481213, 0.004368355, 0.008008381,0.014505686, 
                      0.024617632)
                    ,each = (12*5))
  
  #Takes to long to calibrate if we use the whole length (we do not need the whole lenght to calibrate)
  if (calibrate) {
    n <- 40
  }else{
    n <- 121
  }
  
  #The rate of getting an AE/ILD and being discontinued.
  Ft_AE <- c()
  Ft_ILD <- c()
  
  Ht_AE <- c()
  Ht_ILD <- c()
  
  r_pf2AE <- c()
  r_pf2ILD <- c()
  for (t in 1:(n + 1)) {
    Ft_AE  <- append(Ft_AE, alpha_AE*(1 - exp(-gamma_AE*t)))
    Ft_ILD <- append(Ft_ILD, alpha_ILD*(1 - exp(-gamma_ILD*t)))
    
    Ht_AE  <- append(Ht_AE, -log(1 - Ft_AE[t]))
    Ht_ILD <- append(Ht_ILD, -log(1 - Ft_ILD[t]))
    
    if (t > 1) {
      r_pf2AE  <- append(r_pf2AE, Ht_AE[t] - Ht_AE[t - 1])
      r_pf2ILD <- append(r_pf2ILD, Ht_ILD[t] - Ht_ILD[t - 1])
    }
  }
  
  
  #This are the rates from the seven state run
  r_pfAE2progAE  <- 0.15362337
  r_progAE2death <- 0.07171597
  
  r_pfILD2progILD    <- r_pfAE2progAE
  r_pf2progILD_death <- r_progAE2death
  
  
  #Stores the transition matrices
  P <- list()
  
  for (t in 1:n) {
    #Find the rate matrix
    min_pf <- r_OC_death[t] + r_pf2prog_nodrug + r_pf2prog_drug + 
      r_pf2ILD[t] + r_pf2AE[t]
    G <- matrix(c(-min_pf, r_pf2AE[t], r_pf2ILD[t], r_pf2prog_drug, r_pf2prog_nodrug, 0, 0, r_OC_death[t],
                        0, -(r_OC_death[t]+r_pfAE2progAE), 0, 0, 0, r_pfAE2progAE, 0, r_OC_death[t],
                        0, 0, -(r_OC_death[t]+r_pfILD2progILD), 0, 0, 0, r_pfILD2progILD, r_OC_death[t],
                        0, 0, 0, -(r_OC_death[t] + r_prog_drug2prog_nodrug), r_prog_drug2prog_nodrug, 0, 0, (r_OC_death[t]),
                        0, 0, 0, 0, -(r_OC_death[t] + r_prog_nodrug2death), 0, 0, (r_OC_death[t] + r_prog_nodrug2death),
                        0, 0, 0, 0, 0, -(r_OC_death[t] + r_progAE2death), 0, (r_OC_death[t] + r_progAE2death),
                        0, 0, 0, 0, 0, 0, -(r_OC_death[t] + r_pf2progILD_death), (r_OC_death[t] + r_pf2progILD_death),
                        0, 0, 0, 0, 0, 0, 0, 0), 
                      ncol = 8, nrow = 8,
                      dimnames = list(c("PF", "PF_AE", "PF_ILD", "Prog_drug", "Prog_nodrug", "Prog_AE", "Prog_ILD", "Death"),
                                      c("PF", "PF_AE", "PF_ILD", "Prog_drug", "Prog_nodrug", "Prog_AE", "Prog_ILD", "Death")),
                      byrow = T
    )
    #Transform the rate matrix to a probability matrix
    A <- expm(G)
    
    P[[t]] <- A
  }
  
  return(P)
}




#This function is for evaluation how good the estimate if r_prog_nodrug2death is. Its being used in the optimization function later
eval_variables <- function(r_prog_nodrug2death, 
                           proportion_r_pf2prog_drug,
                           seven_state_var,
                           drug_sequence = "Chemo-Chemo"){
  
  if (drug_sequence == "Chemo-Chemo") {
    median_os <- 16.8
    A_list <- transition_matrices_chemo2chemo(r_prog_nodrug2death, 
                                              proportion_r_pf2prog_drug,
                                              seven_state_var)
  }
  
  
  #Vectors for storing the estimates of the KM curves
  km_os_model <- c()
  
  #number of cycles
  n = 30
  
  #Start states
  s0 <- c(1,0,0,0,0,0,0,0)
  
  #Add the first values
  survival <- s0[1] + s0[4] + s0[5]
  km_os_model <- append(km_os_model, survival)
  
  #Variable for storing the median overall survival
  model_median_os <- 0
  
  #Do the simulations
  for (t in 1:n) {
    
    #Find transition matrix
    A <- A_list[[t]]
    #Do one transition
    s1 <- s0 %*% A
    
    #Add the new value to the km curve
    survival <- km_os_model[t]*(s1[1]+s1[4]+s1[5])/(s1[1]+s1[4]+s1[5]+s0[1]*A[1,8]+s0[4]*A[4,8]+s0[5]*A[5,8])
    km_os_model <- append(km_os_model, survival)
    
    #Find the median os and stops the for loop once this is found
    if ((km_os_model[t + 1] < 0.5) && (model_median_os == 0)) {
      model_median_os <- t
      break
    }
    
    # Update the current state
    s0 <- s1
  }
  
  #Calculate the likelihood
  v_target <- c(median_os)
  v_target_sd <- c(1)
  v_output <- c(model_median_os)
  
  nllik <- -sum(dnorm(x = v_target, mean = v_output, sd = v_target_sd, 
                      log = TRUE))
  
  return(nllik)
}


eval_variables2 <- function(x){
  return(eval_variables(x,
                        proportion_r_pf2prog_drug,
                        seven_state_var,))
}
#Find the optimal values (Here we could do a wider search for starting values since we might have missed the global optimum.)
fit_out <- optimize(eval_variables2, c(0.01, 0.7))

#Extract the optimal value
opt_var <- fit_out$minimum
r_prog_nodrug2death <- opt_var

opt_var
#0.1107312







tdxd_chemo_var <- function(r_prog_nodrug2death, seven_state_var, proportion_r_pf2prog_drug){
  r_pf2prog_drug <- seven_state_var[1]*seven_state_var[3]*proportion_r_pf2prog_drug
  r_pf2prog_nodrug <- seven_state_var[1]*seven_state_var[3]*(1-proportion_r_pf2prog_drug)

  #Variables for the AEs in Tdxd
  gamma_AE <- seven_state_var[6]
  alpha_AE <- seven_state_var[7]
  
  #Variables for the ILDs in Tdxd
  gamma_ILD <- seven_state_var[8]
  alpha_ILD <- seven_state_var[9]
  
  return(c(r_pf2prog_drug, r_pf2prog_nodrug, r_prog_nodrug2death, gamma_AE, alpha_AE, gamma_ILD, alpha_ILD))
}


chemo_tdxd_var <- function(r_prog_nodrug2death, seven_state_var, proportion_r_pf2prog_drug){
  r_pf2prog_drug   <- seven_state_var[1]*proportion_r_pf2prog_drug
  r_pf2prog_nodrug <- seven_state_var[1]*(1 - proportion_r_pf2prog_drug)
  
  #Variables for the AEs in Chemo
  gamma_AE <- seven_state_var[4]
  alpha_AE <- seven_state_var[5]
  
  #Variables for the ILDs in Chemo
  gamma_ILD <- 0
  alpha_ILD <- 0
  
  return(c(r_pf2prog_drug, r_pf2prog_nodrug, r_prog_nodrug2death, gamma_AE, alpha_AE, gamma_ILD, alpha_ILD))
}


tdxd_sg_var <- function(r_prog_nodrug2death, seven_state_var, proportion_r_pf2prog_drug){
  r_pf2prog_drug   <- seven_state_var[1]*seven_state_var[3]*proportion_r_pf2prog_drug
  r_pf2prog_nodrug <- seven_state_var[1]*seven_state_var[3]*(1 - proportion_r_pf2prog_drug)
  
  #Variables for the AEs in Tdxd
  gamma_AE <- seven_state_var[6]
  alpha_AE <- seven_state_var[7]
  
  #Variables for the ILDs in Tdxd
  gamma_ILD <- seven_state_var[8]
  alpha_ILD <- seven_state_var[9]
  
  return(c(r_pf2prog_drug, r_pf2prog_nodrug, r_prog_nodrug2death, gamma_AE, alpha_AE, gamma_ILD, alpha_ILD))
}






#Function for finding the transition matrices such that they can be imported to the CEA file
#Takes as its input the optimal values for the eval_variables function and outputs a list of transition matrices
transition_matrices <- function(r_prog_drug2prog_nodrug, 
                                proportion_r_pf2prog_drug,
                                seven_state_var,
                                calibrate=T, 
                                firstline_2line="TDxD-Chemo"){
  
  

  if(firstline_2line=="TDxD-Chemo"){
    variables <- tdxd_chemo_var(r_prog_nodrug2death, seven_state_var, proportion_r_pf2prog_drug)
  }
  if(firstline_2line=="Chemo-TDxD"){
    variables <- chemo_tdxd_var(r_prog_nodrug2death, seven_state_var, proportion_r_pf2prog_drug)
  }
  if(firstline_2line=="TDxD-SG"){
    variables <- tdxd_sg_var(r_prog_nodrug2death, seven_state_var, proportion_r_pf2prog_drug)
  }
  
  r_pf2prog_drug <- variables[1]
  r_pf2prog_nodrug <- variables[2]
  r_prog_nodrug2death <- variables[3]
  gamma_AE <- variables[4]
  alpha_AE <- variables[5]
  gamma_ILD <- variables[6]
  alpha_ILD <- variables[7]
  
  
  #Rate of dying from other causes (background mortality)
  r_OC_death <- rep(c(0.000455117,0.000667712,0.000939403,0.001487732,0.002481213, 0.004368355, 0.008008381,0.014505686, 0.024617632)
                    ,each=(12*5))
  
  #Takes to long to calibrate if we use the whole length (we do not need the whole lenght to calibrate)
  if(calibrate){
    n <- 40
  }else{
    n <- 121
  }
  #The rate of getting an AE and being discontinued. (Needs updating)
  Ft_AE <- c()
  Ft_ILD <- c()
  
  Ht_AE <- c()
  Ht_ILD <- c()
  
  r_pf2AE <- c()
  r_pf2ILD <- c()
  for(t in 1:(n+1)){
    Ft_AE <- append(Ft_AE, alpha_AE*(1-exp(-gamma_AE*t)))
    Ft_ILD <- append(Ft_ILD, alpha_ILD*(1-exp(-gamma_ILD*t)))
    
    Ht_AE <- append(Ht_AE, -log(1-Ft_AE[t]))
    Ht_ILD <- append(Ht_ILD, -log(1-Ft_ILD[t]))
    
    if(t > 1){
      r_pf2AE <- append(r_pf2AE, Ht_AE[t]-Ht_AE[t-1])
      r_pf2ILD <- append(r_pf2ILD, Ht_ILD[t]-Ht_ILD[t-1])
    }
  }
  
  #This are the rates from chemo ???
  r_pfAE2progAE <- 0.15362337
  r_progAE2death <- 0.07171597
  
  r_pfILD2progILD <- r_pfAE2progAE
  r_pf2progILD_death <- r_progAE2death
  
  
  P <- list()
  
  for(t in 1:n){
    min_pf <- r_OC_death[t]+r_pf2prog_nodrug+r_pf2prog_drug+r_pf2ILD[t]+r_pf2AE[t]
    G <- matrix(c(-min_pf, r_pf2AE[t], r_pf2ILD[t], r_pf2prog_drug, r_pf2prog_nodrug, 0, 0, r_OC_death[t],
                  0, -(r_OC_death[t]+r_pfAE2progAE), 0, 0, 0, r_pfAE2progAE, 0, r_OC_death[t],
                  0, 0, -(r_OC_death[t]+r_pfILD2progILD), 0, 0, 0, r_pfILD2progILD, r_OC_death[t],
                  0, 0, 0, -(r_OC_death[t] + r_prog_drug2prog_nodrug), r_prog_drug2prog_nodrug, 0, 0, (r_OC_death[t]),
                  0, 0, 0, 0, -(r_OC_death[t] + r_prog_nodrug2death), 0, 0, (r_OC_death[t] + r_prog_nodrug2death),
                  0, 0, 0, 0, 0, -(r_OC_death[t] + r_progAE2death), 0, (r_OC_death[t] + r_progAE2death),
                  0, 0, 0, 0, 0, 0, -(r_OC_death[t] + r_pf2progILD_death), (r_OC_death[t] + r_pf2progILD_death),
                  0, 0, 0, 0, 0, 0, 0, 0), 
                ncol = 8, nrow = 8,
                dimnames = list(c("PF", "PF_AE", "PF_ILD", "Prog_drug", "Prog_nodrug", "Prog_AE", "Prog_ILD", "Death"),
                                c("PF", "PF_AE", "PF_ILD", "Prog_drug", "Prog_nodrug", "Prog_AE", "Prog_ILD", "Death")),
                byrow = T
    )
    
    A <- expm(G)
    
    P[[t]] <- A
  }
  
  return(P)
}



#This function is used to calculate the rate of second line to palliative care.
calc_r_prog_drug2prog_nodrug <- function(r_prog_drug2prog_nodrug, median_pf=3.7){
  #Rate of dying from other causes (background mortality)
  r_OC_death <- rep(c(0.000455117,0.000667712,0.000939403,0.001487732,0.002481213, 0.004368355, 0.008008381,0.014505686, 0.024617632)
                    ,each=(12*5))
  P <- list()
  n<-40
  for(t in 1:n){
    G <- matrix(c(-(r_prog_drug2prog_nodrug+r_OC_death[t]), r_prog_drug2prog_nodrug, r_OC_death[t],
                  0, -0.04853368, 0.04853368,
                  0, 0, 0), 
                ncol = 3, nrow = 3,
                dimnames = list(c("PF", "P", "D"),
                                c("PF", "P", "D")),
                byrow = T
    )
    
    A <- expm(G)
    P[[t]] <- A
  }
  
  s0 <- c(1,0,0)
  
  model_pf <- 0
  km_pf_model <- c()
  km_pf_model <- append(km_pf_model, s0[1])
  
  for(t in 1:n){
    
    A <- P[[t]]
    
    s1 <- s0 %*% A
    
    km_pf_model <- append(km_pf_model, km_pf_model[t]*s1[1]/(s1[1]+s0[1]*A[1,2]+s0[1]*A[1,3]))
    
    if((km_pf_model[t+1]<0.5) && (model_pf == 0)){
      model_pf <- t
    }
    s0 <- s1
  }
  
  #Calculate the likelihood
  v_target <- c(median_pf)
  v_target_sd <- c(1)
  v_output <- c(model_pf)
  nllik <- -sum(dnorm(x = v_target, mean = v_output, sd = v_target_sd, log = TRUE))
  return(nllik)
}

eval_chemo_tdxd <- function(r_prog_drug2prog_nodrug, median_pf=9.9){
  calc_r_prog_drug2prog_nodrug(r_prog_drug2prog_nodrug, median_pf)
}
eval_tdxd_chemo <- function(r_prog_drug2prog_nodrug, median_pf=4){
  calc_r_prog_drug2prog_nodrug(r_prog_drug2prog_nodrug, median_pf)
}
eval_tdxd_sg <- function(r_prog_drug2prog_nodrug, median_pf=7){
  calc_r_prog_drug2prog_nodrug(r_prog_drug2prog_nodrug, median_pf)
}

fit_out_tdxd_chemo <- optimize(eval_tdxd_chemo, c(0.01, 0.5))
opt_var_tdxd_chemo <- fit_out_tdxd_chemo$minimum

fit_out_chemo_tdxd <- optimize(eval_chemo_tdxd, c(0.01, 0.5))
opt_var_chemo_tdxd <- fit_out_chemo_tdxd$minimum

fit_out_tdxd_sg <- optimize(eval_tdxd_sg, c(0.01, 0.5))
opt_var_tdxd_sg <- fit_out_tdxd_sg$minimum





#Find the optimal transition matrices (These lists will be send to the CEA_8state.R file)
tm_chemo_chemo <- transition_matrices_chemo2chemo(r_prog_nodrug2death, 
                                                  proportion_r_pf2prog_drug,
                                                  seven_state_var,
                                                  calibrate = F)
tm_tdxd_chemo  <- transition_matrices(opt_var_tdxd_chemo, 
                                      proportion_r_pf2prog_drug,
                                      seven_state_var,
                                      calibrate = F)
tm_chemo_tdxd  <- transition_matrices(opt_var_chemo_tdxd, 
                                      proportion_r_pf2prog_drug,
                                      seven_state_var,
                                      calibrate = F, 
                                      firstline_2line = "Chemo-TDxD")
tm_tdxd_sg <- transition_matrices(opt_var_tdxd_sg, 
                                  proportion_r_pf2prog_drug,
                                  seven_state_var,
                                  calibrate = F, 
                                  firstline_2line = "TDxD-SG")















#DO NOT NEED TO LOOK AT THE CODE UNDER HERE!!!



#Plot the correct KM-curves vs the cureves from our model. This part of the file is used for validation

#The transitions from rate to probability is found in https://journals.sagepub.com/doi/pdf/10.1177/0272989X05282637 Figure 3
plot_function <- function(tm_tdxd_chemo){
  
  
  A_list <- tm_tdxd_chemo

  
  #The extracted Kaplan-Meier values
  km_os_tdxd <- c(1, 0.992, 0.985, 0.967, 0.958, 0.939, 0.928, 0.899, 0.87, 0.857, 0.829, 0.815, 0.793, 0.774)#, 0.742, 0.731, 0.688, 0.661, 0.623, 0.586, 0.566, 0.539, 0.507, 0.500, 0.484)
  km_pf_tdxd <- c(1, 0.988, 0.896, 0.82, 0.809, 0.765, 0.681, 0.633, 0.592, 0.554, 0.491, 0.457, 0.422, 0.386)#, 0.371, 0.363, 0.339, 0.31, 0.295, 0.274, 0.259) 
  
  #Vectors for storing the estimates of the KM curves
  km_os_tdxd_model <- c()
  km_pf_tdxd_model <- c() 
  
  
  n = max(length(km_os_tdxd), length(km_pf_tdxd))
  
  #Start states
  s0_tdxd <- c(1,0,0,0,0,0,0,0)
  
  survival_tdxd <- s0_tdxd[1]+s0_tdxd[4]+s0_tdxd[5]
  km_os_tdxd_model <- append(km_os_tdxd_model, survival_tdxd)
  km_pf_tdxd_model <- append(km_pf_tdxd_model, s0_tdxd[1])
  
  ae_test_tdxd <- c(0)
  ae_test_ild <- c(0)
  
  
  for(t in 1:n){
    A_tdxd <- A_list[[t]]
    s1_tdxd <- s0_tdxd %*% A_tdxd
    
    
    survival_tdxd <- km_os_tdxd_model[t]*(s1_tdxd[1]+s1_tdxd[4]+s1_tdxd[5])/(s1_tdxd[1]+s1_tdxd[4]+s1_tdxd[5]+s0_tdxd[1]*A_tdxd[1,8]+s0_tdxd[4]*A_tdxd[4,8]+s0_tdxd[5]*A_tdxd[5,8])
    km_os_tdxd_model <- append(km_os_tdxd_model, survival_tdxd)
    km_pf_tdxd_model <- append(km_pf_tdxd_model, km_pf_tdxd_model[t]*(s1_tdxd[1])/(s1_tdxd[1]+s0_tdxd[1]*A_tdxd[1,4]+s0_tdxd[1]*A_tdxd[1,5]+s0_tdxd[1]*A_tdxd[1,8]))
    
    ae_test_tdxd <- c(ae_test_tdxd, ae_test_tdxd[t]+s0_tdxd[1]*A_tdxd[1,2])
    ae_test_ild <- c(ae_test_ild, ae_test_ild[t]+s0_tdxd[1]*A_tdxd[1,3])
    
    s0_tdxd <- s1_tdxd
    
  }
  
  
  
  m = min(length(km_os_tdxd), length(km_pf_tdxd))
  idx <- 1:m
  df <- data.frame(idx, km_os_tdxd_model[1:m])
  
  
  
  plot2 <- ggplot(df, aes(x=idx)) + 
    geom_line(aes(y = km_os_tdxd_model[1:m]), color = "red") + 
    geom_line(aes(y = km_os_tdxd[1:m]), color="blue", linetype="twodash") +
    ylim(0, 1) +
    ylab("Probability")+
    xlab("Month")+
    ggtitle("OS T-DxD")+
    scale_x_continuous(breaks = round(seq(0, 18, by = 3),1)) + 
    theme(axis.text = element_text(size = 12),
          axis.title = element_text(size = 14))
  
  plot4 <- ggplot(df, aes(x = idx)) + 
    geom_line(aes(y = km_pf_tdxd_model[1:m]), color = "red") + 
    geom_line(aes(y = km_pf_tdxd[1:m]), color = "blue", linetype = "twodash") +
    ylim(0, 1) +
    ylab("Probability") +
    xlab("Month") +
    ggtitle("PFS T-DxD") +
    scale_x_continuous(breaks = round(seq(0, 18, by = 3),1)) + 
    theme(axis.text = element_text(size = 12),
          axis.title = element_text(size = 14))
  
  print(ae_test_tdxd)
  print(ae_test_ild)
  plot(0:(length(ae_test_tdxd)-1), ae_test_tdxd, type = "l")
  plot(0:(length(ae_test_tdxd)-1), ae_test_ild, type = "l")
  
  # ggarrange(plot2, plot4,
  #           ncol = 1, nrow = 2, common.legend = TRUE,legend="bottom") 
  
}


plot_function(tm_tdxd_chemo)























