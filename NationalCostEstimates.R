#This is the script for generating the cost estimates using our model for all payers in the US
#This script is written by James C Dickerson, MD, MS. Inputs can be found in the manuscript appendix
library(writexl)

options(scipen=999)

#Note! This does not address market diffusion of new technology. It is a very rough estimate 

#Number of patients with mBC incidence 
n <- 54394 

HER2low_low_est <- 0.457
HER2low_base_est <- 0.557
HER2low_high_est <- 0.657
HER2_ultra_low_est <- 0.65

#Cost of Care (using our model). The low and high estimates represent two SDs from the mean of the PSA simulation 
#For now, we have just used the mean values for the cost of each strategy given we are varying the inputs on HER2 low numbers and lines of therapy estimates
#This is because there is already tremendous imprecision in these estimates 
#c_chemo_low <- 
c_chemo_base <- 180252
#c_chemo_high <- 

#c_tdxd_first_low <- 
c_tdxd_first_base <- 279810
#c_tdxd_first_high <- 

#c_tdxd_second_low <- 
c_tdxd_second_base <- 243834
#c_tdxd_second_high <- 

#c_ADCs_low <- 
c_ADCs_base <- 310036
#c_ADCs_high <- 
  
#Number receiving third line of therapy
third_low <- 0.35
third_base <- 0.41
third_high <- 0.55


#Create the DF that is the table in the paper 
low_est <- n*HER2low_low_est*third_low
med_est <- n*HER2low_base_est*third_base
high_est <- n*HER2low_high_est*third_high
ultra_est <- n*HER2_ultra_low_est*third_base

t_chemo <- c_tdxd_first_base
chemo_t <- c_tdxd_second_base
adcs <- c_ADCs_base
chemo_chemo <- c_chemo_base


ests <- c(low_est,
          med_est,
          high_est,
          ultra_est)

strats <- c(chemo_chemo,
            chemo_t,
            t_chemo,
            adcs)



final_df <- matrix(data = NA, nrow = 4, ncol = 4)
rownames(final_df) <- c("chemo_chemo",
                        "chemo_t",
                        "t_chemo",
                        "adcs")
colnames(final_df) <- c("low_est",
                        "med_est",
                        "high_est",
                        "ultra_est")

for (i in 1:4) {
  for (j in 1:4) {
    final_df[i, j] <- strats[i] * ests[j]
  }
}

final_df <- as.data.frame(final_df)

output_file <- "data_/final_df.xlsx"

# Save final_df as an Excel file
write_xlsx(final_df, path = output_file)

