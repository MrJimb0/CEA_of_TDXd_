#This code was created by Marcus Moen, MS and edited by James Dickerson, MD, MS and Wesley Suen
#Oversight from Fernando Alarid-Escudero, Ph.D.
options(scipen=999)

library(ggplot2)
library(ggpubr)
library(dampack)
library(grid)
library(gridExtra)
library(dplyr)
library(patchwork)  

source("finding_transition_probabilities_7state.R")
source("finding_transition_probabilities_8state.R")

#Main Graphs for the Manuscript from DAMPACK

df_psa_lng <- read.csv("data/base_case_output.csv")
# Transform long to wide the dataset with group as columns
df_psa <- reshape(df_psa_lng[, -1], idvar = "sim", timevar = "group", direction = "wide")
# v_names_str <- df_psa_lng[1:4, "group"]
v_names_str <- c("Chemo → Chemo",
                 "Chemo → TDXd", 
                 "TDXd → Chemo", 
                 "TDXd → SG")
# expression(paste("Chemo", ""%right%"", "Chemo", sep = ""))
df_costs <- df_psa[, c(2, 5, 8, 11)]
df_qalys <- df_psa[, c(3, 6, 9, 12)]
df_lys   <- df_psa[, c(4, 7, 10, 13)]

## Visualize PSA results for CEA ----
### Create PSA object ----
l_psa <- dampack::make_psa_obj(cost          = df_costs, 
                               effectiveness = df_qalys, 
                               strategies    = v_names_str)
l_psa$strategies <- v_names_str
colnames(l_psa$effectiveness) <- v_names_str
colnames(l_psa$cost) <- v_names_str

#* Vector with willingness-to-pay (WTP) thresholds.
v_wtp <- seq(0, 700000, by = 5000)

### Cost-Effectiveness Scatter plot ----
txtsize <- 16
gg_scattter <- plot(l_psa, txtsize = txtsize) +
  ggthemes::scale_color_colorblind() +
  ggthemes::scale_fill_colorblind() +
  scale_y_continuous("Cost (Thousands of $)", 
                     n.breaks = 6,
                     labels = function(x) x/1000) +
  xlab("Effectiveness (QALYs)") +
  guides(col = guide_legend(nrow = 2)) +
  theme(legend.position = "bottom")
gg_scattter
ggsave(gg_scattter, filename = "figs/cea_scatter.png", width = 8, height = 6, dpi = 300)

### Incremental cost-effectiveness ratios (ICERs) with probabilistic output ----
#* Compute expected costs and effects for each strategy from the PSA
df_out_ce_psa <- summary(l_psa)
df_cea_psa <- dampack::calculate_icers(cost       = df_out_ce_psa$meanCost, 
                                       effect     = df_out_ce_psa$meanEffect,
                                       strategies = df_out_ce_psa$Strategy)
df_cea_psa

### Plot cost-effectiveness frontier with probabilistic output ----
plot(df_cea_psa, label = "all", txtsize = txtsize) +
  theme(legend.position = c(0.8, 0.2))

### Cost-effectiveness acceptability curves (CEACs); CEAF and other DAMPACK code later down with graphs that were not in the final draft
ceac_obj <- dampack::ceac(wtp = v_wtp, psa = l_psa)
#* Regions of highest probability of cost-effectiveness for each strategy
summary(ceac_obj)
#* CEAC & CEAF plot
gg_ceac <- plot(ceac_obj, 
                txtsize = txtsize, xlim = c(0, NA), n_x_ticks = 16) +
  xlab("Cost-Effectiveness Threshold (Thousands $/QALY)") +
  ylab("Probability of being cost-effective") +
  ggthemes::scale_color_colorblind() +
  ggthemes::scale_fill_colorblind() +
  theme(legend.position = c(0.2, 0.48))
gg_ceac
ggsave(gg_ceac, filename = "figs/ceac.png", width = 8, height = 6, dpi = 300)



#Plot efficiency frontiers
tdxd_icers <- calculate_icers(cost = df_res$DiscountedCost, 
                              effect = df_res$DiscountedQALY, 
                              strategies = df_res$Strategy)

plot(tdxd_icers) + 
  scale_x_continuous(n.breaks = 10) + 
  scale_y_continuous(n.breaks = 10) + 
  theme(legend.position = c(0.3, 0.7))



#End of DAMPACK 
#Start of the code for other graphs

#This graph shows the evolution of the patients in the different treatment sequences (appendix graph)
source("CEA_8state.R")
#Plot evolution of patients
plot_evolution <- function(df_plot, title){
  p <- ggplot(df_plot, aes(x = cycle, y = value, group = state, color = state)) +
    geom_line() + 
    ggtitle(title) +
    xlab("Months") + ylab("Probability") + 
    theme_minimal() +
    theme(
      plot.title = element_text(color="black", size=14, face="bold.italic", hjust = 0.5), # Center title
      axis.title.x = element_text(color="black", size=14, face="bold"),
      axis.title.y = element_text(color="black", size=14, face="bold"),
      axis.line = element_line(color = "black"),  # Set color of axis lines
      legend.text = element_text(family = "Arial"),  # Set legend font family
      legend.title = element_text(family = "Arial", face = "bold")  # Set legend title font family
    ) +
    scale_x_continuous(breaks = seq(0, max(df_plot$cycle), by = 6), limits = c(0,90)) +  # Adjust x-axis breaks
    scale_color_manual(values = c("Dead" = "red", "Progressed_drug" = "gold", "Progressed_nodrug" = "green", "ProgressedAE"="cyan", "ProgressedILD"= "blue", "ProgressionFree"= "purple", "ProgressionFreeAE"= "magenta", "ProgressionFreeILD"="pink"),
                       labels = c("Dead", "Progressed & candidate for other drug", "Progressed & no further drug", "Progressed after AE", "Progressed after ILD", "Progression-Free", "Progression-Free after AE", "Progression-Free after ILD"),
                       name = "State")  # Update legend title
  
  return(p)
}

# Clear plotting region
dev.off()

# Example calls to the function
plot5 <- plot_evolution(df_plot_tdxd_chemo, "Evolution of T-DXd → Chemo Cohort")
plot6 <- plot_evolution(df_plot_chemo_chemo, "Evolution of Chemo → Chemo Cohort")
plot7 <- plot_evolution(df_plot_chemo_tdxd, "Evolution of Chemo → T-DXd Cohort")
plot8 <- plot_evolution(df_plot_tdxd_sg, "Evolution of T-DXd → SG Cohort")

# Combine plots vertically
# Just main comparison of chemo-chemo and tdxd-chemo is in the appendix of the manuscript
combined_plots <- arrangeGrob(plot5, plot6, ncol = 1)
#combined_plots <- arrangeGrob(plot7, plot8, ncol = 1)
print(combined_plots)
grid.draw(combined_plots)
ggsave(combined_plots, filename = "figs/evolution_plots.png", width = 8, height = 6, dpi = 300)



#NEXT GRAPH
#NEXT GRAPH
#NEXT GRAPH
#NEXT GRAPH
#NEXT GRAPH

#This graph shows the Kaplan-Meier curves for the different treatment sequences

#Store the transition matrices from the "finding_transition_probabilities_8state.R" file
#Store the transition matrices from the "finding_transition_probabilities_8state.R" file
A_tdxd_chemo <- tm_tdxd_chemo
A_chemo_chemo <- tm_chemo_chemo
A_chemo_tdxd <- tm_chemo_tdxd
A_tdxd_sg <- tm_tdxd_sg

input_var <- list(A_tdxd_chemo,
                  A_chemo_chemo,
                  A_chemo_tdxd,
                  A_tdxd_sg)

plot_function <- function(input_var){
  A_tdxd_chemo <- input_var[[1]]
  A_chemo_chemo <- input_var[[2]]
  A_chemo_tdxd <- input_var[[3]]
  A_tdxd_sg <- input_var[[4]]
  
  
  #The extracted Kaplan-Meier values
  km_os_chemo_chemo <- c(1, 0.986, 0.981, 0.957, 0.94, 0.927, 0.884, 0.841, 0.793, 0.738, 0.712, 0.683, 0.669, 0.634, 0.610, 0.566, 0.520, 0.484, 0.459, 0.460, 0.457, 0.419, 0.375)
  km_pf_chemo_chemo <- c(1, 0.983, 0.756, 0.62, 0.603, 0.501, 0.439, 0.398, 0.34, 0.277, 0.268, 0.226, 0.217, 0.185) 
  km_os_tdxd_chemo <- c(1, 0.992, 0.985, 0.967, 0.958, 0.939, 0.928, 0.899, 0.87, 0.857, 0.829, 0.815, 0.793, 0.774, 0.742, 0.731, 0.688, 0.661, 0.623, 0.586, 0.566, 0.539, 0.507, 0.500, 0.484)
  km_pf_tdxd_chemo <- c(1, 0.988, 0.896, 0.82, 0.809, 0.765, 0.681, 0.633, 0.592, 0.554, 0.491, 0.457, 0.422, 0.386, 0.371, 0.363, 0.339, 0.31, 0.295, 0.274, 0.259) 
  
  #Vectors for storing the estimates of the KM curves
  km_os_chemo_chemo_model <- c()
  km_pf_chemo_chemo_model <- c()
  
  km_os_tdxd_chemo_model <- c()
  km_pf_tdxd_chemo_model <- c() 
  
  km_os_chemo_tdxd_model <- c()
  km_pf_chemo_tdxd_model <- c() 
  
  km_os_tdxd_sg_model <- c()
  km_pf_tdxd_sg_model <- c() 
  
  
  n <- 24
  #n = max(length(km_os_chemo_chemo), length(km_pf_chemo_chemo), length(km_os_tdxd_chemo), length(km_pf_tdxd_chemo))
  
  #The initial state
  s0_chemo_chemo <- matrix(c(1,0,0,0,0,0,0,0), ncol = 8)
  s0_tdxd_chemo <- matrix(c(1,0,0,0,0,0,0,0), ncol = 8)
  s0_chemo_tdxd <- matrix(c(1,0,0,0,0,0,0,0), ncol = 8)
  s0_tdxd_sg <- matrix(c(1,0,0,0,0,0,0,0), ncol = 8)
  
  survival_chemo_chemo <- s0_chemo_chemo[1]+s0_chemo_chemo[4]+s0_chemo_chemo[5]
  km_os_chemo_chemo_model <- append(km_os_chemo_chemo_model, survival_chemo_chemo)
  km_pf_chemo_chemo_model <- append(km_pf_chemo_chemo_model, s0_chemo_chemo[1])
  
  survival_tdxd_chemo <- s0_tdxd_chemo[1]+s0_tdxd_chemo[4]+s0_tdxd_chemo[5]
  km_os_tdxd_chemo_model <- append(km_os_tdxd_chemo_model, survival_tdxd_chemo)
  km_pf_tdxd_chemo_model <- append(km_pf_tdxd_chemo_model, s0_tdxd_chemo[1])
  
  survival_chemo_tdxd <- s0_chemo_tdxd[1]+s0_chemo_tdxd[4]+s0_chemo_tdxd[5]
  km_os_chemo_tdxd_model <- append(km_os_chemo_tdxd_model, survival_chemo_tdxd)
  km_pf_chemo_tdxd_model <- append(km_pf_chemo_tdxd_model, s0_chemo_tdxd[1])
  
  survival_tdxd_sg <- s0_tdxd_sg[1]+s0_tdxd_sg[4]+s0_tdxd_sg[5]
  km_os_tdxd_sg_model <- append(km_os_tdxd_sg_model, survival_tdxd_sg)
  km_pf_tdxd_sg_model <- append(km_pf_tdxd_sg_model, s0_tdxd_sg[1])
  
  #ae_test_chemo <- c(0)
  #ae_test_tdxd <- c(0)
  #ae_test_ild <- c(0)
  
  
  for(t in 1:n){ # t <- 1
    
    s1_chemo_chemo <- s0_chemo_chemo %*% A_chemo_chemo[[t]]
    s1_tdxd_chemo <- s0_tdxd_chemo %*% A_tdxd_chemo[[t]]
    s1_chemo_tdxd <- s0_chemo_tdxd %*% A_chemo_tdxd[[t]]
    s1_tdxd_sg <- s0_tdxd_sg %*% A_tdxd_sg[[t]]
    
    survival_chemo_chemo <- km_os_chemo_chemo_model[t] * (s1_chemo_chemo[1] + s1_chemo_chemo[4] + s1_chemo_chemo[5]) / 
      (s1_chemo_chemo[1] + s1_chemo_chemo[4] + s1_chemo_chemo[5] + s0_chemo_chemo[1] * A_chemo_chemo[[t]][1, 8] + 
         s0_chemo_chemo[4] * A_chemo_chemo[[t]][4, 8] + s0_chemo_chemo[5] * A_chemo_chemo[[t]][5, 8])
    km_os_chemo_chemo_model <- append(km_os_chemo_chemo_model, survival_chemo_chemo)
    km_pf_chemo_chemo_model <- append(km_pf_chemo_chemo_model, km_pf_chemo_chemo_model[t] * s1_chemo_chemo[1] / 
                                        (s1_chemo_chemo[1] + s0_chemo_chemo[1] * A_chemo_chemo[[t]][1, 4] + 
                                           s0_chemo_chemo[1] * A_chemo_chemo[[t]][1, 5] + s0_chemo_chemo[1] * A_chemo_chemo[[t]][1, 8]))
    
    survival_tdxd_chemo <- km_os_tdxd_chemo_model[t] * (s1_tdxd_chemo[1] + s1_tdxd_chemo[4] + s1_tdxd_chemo[5]) / 
      (s1_tdxd_chemo[1] + s1_tdxd_chemo[4] + s1_tdxd_chemo[5] + s0_tdxd_chemo[1] * A_tdxd_chemo[[t]][1, 8] + 
         s0_tdxd_chemo[4] * A_tdxd_chemo[[t]][4, 8] + s0_tdxd_chemo[5] * A_tdxd_chemo[[t]][5, 8])
    km_os_tdxd_chemo_model <- append(km_os_tdxd_chemo_model, survival_tdxd_chemo)
    km_pf_tdxd_chemo_model <- append(km_pf_tdxd_chemo_model, km_pf_tdxd_chemo_model[t] * s1_tdxd_chemo[1] / 
                                       (s1_tdxd_chemo[1] + s0_tdxd_chemo[1] * A_tdxd_chemo[[t]][1, 4] + 
                                          s0_tdxd_chemo[1] * A_tdxd_chemo[[t]][1, 5] + s0_tdxd_chemo[1] * A_tdxd_chemo[[t]][1, 8]))
    
    survival_chemo_tdxd <- km_os_chemo_tdxd_model[t] * (s1_chemo_tdxd[1] + s1_chemo_tdxd[4] + s1_chemo_tdxd[5]) / 
      (s1_chemo_tdxd[1] + s1_chemo_tdxd[4] + s1_chemo_tdxd[5] + s0_chemo_tdxd[1] * A_chemo_tdxd[[t]][1, 8] + 
         s0_chemo_tdxd[4] * A_chemo_tdxd[[t]][4, 8] + s0_chemo_tdxd[5] * A_chemo_tdxd[[t]][5, 8])
    km_os_chemo_tdxd_model <- append(km_os_chemo_tdxd_model, survival_chemo_tdxd)
    km_pf_chemo_tdxd_model <- append(km_pf_chemo_tdxd_model, km_pf_chemo_tdxd_model[t] * s1_chemo_tdxd[1] / 
                                       (s1_chemo_tdxd[1] + s0_chemo_tdxd[1] * A_chemo_tdxd[[t]][1, 4] + 
                                          s0_chemo_tdxd[1] * A_chemo_tdxd[[t]][1, 5] + s0_chemo_tdxd[1] * A_chemo_tdxd[[t]][1, 8]))
    
    survival_tdxd_sg <- km_os_tdxd_sg_model[t] * (s1_tdxd_sg[1] + s1_tdxd_sg[4] + s1_tdxd_sg[5]) / 
      (s1_tdxd_sg[1] + s1_tdxd_sg[4] + s1_tdxd_sg[5] + s0_tdxd_sg[1] * A_tdxd_sg[[t]][1, 8] + 
         s0_tdxd_sg[4] * A_tdxd_sg[[t]][4, 8] + s0_tdxd_sg[5] * A_tdxd_sg[[t]][5, 8])
    km_os_tdxd_sg_model <- append(km_os_tdxd_sg_model, survival_tdxd_sg)
    km_pf_tdxd_sg_model <- append(km_pf_tdxd_sg_model, km_pf_tdxd_sg_model[t] * s1_tdxd_sg[1] / 
                                    (s1_tdxd_sg[1] + s0_tdxd_sg[1] * A_tdxd_sg[[t]][1, 4] + 
                                       s0_tdxd_sg[1] * A_tdxd_sg[[t]][1, 5] + s0_tdxd_sg[1] * A_tdxd_sg[[t]][1, 8]))
    
    s0_chemo_chemo <- s1_chemo_chemo
    s0_tdxd_chemo <- s1_tdxd_chemo
    s0_chemo_tdxd <- s1_chemo_tdxd
    s0_tdxd_sg <- s1_tdxd_sg
    
  }
  
  m = 24
  #m = min(length(km_os_chemo_chemo), length(km_pf_chemo_chemo), length(km_os_tdxd_chemo), length(km_pf_tdxd_chemo))
  idx <- 1:m
  
  # Create the data frame for plotting
  df <- data.frame(
    idx = idx,
    km_os_chemo_chemo_model = km_os_chemo_chemo_model[1:m],
    km_os_tdxd_chemo_model = km_os_tdxd_chemo_model[1:m],
    km_os_chemo_chemo = km_os_chemo_chemo[1:m],
    km_os_tdxd_chemo = km_os_tdxd_chemo[1:m],
    km_os_chemo_tdxd_model = km_os_chemo_tdxd_model[1:m],
    km_os_tdxd_sg_model = km_os_tdxd_sg_model[1:m]
  )
  
  # Convert the y-values to percentages
  df$km_os_chemo_chemo_model <- df$km_os_chemo_chemo_model * 100
  df$km_os_tdxd_chemo_model <- df$km_os_tdxd_chemo_model * 100
  df$km_os_chemo_chemo <- df$km_os_chemo_chemo * 100
  df$km_os_tdxd_chemo <- df$km_os_tdxd_chemo * 100
  df$km_os_chemo_tdxd_model <- df$km_os_chemo_tdxd_model * 100
  df$km_os_tdxd_sg_model <- df$km_os_tdxd_sg_model * 100  
  
  df_plot1 <- rbind(data.frame(Source = "Modeled",
                               Strategy = "Chemo → Chemo", 
                               Time = 0:24,
                               Survival = km_os_chemo_chemo_model[1:25]),
                    data.frame(Source = "Observed (Kaplan-Meier)",
                               Strategy = "Chemo → Chemo", 
                               Time = 0:22,
                               Survival = km_os_chemo_chemo[1:23]),
                    data.frame(Source = "Modeled",
                               Strategy = "Chemo → T-DXd", 
                               Time = 0:24,
                               Survival = km_os_chemo_tdxd_model[1:25]),
                    data.frame(Source = "Modeled",
                               Strategy = "T-DXd → Chemo", 
                               Time = 0:24,
                               Survival = km_os_tdxd_chemo_model[1:25]),                    
                    data.frame(Source = "Observed (Kaplan-Meier)", 
                               Strategy = "T-DXd → Chemo", 
                               Time = 0:23,
                               Survival = km_os_tdxd_chemo[1:24]),
                    data.frame(Source = "Modeled",
                               Strategy = "T-DXd → SG", 
                               Time = 0:24,
                               Survival = km_os_tdxd_sg_model[1:25]))
  
  plot1_alt <- ggplot(df_plot1, aes(x = Time, y = Survival, color = Strategy, linetype = Source)) +
    geom_line(size = 1.2) +
    ggthemes::scale_color_colorblind() +
    scale_y_continuous(breaks = seq(0, 1, by = 0.1), limits = c(0, 1), expand = c(0, 0)) +
    ylab("Overall Survival Probability") +
    xlab("Months") +
    ggtitle("Modeled Overall Survival for the Treatment Sequences") +
    scale_x_continuous(breaks = round(seq(0, 24, by = 3), 1), limits = c(0, 24), expand = c(0, 0)) +
    theme_classic(base_size = 16) + 
    theme(
      text = element_text(family = "Arial"),
      axis.text = element_text(size = 12, face = "bold"),
      axis.title = element_text(size = 12, face = "bold"),
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      panel.grid.major.y = element_blank(),
      panel.grid.minor.y = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.x = element_blank(),
      legend.position = "right",  # Place legend at the bottom of the plot
      legend.text = element_text(family = "Arial", size = 12),  # Set legend text font and size
      legend.title = element_text(family = "Arial", size = 14, face = "bold")  # Set legend title font and size
    )
  plot1_alt

  ggsave(plot1_alt, filename = "figs/overall_survival_plot.png", width = 8, height = 6, dpi = 300)
  
}

opt_var <- list(A_tdxd_chemo,
                A_chemo_chemo,
                A_chemo_tdxd,
                A_tdxd_sg)

plot_function(opt_var)

#NEXT GRAPH
#NEXT GRAPH
#NEXT GRAPH
#NEXT GRAPH

#One Way Sensitivity Analysis of Cost (appendix)
#9/1/24 updated the QALY and cost values to correct ones
                     

one_way_sensitivity_tdxd_price <- function(df, dr_v){
  
  cost <- seq(1000, 15000, 500)
  
  n_sims = length(cost)
  
  icer <- c()
  
  tdxd_chemo_qaly <- c(qaly_pf = 0.65, 
                       qaly_p_drug = 0.54, 
                       qaly_p_nodrug = 0.40, 
                       qaly_pfAE = 0.547, 
                       qaly_pAE = 0.54, 
                       qaly_pfILD = 0.547, 
                       qaly_pILD = 0.54, 
                       decrement_qaly_ae = 0.05604)
  chemo_tdxd_qaly <- c(qaly_pf = 0.65, 
                       qaly_p_drug = 0.54, 
                       qaly_p_nodrug = 0.40, 
                       qaly_pfAE = 0.547, 
                       qaly_pAE = 0.54, 
                       qaly_pfILD = 0.547, 
                       qaly_pILD = 0.54, 
                       decrement_qaly_ae = 0.0714327)
  chemo_chemo_qaly <- c(qaly_pf = 0.65, 
                        qaly_p_drug = 0.54, 
                        qaly_p_nodrug = 0.40, 
                        qaly_pfAE = 0.547, 
                        qaly_pAE = 0.54, 
                        qaly_pfILD = 0.547, 
                        qaly_pILD = 0.54, 
                        decrement_qaly_ae = 0.0714327)
  
  df_tdxd_chemo = df[[1]]
  df_chemo_tdxd = df[[2]]
  
  # Initialize empty data frames to store results
  df_results <- data.frame()
  
  for(i in cost){
    
    tdxd_chemo_cost <- c(cost_pf = i, 
                         cost_p_drug = 7203.56, 
                         cost_p_nodrug = 10749.48, 
                         cost_pfAE = 5093.37, 
                         cost_pAE = 10749.48, 
                         cost_pfILD = 5093.37, 
                         cost_pILD = 10749.48,
                         additional_cost = 7461.73)
    
    chemo_tdxd_cost <- c(cost_pf = 7203.56, 
                         cost_p_drug = i, 
                         cost_p_nodrug = 10749.48, 
                         cost_pfAE = 5093.37, 
                         cost_pAE = 10749.48, 
                         cost_pfILD = 5093.37, 
                         cost_pILD = 10749.48,
                         additional_cost = 10349.86)
  
      chemo_chemo_cost <- c(cost_pf = 7203.56, 
                          cost_p_drug = 5093.37, 
                          cost_p_nodrug = 10749.48, 
                          cost_pfAE = 5093.37, 
                          cost_pAE = 10749.48, 
                          cost_pfILD = 5093.37, 
                          cost_pILD = 10749.48,
                          additional_cost = 10349.86)
    
    
      # Calculate summary data for each comparison
      res_vec_tdxd_chemo <- calc_summary_data(df = df_tdxd_chemo, 
                                              cost_data = tdxd_chemo_cost, 
                                              qaly_data = tdxd_chemo_qaly, 
                                              dr_v = dr_v)
      
      res_vec_chemo_tdxd <- calc_summary_data(df = df_chemo_tdxd, 
                                              cost_data = chemo_tdxd_cost, 
                                              qaly_data = chemo_tdxd_qaly, 
                                              dr_v = dr_v)
      
      res_vec_chemo_chemo <- calc_summary_data(df = df_chemo_tdxd, 
                                               cost_data = chemo_chemo_cost, 
                                               qaly_data = chemo_chemo_qaly, 
                                               dr_v = dr_v)
      
      # Calculate ICER for each scenario and store in the results vector
      icer_tdxd_chemo <- (res_vec_tdxd_chemo[2] - res_vec_chemo_tdxd[2]) / (res_vec_tdxd_chemo[4] - res_vec_chemo_tdxd[4])
      icer_tdxd_chemo_chemo <- (res_vec_tdxd_chemo[2] - res_vec_chemo_chemo[2]) / (res_vec_tdxd_chemo[4] - res_vec_chemo_chemo[4])
      
      icer <- c(icer, icer_tdxd_chemo, icer_tdxd_chemo_chemo)
      
      # Prepare data frame to store results
      df_results <- rbind(df_results,
                          data.frame(Willingness2Pay = c(icer_tdxd_chemo, icer_tdxd_chemo_chemo),
                                     Cost = rep(i, 2),
                                     Comparison = c('T-DXd → chemo vs chemo → T-DXd', 'T-DXd → chemo vs. chemo → chemo')))
      
  }
  
  m2 <- (df_results$Cost[3]-df_results$Cost[1])/(df_results$Willingness2Pay[3]-df_results$Willingness2Pay[1])
  b2 <- df_results$Cost[3]-m2*df_results$Willingness2Pay[3]
  cost2 <- 150000*m2+b2
  
  m1 <- (df_results$Cost[4]-df_results$Cost[2])/(df_results$Willingness2Pay[4]-df_results$Willingness2Pay[2])
  b1 <- df_results$Cost[4]-m1*df_results$Willingness2Pay[4]
  cost1 <- 150000*m1+b1
  
  p_reduction1 <- 1-cost1/14300.690
  p_reduction2 <- 1-cost2/14300.690
  
  # Plotting the results
  dot_df <- data.frame(
    Willingness2Pay = 150000,
    Cost = c(cost1, cost2),  # adjust these values to match the y-values of the lines
    Comparison = c("Comparison 1", "Comparison 2"))
  
  ggplot(data = df_results, aes(x = Cost, y = Willingness2Pay, color = Comparison, group = Comparison)) +
    annotate("text", y = 50000, x = 14300.690-522.43, label = "T-DXd Actual Monthly Cost", family = "Arial", size = 3) +
    annotate("text", y = 175000, x = cost2, label = paste0(round(cost2, 1), " USD"), family = "Arial", size = 3, vjust = 0.5) +
    annotate("text", y = 175000, x = cost1, label = paste0(round(cost1, 1), " USD"), family = "Arial", size = 3, vjust = 0.5) +
    geom_line(aes(linetype = Comparison, show.legend = FALSE)) +
    geom_point(aes(shape = Comparison), size = 3, show.legend = FALSE) +  # fixed size to 3
    geom_vline(xintercept = 14113.690, linetype = "dashed", color = "black") + 
    geom_hline(yintercept = 150000, linetype = "dotted", color = "black") + 
    geom_point(data = dot_df, aes(x = Cost, y = Willingness2Pay), shape = 23, fill = "green", size = 5, show.legend = FALSE) + 
    scale_color_manual(values = c("black", "black", "blue", "red")) + 
    scale_linetype_manual(values = c("solid", "dashed", "dotdash", "longdash")) + 
    scale_shape_manual(values = c(16, 17, 18, 19)) + 
    labs(x = "T-DXd Monthly Dose Cost ($)", y = "Cost-Effectiveness ($)", color = "Compared Strategies", linetype = "Compared Strategies", shape = "Compared Strategies") +
    scale_y_continuous(limits = c(0, 400000), breaks = seq(0, 500000, 40000)) +
    scale_x_continuous(limits = c(5000, 15000), breaks = seq(5000, 15000, 2500)) +
    theme_minimal() +
    ggtitle("One-way Sensitivity Analysis of T-DXd Cost") +
    theme(axis.line = element_line(size = 1, color = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          legend.position = "none",
          axis.text = element_text(family = "Arial"),
          axis.title = element_text(family = "Arial"),
          plot.title = element_text(hjust = 0.5, size = 18, face = "bold", family = "Arial"))
}

#Appendix Plot
df_one_way = list(df_tdxd_chemo, df_chemo_tdxd, df_tdxd_sg)
one_way_sensitivity_tdxd_price(df_one_way, dr_v = df_list_tdxd_chemo[[3]])
































#Plots that are not used in the final version of the manuscript

### Expected Loss Curves (ELCs) ----
elc_obj <- dampack::calc_exp_loss(wtp = v_wtp, psa = l_psa)
elc_obj

#* ELC plot
gg_elc <- plot(elc_obj, log_y = FALSE, 
               txtsize = txtsize, xlim = c(0, NA), n_x_ticks = 14,
               col = "full") +
  ggthemes::scale_color_colorblind() +
  ggthemes::scale_fill_colorblind() +
  # geom_point(aes(shape = as.name("Strategy"))) +
  scale_y_continuous("Expected Loss (Thousand $)", 
                     n.breaks = 10,
                     labels = function(x) x/1000) +
  theme(legend.position = c(0.4, 0.7),)
gg_elc

### Expected value of perfect information (EVPI) ----
#* Function included in "R/Functions.R". The latest version can be found in `dampack` package
evpi <- dampack::calc_evpi(wtp = v_wtp, psa = l_psa)
#* EVPI plot
gg_evpi <- plot(evpi, effect_units = "QALY", 
                txtsize = txtsize, xlim = c(0, NA), n_x_ticks = 14) +
  scale_y_continuous("EVPI (Thousand $)", 
                     n.breaks = 10,
                     labels = function(x) x/1000)
gg_evpi

### Combine all figures into one ----
patched_cea <- (gg_scattter +  gg_ceac + plot_layout(guides = "keep"))/(gg_elc + gg_evpi)
gg_psa_plots <- patched_cea + 
  plot_annotation(tag_levels = 'A')
gg_psa_plots

# Pairwise comparisons ----
## Chemo-Chemo vs Chemo-TDxD ----
l_psa_chemochemo_chemotdxd <- dampack::make_psa_obj(cost          = df_costs[, c(1, 2)],
                                                    effectiveness = df_qalys[, c(1, 2)],
                                                    strategies    = v_names_str[c(1, 2)])
### Incremental cost-effectiveness ratios (ICERs) with probabilistic output ----
#* Compute expected costs and effects for each strategy from the PSA
df_out_ce_psa_chemochemo_chemotdxd <- summary(l_psa_chemochemo_chemotdxd)
df_cea_psa_chemochemo_chemotdxd <- dampack::calculate_icers(cost       = df_out_ce_psa_chemochemo_chemotdxd$meanCost, 
                                                            effect     = df_out_ce_psa_chemochemo_chemotdxd$meanEffect,
                                                            strategies = df_out_ce_psa_chemochemo_chemotdxd$Strategy)
df_cea_psa_chemochemo_chemotdxd
### Plot cost-effectiveness frontier with probabilistic output ----
plot(df_cea_psa_chemochemo_chemotdxd, label = "all", txtsize = txtsize) +
  theme(legend.position = c(0.8, 0.2))

### Cost-effectiveness acceptability curves (CEACs) and frontier (CEAF) ---
ceac_obj_chemochemo_chemotdxd <- dampack::ceac(wtp = v_wtp, psa = l_psa_chemochemo_chemotdxd)
#* Regions of highest probability of cost-effectiveness for each strategy
summary(ceac_obj_chemochemo_chemotdxd)
#* CEAC & CEAF plot
gg_ceac_chemochemo_chemotdxd <- plot(ceac_obj_chemochemo_chemotdxd, txtsize = txtsize, xlim = c(0, NA), n_x_ticks = 14) +
  ggthemes::scale_color_colorblind() +
  ggthemes::scale_fill_colorblind() +
  theme(legend.position = c(0.8, 0.48))
gg_ceac_chemochemo_chemotdxd

## Chemo-Chemo vs TDxd-Chemo ----
l_psa_chemotdxd_chemochemo <- dampack::make_psa_obj(cost          = df_costs[, c(1, 3)],
                                                    effectiveness = df_qalys[, c(1, 3)],
                                                    strategies    = v_names_str[c(1, 3)])




#Old version of the PSA plot from Wes 
#source("PSA_8state.R")
#CEAC plot
#plot_ceac_curve <- function(df_psa_res, group_name1, group_name2, group_name3, group_name4){
#  X<-split(df_psa_res, df_psa_res$group)
#  Y = X$`TDxD-Chemo`
#  
#  
#  x = seq(50000, 400000, 10000)
#  y = c()
#  for(i in x){
#    y = append(y, sum(Y$'ICER'<i)/length(Y$'ICER'))
#  }

#  df <- data.frame(threshold = x,
#                   probability = y)
#  print(df)
#  
#  ggplot(data=df, aes(x=threshold, y=probability, group=1)) +
#    geom_line()+
#    geom_point()
#}

#plot_ceac_curve(res[[1]], "Chemo-Chemo","Chemo-TDxd", "TDxD-Chemo", "TDxD-SG")

#This plot function is used to plot the Kaplan-Meier curves for the different treatment sequences; currently just displaying OS/PSF for T-DXd (no overlay)
source("finding_transition_probabilities_8state.R")

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
    ggtitle("OS T-DXd")+
    scale_x_continuous(breaks = round(seq(0, 18, by = 3),1)) + 
    theme(axis.text = element_text(size = 12),
          axis.title = element_text(size = 14))
  
  plot4 <- ggplot(df, aes(x = idx)) + 
    geom_line(aes(y = km_pf_tdxd_model[1:m]), color = "red") + 
    geom_line(aes(y = km_pf_tdxd[1:m]), color = "blue", linetype = "twodash") +
    ylim(0, 1) +
    ylab("Probability") +
    xlab("Month") +
    ggtitle("PFS T-DXd") +
    scale_x_continuous(breaks = round(seq(0, 18, by = 3),1)) + 
    theme(axis.text = element_text(size = 12),
          axis.title = element_text(size = 14))
  
  print(ae_test_tdxd)
  print(ae_test_ild)
  plot(0:(length(ae_test_tdxd)-1), ae_test_tdxd, type = "l")
  plot(0:(length(ae_test_tdxd)-1), ae_test_ild, type = "l")
  
  ggarrange(plot2, plot4,
            ncol = 1, nrow = 2, common.legend = TRUE,legend="bottom") 
  
}

plot_function(tm_tdxd_chemo)





