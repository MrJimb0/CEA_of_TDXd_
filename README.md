# Cost-Effectiveness Analysis of T-DXd in Late-Stage HER2-Low Breast Cancer

This repository contains the code and data for a cost-effectiveness analysis (CEA) of Trastuzumab Deruxtecan (T-DXd) in patients with late-stage HER2-low breast cancer. The study employs Markov models to evaluate the economic value of T-DXd therapy compared to standard treatments.

Link to Article: https://ascopubs.org/doi/10.1200/JCO-24-01960 

## Project Overview

The analysis focuses on assessing the cost-effectiveness of T-DXd by modeling patient transitions through various health states. Two distinct models are utilized:

- 7-State Model: Used to derive specific transition probabilities.
- 8-State and CEA Model: Serves as the primary model for the base case analysis.

Probabilistic Sensitivity Analysis (PSA) is conducted to account for uncertainty in model parameters.

## Getting Started
### Prerequisites
- **R** (Version 4.0 or higher)
- Necessary packages:
install.packages(c("dampack", "ggplot2", "dplyr", "readr", "reshape2"))

## Running the Analysis
1. **Clone the Repository**
2. **Open the Project**: Launch RStudio and open the CEA_of_TDXd_.Rproj file.
3. **Execute Scripts**:
  - Run finding_transition_probabilities_7state.R to compute transition probabilities for the 7-state model.
  - Run finding_transition_probabilities_8state.R for the 8-state model.
  - Execute CEA_8state.R to perform the base case cost-effectiveness analysis.
  - Use PSA_8state.R to conduct the probabilistic sensitivity analysis.
  - Generate visualizations by running plot_file_8state.R.

## Notes
- The 7-state model is primarily used to derive specific transition rates, which are then integrated into the 8-state model for comprehensive analysis.
- The TNBC/ directory contains data pertinent to Triple-Negative Breast Cancer, which may be used for comparative purposes or extended analyses.
- The data folder contains the csv output from the main PSA used to create the graphs for the manuscript.
