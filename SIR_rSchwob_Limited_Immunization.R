
library(readr)
library(ggplot2)
library(tidyverse)
library(reshape2)
library(plotly)
library(RColorBrewer)
library(deSolve)

### ### ### Create an SIR function
sir <- function(time, state, parameters) {
  with(as.list(c(state, parameters)), {
    
    I = I1 + I2
    
    dS1 = - beta * S1 * I + gamma * I2
    dI1 = beta * S1 * I - gamma * I1
    
    dR = gamma * I1
    
    dS2 = - beta * S2 * I
    dI2 = beta * S2 * I - gamma * I2
    
    return(list(c(dS1, dI1, dR, dS2, dI2)))
  })
}

### ### ### Create a function to run the SIR model
run_sir <- function(s2_rate = 0, b = 1.5, g = 0.5, t = seq(0, 70, by = 0.01)) {
  ## s2_rate = Amount of people that only get immune after 2 infections

  ## Proportion in each compartment
  i_tot_zero <- 1e-6
  i2_zero <- i_tot_zero * s2_rate
  i1_zero <- i_tot_zero - i2_zero
  s_tot_zero <- 1 - i_tot_zero
  s2_zero <- s_tot_zero * s2_rate
  s1_zero <- s_tot_zero - s2_zero
  r_zero <- 0
  
  init <- c(S1 = s1_zero, I1 = i1_zero, R = r_zero, S2 = s2_zero, I2 = i2_zero)
  
  ## beta: infection parameter; gamma: recovery parameter
  parameters <- c(beta = b, gamma = g)
  
  ## Time frame
  times <- t
  
  ### ### ### Calculate S, I and R by time
  ## Solve using ode (General Solver for Ordinary Differential Equations)
  SIR_out <- ode(y = init, times = times, func = sir, parms = parameters)
  
  ## Change to data frame, name columns appropriately
  SIR_df <- as.data.frame(SIR_out)
  SIR_df$S <- SIR_df$S1 + SIR_df$S2
  SIR_df$I <- SIR_df$I1 + SIR_df$I2
  
  ### ### ### Calculate important values
  I_max <- max(SIR_df$I)
  HIT <- 1 - SIR_df$S[which.max(SIR_df$I)]
  t_HIT <- SIR_df$time[which.max(SIR_df$I)]
  S_end <- tail(SIR_df$S, n = 1)
  R_end <- tail(SIR_df$R, n = 1)
  
  ### ### ### Output
  names(SIR_df) = c("Time","Susceptible 1","Infected 1","Recovered", "Susceptible 2","Infected 2", "Susceptible Total", "Infected Total")
  return(list("I_max" = I_max, "HIT" = HIT, "t_HIT" = t_HIT, "S_end" = S_end, "R_end" = R_end, "SIR_df" = SIR_df))
}

### ### ### Create a function to plot the results of the SIR model
plot_sir <- function(SIR_df, detailed = FALSE, pdf = "") {
  ## Make linear data frame (with additional column for measure)
  if (detailed){
    SIR_lin = melt(SIR_df,id = "Time", measure = c("Susceptible 1","Infected 1","Recovered", "Susceptible 2","Infected 2"))
  }
  else{
    SIR_lin <- melt(SIR_df,id = "Time", measure = c("Susceptible Total","Infected Total","Recovered"))
  }
  
  names(SIR_lin) <- c("Time", "Compartment", "Value")
  
  ### ### ### Plot
  ## Set colours
  if (detailed){
    colours <- c("#999999", "#E69F00", "#56B4E9", "#999900", "#E69F99")
  }
  else{
    colours <- c("#999999", "#E69F00", "#56B4E9")
  }
  
  ## Make a plot
  pp1 = ggplot(SIR_lin) +
    geom_line(aes(x = Time, y = Value, color = Compartment), size = 1.2) +
    theme_minimal() +
    xlab ("Time") +
    ylab("Proportion of Population") +
    theme_classic() + 
    theme(text = element_text(size = 20)) +
    ylim(0,1) + 
    scale_color_manual(values = colours)

  ## Show graph or save as pdf
  if (pdf != ""){
    if (detailed){ ggsave(paste0(pdf, "_(detailed).pdf"), height = 4 , width = 7) }
    else{ ggsave(paste0(pdf, ".pdf"), height = 4 , width = 7) }
  }
  else{
    show(pp1)
  }
}

### ### ### Run the model with multiple s2 rates and compare the output (for a given beta)
compare_s2_rates <- function(s2_rates, beta, gamma, t = seq(0, 70, by = 0.01), csv = FALSE){
  num_runs <- length(s2_rates)
  results_df <- data.frame("s2_rate" = c(num_runs), "beta" = c(num_runs), "gamma" = c(num_runs), "I_max" = c(num_runs), "HIT" = c(num_runs), "t_HIT" = c(num_runs), "S_end" = c(num_runs), "R_end" = c(num_runs))
  
  cur_row <- 1
  for (s2_rate in s2_rates){
    results <- run_sir(s2_rate, beta, gamma, t)
    results_df[cur_row,] <- c(s2_rate, beta, gamma, results$I_max, results$HIT, results$t_HIT, results$S_end, results$R_end)
    cur_row <- cur_row + 1
  }
  
  if(csv == TRUE){write.csv2(results_df, file = paste0("beta_", beta, "_gamma_", gamma, ".csv"))}

  return(results_df)
}

### ### ### Functions to plot and compare important values of the model for different values of S2_rate 

## Maximum Infected reached
plot_I_max_onePlot <- function(s2_rates, betas, gamma, t = seq(0, 70, by = 0.01), csv = FALSE){
  pdf("I_max_by_s2_rates_for_different_betas_onePlot.pdf", 5, 7)
  plot(NA, type = "b", xlab = "S2 rate", ylab = "Maximum Infected", xlim=range(0, 1), ylim = c(0, 0.9))
  i <- 1
  for (beta in betas){
    results_df <- compare_s2_rates(s2_rates, beta, gamma, t, csv = csv)
    lines(results_df$s2_rate, results_df$I_max, type = "b")
    text(0.8, c(0.068, 0.16, 0.31, 0.61, 0.8)[i], label = paste0("beta = ", beta), cex = 0.8, pos = c(3, 3, 3, 3, 3)[i])
    
    i <- i + 1
  }
  dev.off()
}

## Time until HIT is reached
plot_t_HIT_onePlot <- function(s2_rates, betas, gamma, t = seq(0, 70, by = 0.01), csv = FALSE){
  pdf("t_HIT_by_s2_rates_for_different_betas_onePlot.pdf", 5, 7)
  plot(NA, type = "b", xlab = "S2 rate", ylab = "Time to HIT", xlim=range(0, 1), ylim = c(2, 50))
  i <- 1
  for (beta in betas){
    results_df <- compare_s2_rates(s2_rates, beta, gamma, t, csv = csv)
    lines(results_df$s2_rate, results_df$t_HIT, type = "b")
    text(0.5, results_df$t_HIT[which(results_df$s2_rate == 0.5)] * 0.99, label = paste0("beta = ", beta), cex = 0.8, pos = c(1, 1, 1, 1, 1)[i])
    
    i <- i + 1
  }
  dev.off()
}

## Final S value

plot_S_end_onePlot <- function(s2_rates, betas, gamma, t = seq(0, 70, by = 0.01), csv = FALSE){
  pdf("S_end_by_s2_rates_for_different_betas_onePlot.pdf", 5, 7)
  plot(NA, type = "b", xlab = "S2 rate", ylab = "Final Susceptible", xlim=range(0, 1), ylim = c(0, 0.7))
  i <- 1
  for (beta in betas){
    results_df <- compare_s2_rates(s2_rates, beta, gamma, t, csv = csv)
    s2_peak <- results_df$s2_rate[which.max(results_df$S_end)]
    lines(results_df$s2_rate, results_df$S_end, type = "b")
    text(s2_peak, max(results_df$S_end), label = paste0("beta = ", beta), cex = 0.8, pos = c(3, 3, 3, 3, 1)[i])
    
    i <- i + 1
  }
  dev.off()
}

plot_S_end_separate <- function(s2_rates, betas, gamma, t = seq(0, 70, by = 0.01), csv = FALSE){
  pdf("S_end_by_s2_rates_for_different_betas_separate.pdf", 12, 4)
  par(mfrow = c(1, 3))
  i <- 1
  for (beta in betas){
    results_df <- compare_s2_rates(s2_rates, beta, gamma, t, csv = csv)
    s2_peak <- results_df$s2_rate[which.max(results_df$S_end)]
    R0 <- beta / gamma
    plot(results_df$s2_rate, results_df$S_end, type = "b", xlab = "S2 rate", ylab = "Final Susceptible", main = paste0("beta = ", beta, ", gamma = ", gamma, " (R0 = ", R0, ")"))
    text(0.4, min(results_df$S_end) + (max(results_df$S_end)-min(results_df$S_end))/2, label = paste0("S2 peak: ", s2_peak))
    i <- i + 1
  }
  dev.off()
}



### ### ### Set parameters and run the model
gamma <- 1

## Create SIR plots (compartments by time) for betas and S2_rates, both with and without detailed depiction of S1 and S2
s2_rates <- seq(0, 1, by = 0.5) # seq(0, 1, by = 0.25)
betas <- 1.5 # c(1.25, 1.5, 2, 3, 4, 8)
for (s2_rate in s2_rates){
  for (beta in betas){
    results <- run_sir(s2_rate, beta, gamma)
    plot_sir(results$SIR_df, detailed = TRUE, pdf = paste0(s2_rate, "_", beta))
    plot_sir(results$SIR_df, detailed = FALSE, pdf = paste0(s2_rate, "_", beta))
  }
}

## Calculate and plot important values for different values for beta and S2_rate
s2_rates <- seq(0, 1, by = 0.05)
betas <- c(1.25, 1.5, 2, 4, 8)
times = seq(0, 70, by = 0.01)

plot_I_max_onePlot(s2_rates, betas, gamma, times, csv = FALSE)
plot_t_HIT_onePlot(s2_rates, betas, gamma, times, csv = FALSE)

times = seq(0, 700, by = 0.01) # To ensure that always the final value for S is calculated, end time is increased
plot_S_end_onePlot(s2_rates, betas, gamma, times, csv = FALSE)

betas <- c(1.25, 2, 4) # Compare only 3 betas in separate plots
plot_S_end_separate(s2_rates, betas, gamma, times, csv = TRUE) # Also create csv file with all important values
