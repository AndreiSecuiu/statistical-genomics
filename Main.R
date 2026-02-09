########## Description ##########

#' This file contains the main code used for the simulation study in which 
#' Simulated Annealing (SA) is compared to structure MCMC and Gradient Ascent

########## Preamble ##########

# Libraries
library(stats)
library(ROC)
library(limma)

library(ggplot2)
library(reshape2) 
library(dplyr)

# Import custom functions (Marco's and mine)

# CHANGE THE PATH TO YOUR LOCAL WORKING DIRECTORY
folder_path = "C:\\Users\\andre\\Python Notebooks\\statistical-genomics"
setwd(folder_path)

source("Functions.R")

set.seed(42)

# - Global parameters

true_Net = make_true_Net()        # The gold standard incidence matrix 

# Concerning the running of functions

N_repl_per_m = 20         # How many data sets to generate per 'm' in m_list    
m_list = c(20, 50, 100)   # How many samples 'm' to generate                   
nr_restarts = 6           # How many restarts to use per data set                      

N_MCMC = 10^5            # Number MCMC iterations 
N_SA   = 2 * 10^3        # Number SA   iterations (default value)
N_GA   = 10^2            # Maximum number of iterations for GA

thin_factor = 100        # Thinning factor (reduce auto-correlation)
frac_burn = 0.5          # Fraction of the samples to eliminate, for burn in
noise_var = 1            # Noise parameter used during data generation

# The methods to be investigated (incl. temp. schedules)

method_names = c("strMCMC", "GA", "SA (exp.2)", "SA (exp.4)", "SA (pow.2)",
                 "SA (pow.4)")
nr_methods = length(method_names)

########## Custom functions ##########  

#' Global parameters that affect the functions below (e.g. the number of MCMC 
#' iterations) are *hard coded*, but they can be easily changed in the preamble

# ---------- For running the main algorithms ----------

run_strMCMC = function(data, seed=NULL){
  #' Wrapper function for running strMCMC with the pre-defined parameters
  #' (fixed globally) and eliminating the burn-in samples
  
  if(!is.null(seed)){set.seed(seed)}
  
  OUT = strMCMC(data, matrix(data = 0, nrow = 11, ncol = 11), N_MCMC, thin_factor) 
  
  # Get rid of burn-in samples
  n = dim(OUT[[1]])[3]
  n_burn = ceiling(frac_burn * n)
  
  OUT[[2]] = OUT[[2]][(n_burn:n)]    # BGe of MCMC samples (for trace plots)
  OUT[[1]] = OUT[[1]][,,(n_burn:n)]  # List of all incidence matrices
  
  return(OUT)
}

run_strGA = function(data, inc_mat_start = NULL, N_iter = N_GA,
                     restarts = 1, seed=NULL){
  #' Wrapper function for running strGA with the pre-defined parameters
  #' (fixed globally). The function behaves differently depending on the hyperparameters:
  #' 
  #'   ~ if inc_mat_start is provided (= the incidence matrix for the start of 
  #'   the investigation), then the function will output the result of strGA for
  #'   that starting point
  #'   ~ if inc_mat_start is NULL, then the function outputs a list with the outputs
  #'   of "restarts" straGA runs, where the starting points are *randomly* initialized
  
  if(!is.null(inc_mat_start)){
    OUT = strGA(data, inc_mat_start, N_iter)
  } else{
    
    # Initialize the seed, if given
    if(!is.null(seed)){set.seed(seed)}
    
    # Run strGA repeatedly
    OUT = list()
    
    for(i in 1:restarts){
      inc_mat_start = GENERATE_RANDOM_DAG(11) 
      OUT[[i]] = strGA(data, inc_mat_start, N_iter)
    }
  }
  
  return(OUT)
}

run_strSA = function(data, temperature_schedule, 
                     restarts = 1, seed=NULL){
  #' Wrapper function for running strSA with the pre-defined parameters
  #' (fixed globally). The starting DAG is always the empty graph.
  #' 
  #' The function outputs a list with the outputs of "restarts" strSA runs
  
  if(!is.null(seed)){set.seed(seed)}
  
  OUT = list()
  N_iter = length(temperature_schedule)
  
  for(i in 1:restarts){
    OUT[[i]] = strSA(data, matrix(data = 0, nrow = 11, ncol = 11), N_iter, temperature_schedule)
  }
  
  return(OUT)
}

# ---------- For calculating & comparing AUROC values ----------

get_AUROC = function(data, temperature_schedules, seed=NULL){  
  #' For one data set, this function computes the AUROC values for all methods
  #' investigated by this project
  #' 
  #' @section Inputs:
  #'   
  #'   @param data A matrix of dimensions (n, m), where n = 11 is the number of 
  #'   nodes and m is the number of data points
  #'   
  #'   @param temperature_schedules List of arrays. Each array is the temperature
  #'   schedule to be used for SA, at the corresponding method. In order: Exp. 2,
  #'   Exp. 4, Pow. 2, Pow. 4
  #'   
  #'   @param seed A positive integer which provides the seed to be initialized.
  #'   If not provided, then no seed is initialized
  #'   
  #' @section Output:
  #' 
  #'   An array of length 6 containing AUROC values for the following routines:
  #'   
  #'   ~ strMCMC (the baseline)
  #'   ~ GA
  #'   ~ SA: one output for each temperature schedule
  
  if(!is.null(seed)){set.seed(seed)}
  
  # Initialize list of AUROC to be returned
  OUT = array(dim = nr_methods)
  
  # ------ strMCMC ------

  OUT_MCMC = run_strMCMC(data)   
  n = dim(OUT_MCMC[[1]])[3]      # All effective MCMC samples left
  
  # Empirical marginal edge posterior
  SCORES = apply(OUT_MCMC[[1]], 1:2, sum)/n
  OUT[1] = compute_AUROC(SCORES, true_Net)
  
  # ------ GA ------
  
  OUT_GA = run_strGA(data, restarts = nr_restarts)
  
  # Among all local maxima, find the largest
  
  log_BGe_GA = array(dim=nr_restarts)

  for(i in 1:nr_restarts){
    temp = OUT_GA[[i]]$log_bge_trace
    log_BGe_GA[i] = temp[length(temp)]
  }
  
  max_idx = which.max(log_BGe_GA)
  SCORES = OUT_GA[[max_idx]]$incidence_best
  OUT[2] = compute_AUROC(SCORES, true_Net)
  
  # ------ SA ------
  
  # Iterate over the list of temperature schedules
  # Take the DAG with the maximum score among "nr_restarts" SA runs
  
  for(i in 1:length(temperature_schedules)){
    log_BGe_SA = array(dim=nr_restarts)
    temperature_schedule = temperature_schedules[[i]]
    
    OUT_SA = run_strSA(data, temperature_schedule, restarts = nr_restarts)
    
    for(j in 1:nr_restarts){
      temp = OUT_SA[[j]]$log_bge_trace
      log_BGe_SA[j] = temp[length(temp)]
    }
    
    # Take the DAG with the best score
    
    max_idx = which.max(log_BGe_SA)
    SCORES = OUT_SA[[max_idx]]$incidence_final
    OUT[i+2] = compute_AUROC(SCORES, true_Net)
  }
  
  return(OUT)
}

compare_AUROC = function(cat_1, cat_2, AUROC_matrix, mean_std_table){
  #' This function performs the comparison of vectors of AUROC replications in a
  #' statistically informed manner. Following Marco's suggestion:
  #' 
  #' In the t-test confidence interval plot, see if the CI overlap:
  #'   - no : methods are significantly different
  #'   - yes: look at the difference among replications
  #'     and statistically test that the mean is zero (rely on the (un)paired t-test)
  #'       - reject: methods are significantly different
  #'       - accept: no evidence of them being significantly different
  #' 
  #' *Note*: The paired test should be used if the 'm's of the two categories 
  #' coincide because AUROC values were computed based on the same data set. 
  #' Otherwise there is no natural pairing
  #' 
  #'
  #' @section Inputs: 
  #' 
  #' @param cat_1, cat_2: array of length 2. The first entry represents 'm' (i.e. the 
  #' index of the desired 'm' in m_list), the second represents the method used to
  #' obtain AUROC (i.e. the index of the desired method in method_names)
  #' 
  #' @param AUROC_matrix, mean_std_table: same as in the main file
  #' 
  #' @section Output: 
  #' 
  #' If the CI do not overlap, return -1. Else return the p-value of the relevant
  #' t-test
  
  # - Unpack the replicated values & CI
  
  vec_1 = AUROC_matrix[,cat_1[1], cat_1[2]]
  vec_2 = AUROC_matrix[,cat_2[1], cat_2[2]]
  
  CI_1 = c(mean_std_table[cat_1[1],cat_1[2],1] - mean_std_table[cat_1[1],cat_1[2],3], 
           mean_std_table[cat_1[1],cat_1[2],1] + mean_std_table[cat_1[1],cat_1[2],3])
  CI_2 = c(mean_std_table[cat_2[1],cat_2[2],1] - mean_std_table[cat_2[1],cat_2[2],3], 
           mean_std_table[cat_2[1],cat_2[2],1] + mean_std_table[cat_2[1],cat_2[2],3])
  
  # - Do the comparisons
  
  # 1. Do the CI overlap?
  
  if( (CI_1[2] < CI_2[1]) | (CI_2[2] < CI_1[1])){return(-1)} else{
    # 2. Should the samples be paired?
    
    if(cat_1[1] == cat_2[1]){
      test_OUT = t.test(vec_1, vec_2, paired = TRUE)
    } else{
      test_OUT = t.test(vec_1, vec_2, paired = FALSE)
    }
    
    return(test_OUT$p.value)
  }
}

########## Preliminary investigations ##########

# Default values for preliminary investigations:

temperature_schedule_list = cbind((1 + 10^(-2))^((-1)* 1:N_SA), 
                                  (1 + 10^(-4))^((-1)* 1:N_SA),
                                  (1 - (1:N_SA)/(N_SA + 1))^2,
                                  (1 - (1:N_SA)/(N_SA + 1))^4) 
colnames(temperature_schedule_list) = c("Exp. (2)", "Exp. (4)", 
                                        "Pow. (2)", "Pow. (4)")

# --- Convergence diagnostics of MCMC iterates with the current settings ---

#' Investigate whether MCMC iterations converge to the stationary distribution
#'   ~ trace plots of log_BGe scores for repeated runs
#'   ~ scatter plots of marginal edge posterior probabilities

set.seed(42)

idx_m = 3  # Go through all "m" in m_list
data_diagnostic = make_test_Data(m_list[idx_m], noise_var)

# Run 3 simulations (time: 2min30sec)
init_mat = matrix(data = 0, nrow = 11, ncol = 11)  # Starting from the empty graph

MCMC_1 = strMCMC(data_diagnostic, init_mat, N_MCMC, thin_factor)
MCMC_2 = strMCMC(data_diagnostic, init_mat, N_MCMC, thin_factor)
MCMC_3 = strMCMC(data_diagnostic, init_mat, N_MCMC, thin_factor)

# - i) Trace plot diagnostics

BGe_1 = unlist(MCMC_1[[2]])
BGe_2 = unlist(MCMC_2[[2]])
BGe_3 = unlist(MCMC_3[[2]])

n_sam = length(BGe_1)
idx = seq(1, n_sam, 1)
burn_in_idx = ceiling(frac_burn * n_sam)

plot(idx[burn_in_idx:n_sam], BGe_1[burn_in_idx:n_sam], col = "blue", type = "l",
     xlab = "Index", ylab = "log(BGe)", 
     main = paste0("Trace plot diagnostic\nm = ", m_list[idx_m]))
lines(idx[burn_in_idx:n_sam], BGe_2[burn_in_idx:n_sam], col = "red")
lines(idx[burn_in_idx:n_sam], BGe_3[burn_in_idx:n_sam], col = "green")

# - ii) Scatter plot of marginal edge posteriors

# Find the CPDAGs
marg_edge_1 = cpdag_list(MCMC_1[[1]], burn_in_idx)
marg_edge_2 = cpdag_list(MCMC_2[[1]], burn_in_idx)
marg_edge_3 = cpdag_list(MCMC_3[[1]], burn_in_idx)

# Flatten and discard diagonal entries
marg_edge_flat_1 = marg_edge_1[[3]]
marg_edge_flat_1 = marg_edge_flat_1[row(marg_edge_flat_1) != col(marg_edge_flat_1)]
marg_edge_flat_2 = marg_edge_2[[3]]
marg_edge_flat_2 = marg_edge_flat_2[row(marg_edge_flat_2) != col(marg_edge_flat_2)]
marg_edge_flat_3 = marg_edge_3[[3]]
marg_edge_flat_3 = marg_edge_flat_3[row(marg_edge_flat_3) != col(marg_edge_flat_3)]

# Plot
plot(marg_edge_flat_1, marg_edge_flat_2, pch=4, xlab="Run 1", ylab="Run 2", 
     main = "1 versus 2")
lines(0:1, 0:1, type="l", col = "red")

plot(marg_edge_flat_2, marg_edge_flat_3, pch=4, xlab="Run 2", ylab="Run 3", 
     main = "2 versus 3")
lines(0:1, 0:1, type="l", col = "red")

plot(marg_edge_flat_1, marg_edge_flat_3, pch=4, xlab="Run 1", ylab="Run 3", 
     main = "1 versus 3")
lines(0:1, 0:1, type="l", col = "red")

# --- How many iterations does it take, on average, for GA to converge? ---

#' Across many repetitions, investigate how many steps it takes the GA algorithm
#' to converge to a local minimum. Repeat for a few data sets, to investigate
#' robustness of the conclusions

set.seed(42)

idx_m = 3          # Go through all "m" in m_list
data_diagnostic = make_test_Data(m_list[idx_m], noise_var)

repeat_test = 30   # How many times should GA be restarted, for the test?

OUT = run_strGA(data_diagnostic, restarts=repeat_test)
iterates_until_convergence = array(dim=repeat_test)
successful_runs = array(dim=repeat_test)

for(i in 1:repeat_test){
  iterates_until_convergence[i] = length(OUT[[i]]$log_bge_trace)
  successful_runs[i] = OUT[[i]]$optimum_found
}

# Show how many iterations actually converged
cat("Runs that converged within ", N_GA, " iterations: ", 
    sum(successful_runs), "/", repeat_test, "\n")

# Display the distribution via a simple boxplot
boxplot(iterates_until_convergence[successful_runs], 
        ylab="Iter. until conv.")  # Only for successful runs!

min(iterates_until_convergence[successful_runs])
quantile(iterates_until_convergence[successful_runs], 0.25)
quantile(iterates_until_convergence[successful_runs], 0.50)
quantile(iterates_until_convergence[successful_runs], 0.75)
max(iterates_until_convergence[successful_runs])

# --- Evolution of SA iterates with time ---

#' Track the BGe scores for multiple SA iterates across iterations. Behaviour?
#' Track the changes across different temperature schedules too

#' For the given temperature schedules, we note:
#'   ~ whether convergence is achieved within the given number of iterations 
#'   (and if so, where)
#'   ~ t-test CI for the log_BGe values at the end of each routine
#'   ~ repeat for other datasets/m to see if the conclusions hold

# - Visualize the temperature schedules

plot_mat = melt(log(temperature_schedule_list))
colnames(plot_mat) = c("Time", "Schedule", "Temperature")

ggplot(plot_mat, aes(x = Time, y = Temperature, color = Schedule)) +
  geom_line(linewidth = 1) +
  labs(title = "Temperature schedules",
       x = "Iteration t",
       y = "log(Temperature)") +
  theme_minimal() 

# - The convergence investigations

set.seed(42)

repeat_test = 5    # How many times should SA be restarted, for the test?
idx_m = 2          # Vary m over this

data_diagnostic = make_test_Data(m_list[idx_m], noise_var)


# Save all traces in one compact array
log_BGe_traces = array(dim = c(N_SA, repeat_test, ncol(temperature_schedule_list)))

for(idx_temp in 1:ncol(temperature_schedule_list)){
  OUT = run_strSA(data_diagnostic, temperature_schedule_list[,idx_temp], restarts=repeat_test)
  
  for(i in 1:repeat_test){log_BGe_traces[,i,idx_temp] = OUT[[i]]$log_bge_trace}
}

# Plot

dimnames(log_BGe_traces) = list(
  Iteration = 1:N_SA,
  Run = 1:repeat_test, 
  Schedule = colnames(temperature_schedule_list) 
)

plot_data = melt(log_BGe_traces, value.name = "log_BGe")

ggplot(plot_data, aes(x = Iteration, y = log_BGe, color = Schedule)) + 
  geom_line(aes(group = interaction(Run, Schedule)), alpha = 0.8, linewidth = 0.7) + 
  scale_y_continuous(limits = c(-500, -450)) +
  labs(
    title = paste("Trace plots for SA repeated runs distinguished by temperature schedules\nm = ", m_list[idx_m]),
    x = "Iteration t",
    y = "log BGe",
    color = "Schedule"
  ) +
  facet_wrap(~ Schedule, scales = "fixed") +
  theme_minimal()

# - t-test CI for log_BGe per method

alpha = 0.05
means = apply(log_BGe_traces[N_SA,,], 2, mean)
widths = apply(log_BGe_traces[N_SA,,], 2, sd) * qt(1 - alpha/2, repeat_test - 1) / sqrt(repeat_test)
max_sample = apply(log_BGe_traces[N_SA,,], 2, max)

ci_table = rbind(Mean = means, 
                 CI_width = widths,
                 max = max_sample) 
ci_table

# --- How many restarts should one make for SA/GA? ---

#' One way to investigate the quality of methods is via the probability of a new
#' estimate to be larger than all previous samples (analytic approximation available). 
#' We may also wish to determine the expected improvement of the score

# Initialize relevant variables

repeat_test_GA = 31    # How many times should GA/SA be restarted, for the test?
repeat_test_SA = 31    

# - Compute the log_BGe for repeated routines

if(FALSE){  # TRUE = Compute the log scores for GA/SA again. FALSE = use info from the saved file
  log_BGe_GA = array(dim=c(length(m_list), repeat_test_GA))
  log_BGe_SA = array(dim=c(length(m_list), repeat_test_SA, ncol(temperature_schedule_list)))
  
  set.seed(42)  # Set the seed for reproducibility
  
  for(idx_m in 1:length(m_list)){
    # Generate one data set for the new m
    data_diagnostic = make_test_Data(m_list[idx_m], noise_var)
    
    # Run GA iterations & save log_BGe for end graphs
    OUT_GA = run_strGA(data_diagnostic, restarts=repeat_test_GA)
    
    for(i in 1:repeat_test_GA){
      temp = OUT_GA[[i]]$log_bge_trace
      log_BGe_GA[idx_m, i] = temp[length(temp)]
    }
    
    # Repeat for SA
    for(idx_temp in 1:ncol(temperature_schedule_list)){
      OUT_SA = run_strSA(data_diagnostic, temperature_schedule_list[, idx_temp],
                         restarts=repeat_test_SA)
      
      for(i in 1:repeat_test_SA){
        temp = OUT_SA[[i]]$log_bge_trace
        log_BGe_SA[idx_m, i, idx_temp] = temp[length(temp)]
      }
    }
  }
  
  # Save array for future replotting
  saveRDS(log_BGe_GA, "log_BGe_GA_improvAnal.rds")
  saveRDS(log_BGe_SA, "log_BGe_SA_improvAnal.rds")

} else{
  log_BGe_GA = readRDS("log_BGe_GA_improvAnal.rds")
  log_BGe_SA = readRDS("log_BGe_SA_improvAnal.rds")
}

# - Do the expected improvement analysis

exp_imp_GA = array(dim=c(length(m_list), repeat_test_GA-1))  # Save the expected improvements per m
exp_max_GA = array(dim=length(m_list))                       # Save the expected max FOR K = 1

exp_imp_SA = array(dim=c(length(m_list), repeat_test_GA-1, ncol(temperature_schedule_list)))  
exp_max_SA = array(dim=c(length(m_list),  ncol(temperature_schedule_list)))                     

for(idx_m in 1:length(m_list)){
  # - GA

  OUT = restarts_benefit_analysis(log_BGe_GA[idx_m,])
  
  exp_max_GA[idx_m] = OUT$exp_max[1]  # Save expected max for one sample
  exp_imp_GA[idx_m,] = OUT$exp_imp    # Save the improvements 

  # - SA
  
  for(idx_temp in 1:ncol(temperature_schedule_list)){
    OUT = restarts_benefit_analysis(log_BGe_SA[idx_m,,idx_temp])
    
    exp_max_SA[idx_m, idx_temp] = OUT$exp_max[1]
    exp_imp_SA[idx_m,,idx_temp] = OUT$exp_imp
  }
}

# - Print expected maxima for K = 1

cat("Expected maximum (K = 1) - GA:\n")
cat(exp_max_GA, "\n\n")

cat("Expected maximum (K = 1) - SA:\n")
for(idx_temp in 1:ncol(temperature_schedule_list)){
  cat("  -> ", sa_schedule_names[idx_temp], ":", exp_max_SA[,idx_temp], "\n")
}

# - Plot

# i) Prepare the data for plotting

sa_schedule_names = c("Exp. (2)", "Exp. (4)", "Pow. (2)", "Pow. (4)")

df_ga = melt(exp_imp_GA)
colnames(df_ga) = c("m_idx", "repeats", "value")

df_ga = df_ga %>%
  mutate(
    m = factor(m_list[m_idx]),     
    method_type = "GA",            
    schedule = ""                
  )

df_sa = melt(exp_imp_SA)
colnames(df_sa) = c("m_idx", "repeats", "sched_idx", "value")

df_sa = df_sa %>%
  mutate(
    m = factor(m_list[m_idx]),
    method_type = "SA",
    schedule = factor(sa_schedule_names[sched_idx], levels = sa_schedule_names)
  )

# Combine data
plot_data = bind_rows(
  df_ga %>% select(repeats, value, m, method_type, schedule),
  df_sa %>% select(repeats, value, m, method_type, schedule)
)

# Ensure 'm' factor levels are ordered correctly for the legend
plot_data$m = factor(plot_data$m, levels = sort(m_list))

# ii) Actual plotting

ggplot(plot_data, aes(x = repeats, y = value)) +
  geom_step(aes(color = schedule, linetype = m), linewidth = 1) +
  #facet_grid(. ~ method_type) + 
  facet_wrap(~ method_type + schedule, ncol = 2, scales = "fixed") +
  scale_x_continuous(
    name = "Number of repeats",
    limits = c(1, 10),
    breaks = seq(1, max(plot_data$repeats), 1), 
    labels = seq(1, max(plot_data$repeats), 1),
    expand = c(0, 0)                            
  ) +
  scale_y_continuous(name = "max log BGe. diff.") +
  scale_color_manual(
    name = "Schedule",
    values = c("GA" = "black", 
               "Exp. (2)" = "#E41A1C", 
               "Exp. (4)" = "#377EB8", 
               "Pow. (2)" = "#4DAF4A", 
               "Pow. (4)" = "#984EA3")
  ) +
  scale_linetype_discrete(name = "Data Set Size (m)") +
  ggtitle("Expected improvement of score via repeated fitting") +
  theme_bw() +
  theme(
    strip.background = element_rect(fill = "grey90"), 
    strip.text = element_text(size = 12, face = "bold"), 
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    legend.position = "bottom", 
    legend.box = "vertical",    
    panel.grid.minor.x = element_blank() 
  )

########## Compute AUROC for all methods ##########

#' Save the obtained AUROC values in a separate R file. 
#' Then read the file so we don't have to re-run everything

# - Define list of temperature schedules (based on earlier findings)

#temperature_schedules = split(temperature_schedule_list, col(temperature_schedule_list))

# Scale of the number of iterations for Pow. scale
N_SA_pow = 5 * 10^3  
# Fractions found empirically (in the preliminary investigation)
iter_pow = c(ceiling(N_SA_pow * 0.9), ceiling(N_SA_pow * 0.65))

temperature_schedules = list(
  (1 + 10^(-2))^((-1)* 1:1000),
  (1 + 10^(-4))^((-1)* 1:46000),
  ((1 - (1:N_SA_pow)/(N_SA_pow + 1))^2)[1:iter_pow[1]],
  ((1 - (1:N_SA_pow)/(N_SA_pow + 1))^4)[1:iter_pow[2]]
)

# - The actual computation

# One full run, per m, per data set: 3min16sec
# Estimated time for full run: 11760 sec = 3.3h

saved_file = TRUE

if(saved_file){AUROC_matrix = readRDS("AUROC_matrix.rds")} else{
  set.seed(42)
  
  # Init. final list of AUROC
  # Shape: (Replications per m) x (Number of m) x (Number AUROC methods)
  AUROC_matrix = array(dim=c(N_repl_per_m, length(m_list), length(method_names)))
  
  for(i in 1:length(m_list)){  # Loop over m
    m = m_list[i]
    
    for(j in 1:N_repl_per_m){    # Repeat 'N_repl_per_m' times
      data = make_test_Data(m, noise_var)
      
      # Get all desired AUROCs & store
      OUT = get_AUROC(data, temperature_schedules)
      AUROC_matrix[j, i, ] = OUT
    }
  }
  
  saveRDS(AUROC_matrix, "AUROC_matrix.rds")
}

# --- Table of means, standard deviations & t-test CI

mean_std_table = array(dim = c(dim(AUROC_matrix)[2], dim(AUROC_matrix)[3], 3))
alpha = 0.05

# mean_std_table[i, j, 1] stores the mean of AUROC_matrix[,i,j]
# mean_std_table[i, j, 2] stores the standard deviation of AUROC_matrix[,i,j]
# mean_std_table[i, j, 3] stores the 0.5*length of the t-test CI of AUROC_matrix[,i,j]

mean_std_table[,,1] = apply(AUROC_matrix, c(2, 3), mean)
mean_std_table[,,2] = apply(AUROC_matrix, c(2, 3), sd)
mean_std_table[,,3] = mean_std_table[,,2] * qt(1 - alpha/2, N_repl_per_m - 1) / sqrt(N_repl_per_m)


########## Systematic comparison of AUROC ##########

# --- 1. Box plots ---

# For displaying in the plots
m_list_char = as.character(m_list)

# Set dimension names for clarity during melting 
dimnames(AUROC_matrix) = list(
  NULL, 
  m_list_char, 
  method_names
)

# - Reshape the 3D array into a long data frame ---
AUROC_df = melt(AUROC_matrix, 
                varnames = c("Replication", "m_value", "Method"), 
                value.name = "AUROC")

# - Ensure factors are set correctly for plotting order
AUROC_df$m_value = factor(AUROC_df$m_value, levels = m_list_char)
AUROC_df$Method = factor(AUROC_df$Method, levels = method_names)

# - Create the plot

p = ggplot(AUROC_df, aes(x = Method, y = AUROC, fill = Method)) +
  geom_boxplot(outlier.shape = NA) + 
  geom_point(position = position_jitter(width = 0.15), 
             alpha = 0.4, size = 1.5, shape = 16, color = "black") + 
  facet_wrap(~ m_value, 
             nrow = 1, 
             labeller = labeller(m_value = function(x) paste("m =", x))) +
  labs(title = "Boxplot representation of AUROC values for all methods\nSeparated by 'm'",
       y = "AUROC Value",
       x = "Method") +
  scale_y_continuous(limits = c(0.5, 1)) +
  scale_fill_brewer(palette = "Set1") +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    strip.text = element_text(face = "bold", size = 12), # Subplot titles 
    axis.text.x = element_text(angle = 45, hjust = 1),   # Rotated x-axis labels
    legend.position = "none"                             # Hide legend
  )

print(p)

# --- 2. Error-bar plot ---

# - Prepare data frame for plotting

Mean_matrix = mean_std_table[,,1]
SD_matrix = mean_std_table[,,2]

# Ensure column names for melting are the method names
colnames(Mean_matrix) = method_names
colnames(SD_matrix) = method_names

# Add row names (m_list values) to ensure they are captured during melting
rownames(Mean_matrix) = m_list_char
rownames(SD_matrix) = m_list_char

mean_df = melt(Mean_matrix, varnames = c("m_value", "Method"), value.name = "Mean_AUROC")
sd_df = melt(SD_matrix, varnames = c("m_value", "Method"), value.name = "SD_AUROC")

# Combine all into the final data frame
errorbar_df = data.frame(
  m_value = mean_df$m_value,
  Method = mean_df$Method,
  Mean_AUROC = mean_df$Mean_AUROC,
  SD_AUROC = sd_df$SD_AUROC # melt() preserves order -> correctly merged
)

# Ensure the below are factors 
errorbar_df$m_value_factor = factor(errorbar_df$m_value, levels = m_list_char)
errorbar_df$Method = factor(errorbar_df$Method, levels = method_names)

# - Actual plot

# For jitter on horizontal axis
pd = position_dodge(width = 0.1)

p_errorbar = ggplot(errorbar_df, aes(x = m_value_factor, y = Mean_AUROC, group = Method, color = Method)) +
  geom_line(linewidth = 0.8, position = pd) +
  geom_point(size = 3, position = pd) +
  geom_errorbar(aes(ymin = Mean_AUROC - SD_AUROC, 
                    ymax = Mean_AUROC + SD_AUROC), 
                width = 0.1, linewidth = 1, 
                position = pd) + 
  labs(title = "Errorbar plot of mean and std versus 'm'\n for each method",
       y = "AUROC",
       x = "No. sampl. 'm'",
       color = "Method") +
  scale_y_continuous(limits = c(0.5, 1)) +
  scale_color_brewer(palette = "Set1") +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    legend.position = "bottom"
  )

print(p_errorbar)

# --- 3. t-test confidence intervals ---

# - Prepare data frame for plotting

Mean_matrix = mean_std_table[,,1]
T_matrix = mean_std_table[,,3]

# Ensure column names for melting are the method names
colnames(Mean_matrix) = method_names
colnames(T_matrix) = method_names

# Add row names (m_list values) to ensure they are captured during melting
rownames(Mean_matrix) = m_list_char
rownames(T_matrix) = m_list_char

mean_df = melt(Mean_matrix, varnames = c("m_value", "Method"), value.name = "Mean_AUROC")
t_df = melt(T_matrix, varnames = c("m_value", "Method"), value.name = "T_AUROC")

# Combine all into the final data frame
errorbar_df = data.frame(
  m_value = mean_df$m_value,
  Method = mean_df$Method,
  Mean_AUROC = mean_df$Mean_AUROC,
  T_AUROC = t_df$T_AUROC 
)

# Ensure the below are factors 
errorbar_df$m_value_factor = factor(errorbar_df$m_value, levels = m_list_char)
errorbar_df$Method = factor(errorbar_df$Method, levels = method_names)

# - Actual plot

# For jitter on horizontal axis
pd = position_dodge(width = 0.2)

p_CI = ggplot(errorbar_df, aes(x = m_value_factor, y = Mean_AUROC, group = Method, color = Method)) +
  geom_line(linewidth = 0.8, position = pd) +
  geom_point(size = 3, position = pd) +
  geom_errorbar(aes(ymin = Mean_AUROC - T_AUROC, 
                    ymax = Mean_AUROC + T_AUROC), 
                width = 0.1, linewidth = 1, 
                position = pd) + 
  labs(title = "t-test CI plot of AUROC values versus 'm'\n for each method",
       y = "AUROC",
       x = "No. sampl. 'm'",
       color = "Method") +
  scale_y_continuous(limits = c(0.5, 1)) +
  scale_color_brewer(palette = "Set1") +
  theme_bw() +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
    legend.position = "bottom"
  )

print(p_CI)

# --- 4. Statistical testing ---

# - Perform pairwise comparisons of all method's AUROC using compare_AUROC()

N_m = dim(AUROC_matrix)[2]
N_meth = dim(AUROC_matrix)[3]

# Create a data frame of all category coordinates (i, j)
category_coords = expand.grid(m_idx = 1:N_m, meth_idx = 1:N_meth)
total_categories = nrow(category_coords) # N_m * N_meth

# Identify all unordered pairs of category indices
pair_indices = combn(total_categories, 2)
num_comparisons = ncol(pair_indices)

# Compare all categories using compare_AUROC() and apply()
comparison_results = apply(pair_indices, 2, function(col_indices){
  # Unpack the relevant indices
  cat_1 = unlist(category_coords[col_indices[1], ]) 
  cat_2 = unlist(category_coords[col_indices[2], ])
  
  result = compare_AUROC(cat_1, cat_2, AUROC_matrix, mean_std_table)
  
  return(result)
})

# - Display the results in a more comfortable format

comparison_df = data.frame(Result = comparison_results)

comparison_df$m_idx_1 = category_coords[pair_indices[1,], "m_idx"]
comparison_df$meth_idx_1 = category_coords[pair_indices[1,], "meth_idx"]
comparison_df$m_idx_2 = category_coords[pair_indices[2,], "m_idx"]
comparison_df$meth_idx_2 = category_coords[pair_indices[2,], "meth_idx"]

comparison_df$m_1 = m_list[comparison_df$m_idx_1]
comparison_df$meth_1 = method_names[comparison_df$meth_idx_1]
comparison_df$m_2 = m_list[comparison_df$m_idx_2]
comparison_df$meth_2 = method_names[comparison_df$meth_idx_2]

# All entries
comparison_df

# Some entries

avoid_displaying = c("m_idx_1", "m_idx_2", "meth_idx_1", "meth_idx_2")

mask = (comparison_df["m_idx_1"] == comparison_df["m_idx_2"]) & (comparison_df["m_2"] == 20)  
#mask = (comparison_df["meth_idx_1"] == comparison_df["meth_idx_2"]) & (comparison_df["meth_1"] == "SA (exp.2)")

comparison_df[mask, !names(comparison_df) %in% avoid_displaying]

# - Interpret the output

# If we use the p-values for statistical tests, need to put an appropriate 
# level. Implement the Bonferroni correction (conservative)

alpha_bonferroni = alpha / ( num_comparisons )

print(paste("Number of tests performed: ", num_comparisons))
print(paste("Significance (via Bonferroni) to prevent p-hacking: ", alpha_bonferroni))