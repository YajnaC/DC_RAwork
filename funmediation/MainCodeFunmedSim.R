library(mgcv);
library(refund);
library(boot);
source("simulate_functional_mediation_example.R");
source("funreg_mediation.R");
source("tvem.R");
source("select_tvem.R");
source("plot_funreg_mediation.R");
source("print_funreg_mediation.R");
source("plot_tvem.R");
source("print_tvem.R");
set.seed(the_seed);
nboot <- 199;
answers <- NULL;
start_time <- Sys.time(); 
for (this_sim in 1:nsim) { 
  simulated_data <- simulate_functional_mediation_example(nsub=nsub,
                                                          simulate_binary_Y=TRUE);
  true_beta <- simulated_data$true_beta;
  true_gamma_int <- simulated_data$true_gamma_int;
  true_gamma_X <- simulated_data$true_gamma_X;
  true_alpha_int <- simulated_data$true_alpha_int;
  true_alpha_X <- simulated_data$true_alpha_X;
  true_alpha_M <- simulated_data$true_alpha_M;
  long_data <- simulated_data$dataset;
  fun_med_results <- funreg_mediation(data=long_data,
                                      treatment=X,
                                      mediator=M,
                                      outcome=Y,
                                      id=subject_id,
                                      time=t,
                                      logistic=TRUE,
                                      nboot=nboot);
  
  tvem_model_summary <- summary(fun_med_results$original_results$tvem_XM_details$back_end_model);
  gamma_int_residuals <- fun_med_results$original_results$gamma_int_estimate - true_gamma_int;
  gamma_int_se <- fun_med_results$original_results$gamma_int_se;
  gamma_X_residuals <- fun_med_results$original_results$gamma_X_estimate - true_gamma_X;
  gamma_X_se <- fun_med_results$original_results$gamma_X_se;
  alpha_M_residuals <- fun_med_results$original_results$alpha_M_estimate - true_alpha_M;
  alpha_M_se <- fun_med_results$original_results$alpha_M_se;
  norm_cover <- (fun_med_results$bootstrap_results$beta_boot_norm_lower < true_beta) & 
    (fun_med_results$bootstrap_results$beta_boot_norm_upper > true_beta);
  basic_cover <- (fun_med_results$bootstrap_results$beta_boot_basic_lower < true_beta) & 
    (fun_med_results$bootstrap_results$beta_boot_basic_upper > true_beta);
  perc_cover <- (fun_med_results$bootstrap_results$beta_boot_perc_lower < true_beta) & 
    (fun_med_results$bootstrap_results$beta_boot_perc_upper > true_beta);
  norm_power <- (sign(fun_med_results$bootstrap_results$beta_boot_norm_lower)==
                   sign(fun_med_results$bootstrap_results$beta_boot_norm_upper));
  basic_power <- (sign(fun_med_results$bootstrap_results$beta_boot_basic_lower)==
                    sign(fun_med_results$bootstrap_results$beta_boot_basic_upper));
  perc_power <- (sign(fun_med_results$bootstrap_results$beta_boot_perc_lower)==
                   sign(fun_med_results$bootstrap_results$beta_boot_perc_upper));
  current_time <- Sys.time();
  answers <- rbind(answers, 
                   # Create structures to hold results:  
                   # ... for the effect of X on M: 
                   c(nsub=nsub, 
                     this_sim=this_sim,
                     gamma_int_estimate_bias =  mean(gamma_int_residuals), 
                     gamma_int_estimate_mse =  mean(gamma_int_residuals^2),
                     gamma_int_estimate_mean_std_err = mean(gamma_int_se),
                     gamma_int_estimate_coverage =  mean(abs(gamma_int_residuals)<1.96*gamma_int_se), 
                     gamma_int_estimate_family_coverage =  all(abs(gamma_int_residuals)<1.96*gamma_int_se), 
                     gamma_int_estimate_pvalue_overall =  as.numeric(tvem_model_summary$p.coeff["(Intercept)"]), 
                     gamma_int_estimate_pvalue_varying =  tvem_model_summary$s.table["s(time)","p-value"],
                     gamma_X_estimate_bias =  mean(gamma_X_residuals), 
                     gamma_X_estimate_mse =  mean(gamma_X_residuals^2),
                     gamma_X_estimate_mean_std_err = mean(gamma_X_se),
                     gamma_X_estimate_coverage =  mean(abs(gamma_X_residuals)<1.96*gamma_X_se), 
                     gamma_X_estimate_family_coverage =  all(abs(gamma_X_residuals)<1.96*gamma_X_se), 
                     gamma_X_estimate_pvalue_overall =  as.numeric(tvem_model_summary$p.coeff["(Intercept)"]),   
                     gamma_X_estimate_pvalue_varying =  tvem_model_summary$s.table["s(time):treatment","p-value"],  
                     # ... for the joint effect of X and M on Y: 
                     alpha_int_bias =  fun_med_results$original_results$alpha_int_estimate - true_alpha_int,
                     alpha_int_mse = (fun_med_results$original_results$alpha_int_estimate - true_alpha_int)^2,
                     alpha_int_se =  fun_med_results$original_results$alpha_int_se,
                     alpha_int_coverage = abs(fun_med_results$original_results$alpha_int_estimate - true_alpha_int) < 1.96*fun_med_results$original_results$alpha_int_se,
                     alpha_X_bias =  fun_med_results$original_results$alpha_X_estimate - true_alpha_X,
                     alpha_X_mse = (fun_med_results$original_results$alpha_X_estimate - true_alpha_X)^2,
                     alpha_X_se =  fun_med_results$original_results$alpha_X_se,
                     alpha_X_coverage = abs(fun_med_results$original_results$alpha_X_estimate - true_alpha_X) < 1.96*fun_med_results$original_results$alpha_X_se,
                     alpha_M_estimate_bias =  mean(alpha_M_residuals), 
                     alpha_M_estimate_mse =  mean(alpha_M_residuals^2), 
                     alpha_M_estimate_coverage =  mean(abs(alpha_M_residuals)<1.96*alpha_M_se), 
                     alpha_M_estimate_family_coverage =  all(abs(alpha_M_residuals)<1.96*alpha_M_se), 
                     alpha_M_pvalue =   fun_med_results$original_results$alpha_M_pvalue,  
                     # ... for the direct effect of X on Y: 
                     delta_int_estimate =  fun_med_results$original_results$delta_int_estimate, 
                     delta_int_se =  fun_med_results$original_results$delta_int_se, 
                     delta_X_estimate =  fun_med_results$original_results$delta_X_estimate, 
                     delta_X_se =  fun_med_results$original_results$delta_X_se,  
                     # ... for the mediation effect of X on Y through M:  
                     beta_bias = fun_med_results$bootstrap_results$beta_boot_estimate - true_beta,
                     beta_mse = (fun_med_results$bootstrap_results$beta_boot_estimate - true_beta)^2,
                     beta_boot_se = fun_med_results$bootstrap_results$beta_boot_se, 
                     beta_boot_norm_lower = fun_med_results$bootstrap_results$beta_boot_norm_lower, 
                     beta_boot_norm_upper = fun_med_results$bootstrap_results$beta_boot_norm_upper, 
                     beta_boot_norm_coverage = norm_cover,
                     beta_boot_norm_power = norm_power,
                     beta_boot_basic_lower = fun_med_results$bootstrap_results$beta_boot_basic_lower, 
                     beta_boot_basic_upper = fun_med_results$bootstrap_results$beta_boot_basic_upper, 
                     beta_boot_basic_coverage = basic_cover, 
                     beta_boot_basic_power = basic_power,
                     beta_boot_perc_lower = fun_med_results$bootstrap_results$beta_boot_perc_lower,
                     beta_boot_perc_upper = fun_med_results$bootstrap_results$beta_boot_perc_upper,
                     beta_boot_perc_coverage = perc_cover,
                     beta_boot_perc_power = perc_power
                   ));
  save.image("working.rdata");
}
finish_time <- Sys.time();
print(c(nsim, nboot));
time_required <- difftime(finish_time,start_time);
print(time_required);
mean_answers <- apply(answers,2,mean);
print(round(mean_answers,3));
save.image(paste("answers",nsub,"_",the_seed,".rdata",sep=""));

