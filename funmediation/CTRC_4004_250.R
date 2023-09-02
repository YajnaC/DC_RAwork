library(refund)
seed_values = 4
nsub_values = c(250)
for (seed in 1001*seed_values) {
  for (nsub in nsub_values) {
    file_name = paste("run_sim_",seed,"_",nsub,".R",sep="")
    source(file_name)
    
    results <- list(answers=answers,
                    fun_med_results=fun_med_results,
                    long_data=long_data, 
                    simulated_data=simulated_data,
                    tvem_model_summary=tvem_model_summary)
    
    results_seed_nsub = paste("results_",seed,"_",nsub,".rda",sep = "")
    save(results, file = results_seed_nsub)
    
    write(paste("the_seed =",seed,";"),file=file_name,append=FALSE)
    write(paste("nsub =",nsub,";"),file=file_name,append=TRUE)
    write("nsim = 1000;",file=file_name,append=TRUE)
    write("source('MainCodeFunmedSim.R');",file=file_name,append=TRUE)
  }
}
