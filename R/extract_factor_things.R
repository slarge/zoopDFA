#extract factor loadings/spatial predictions

library(TMB)
library(VAST)
library(FishStatsUtils)
library(dplyr)
library(tidyverse)
library(ggrepel)

#load fit object


Report = fit$Report
ParHat = fit$ParHat
Data = fit$data_list
Obj = fit$tmb_list$Obj
SD = fit$parameter_estimates$SD
year_labels = fit$year_labels
category_names = spp_list
RotationMethod = "Varimax"
mapdetails_list = mapdetails_list=results$map_list
Dim_year = NULL
Dim_species = NULL
projargs = fit$extrapolation_list$projargs
plotdir=paste0(moddir,"/")
land_color = "grey"
zlim = NA
testcutoff = 0.1



# Extract Options and Options_vec (depends upon version)
if( "Options_list" %in% names(Data) ){
  Options_vec = Data$Options_list$Options_vec
  Options = Data$Options_list$Options
}


# Dimensions for plotting
Dim = function( num ) c(ceiling(sqrt(num)), ceiling(num/ceiling(sqrt(num))) )
Dim_year = Dim(length(year_labels))
Dim_species = Dim(length(category_names))

# Extract loadings matrices (more numerically stable than extracting covariances, and then re-creating Cholesky)
Psi2prime_list = Psiprime_list = Psi2prime_SE_list = Lprime_SE_list = Hinv_list = Lprime_list = L_list = vector("list", length=1)    # Add names at end so that NULL doesn't interfere

# Loop through
i<-1
  
  # Variable names
  Par_name = "EpsilonTime1"
  Lpar_name = "Ltime_epsilon1_z"
  
  # Backwards compatible loading of variables and names
  Var_name = "Epsiloninput1_sff" 
  Var2_name = "Epsiloninput1_gff"
  L_name = "Ltime_epsilon1_tf" 

  # Continue if component is included

    # Get loadings matrix
      L_list[[i]] = Report[[L_name]]

    # Get covariance
     Psi_sjt = ParHat[[Var_name]]
     Psi_gjt = Report[[Var2_name]]
    ## the betas and EpsilonTime are transposed compared to others so fix that here
      Psi_sjt = aperm( Psi_sjt, c(1,3,2) )
      Psi_gjt = aperm( Psi_gjt, c(1,3,2) )
    
    
    logkappa = unlist(ParHat['logkappa1'])
      tau = 1 / (exp(logkappa) * sqrt(4*pi));

    # Rotate stuff
    Var_rot = rotate_factors( L_pj=L_list[[i]], Psi=Psi_sjt/tau, RotationMethod=RotationMethod, testcutoff=testcutoff, quiet=TRUE )
    Var_rot$Psi_rot = aperm( Var_rot$Psi_rot, c(1,3,2) )
  
    Report_tmp = list("D_xct"=Var_rot$Psi_rot, "Epsilon1_sct"=Var_rot$Psi_rot)
    Lprime_list[[i]] = Var_rot$L_pj_rot
    Psiprime_list[[i]] = Var_rot$Psi_rot
    Hinv_list[[i]] = Var_rot$Hinv
    

    # Extract SEs if available
    # Could also edit to extract L_SE_list and Psi2_SE_list
    if( !missing(Obj) && class(SD)=="sdreport" ){
      L_cfr = sample_variable( Sdreport=SD, Obj=Obj, variable_name=L_name, n_samples=100, sample_fixed=TRUE, seed=123456 )
      Psi_gjtr = sample_variable( Sdreport=SD, Obj=Obj, variable_name=Var2_name, n_samples=100, sample_fixed=TRUE, seed=123456 )
      if( Par_name %in% c("EpsilonTime1","EpsilonTime2") ){
        Psi_gjtr = aperm( Psi_gjtr, c(1,3,2,4) )
      }
      Lprime_cfr = array(NA, dim=dim(L_cfr) )
      Psiprime_gjtr = array(NA, dim=dim(Psi_gjtr) )
      for( rI in 1:100 ){
        tmplist = rotate_factors( L_pj=array(L_cfr[,,rI],dim=dim(L_cfr)[1:2]),
                                  Psi_sjt=array(Psi_gjtr[,,,rI],dim=dim(Psi_gjtr)[1:3])/tau, RotationMethod=RotationMethod, quiet=TRUE,testcutoff = testcutoff )
        Lprime_cfr[,,rI] = tmplist$L_pj_rot
        Psiprime_gjtr[,,,rI] = tmplist$Psi_rot
      }
      if( Par_name %in% c("EpsilonTime1","EpsilonTime2") ){
        Psiprime_gjtr = aperm( Psiprime_gjtr, c(1,3,2,4) )
      }
      Lprime_SE_list[[i]] = apply(Lprime_cfr, MARGIN=1:2, FUN=sd)
      if( nrow(Lprime_SE_list[[i]])==length(category_names) ){
        rownames(Lprime_SE_list[[i]]) = category_names
      }
      Psi2prime_SE_list[[i]] = apply(Psiprime_gjtr, MARGIN=1:3, FUN=sd)
    }
    
    # Extract projected factors is available
      Var2_rot = rotate_factors( L_pj=L_list[[i]], Psi=Psi_gjt/tau, RotationMethod=RotationMethod, testcutoff=testcutoff, quiet=TRUE )
        Var2_rot$Psi_rot = aperm( Var2_rot$Psi_rot, c(1,3,2) )
      
      Report2_tmp = list("D_xct"=Var2_rot$Psi_rot, "Epsilon1_sct"=Var2_rot$Psi_rot)
      Psi2prime_list[[i]] = Var2_rot$Psi_rot
    

      dir<-"~/Desktop/indicators/indicators_working_output/zoop_output"
      
      write.table(Lprime_list, file = paste0(dir, "/eps_t_loadings.txt"), row.names = F)
      write.table(Lprime_SE_list, file = paste0(dir, "/eps_t_SE.txt"), row.names = F)
      
            saveRDS(Report2_tmp, file = paste0(dir, "/spat_eps_t.txt") )
      