library(tidyverse)
library(rEDM)

### Analysis wrappers

make_embedded_dfr <- function(time_series,lib,E,tp=1){
  map_dfr(1:NROW(lib),function(i_lib){
    map_dfc(0:(E-1),function(j_lag){
      x_j_lag <- data.frame(x=lag(time_series[lib[i_lib,1]:lib[i_lib,2]],j_lag))
      names(x_j_lag) <- paste0("x_minus_",j_lag)
      return(x_j_lag)
    }) %>%
      bind_cols(data.frame(x_Tp=lead(time_series[lib[i_lib,1]:lib[i_lib,2]],tp)))
  })
}

do_ccm_on_embedded <- function(block_embedded_ccm,lib_sizes,num_samples,RNGseed=8675309){
  
  set.seed(RNGseed)
  lib_sizes <- lib_sizes[lib_sizes <= NROW(block_embedded_ccm)]
  vars <- names(block_embedded_ccm[-1])
  E <- length(vars)-1
  
  map_dfr(lib_sizes,function(L_i){
    map_dfr(1:num_samples,function(j){
      
      I_ij <- sample(1:NROW(block_embedded_ccm),L_i,replace=F)
      
      block_embedded_ccm_ij <- block_embedded_ccm[I_ij,]
      # lib_ij <- paste(rep(lib_ij,each=2))
      
      ccm_out <- Simplex(dataFrame=block_embedded_ccm_ij,lib=paste("1",L_i),pred=paste("1",L_i),
                         E=E,
                         target=vars[E+1],Tp=0,
                         embedded=TRUE,
                         columns = paste(vars[1:E],collapse=" "))
      
      
      # stats <- ComputeError(model_out$Predictions,model_out$Observations)
      # names(stats) <- c("mae","rho","rmse")
      
      stats <- compute_stats(ccm_out$Predictions,ccm_out$Observations)
      
      return(c(lib_size=L_i,index=j,stats))
    })
  })
  
}




make_embedded_ccm_dfr <- function(x,y,lib,E,tp){
  map_dfr(1:NROW(lib),function(i_lib){
    map_dfc(0:(E-1),function(j_lag){
      x_j_lag <- data.frame(x=lag(x[lib[i_lib,1]:lib[i_lib,2]],j_lag))
      names(x_j_lag) <- paste0("x_minus_",j_lag)
      return(x_j_lag)
    }) %>%
      bind_cols(data.frame(y_Tp=leadlag(y[lib[i_lib,1]:lib[i_lib,2]],tp)))
  })
}

## Do pairwise CCM for all columns in a block (originally of gene expression data)
compute_ccm_on_block <- function( block,
                                  target,
                                  predictor,
                                  lib = c(1,NROW(block)),
                                  lib_sizes = seq(20,100,by=10),
                                  num_samples = 100,
                                  max_E = 8,
                                  results_file = NULL)
{
  lib <- matrix(lib,ncol=2)
  
    ## Inside here is really what you need
    
    ccm_tp_1 <- map_df(1:8, function(E_j) {
      
      block_ij <- make_embedded_ccm_dfr(block[[predictor]],block[[target]],tp=-1,E=E_j,lib=lib)
      block_ij <- block_ij %>%
        filter(complete.cases(.)) %>%
        mutate(time=row_number()) %>%
        select(time,everything())
      
      vars_E_j <- names(block_ij)[-1]
      n_j <- NROW(block_ij)
      
      ccm_out <- Simplex(dataFrame=block_ij,lib=paste("1",n_j),pred=paste("1",n_j),
                         E=E_j,
                         target=vars_E_j[E_j+1],Tp=0,
                         embedded=TRUE,
                         columns = paste(vars_E_j[1:E_j],collapse=" "))
      stats <- compute_stats(ccm_out$Predictions,ccm_out$Observations)
      
      return(c(E=E_j,stats))
    })
    
    best_E <- ccm_tp_1$E[which.max(ccm_tp_1$rho)]
    
    if(length(best_E) != 1) # compute CCM using tp = 0
    { # if invalid, return NA results
      ccm_tp_0 <- data.frame(E = NA,
                             # tau = NA, tp = NA,num_neighbors = NA, 
                             num_pred = NA,
                             lib_column = predictor, target_column = target,
                             # lib_size = NA, num_pred = NA, 
                             rho = NA, mae = NA, rmse = NA,
                             p_val = NA)
      return(ccm_tp_0)
    }
    
    block_i_E_star <- make_embedded_ccm_dfr(block[[predictor]],block[[target]],tp=0,E=best_E,lib=lib)
    block_i_E_star <- block_i_E_star %>%
      filter(complete.cases(.)) %>%
      mutate(time=row_number()) %>%
      select(time,everything())
    
    vars_E_star <- names(block_i_E_star)[-1]
    n_E_star <- NROW(block_i_E_star)
    
    ccm_out <- do_ccm_on_embedded(block_i_E_star,lib_sizes,num_samples)
    
    # ccm_out <- Simplex(dataFrame=block_i_E_star,lib=paste("1",n_E_star),pred=paste("1",n_E_star),
    #                    E=best_E,
    #                    target=vars_E_star[best_E+1],Tp=0,
    #                    embedded=TRUE,
    #                    columns = paste(vars_E_star[1:best_E],collapse=" "))
    # stats <- compute_stats(ccm_out$Predictions,ccm_out$Observations)
    
    ccm_out <- ccm_out %>%
      mutate(E = best_E,lib_column = predictor, target_column = target) %>%
      select(lib_column,target_column,E,everything())

    return(ccm_out)
    
  }
