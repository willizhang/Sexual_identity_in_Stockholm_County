# title: "Helper function"
# author: Guoqiang Zhang
# email: guoqiang.zhang@ki.se


##### Function for Weighted Regression Models #####

fit_models <- function( formula, design ) {
  
  models <- list()
  
  # Poisson regression
  models$poisson <- tryCatch(
    svyglm( formula, design = design, family = quasipoisson( link = "log" ) ),
    error = function( e ) {
      message( e$message )
      return ( NULL )
      }
  )
  
  # log-binomial regression
  models$log_binomial <- tryCatch(
    svyglm( formula, design = design, family = quasibinomial( link = "log" ) ),
    error = function( e ) {
      message( e$message )
      return ( NULL )
    }
  )
  
  # logistic regression
  models$logistic <- tryCatch(
    svyglm( formula, design = design, family = quasibinomial( link = "logit" ) ),
    error = function( e ) {
      message( e$message )
      return ( NULL )
    }
  )
  
  return( models )
}


##### Function to Extract Coefficients and Confidence Intervals (on log scale) #####

extract_results_coef <- function( all_models_list, exposure, outcome, model_type, variables, design ) {
  
  model <- all_models_list[[ exposure ]][[ outcome ]][[ model_type ]]
  
  # check if model is available (i.e., converged)
  if ( !is.null( model ) ) {
    
    results <- cbind( coef( model )[ variables ], 
                      confint( model, ddf = degf( design ) )[ variables, , drop = FALSE ] )
    
    return( results )
    } else {
    return( NULL )
    }
  }


##### Function to Extract Proportion Ratios and Confidence Intervals for Each Model Across Exposures #####

extract_results_pr <- function( model_base, model_type, identities, exposure ) {
  results_list <- list()
  
  for (identity in identities) {

    model_name <- paste0( "fml_", identity, "_", exposure, "_", model_type )
    model <- model_base[[ model_name ]]
    
    if ( !is.null( model ) ) {
      exp_model <- exp( model )
      temp_df <- as.data.frame( exp_model )
      
      colnames( temp_df ) <- c( paste0( "Point_estimate_", identity, "_2021" ),
                                paste0( "Lower_ci_", identity, "_2021" ),
                                paste0( "Upper_ci_", identity, "_2021" ) )
      
      results_list[[ identity ]] <- temp_df
    }
  }
  
  names( results_list ) <- NULL
  
  combined_df <- do.call( cbind, results_list )
  
  return( combined_df )
}
