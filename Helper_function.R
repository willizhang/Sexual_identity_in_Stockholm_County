# title: "Helper function"
# author: Guoqiang Zhang
# email: guoqiang.zhang@ki.se


##### Function for Weighted Regression Models #####

### crude analysis ###
fit_model_crude <- function( formula, design ) {
  
  # Poisson regression
 tryCatch(
    svyglm( formula, design = design, family = quasipoisson( link = "log" ) ),
    error = function( e ) {
      message( e$message )
      return ( NULL )
    }
  )
}


### adjusted analysis ###
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

### crude analysis ###
extract_results_coef_crude <- function( all_models_list, exposure, outcome, variables ) {
  
  model <- all_models_list[[ exposure ]][[ outcome ]]
  
  # check if model is available (i.e., converged)
  if ( !is.null( model ) ) {
    
    results <- cbind( coef( model )[ variables ], 
                      confint( model, ddf = degf( model$survey.design ) )[ variables, , drop = FALSE ] )
    
    return( results )
  } else {
    return( NULL )
  }
}


### adjusted analysis ###
extract_results_coef <- function( all_models_list, exposure, outcome, model_type, variables ) {
  
  model <- all_models_list[[ exposure ]][[ outcome ]][[ model_type ]]
  
  # check if model is available (i.e., converged)
  if ( !is.null( model ) ) {
    
    results <- cbind( coef( model )[ variables ], 
                      confint( model, ddf = degf( model$survey.design ) )[ variables, , drop = FALSE ] )
    
    return( results )
    } else {
    return( NULL )
    }
  }


##### Function to Extract Proportion Ratios and Confidence Intervals for Each Model Across Exposures #####

### crude analysis ###
extract_results_pr_crude <- function( model_base, identities, exposure ) {
  results_list <- list()
  
  for ( identity in identities ) {
    
    model_name <- paste0( "fml_", identity, "_", exposure, "_crude" )
    model <- model_base[[ model_name ]]
    
    if ( !is.null( model ) ) {
      exp_model <- exp( model )
      temp_df <- as.data.frame( exp_model )
      
      colnames( temp_df ) <- c( paste0( identity, "_point_estimate_crude" ),
                                paste0( identity, "_lower_ci_crude" ),
                                paste0( identity, "_upper_ci_crude" ) )
      
      results_list[[ identity ]] <- temp_df
    }
  }
  
  names( results_list ) <- NULL
  
  combined_df <- tibble::rownames_to_column(
    do.call( cbind, results_list ),
    var = "Exposure" )
  
  return( combined_df )
}


### adjusted analysis ###
extract_results_pr <- function( model_base, model_type, identities, exposure ) {
  results_list <- list()
  
  for ( identity in identities ) {

    model_name <- paste0( "fml_", identity, "_", exposure, "_", model_type )
    model <- model_base[[ model_name ]]
    
    if ( !is.null( model ) ) {
      exp_model <- exp( model )
      temp_df <- as.data.frame( exp_model )
      
      colnames( temp_df ) <- c( paste0( identity, "_point_estimate" ),
                                paste0( identity, "_lower_ci" ),
                                paste0( identity, "_upper_ci" ) )
      
      results_list[[ identity ]] <- temp_df
    }
  }
  
  names( results_list ) <- NULL
  
  combined_df <- tibble::rownames_to_column(
    do.call( cbind, results_list ),
    var = "Exposure" )
  
  return( combined_df )
}
