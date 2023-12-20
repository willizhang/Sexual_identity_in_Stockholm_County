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


### sexual identity fluidity ###
# for data SPHC-F 2010-2021
extract_results_coef_fluidity <- function( all_models_list, exposure, model_type, variables ) {
  
  model <- all_models_list[[ exposure ]][[ model_type ]]
  
  # check if model is available
  if ( !is.null( model ) ) {
    
    results <- cbind( coef( model )[ variables ], 
                      coefci( model, vcov = sandwich )[ variables, , drop = FALSE ] )
    
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



##### Function to Extract Risk Ratios and Confidence Intervals for Each Model Across Exposures #####

extract_rr_fluidity <- function( model_results ) {
  
  results_list <- list()
  
  for ( model_name in names( model_results ) ) {
    
    exp_results <- exp( model_results[[ model_name ]] )
    
    temp_df <- as.data.frame( exp_results )
    
    type <- case_when(
      grepl( "crude", model_name )   ~ "crude",
      grepl( "model_1", model_name ) ~ "model_1",
      grepl( "model_2", model_name ) ~ "model_2" )
    
    colnames( temp_df ) <- c( paste0( "point_estimate_", type ),
                              paste0( "lower_ci_", type ),
                              paste0( "upper_ci_", type ) )
    
    results_list[[ model_name ]] <- temp_df
  }
  
  names( results_list ) <- NULL
  
  combined_df <- tibble::rownames_to_column(
    do.call( cbind, results_list ),
    var = "Exposure" )
  
  return( combined_df )
}