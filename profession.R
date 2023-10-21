# profession
table( d_2021$SSYK_kl, useNA = "always" )
d_2021$profession <- factor( ifelse( d_2021$SSYK_kl == "Yrken inom byggverksamhet och tillverkning" |
                                       d_2021$SSYK_kl == "Yrken inom lantbruk, trädgård, skogsbruk och fiske" |
                                       d_2021$SSYK_kl == "Yrken inom maskinell tillverkning och transport m.m.",
                                     "Manual and field trades",
                                     ifelse( d_2021$SSYK_kl == "Service-, omsorgs- och försäljningsyrken" |
                                               d_2021$SSYK_kl == "Yrken inom administration och kundtjänst" |
                                               d_2021$SSYK_kl == "Yrken med krav på kortare utbildning eller introduktion",
                                             "Service and support",
                                             ifelse( d_2021$SSYK_kl == "Yrken med krav på fördjupad högskolekompetens" |
                                                       d_2021$SSYK_kl == "Yrken med krav på högskolekompetens eller motsvarande" |
                                                       d_2021$SSYK_kl == "Chefsyrken",
                                                     "Expertise and leadership", NA ) 
                                     )
),
levels = c( "Manual and field trades", "Service and support", "Expertise and leadership" ) 
)
table( d_2021$profession, useNA = "always" )


##### 3.1.8. By profession
```{r}
survey_design_cc_profession <- twophase( id = list( ~ 1, ~ 1 ), 
                                         strata = list( ~ sampling_strata_region, ~ sampling_strata_region ), 
                                         method = "simple", 
                                         weights = list( ~ design_weight, NULL ), 
                                         subset = ~ !is.na( sexual_identity ) & !is.na( profession ), 
                                         data = d_2021_complete_cc 
)

survey_design_cc_profession <- survey::calibrate( survey_design_cc_profession, 
                                                  formula = ~ sampling_strata_region,
                                                  phase = 2 
)

list_of_df <- list()

for ( cat in categories ) {
  prop_cc_profession <- svyby( formula = as.formula( paste0( "~ I( sexual_identity == '", cat, "' )" ) ),
                               by = ~ profession,
                               design = survey_design_cc_profession,
                               FUN = svyciprop,
                               vartype = "ci",
                               method = "beta" )
  
  colnames( prop_cc_profession ) <- c( "subgroup", paste0( cat, "_point_estimate_2021" ), paste0( cat, "_lower_ci_2021" ), paste0( cat, "_upper_ci_2021" ) )
  
  list_of_df[[cat]] <- prop_cc_profession
}

prop_cc_profession <- Reduce( function( df1, df2 ) { 
  merge( df1, df2, by = "subgroup" )
}, 
list_of_df )

prop_cc_profession <- left_join( prop_cc_profession,
                                 d_2021_complete[ !is.na( d_2021_complete$sexual_identity ) &
                                                    !is.na( d_2021_complete$profession ), ] %>% 
                                   group_by( subgroup = profession ) %>% 
                                   summarise( sample_size_2021 = n() ),
                                 by = "subgroup" ) %>%
  mutate( sample_size_2021 = prettyNum( sample_size_2021, big.mark = ",", preserve.width = "none" ) )
```
