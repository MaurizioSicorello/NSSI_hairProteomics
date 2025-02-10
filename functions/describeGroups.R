
# helper function which returns means, SDs and a t-test for var by groups
describeGroups_help <- function(var, groups){
  
  l <- levels(as.factor(groups))
  d <- describeBy(var, group = groups, mat = T)
  
  
  # Descriptives for group 1
  g1 = l[1]
  m1 = round(d[1, 5], 2)
  sd1 = round(d[1, 6], 2)
  n1 = length(var[groups == l[1]])
  
  cat("Group: ", g1, "\n Mean: ", m1, "\n SD: ", sd1)
  cat("\n ----------------- \n")
  
  
  # Descriptives for group 2
  g2 = l[2]
  m2 = round(d[2, 5], 2)
  sd2 = round(d[2, 6], 2)
  n2 = length(var[groups == l[2]])
  
  cat("Group: ", g2, "\n Mean: ", m2, "\n SD: ", sd2)
  cat("\n ----------------- \n")
  
  
  # T-test
  t = t.test(var[groups == l[1]], var[groups == l[2]])
  s = round(t$statistic, 2)
  df = round(t$parameter, 2)
  p = round(t$p.value, 3)
  
  if(p <= 0.001){sign = "***"}else{
    if(p <= 0.01){sign = "**"}else{
      if(p <= 0.05){sign = "*"}else{
        if(p > 0.05){sign = "n.s."}
      }
    }
  }

  SD_pooled = sqrt((n1*sd1^2 + n2*sd2^2)/(n1+n2))
  cohensd = round((m1-m2) / SD_pooled, 2)


  cat("t-value: ", s, "\n df: ", df, "\n p-value: ", p, sign, "\n Cohens d: ", cohensd, "\n\n\n")

  # Save results as a formatted row vector
  result_row <- c(
    Variable = deparse(substitute(var)),
    `Clinical Group (M ± SD)` = sprintf("%.2f ± %.2f", m1, sd1),
    `Control Group (M ± SD)` = sprintf("%.2f ± %.2f", m2, sd2),
    `t(df)` = sprintf("%.2f(%d)", s, as.integer(df)),
    `p-value` = sprintf("%.3f", p),
    `Cohen’s d` = sprintf("%.2f", cohensd)
  )

  return(result_row)


}




# # Applies helper function to a list of variables (listVars) in a dataframe (data)
# describeGroups <- function(data,  groups, listVars = NULL){
#   
#   for(i in 1:length(listVars)){
#     
#     cat("-----------------", listVars[i], "-----------------\n\n")
#     describeGroups_help(var = data[, listVars[i]], groups = data[, groups])
#     
#     
#   }
# }


describeGroups <- function(data, groups, listVars = NULL) {
  
  results_list <- list()  # Initialize an empty list to store results
  
  for (i in seq_along(listVars)) {
    
    cat("-----------------", listVars[i], "-----------------\n\n")
    
    # Call the helper function and store results
    completeInd <- !is.na(data[[listVars[i]]])
    result_row <- describeGroups_help(var = data[[listVars[i]]][completeInd], groups = data[[groups]][completeInd])
    
    # Append result as a row to the results list
    results_list[[i]] <- result_row
  }
  
  # Combine list into a dataframe
  results_df <- do.call(rbind, results_list) |> as.data.frame()
  
  # Return results as a dataframe
  return(results_df)
}




formatMeanSD <- function(var) {
  
  # Compute mean and SD, rounding to 2 decimal places
  m <- round(mean(var, na.rm = TRUE), 2)
  sd <- round(sd(var, na.rm = TRUE), 2)
  
  # Format as "M ± SD"
  mean_sd_str <- sprintf("%.2f ± %.2f", m, sd)
  
  # Return a formatted row with placeholders for other columns
  result_row <- c(
    Variable = deparse(substitute(var)),
    `Clinical Group (M ± SD)` = mean_sd_str,
    `Control Group (M ± SD)` = "",
    `t(df)` = "",
    `p-value` = "",
    `Cohen’s d` = ""
  )
  
  return(result_row)
}

