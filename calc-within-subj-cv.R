# Reference: https://rowannicholls.github.io/R/statistics/agreement/coefficient_of_variation.html
# The results from here match with the MedCalc program results + me following their instructions: https://www.medcalc.org/manual/cvfromduplicates.php

# Set up new function 
calculate_within_subject_cv <- function(stats, method) {
  if (method == "root_mean") {
    # Original function's method of calculating 
    n <- stats$based.on
    sum_squares <- sum((stats$diffs / stats$means)^2, na.rm = T)
    cv_ws <- sqrt(sum_squares / (2 * n))
    
    # KH way of calculating to match the naming scheme on Martin Bland's York site (sanity check)
    ## Calculate the within-subject variance for the natural scale values. 
    ## (Within-subject variance is given by difference squared over 2 when we have pairs of subjects.)
    s2 <- (stats$diffs)^2/2
    ## Calculate subject mean and s squared / mean squared, i.e. CV squared.
    m <- stats$means
    s2m2 <- s2/m^2
    ## Calculate mean of s squared / mean squared.
    mean_s2m2 <- mean(s2m2)
    ## The within-subject CV is the square root of the mean of s squared / mean squared:
    cv_ws_KH <- sqrt(mean_s2m2)
    
    # Check if their calc matches my calc
    if (cv_ws == cv_ws_KH) {
      print("The within-subject CV estimated by KH and the original function match.")
    }
    
    # Calc CIs according to Martin Bland York site 
    ## For the root mean square method, this is very direct. 
    ## We have the mean of the squared CV == mean_s2m2 (i.e., the term right before you take the square root)... 
    ## ...so we use the usual confidence interval for a mean on this then take the square root of the sample size.
    
    ## i.e., the standard error is the standard deviation of the CVs divided by the square root of the sample size
    SE_term <- sd(s2m2)/sqrt(n)
    
    ## Then, the 95% confidence interval for the squared CV can be found by the mean minus or plus 1.96 standard errors. 
    ## HOWEVER, if the sample is small we should use the t distribution INSTEAD of the z-distribution. 
    ## (Bland only used 1.96 because he was using a z-distribution because his simulation was with 100 subjects)
    ## (However, Bland also notes: the squared CVs are unlikely to be Normal, so the CI will still be very approximate.)
    ### Find the t_critical value for a 95% confidence interval using the t-distribution
    ### 95% CI: 1-(1-.95)/2 = 0.975
    #crit_val = qt(0.975, (n-1))
    ### 90% CI: 1-(1-.90)/2 = 0.95
    crit_val = qt(0.95, (n-1))
    ### 83.4% CI: 1-(1-.834)/2 = 0.917
    #crit_val = qt(0.917, (n-1))
    
    ### THEN: do mean +/- t_critical*SE  
    low_bound_sq <- mean_s2m2 - crit_val*SE_term
    high_bound_sq <- mean_s2m2 + crit_val*SE_term
    
    ## The square roots of these limits give the XX% confidence interval for the CV.
    lower_CI = sqrt(low_bound_sq)
    upper_CI = sqrt(high_bound_sq)
    
    # Return CV calc and CIs 
    results_list <- list("Root_Mean_CV" = cv_ws*100, "Lower_CI" = lower_CI*100, "Upper_CI" = upper_CI*100);
    return(results_list);
    
  } else if (method == "logarithmic") {
    n <- stats$based.on
    sl <- sum((log(stats$group$group1) - log(stats$group$group2))^2, na.rm = T)
    cv_ws <- exp(sqrt(sl / (2 * n))) - 1
    
    # KH way of calculating to match the naming scheme on Martin Bland's York site (sanity check)
    ## First we log transform.
    lx=log(stats$group$group1)
    ly=log(stats$group$group2)
    ## Calculate the within-subject variance for the log values.
    s21=(lx-ly)^2/2
    ## The within-subject standard deviation on the log scale is the square root of the mean within-subject variance. 
    sw=sqrt(mean(s21))
    ## The CV is the antilog (exponent since we are using natural logarithms) minus one.
    cv_ws_KH=exp(sw)-1
    
    # Check if their calc matches my calc
    if (cv_ws == cv_ws_KH) {
      print("The within-subject CV estimated by KH and the original function match.")
    }
    
    # Calc CIs for the log method according to Bland site 
    ## For the log method, we can find a confidence interval for the within-subject standard deviation on the log scale. 
    ## The standard error is sw/root(2n(m-1)), where: 
    ##  sw is the within-subject standard deviation 
    ##  n is the number of subjects
    ##  m is the number of observations per subject.
    SE_forCI=sw/sqrt(2*n*(2-1))
    
    # The 95% confidence interval is sw - [critical_val]*SE_forCI to sw + [critical_val]*SE_forCI
    ## To account for small sample sizes, we'll use the t-distribution
    ## Get t-critical value for 95% CI (1 - 0.05 / 2 = 0.975 quantile from the t-distribution with DF=n-1)
    crit_val = qt(0.975, (n-1))
    
    # Get 95% CI 
    low_bound_pre = sw - crit_val*SE_forCI
    upper_bound_pre = sw + crit_val*SE_forCI
    
    # Finally, we antilog these limits and subtract one to give confidence limits for the CV.
    ## These are slightly narrower than the root mean square confidence limits, but very similar.
    lower_CI=exp(low_bound_pre)-1
    upper_CI=exp(upper_bound_pre)-1
    
    # Return CV calc and CIs 
    results_list <- list("Log_CV" = cv_ws*100, "Lower_CI" = lower_CI*100, "Upper_CI" = upper_CI*100);
    return(results_list);
    
  } else if (method == "whole_dataset") {
    n <- stats$based.on
    sd <- sqrt(sum(stats$diffs^2) / (2 * n))
    cv_ws <- sd / mean(stats$means)
    
    # Return CV (no CI)
    results_list <- list("Whole_Dataset_CV" = cv_ws*100);
    return(results_list);
  }
  
}

# Example use: 
## stats <- bland.altman.stats(x, y)
## cv_ws <- calculate_within_subject_cv(stats, "whole_dataset")
