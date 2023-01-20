###############################################################################################################
# R function to calculate within-subject CVs and confidence intervals 
#   This function calculates within-subject CVs and CIs, according to Martin Bland's method, 
#   which is described in detail here: https://www-users.york.ac.uk/~mb55/meas/cv.htm
#
#   Rowan Nicholls provides another great reference on within-subject CV calculation: 
#   https://rowannicholls.github.io/R/statistics/agreement/coefficient_of_variation.html
# 
#   CVs and CIs from this function match with those calculated by the MedCalc calculator: 
#   https://www.medcalc.org/manual/cvfromduplicates.php
#
#   You may select one of three methods for within-subject CV calculation: 
#     root_mean, logarithmic, and whole_dataset
#   See Bland site for details on each. Note that the whole_dataset method is *not* recommended. 
#
#   You must enter the preferred confidence interval, for instance: 
#     0.95, 0.90, 0.85
# 
#   Example usage to calculate the within-subject CV using each method
#     library(BlandAltmanLeh) # required package
#     t1 <- c(10,15,25,30,22,14) # time point 1 
#     t2 <- c(11,14.5,22.5,31,21,15) # time point 2 
#     stats <- bland.altman.stats(t1, t2) # run bland.altman.stats from the BlandAltmanLeh library
#     calculate_within_subject_cv(stats, "root_mean", 0.90) # use root mean method and 90% CIs
#     calculate_within_subject_cv(stats, "logarithmic", 0.90) # use logarithmic method and 90% CIs
#     calculate_within_subject_cv(stats, "whole_dataset") # use whole_group method (CI calc not possible)
#
###############################################################################################################

calculate_within_subject_cv <- function(stats, method, CI) {
  
  ### ROOT MEAN METHOD ### 
  
  if (method == "root_mean") {
    # Calculate the within-subject variance for the natural scale values. 
    ## Within-subject variance is given by difference squared over 2 when we have pairs of subjects.
    ## (Match the variable naming scheme on Martin Bland's site.)
    s2 <- (stats$diffs)^2/2
  
    # Calculate subject mean and s squared / mean squared, i.e. CV squared.
    m <- stats$means
    s2m2 <- s2/m^2
    
    # Calculate mean of s squared / mean squared.
    mean_s2m2 <- mean(s2m2)
    
    # The within-subject CV is the square root of the mean of s squared / mean squared:
    cv_ws <- sqrt(mean_s2m2)
    
    # Calculate confidence intervals.
    ## For the root mean square method, this is very direct. 
    ## We have the mean of the squared CV == mean_s2m2 (i.e., the term right before you take the square root)... 
    ## ...so we use the usual confidence interval for a mean on this then take the square root of the sample size.
    ## i.e., the standard error is the standard deviation of the CVs divided by the square root of the sample size.
    n <- stats$based.on
    SE_term <- sd(s2m2)/sqrt(n)
    
    # Then, the 95% confidence interval for the squared CV can be found by the mean minus or plus 1.96 standard errors. 
    ## HOWEVER, if the sample is small we should use the t-distribution INSTEAD of the z-distribution. 
    ## (Bland used 1.96 because he was using a z-distribution because his simulation was with 100 subjects)
    ## (However, Bland also notes: the squared CVs are unlikely to be Normal, so the CI will still be very approximate.)
    
    # Find the t_critical value for the CI using the t-distribution; for example: 
    ## 95% CI prob: 1-(1-.95)/2 = 0.975
    ## 90% CI prob: 1-(1-.90)/2 = 0.950
    prob <- 1-(1-CI)/2
    
    ## 95% CI crit_val with n=20 subjects: 2.093024
    ## 90% CI crit_val with n=20 subjects: 1.729133
    crit_val = qt(prob, (n-1))

    ## Next, calculate mean +/- t_critical*SE  
    low_bound_sq <- mean_s2m2 - crit_val*SE_term
    high_bound_sq <- mean_s2m2 + crit_val*SE_term
    
    ## The square roots of these limits give the confidence interval for the CV.
    lower_CI = sqrt(low_bound_sq)
    upper_CI = sqrt(high_bound_sq)
    
    # Return CV calc and CIs 
    results_list <- list("Root_Mean_CV" = cv_ws*100, "Lower_CI" = lower_CI*100, "Upper_CI" = upper_CI*100);
    return(results_list);
  
  ### LOGARITHMIC METHOD ### 
    
  } else if (method == "logarithmic") {
    # First we log transform.
    ## (Match the variable naming scheme on Martin Bland's site.)
    lx=log(stats$group$group1)
    ly=log(stats$group$group2)
    
    ## Calculate the within-subject variance for the log values.
    s21=(lx-ly)^2/2
    
    ## The within-subject standard deviation on the log scale is the square root of the mean within-subject variance. 
    sw=sqrt(mean(s21))
    
    ## The CV is the antilog (exponent since we are using natural logarithms) minus one.
    cv_ws=exp(sw)-1
    
    # Calc CIs for the log method.
    ## For the log method, we can find a confidence interval for the within-subject standard deviation on the log scale. 
    ## The standard error is sw/root(2n(m-1)), where: 
    ##  sw is the within-subject standard deviation 
    ##  n is the number of subjects
    ##  m is the number of observations per subject.
    n <- stats$based.on
    SE_forCI=sw/sqrt(2*n*(2-1))
    
    ## The 95% confidence interval is sw - [critical_val]*SE_forCI to sw + [critical_val]*SE_forCI
    ## To account for small sample sizes, we'll use the t-distribution
    ## Find the t_critical value for the CI using the t-distribution; for example: 
    ## 95% CI prob: 1-(1-.95)/2 = 0.975
    ## 90% CI prob: 1-(1-.90)/2 = 0.950
    prob <- 1-(1-CI)/2
    
    ## 95% CI crit_val with n=20 subjects: 2.093024
    ## 90% CI crit_val with n=20 subjects: 1.729133
    crit_val = qt(prob, (n-1))
    
    ## Get the CI 
    low_bound_pre = sw - crit_val*SE_forCI
    upper_bound_pre = sw + crit_val*SE_forCI
    
    ## Finally, we antilog these limits and subtract one to give confidence limits for the CV.
    ## These are slightly narrower than the root mean square confidence limits, but very similar.
    lower_CI=exp(low_bound_pre)-1
    upper_CI=exp(upper_bound_pre)-1
    
    # Return CV calc and CIs 
    results_list <- list("Log_CV" = cv_ws*100, "Lower_CI" = lower_CI*100, "Upper_CI" = upper_CI*100);
    return(results_list);

    ### WHOLE DATASET METHOD ### 
        
  } else if (method == "whole_dataset") {
    # Use the 'whole dataset' method; note that this is *NOT* recommended
    n <- stats$based.on
    sd <- sqrt(sum(stats$diffs^2) / (2 * n))
    cv_ws <- sd / mean(stats$means)
    
    # Return CV only; a CI calculation is not possible here 
    results_list <- list("Whole_Dataset_CV" = cv_ws*100);
    return(results_list);
  }
  
}
