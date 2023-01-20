# within-subject-cv
This is an R function to calculate within-subject CVs and confidence intervals 

This function calculates within-subject CVs and CIs, according to Martin Bland's method, which is described in detail here: https://www-users.york.ac.uk/~mb55/meas/cv.htm

Rowan Nicholls provides another great reference on within-subject CV calculation: https://rowannicholls.github.io/R/statistics/agreement/coefficient_of_variation.html

CVs and CIs from this function match with those calculated by the MedCalc calculator: 
  https://www.medcalc.org/manual/cvfromduplicates.php

You may select one of three methods for within-subject CV calculation: 
  root_mean, logarithmic, and whole_dataset.
See Bland site for details on each. Note that the whole_dataset method is *not* recommended. 

You must enter the preferred confidence interval for the root_mean and logarithmic methods, for instance: 
  0.95, 0.90, 0.85