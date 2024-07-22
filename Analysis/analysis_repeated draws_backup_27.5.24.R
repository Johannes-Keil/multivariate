################### LOAD PACKAGES ######################
#Load Packages
packages <- c("tidyverse", "lme4", "lmerTest", "ggplot2", "forestplot", "nlme", "esmpack", "latticeExtra")

# Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())
if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages])
}

# install the 'esmpack' package (version 0.1-17)
if(!require("esmpack"))requireremotes::install_github("wviechtb/esmpack@v0.1-17")

invisible(lapply(packages, library, character.only = TRUE))
invisible(library(esmpack))


######################################################################################################
##################### DATA PREPARATION  ##############################################################
######################################################################################################

#read in data
dat <- read.table("data_prep.dat", header=TRUE, sep="\t", na.strings="", as.is=TRUE) %>% 
  data.frame() %>% 
  mutate(
    #ensure data are numeric
    study = as.numeric(study),
    status = as.numeric(status),
    pa = as.numeric(pa),
    na = as.numeric(na),
    beep = as.numeric(beep)
  )

#################### EXPLORE DATA ##############

#Check missingness
check.nomiss(beep, subjno, data=dat)
check.nomiss(pa, subjno, data = dat)
check.nomiss(na, subjno, data = dat)

#there is missingness in positive & negative affect
#sometimes entire sub-studies have very high missingness, and thus need to be excluded
study_list <- dat$study %>% unique()

missingness <- c()

#what proportion of entries is allowed to be missing in either pa or na for each group/sub-study?
missingness_criterion <- 0.9

for (i in study_list) {
  s <- filter(dat, study == i)
  missing_na <- sum(is.na(s$na)) / length(s$na)
  missing_pa <- sum(is.na(s$pa)) /length(s$pa)
  
  if(missing_na > missingness_criterion | missing_pa > missingness_criterion ) {
    print(i)
    missingness[length(missingness) + 1] <- i
  }
}

#this gives a list of all studies with above-criterion missingness
missingness

#Check time invariance of status variables
check.timeinvar(status, subjno, data=dat, na.rm=FALSE)

#frequency tables for out of range values
#nothing out of the ordinary
table(dat$subjno)
table(dat$study)
table(dat$beep)
table(dat$status)
hist(dat$pa)
hist(dat$na)

######################################################################################################
##################### ANALYSE WITH FIXED SAMPLE SIZE WITH REPEATED SAMPLING ##########################
######################################################################################################

# The previous analysis suffers from the problem that we only draw one sample for sub-study.
# This means that we can show that classification errors occur depending on which IV/DV we choose
# but not how often.
# For this purpose, draw sub-samples for every study for a large number of times.
# This time, do not restrict the number of beeps per patient.

############ OPTIONS AND PARAMETERS ###########  ###############

#how large is each sample?
#smaller samples = lower power and more unstable results (hence, direction-related classification errors become more likely)
sample_size = 20

#how many samples are drawn per study/group combination?
n_samples = 2

#increase limit for iterations until convergence is expected
options <- lmeControl(maxIter = 10000000)

#seed for random-number generator
set.seed(080545)

#show progress messages during the run?
show_messages = TRUE
show_iterations = FALSE

#show error messages during the run?
show_errors = FALSE

#save the resulting data in an rData file?
save_models = FALSE


# DATA EXCLUSION
# looks like studies 2 and 4 don't meet the missingness criterion of having valid data for at
# #least 10% of cases, so exclude these studies (see data checks above)
# #for now, also exclude study 9 due to convergence issues

missingness

#which studies to exclude?
excluded <- c(2, 4, 9)
dat <- dat %>% 
  filter(!(study %in% excluded))

#save a list of studies to iterate over
study_list <- dat$study %>% unique()


############ DEFINE MODEL CONTAINERS #############################

#store data for the first model
simple_pa <- tibble(study = numeric(),                #which study does this dataset correspond to?
                    status = numeric(),               #what's the diagnostic status
                    simulation = numeric(),           #which run/sample is this model based on?
                    m_simple_pa = list(),         #simple random slopes + intercept model with one direction of prediction
)

#store data for the second model
simple_na <- tibble(study = numeric(),                #which study does this dataset correspond to?
                    status = numeric(),               #what's the diagnostic status
                    simulation = numeric(),           #which run/sample is this model based on?
                    m_simple_na = list(),         #simple random slopes + intercept model with one direction of prediction
)

#store data for the third model
centred_pa <- tibble(study = numeric(),                #which study does this dataset correspond to?
                    status = numeric(),               #what's the diagnostic status
                    simulation = numeric(),           #which run/sample is this model based on?
                    m_centred_pa = list(),         #simple random slopes + intercept model with one direction of prediction
)

#store data for the fourth  model
centred_na <- tibble(study = numeric(),                #which study does this dataset correspond to?
                     status = numeric(),               #what's the diagnostic status
                     simulation = numeric(),           #which run/sample is this model based on?
                     m_centred_na = list(),         #simple random slopes + intercept model with one direction of prediction
)

#store data for the fifth model
centred_pa_slopes <- tibble(study = numeric(),                #which study does this dataset correspond to?
                     status = numeric(),               #what's the diagnostic status
                     simulation = numeric(),           #which run/sample is this model based on?
                     m_centred_pa_slopes = list(),         #simple random slopes + intercept model with one direction of prediction
)

#store data for the sixth  model
centred_na_slopes <- tibble(study = numeric(),                #which study does this dataset correspond to?
                            status = numeric(),               #what's the diagnostic status
                            simulation = numeric(),           #which run/sample is this model based on?
                            m_centred_na_slopes = list(),         #simple random slopes + intercept model with one direction of prediction
)

#store data for the seventh  model
multivar <- tibble(study =  numeric(),                #which study does this dataset correspond to?
                            status = numeric(),       #what's the diagnostic status
                            simulation = numeric(),   #which run/sample is this model based on?
                            m_multi = list(),         #simple random slopes + intercept model with one direction of prediction
)




############ EXECUTE ANALYSIS LOOP ###########################################

for (study_no in study_list) {
  
  #iterate across all status groups within this study
  status_list <- filter(dat, study == study_no)$status %>% unique()
  status_list
  
  #currently the set of participants to sample from is determined separately for each iteration
  #this means that, in the final dataset, the values in each row are actually based on separate draws from the data
  #this, of course, is unacceptable, because our conclusions wouldn't be supported.
  #Solution: Return to an algorithm that draws data once then computes ALL the models
  #Runtime will have to suffer.
  
  for (status_no in status_list) {
    
    ############## MODEL 1: Univariate positive_affect ~ negative affect #############
    
    #reset sample number iterator
    sample_no = 1
    
    #iterate across the required number of draws
    #Sidenote: The overall sample size is very low, leading to frequent convergence issues.
    #To improve runtime, this while-loop is repeated individually for each model
    while (sample_no <= n_samples) {
      
      ################ DATA PREPARATION #################
      
      #now sample a sub-set of participants
      participant_list <- dat %>% 
        filter(study == study_no & status == status_no) %>% 
        unique() %>% 
        sample_n(sample_size) %>% 
        select(subjno)
      
      #create an empty dataframe to store results
      this_dat <- data.frame(matrix(nrow = 0, ncol = 6))
      colnames(this_dat) <- c("subjno", "study", "status", "beep", "pa", "na")
      
      #sample the expected participants with replacement
      for (id in participant_list$subjno) {
        
        this_dat <- dat %>% 
          filter(study == 1 & status == 2,
                 subjno == id) %>% 
          rbind(this_dat, .)
      }
      
      #remove the participants with missing data
      #and transform data to match the multivariate model format
      this_dat <- this_dat %>% 
        filter(
          !is.na(pa) & !is.na(na)
        )
      
      ########## MODEL FITTING #####################
      
      #fit a linear model predicting PA from NA
      model <-try(lme(pa ~ na, random =~ na | subjno, data=this_dat, na.action=na.omit, control = options), silent = !show_errors)
      
      #if convergence problems, retry loop
      if (inherits(model, "try-error")) {
        if(show_messages) print("Retry")
        next
      }
      
      #extract parameters
      simple_pa <- simple_pa %>% 
        add_row(study = study_no, status = status_no, simulation = sample_no,
                m_simple_pa = list(model)
        )
      
      if(show_messages)
        paste(toString(sample_no), ". Iteration Successful", sep = "") %>% print()
      
      sample_no = sample_no + 1
      
    }#sample_no
    
    if(show_messages) print("Model 1 Successful")
    
    ############## MODEL 2: Univariate negative_affect ~ positive_affect #############
    
    #reset sample number iterator
    sample_no = 1
    
    #iterate across the required number of draws
    #Sidenote: The overall sample size is very low, leading to frequent convergence issues.
    #To improve runtime, sampling is done individually for each model
    while (sample_no <= n_samples) {
      
      ################ DATA PREPARATION #################
      #now sample a sub-set of participants
      participant_list <- dat %>% 
        filter(study == study_no & status == status_no) %>% 
        unique() %>% 
        sample_n(sample_size) %>% 
        select(subjno)
      
      #create an empty dataframe to store results
      this_dat <- data.frame(matrix(nrow = 0, ncol = 6))
      colnames(this_dat) <- c("subjno", "study", "status", "beep", "pa", "na")
      
      #sample the expected participants with replacement
      for (id in participant_list$subjno) {
        
        this_dat <- dat %>% 
          filter(study == 1 & status == 2,
                 subjno == id) %>% 
          rbind(this_dat, .)
      }
      
      #remove the participants with missing data
      #and transform data to match the multivariate model format
      this_dat <- this_dat %>% 
        filter(
          !is.na(pa) & !is.na(na)
        )
      
      ########## MODEL FITTING #####################
      
      #fit a linear model predicting NA from PA
      model <-try(lme(na ~ pa, random =~ pa | subjno, data=this_dat, na.action=na.omit, control = options), silent = !show_errors)
      
      #if convergence problems, retry loop
      if (inherits(model, "try-error")) {
        if(show_messages) print("Retry")
        next
      }
      
      #extract parameters
      simple_na <- simple_na %>% 
        add_row(study = study_no, status = status_no, simulation = sample_no,
                m_simple_na = list(model)
        )
      
      if(show_messages)
        paste(toString(sample_no), ". Iteration Successful", sep = "") %>% print()
      
      sample_no = sample_no + 1
      
    }#sample_no
    
    if(show_messages) print("Model 2 Successful")
    
    ############## MODEL 3: Centred positive affect ~ negative_affect + NO random-effects #############
    
    #reset sample number iterator
    sample_no = 1
    
    #iterate across the required number of draws
    #Sidenote: The overall sample size is very low, leading to frequent convergence issues.
    #To improve runtime, sampling is done individually for each model
    while (sample_no <= n_samples) {
      
      ################ DATA PREPARATION #################
      #now sample a sub-set of participants
      participant_list <- dat %>% 
        filter(study == study_no & status == status_no) %>% 
        unique() %>% 
        sample_n(sample_size) %>% 
        select(subjno)
      
      #create an empty dataframe to store results
      this_dat <- data.frame(matrix(nrow = 0, ncol = 6))
      colnames(this_dat) <- c("subjno", "study", "status", "beep", "pa", "na")
      
      #sample the expected participants with replacement
      for (id in participant_list$subjno) {
        
        this_dat <- dat %>% 
          filter(study == 1 & status == 2,
                 subjno == id) %>% 
          rbind(this_dat, .)
      }
      
      #remove the participants with missing data
      this_dat <- this_dat %>% 
        filter(
          !is.na(pa) & !is.na(na)
        )
      
      #now add standardized values to the dataset (with respect to the sample means)
      
      this_dat$zpa <- ave(this_dat$pa, this_dat$subjno, FUN = function(x) (x - mean(x, na.rm=TRUE)) / sd(x, na.rm=TRUE))
      this_dat$zna <- ave(this_dat$na, this_dat$subjno, FUN = function(x) (x - mean(x, na.rm=TRUE)) / sd(x, na.rm=TRUE))
      
      
      ########## MODEL FITTING #####################
      
      #fit a linear model predicting NA from PA
      model <-try(lm(zpa ~ 0 + zna, data=this_dat, na.action=na.omit, control = options), silent = !show_errors)
      
      #if convergence problems, retry loop
      if (inherits(model, "try-error")) {
        if(show_messages) print("Retry")
        next
      }
      
      #extract parameters
      centred_pa <- centred_pa %>% 
        add_row(study = study_no, status = status_no, simulation = sample_no,
                m_centred_pa = list(model)
        )
      
      if(show_messages)
        paste(toString(sample_no), ". Iteration Successful", sep = "") %>% print()
      
      sample_no = sample_no + 1
      
    }#sample_no
    
    if(show_messages) print("Model 3 Successful")
    
    ############## MODEL 4: Centred negative_affect ~ positive affect + NO random-effects #############
    
    #reset sample number iterator
    sample_no = 1
    
    #iterate across the required number of draws
    #Sidenote: The overall sample size is very low, leading to frequent convergence issues.
    #To improve runtime, sampling is done individually for each model
    while (sample_no <= n_samples) {
      
      ################ DATA PREPARATION #################
      #now sample a sub-set of participants
      participant_list <- dat %>% 
        filter(study == study_no & status == status_no) %>% 
        unique() %>% 
        sample_n(sample_size) %>% 
        select(subjno)
      
      #create an empty dataframe to store results
      this_dat <- data.frame(matrix(nrow = 0, ncol = 6))
      colnames(this_dat) <- c("subjno", "study", "status", "beep", "pa", "na")
      
      #sample the expected participants with replacement
      for (id in participant_list$subjno) {
        
        this_dat <- dat %>% 
          filter(study == 1 & status == 2,
                 subjno == id) %>% 
          rbind(this_dat, .)
      }
      
      #remove the participants with missing data
      #and transform data to match the multivariate model format
      this_dat <- this_dat %>% 
        filter(
          !is.na(pa) & !is.na(na)
        )
      
      #now add standardized values to the dataset (with respect to the sample means)
      
      this_dat$zpa <- ave(this_dat$pa, this_dat$subjno, FUN = function(x) (x - mean(x, na.rm=TRUE)) / sd(x, na.rm=TRUE))
      this_dat$zna <- ave(this_dat$na, this_dat$subjno, FUN = function(x) (x - mean(x, na.rm=TRUE)) / sd(x, na.rm=TRUE))
      
      
      ########## MODEL FITTING #####################
      
      #fit a linear model predicting NA from PA
      model <-try(lm(zna ~ 0 + zpa, data=this_dat, na.action=na.omit, control = options), silent = !show_errors)
      
      #if convergence problems, retry loop
      if (inherits(model, "try-error")) {
        if(show_messages) print("Retry")
        next
      }
      
      #extract parameters
      centred_na <- centred_na %>% 
        add_row(study = study_no, status = status_no, simulation = sample_no,
                m_centred_na = list(model)
        )
      
      if(show_messages)
        paste(toString(sample_no), ". Iteration Successful", sep = "") %>% print()
      
      sample_no = sample_no + 1
      
    }#sample_no
    
    if(show_messages) print("Model 4 Successful")
    
    
    ############## MODEL 5: Centred positive affect ~ negative_affect + Random Effects #############
    
    #reset sample number iterator
    sample_no = 1
    
    #iterate across the required number of draws
    #Sidenote: The overall sample size is very low, leading to frequent convergence issues.
    #To improve runtime, sampling is done individually for each model
    while (sample_no <= n_samples) {
      
      ################ DATA PREPARATION #################
      #now sample a sub-set of participants
      participant_list <- dat %>% 
        filter(study == study_no & status == status_no) %>% 
        unique() %>% 
        sample_n(sample_size) %>% 
        select(subjno)
      
      #create an empty dataframe to store results
      this_dat <- data.frame(matrix(nrow = 0, ncol = 6))
      colnames(this_dat) <- c("subjno", "study", "status", "beep", "pa", "na")
      
      #sample the expected participants with replacement
      for (id in participant_list$subjno) {
        
        this_dat <- dat %>% 
          filter(study == 1 & status == 2,
                 subjno == id) %>% 
          rbind(this_dat, .)
      }
      
      #remove the participants with missing data
      #and transform data to match the multivariate model format
      this_dat <- this_dat %>% 
        filter(
          !is.na(pa) & !is.na(na)
        )
      
      
      #now add standardized values to the dataset (with respect to the sample means)
      
      this_dat$zpa <- ave(this_dat$pa, this_dat$subjno, FUN = function(x) (x - mean(x, na.rm=TRUE)) / sd(x, na.rm=TRUE))
      this_dat$zna <- ave(this_dat$na, this_dat$subjno, FUN = function(x) (x - mean(x, na.rm=TRUE)) / sd(x, na.rm=TRUE))
      
      
      ########## MODEL FITTING #####################
      
      #fit a linear model predicting NA from PA
      model <-try(lme(zpa ~ 0 + zna, random =~ 0 + zna | subjno, data=this_dat, na.action=na.omit, control = options), silent = !show_errors)
      
      #if convergence problems, retry loop
      if (inherits(model, "try-error")) {
        if(show_messages) print("Retry")
        next
      }
      
      #extract parameters
      centred_pa_slopes <- centred_pa_slopes %>% 
        add_row(study = study_no, status = status_no, simulation = sample_no,
                m_centred_pa_slopes = list(model)
        )
      
      if(show_messages)
        paste(toString(sample_no), ". Iteration Successful", sep = "") %>% print()
      
      sample_no = sample_no + 1
      
    }#sample_no
    
    if(show_messages) print("Model 5 Successful")
    
    ############## MODEL 6: Centred negative_affect ~ positive_affect + Random Effects #############
    
    #reset sample number iterator
    sample_no = 1
    
    #iterate across the required number of draws
    #Sidenote: The overall sample size is very low, leading to frequent convergence issues.
    #To improve runtime, sampling is done individually for each model
    while (sample_no <= n_samples) {
      
      ################ DATA PREPARATION #################
      ################ DATA PREPARATION #################
      #now sample a sub-set of participants
      participant_list <- dat %>% 
        filter(study == study_no & status == status_no) %>% 
        unique() %>% 
        sample_n(sample_size) %>% 
        select(subjno)
      
      #create an empty dataframe to store results
      this_dat <- data.frame(matrix(nrow = 0, ncol = 6))
      colnames(this_dat) <- c("subjno", "study", "status", "beep", "pa", "na")
      
      #sample the expected participants with replacement
      for (id in participant_list$subjno) {
        
        this_dat <- dat %>% 
          filter(study == 1 & status == 2,
                 subjno == id) %>% 
          rbind(this_dat, .)
      }
      
      #remove the participants with missing data
      #and transform data to match the multivariate model format
      this_dat <- this_dat %>% 
        filter(
          !is.na(pa) & !is.na(na)
        )
      
      
      #now add standardized values to the dataset (with respect to the sample means)
      
      this_dat$zpa <- ave(this_dat$pa, this_dat$subjno, FUN = function(x) (x - mean(x, na.rm=TRUE)) / sd(x, na.rm=TRUE))
      this_dat$zna <- ave(this_dat$na, this_dat$subjno, FUN = function(x) (x - mean(x, na.rm=TRUE)) / sd(x, na.rm=TRUE))
      
      
      ########## MODEL FITTING #####################
      
      #fit a linear model predicting NA from PA
      model <-try(lme(zna ~ 0 + zpa, random =~ 0 + zpa | subjno, data=this_dat, na.action=na.omit, control = options), silent = !show_errors)
      
      #if convergence problems, retry loop
      if (inherits(model, "try-error")) {
        if(show_messages) print("Retry")
        next
      }
      
      #extract parameters
      centred_na_slopes <- centred_na_slopes %>% 
        add_row(study = study_no, status = status_no, simulation = sample_no,
                m_centred_na_slopes = list(model)
        )
      
      if(show_messages)
        paste(toString(sample_no), ". Iteration Successful", sep = "") %>% print()
      
      sample_no = sample_no + 1
      
    }#sample_no
    
    if(show_messages) print("Model 6 Successful")
    
    
    ############## MODEL 7: Multivariate Intercept-only model #############
    
    #reset sample number iterator
    sample_no = 1
    
    #iterate across the required number of draws
    #Sidenote: The overall sample size is very low, leading to frequent convergence issues.
    #To improve runtime, sampling is done individually for each model
    while (sample_no <= n_samples) {
      
      ################ DATA PREPARATION #################
      
      #now sample a sub-set of participants
      participant_list <- dat %>% 
        filter(study == study_no & status == status_no) %>% 
        unique() %>% 
        sample_n(sample_size) %>% 
        select(subjno)
      
      #create an empty dataframe to store results
      this_dat <- data.frame(matrix(nrow = 0, ncol = 6))
      colnames(this_dat) <- c("subjno", "study", "status", "beep", "pa", "na")
      
      #sample the expected participants with replacement
      for (id in participant_list$subjno) {
        
        this_dat <- dat %>% 
          filter(study == 1 & status == 2,
                 subjno == id) %>% 
          rbind(this_dat, .)
      }
      
      #remove the participants with missing data
      #and transform data to match the multivariate model format
      this_dat <- this_dat %>% 
        filter(
          !is.na(pa) & !is.na(na)
        ) %>% 
        pivot_longer(
          cols = c("pa", "na"),
          names_to = "affect_name",
          values_to = "affect"
        ) %>% 
        mutate(
          #create dummy variable representing which type of affect is being predicted
          is_pa = ifelse(affect_name == "pa", 1, 0),
          is_na = ifelse(affect_name == "na", 1, 0)
        )
      
      ########## MODEL FITTING #####################
      
      #fit a linear model predicting NA from PA
      model <-try(lme(zna ~ 0 + zpa, random =~ 0 + zpa | subjno, data=this_dat, na.action=na.omit, control = options), silent = !show_errors)
      
      #if convergence problems, retry loop
      if (inherits(model, "try-error")) {
        if(show_messages) print("Retry")
        next
      }
      
      #extract parameters
      centred_na_slopes <- centred_na_slopes %>% 
        add_row(study = study_no, status = status_no, simulation = sample_no,
                m_centred_na_slopes = list(model)
        )
      
      if(show_messages)
        paste(toString(sample_no), ". Iteration Successful", sep = "") %>% print()
      
      sample_no = sample_no + 1
      
    }#sample_no
    
    if(show_messages) print("Model 7 Successful")
    
    
    ############ Get ready for the next run ################
    
    if(show_messages) print("****************Next Status*********************")
    
  }#status_no
  
  if(show_messages) print(paste("################# ", study_no, ". Study Done #####################"))
  
}#study_no


############ POST_PROCESSING ####################

#Unfortunately, accessing the model's content is a bit cumbersome
#accessing the models for checking
#to illustrate
#x <- res$m_simple_na[1][[1]]$coefficients

#putting results in a single tibble
full_models <- simple_pa %>% 
  inner_join(simple_na, by = c("study", "status", "simulation")) %>% 
  inner_join(centred_pa, by = c("study", "status", "simulation")) %>% 
  inner_join(centred_na, by = c("study", "status", "simulation")) %>% 
  inner_join(centred_pa_slopes, by = c("study", "status", "simulation")) %>% 
  inner_join(centred_na_slopes, by = c("study", "status", "simulation"))

#extract coefficients as needed for plotting

#dictionary: p_pa = p-value when predicting pa ~ na in the simple model
#            p_na = p-value when predicting na ~ pa in the simple model
#            b_pa = corresponding beta-value
#            b_na = corresponding beta-value
#            r_uncorr = correlation in the standardized regression model, (from pa ~ na direction, opposite direction should deliver similar results)
#            r_uncorr_lower = lower 95% confidence interval for this correlation
#            r_uncorr_upper = upper 95% confidence interval for this correlation
#            r_multi = correlation from the multivariate for confirmation, should be almost identical to r_uncorr
#            r_corr = correlation after correction with random slopes (from pa ~ na direction, opposite direction should deliver similar results)
#            r_corr_lower = lower 95% confidence interval for this correlation
#            r_corr_upper = upper 95% confidence interval for this correlation

coefficient_columns <- c("p_na", "p_pa","b_na", "b_pa",  "r_uncorr", "r_multi", "r_corr", "r_uncorr_lower", "r_uncorr_upper", "r_corr_lower", "r_corr_upper")
full_models[ , coefficient_columns] <- NA

for(i in 1:nrow(full_models)) {
  
  #p-values for showing the directionality problem
  full_models$p_pa[i] = summary(full_models$m_simple_pa[i][[1]])$tTable[2, 5]
  full_models$p_na[i] = summary(full_models$m_simple_na[i][[1]])$tTable[2, 5]
  
  #coefficients from the univariate non-standardized model
  full_models$b_pa[i] = summary(full_models$m_simple_pa[i][[1]])$coefficients$fixed[2]
  full_models$b_na[i] = summary(full_models$m_simple_na[i][[1]])$coefficients$fixed[2]
  
  #values for the multivariate model
  
  #values for testing the correction method
  full_models$r_uncorr[i] = summary(full_models$m_centred_pa[i][[1]])$coefficients[1]
  full_models$r_uncorr_lower[i] = confint(full_models$m_centred_pa[i][[1]])[1]
  full_models$r_uncorr_upper[i] = confint(full_models$m_centred_pa[i][[1]])[2]
  

  full_models$r_corr[i] = summary(full_models$m_centred_pa_slopes[i][[1]])$coefficients$fixed[1]
  
  #the intervals() function will occasionally cause trouble, in particular if models did not converge correctly
  #(e.g, non-positive definite approximate variance-covariance)
  #As the number of models estimated through this algorithm is too large to examine the underlying problems individually
  #simply set the confidence intervals for these models as NA
  full_models$r_corr_lower[i] = try(intervals(full_models$m_centred_pa_slopes[i][[1]])[[1]][1], silent = TRUE) %>% as.numeric()
  full_models$r_corr_upper[i] = try(intervals(full_models$m_centred_pa_slopes[i][[1]])[[1]][3], silent = TRUE) %>% as.numeric()
}

#for plotting, add some housekeeping variables

full_models <- full_models %>% 
  mutate(
    status_str = ifelse(
      status == 0, "Control",
      ifelse(status == 1, "Psychosis Risk", 
             ifelse(status == 2, "Psychosis", "Depression"))
    ),
    study_str = paste("Study", study)
  )
#optionally save results
if(save_models) saveRDS(full_models, file = "full_models.rData")


######################################################################################################
##################### PLOTTING RESULTS################################################################
######################################################################################################

#### Plot 1: Type I error rate for pa ~ na against Type I error rate for na ~ pa simple models
#### We observe that discrepancies emerge only for psychosis -> Why?
#### Larger effect size for other studies -> Larger power and lower risk of false positives
#### Also: Larger inter-individual differences in Psychosis than the other groups?
#### To check, one could run a simulation, varying (a) overall r and (b) variance in intra-individual correlations.

plot1 <- ggplot(data = full_models, aes(x = p_na, y = p_pa)) +
  facet_wrap(~status_str) + 
  geom_jitter(aes(x = p_na, y = p_pa, shape = study_str),
             show.legend = TRUE, width = 0.001, height = 0.001) +
  scale_shape_manual(values=seq(0,length(study_list))) +
  labs(shape = "Study", color = "Patient Status") +
  geom_hline(yintercept = 0.05) +
  geom_vline(xintercept = 0.05) +
  xlab("Negative Affect as Dependent Variable") +
  ylab("Positive Affect as Dependent Variable") +
  ggtitle("A) P-values by Direction of Prediction") +
  theme_bw() +
  theme(
    strip.background = element_blank()
  )
plot1

#### Plot 2: Plotting the corresponding beta-parameters of the two models against each other
#### We expect an inverse relationship!

mean_pa <- mean(full_models$b_pa, na.action = na.omit)
mean_na <- mean(full_models$b_na, na.action = na.omit)

plot2 <- ggplot(data = full_models, aes(x = b_na, y = b_pa)) +
  geom_jitter(aes(x = b_na, y = b_pa, color = status_str, shape = study_str),
              show.legend = TRUE, width = 0.001, height = 0.001) +
  scale_shape_manual(values=seq(0,length(study_list))) +
  labs(shape = "Study", color = "Patient Status") +
  geom_hline(yintercept = mean_pa) +
  geom_vline(xintercept = mean_na) +
  xlab("Negative Affect as Dependent Variable") +
  ylab("Positive Affect as Dependent Variable") +
  ggtitle("B) Regression Weights by Direction of Prediction") +
  theme_bw()
plot2

### Plot 3: Validation: Plotting the standardized r values against correlation from multivariate model


#### Plot 4: Scatterplot showing change in confidence intervals between correction methods

plot3 <- ggplot(data = full_models) +
  facet_wrap(~status_str) + 
  geom_point(aes(x = r_uncorr, y = r_corr, color = "Point Estimate")) +
  geom_point(aes(x = r_uncorr_lower, y = r_corr_lower, color = "Lower 95% CI")) +
  geom_point(aes(x = r_uncorr_upper, y = r_corr_upper, color = "Upper 95% CI")) +
  geom_abline(intercept = 0, slope = 1) + 
  labs(color = "Correlation") +
  xlab("Uncorrected Correlation") +
  ylab("Corrected Correlation") +
  ggtitle("D) Correlations in Corrected and Uncorrected Model")+
  xlim(c(-1, 0.3)) +
  ylim(c(-1, 0.3)) +
  theme_bw() +
  theme(
    strip.background = element_blank()
  )
plot3



################## Trying to fit a multivariate model ############

#add the necessary dummies to the data

dat_multivar <- dat %>% 
  pivot_longer(
    cols = c("pa", "na"),
    names_to = "affect_name",
    values_to = "affect"
  ) %>% 
  mutate(
    #create dummy variable representing which type of affect is being predicted
    is_pa = ifelse(affect_name == "pa", 1, 0),
    is_na = ifelse(affect_name == "na", 1, 0)
  )
#fitting a simple intercept-only multivariate model
multivar <- lme(fixed = affect ~ -1 + is_pa + is_na,
                  random = ~ -1 + is_pa + is_na|subjno,
                  data=dat_multivar,
                  na.action = na.exclude)

################# FIXING THE SAMPLING PROCEDURE


#now sample a sub-set of participants
participant_list <- dat %>% 
  filter(study == study_no & status == status_no) %>% 
  unique() %>% 
  sample_n(sample_size) %>% 
  select(subjno)

#create an empty dataframe to store results
this_dat <- data.frame(matrix(nrow = 0, ncol = 6))
colnames(this_dat) <- c("subjno", "study", "status", "beep", "pa", "na")

#sample the expected participants with replacement
for (id in participant_list$subjno) {
  
  this_dat <- dat %>% 
    filter(study == 1 & status == 2,
           subjno == id) %>% 
    rbind(this_dat, .)
}
  
#remove the participants with missing data
#and transform data to match the multivariate model format
this_dat <- this_dat %>% 
  filter(
    !is.na(pa) & !is.na(na)
  )

#remove the participants with missing data
#and transform data to match the multivariate model format
this_dat <- this_dat %>% 
  filter(
    !is.na(pa) & !is.na(na)
  ) %>% 
  pivot_longer(
    cols = c("pa", "na"),
    names_to = "affect_name",
    values_to = "affect"
  ) %>% 
  mutate(
    #create dummy variable representing which type of affect is being predicted
    is_pa = ifelse(affect_name == "pa", 1, 0),
    is_na = ifelse(affect_name == "na", 1, 0)
  )


