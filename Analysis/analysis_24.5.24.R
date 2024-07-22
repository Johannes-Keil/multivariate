# install the 'remotes' package
install.packages("remotes")

# install the 'esmpack' package (version 0.1-17)
remotes::install_github("wviechtb/esmpack@v0.1-17")

# install the 'latticeExtra' package
install.packages("latticeExtra")



library(tidyverse)
library(lme4)
library(lmerTest)
library(ggplot2)
library(forestplot)

# load the 'nlme' package
library(nlme)

# load the 'esmpack' package
library(esmpack)

# load the 'lattice' and 'latticeExtra' packages
library(lattice)
library(latticeExtra)

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

#################### GET AN IDEA OF THE DATA ##############

#Check missingness
check.nomiss(beep, subjno, data=dat)
check.nomiss(pa, subjno, data = dat)
check.nomiss(na, subjno, data = dat)

#there is missingness in positive & negative affect
#sometimes entire sub-studies have very high missingness, and thus need to be excluded
study_list <- dat$study %>% unique()

missingness <- c()

for (i in study_list) {
  s <- filter(dat, study == i)
  missing_na <- sum(is.na(s$na)) / length(s$na)
  missing_pa <- sum(is.na(s$pa)) /length(s$pa)
  
  if(missing_na > .9 | missing_pa > .9 ) {
    print(i)
    missingness[length(missingness) + 1] <- i
  }
}

missingness

#looks like studies 2 and 4 don't meet the missingness criterion of data for at
#least 10% of cases, so exclude these studies
#for now, also exclude study 9 due to convergence issues

dat <- dat %>% 
  filter(study != 2 & study != 4 & study != 9)

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

##################### ANALYSE ##################

#set up a container for the results
names <- c("study", "status", "beta_na", "beta_na_ci_low", "beta_na_ci_hi",
           "beta_pa", "beta_pa_ci_low", "beta_pa_ci_hi",
           "t_na", "t_pa", "p_na", "p_pa")

res <- data.frame(matrix(ncol = length(names), nrow = 0))
colnames(res) <- names

#iterate across all unique studies in the dataset
study_list <- dat$study %>% unique()
study_list

#increase limit for iterations until convergence is expected
options <- lmeControl(maxIter = 10000)

for (study_no in study_list) {
  
  #iterate across all status groups within this study
  status_list <- filter(dat, study == study_no)$status %>% unique
  status_list
  
  for (status_no in status_list) {
    
      this_dat <- dat %>% 
        filter(study == study_no & status == status_no)
      
      #fit a linear model predicting PA from NA
      model_pa <-lme(pa ~ na, random =~ na | subjno, data=this_dat, na.action=na.omit, control = options)
      sum_pa <- summary(model_pa)
      int_pa <- intervals(sum_pa)
      
      #fit a multilevel model predicting NA from PA
      model_na <-lme(na ~ pa, random =~ pa | subjno, data=this_dat, na.action=na.omit, control = options)
      sum_na <- summary(model_na)
      int_na <- intervals(sum_na)
      
      #extract parameters
      res <- res %>% 
        add_row(study = study_no, status = status_no, 
                      beta_na = int_na$fixed[2, 2], beta_na_ci_low = int_na$fixed[2, 1], beta_na_ci_hi = int_na$fixed[2, 3],
                      beta_pa = int_pa$fixed[2, 2], beta_pa_ci_low = int_pa$fixed[2, 1], beta_pa_ci_hi = int_pa$fixed[2, 3],
                      t_na = sum_na$tTable[2, 4], t_pa = sum_pa$tTable[2, 4],
                      p_na = sum_na$tTable[2, 5], p_pa = sum_pa$tTable[2, 5]
          )
      }#status_no
}#study_no

#plot data

#the relationship between the beta params should look linear when plotted like this
plot(res$beta_na, 1/res$beta_pa)

#the relationship between t-values should be roughly linear
plot(res$t_na, res$t_pa)

#same for p-values
plot(res$p_na, res$p_pa)


##################### ANALYSE FIXED SAMPLE SIZE ##########################

#barring one study, all p-values are highly significant, suggesting that the direction
#of measurement did not play a huge role here
#This may be an issue of power - most studies had a relatively large sample size
#re-investigate by drawing a sub-set of participants and beeps from each study
#Attention - This section can produce errors, as the models don't always converge

sample_size = 10
beep_n = 20

#set up a container for the results
names <- c("study", "status", "beta_na", "beta_na_ci_low", "beta_na_ci_hi",
           "beta_pa", "beta_pa_ci_low", "beta_pa_ci_hi",
           "t_na", "t_pa", "p_na", "p_pa")

res <- data.frame(matrix(ncol = length(names), nrow = 0))
colnames(res) <- names

#iterate across all unique studies in the dataset
study_list <- dat$study %>% unique()
study_list

#increase limit for iterations until convergence is expected
options <- lmeControl(maxIter = 10000000)

for (study_no in study_list) {
  
  #iterate across all status groups within this study
  status_list <- filter(dat, study == study_no)$status %>% unique
  status_list
  
  for (status_no in status_list) {
      
    this_dat <- dat %>% 
      filter(study == study_no & status == status_no)
    
    #now sample a sub-set of participants
    participant_list <- this_dat %>% 
      unique() %>% 
      sample_n(sample_size) %>% 
      select(subjno)
    
    #select only those participants who are in the list
    #for these participants, consider only the beep_n first beeps
    this_dat <- merge(this_dat, participant_list, by.x = c("subjno")) %>% 
      filter(
        beep <= beep_n & !is.na(pa) & !is.na(na)
      )
    
    #fit a linear model predicting PA from NA
    model_pa <-lme(pa ~ na, random =~ na | subjno, data=this_dat, na.action=na.omit, control = options)
    sum_pa <- summary(model_pa)
    int_pa <- intervals(sum_pa)
    
    #fit a multilevel model predicting NA from PA
    model_na <-lme(na ~ pa, random =~ pa | subjno, data=this_dat, na.action=na.omit, control = options)
    sum_na <- summary(model_na)
    int_na <- intervals(sum_na)
    
    #extract parameters
    res <- res %>% 
      add_row(study = study_no, status = status_no, 
              beta_na = int_na$fixed[2, 2], beta_na_ci_low = int_na$fixed[2, 1], beta_na_ci_hi = int_na$fixed[2, 3],
              beta_pa = int_pa$fixed[2, 2], beta_pa_ci_low = int_pa$fixed[2, 1], beta_pa_ci_hi = int_pa$fixed[2, 3],
              t_na = sum_na$tTable[2, 4], t_pa = sum_pa$tTable[2, 4],
              p_na = sum_na$tTable[2, 5], p_pa = sum_pa$tTable[2, 5]
      )
  }#status_no
}#study_no

#This loop should return the requested data
#however, as only a sub-set of beeps and participants is selected,
#different results may turn up when re-running it.
#some samples won't converge when sample or beep sizes are too low



#########PLOTTING ################

#load example data from the above code, where the respective models converged
#(beep size = 20, sample size = 10)
#some samples are missing, as they did not converge
#Note: I've purposefully picked a run which illustrates the problem of
#diverging p-values!

res <- read_csv(file = "results_full_model.csv")

#the relationship between the beta params should look linear when plotted like this
plot(res$beta_na, 1/res$beta_pa)

#the relationship between t-values should be roughly linear
plot(res$t_na, res$t_pa)

#same for p-values
plot(res$p_na, res$p_pa)
     
#transform example data for plotting

res <- res %>% 
  mutate(
    problem = ifelse(p_na > .05 & p_pa < .05 |p_pa > .05 & p_na < .05, 
                    "Yes", "No"),
    label = ifelse(problem == "Yes", paste("Study", study, "\nStatus", status), "")
  )

plot <- ggplot(data = res, aes(x = p_na, y = p_pa)) +
                  geom_point(aes(x = p_na, y = p_pa, color = problem),
                             show.legend = FALSE) +
                  geom_text(aes(label = label, hjust = -0.1)) +
                  geom_hline(yintercept = 0.05) +
                  geom_vline(xintercept = 0.05) +
                  xlab("Negative Affect as Dependent Variable") +
                  ylab("Positive Affect as Dependent Variable") +
                  theme_bw()
plot



##################### ANALYSIS WITH OUR NEW TRICK ##################

#add standardized values to the dataset


dat$zpa <- ave(dat$pa, dat$subjno, FUN = function(x) (x - mean(x, na.rm=TRUE)) / sd(x, na.rm=TRUE))
dat$zna <- ave(dat$na, dat$subjno, FUN = function(x) (x - mean(x, na.rm=TRUE)) / sd(x, na.rm=TRUE))


#set up a container for the results
names <- c("study", "status", "beta_na", "beta_na_ci_low", "beta_na_ci_hi",
           "beta_pa", "beta_pa_ci_low", "beta_pa_ci_hi",
           "t_na", "t_pa", "p_na", "p_pa", 
           "r", "r_ci_low", "r_ci_hi",
           "r_corr", "r_corr_ci_low", "r_corr_ci_hi"
           )

res <- data.frame(matrix(ncol = length(names), nrow = 0))
colnames(res) <- names

#iterate across all unique studies in the dataset
study_list <- dat$study %>% unique()
study_list

#increase limit for iterations until convergence is expected
options <- lmeControl(maxIter = 10000)

for (study_no in study_list) {
  
  #iterate across all status groups within this study
  status_list <- filter(dat, study == study_no)$status %>% unique
  status_list
  
  for (status_no in status_list) {
    
    this_dat <- dat %>% 
      filter(study == study_no & status == status_no)
    
    #FULL MODEL
    
    #fit a linear model predicting PA from NA
    model_pa <-lme(pa ~ na, random =~ na | subjno, data=this_dat, na.action=na.omit, control = options)
    sum_pa <- summary(model_pa)
    int_pa <- intervals(sum_pa)
    
    #fit a multilevel model predicting NA from PA
    model_na <-lme(na ~ pa, random =~ pa | subjno, data=this_dat, na.action=na.omit, control = options)
    sum_na <- summary(model_na)
    int_na <- intervals(sum_na)
    
    #CENTRED MODEL, no random effects
    
    #both models are computed here,but their results are equivalent
    #barring divergences due to rounding error
    #by default, only values for the first model (na as DV) are saved
    model_na_cen <- lm(zna ~ 0 + zpa, data=this_dat)
    sum_na_cen <- summary(model_na_cen)
    
    point_na_cen <- sum_na_cen$coefficients[1]
    ci_na_cen <- confint(model_na_cen)
    
    model_pa_cen <- lm(zpa ~ 0 + zna, data=this_dat)
    sum_pa_cen <- summary(model_pa_cen)
    
    #Centred Model with random effects
    
    model_na_cen_alt <- lme(zna ~ 0 + zpa, random = ~ 0 + zpa | subjno, data=this_dat, na.action=na.omit)
    sum_na_cen_alt <- summary(model_na_cen_alt)
    point_na_cen_alt <- sum_na_cen_alt$coefficients[1]
    ci_na_cen_alt <- intervals(model_na_cen_alt)
    
    #extract parameters
    res <- res %>% 
      add_row(study = study_no, status = status_no, 
              beta_na = int_na$fixed[2, 2], beta_na_ci_low = int_na$fixed[2, 1], beta_na_ci_hi = int_na$fixed[2, 3],
              beta_pa = int_pa$fixed[2, 2], beta_pa_ci_low = int_pa$fixed[2, 1], beta_pa_ci_hi = int_pa$fixed[2, 3],
              t_na = sum_na$tTable[2, 4], t_pa = sum_pa$tTable[2, 4],
              p_na = sum_na$tTable[2, 5], p_pa = sum_pa$tTable[2, 5],
              r = point_na_cen[1], r_ci_low = ci_na_cen[1], r_ci_hi = ci_na_cen[2],
              r_corr = point_na_cen_alt$fixed[1], r_corr_ci_low = ci_na_cen_alt$fixed[1],r_corr_ci_hi = ci_na_cen_alt$fixed[3]
      )
  }#status_no
}#study_no

#plot data

#the corrected and uncorrected r-values should be almost perfectly linearly related
plot(res$r, res$r_corr)


#plot results

#format
out_data <- res %>% 
  pivot_longer(cols = c("r", "r_corr"), names_to = "form", values_to = "mean")

out_data2 <- res %>%
  pivot_longer(cols = c("r_ci_low", "r_corr_ci_low"), names_to = "form", values_to = "lower") %>% 
  select("lower")

out_data3 <- res %>%
  pivot_longer(cols = c("r_ci_hi", "r_corr_ci_hi"), names_to = "form", values_to = "upper") %>% 
  select("upper")

out_data <- out_data %>%
  cbind(out_data2, out_data3) %>% 
  select(study, status, form, mean, lower, upper) %>% 
  mutate(
    status_str = ifelse(
      status == 0, "Control",
      ifelse(status == 1, "Psychosis Risk", 
             ifelse(status == 2, "Psychosis", "Depression"))
    ),
    this_group = paste("Study", study, status_str)
  )

plot <- ggplot(data = out_data) +
  facet_wrap(~this_group) +
  geom_bar(stat = "identity", aes(x = form, y = mean)) +
  geom_errorbar(aes(x = form, ymin = lower, ymax = upper)) +
  ylab("Mean Correlation") +
  xlab("Sample") +
  theme_bw()
  
plot
