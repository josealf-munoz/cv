# All together
library("data.table")
library("panelr")
library("evd")
library("tidyverse")
library("numDeriv")
library("doParallel")
library("ucminf")
library("optimx")
library("spatstat")
library("Brobdingnag")

load("data_anubis.Rda")
data = data_anubis
data <- data.table(data)
data = data[, list(id, age, choice_, choice_present, lwage, wage_scale, wage_, urban, g_, expe_, more_eighteen, more_twenty, expe_2, age_2, # utilities
                   attrition_, attrit_present, couple_, one_kids, two_kids, distance_6eme_recent, # attrition
                   both_parents_french, one_parent_french,  # variables at 16
                   late_school, occup_father_0, occup_father_1,
                   occup_father_2, occup_mother_0, occup_mother_1,
                   occup_mother_2)]

# # data = data[id==1]
# # imin=imax=1
n=length(unique(data$id))
imin=1
imax=n
data[,id := as.factor(id)]

# data = data[id==4]
# data[,id:=1]



# Fix parameters
discount = 0.95
euler = 0.5772156649


param_schooling = c("urban","more_eighteen", "more_twenty")
nb_schooling = length(param_schooling)

param_work = c("urban","g_", "expe_", "expe_2")
nb_work = length(param_work)

param_home = c("age","age*age")
nb_home = length(param_home)

param_initial = c("both_parents_french", "one_parent_french", "late_school", 
                  "occup_father_1", "occup_father_2", 
                  "occup_mother_1", "occup_mother_2", "constant")
nb_initial = length(param_initial)

types = 2*3
distribution = 1 + 1 
nb_vip = types+distribution

n_params = nb_vip + nb_schooling + nb_work + nb_home + nb_initial

param = c(0.300 ,1.75, -1.06,
          0.300 ,1.75, -1.06,
          0, 4, # distributions parameters
          -1.76,  -1.48, 0.278, # schooling: 18, 20, urban
          0.051, 0.057, -0.002, .075, # wage: g_, expe_, expe^2, urban
          0.08, -0.005, # home: age, age^2
          -0.144,0.0260,-0.317,-0.132,0.078,-0.129,0.007, 1
) 



complete_likelihood = function(param, imin, imax, data){
  
  alpha_1_1 = param[1]
  alpha_2_1 = param[2]
  alpha_3_1 = param[3]
  
  alpha_1_2 = param[4]
  alpha_2_2 = param[5]
  alpha_3_2 = param[6]
  
  sigma = param[7]
  sd_wages = exp(sigma)
  # integration_eta = 1
  integration_eta = as.numeric(gauss.hermite(function(x) exp(x), mu = 0, sd = sd_wages))
  lambda = param[8]
  # lambda = 1
  
  n_0 = (nb_vip) + 1
  n_1 = n_0 + nb_schooling - 1
  beta = param[n_0:n_1] #school
  
  n_2 = n_1 + 1
  n_3 = n_2 + nb_work - 1
  gamma = param[n_2:n_3] #wage
  
  n_4 = n_3 + 1
  n_5 = n_4 + nb_home - 1
  delta = param[n_4:n_5] #home
  
  n_6 = n_5 + 1
  n_7 = n_6 + nb_initial - 1
  omega = param[n_6:n_7] #initial variables
  
  # Initial conditions----------------------------------------------------------------------------------------------
  data[,theta_1_1 := alpha_1_1]
  
  data[,theta_2_1 := alpha_2_1]
  
  data[,theta_3_1 := alpha_3_1]
  
  data[,theta_1_2 := alpha_1_2]
  
  data[,theta_2_2 := alpha_2_2]
  
  data[,theta_3_2 := alpha_3_2]
  
  
  ## Value functions--------------------------------------------------------------------
  
  
  # Initial utilities-------------------------------------------------------------------
  
  # For the value functions I don't include age related variables
  
  data[,utility_school_1 := beta[3]*urban + theta_1_1]
  data[,utility_work_1 := exp(gamma[4]*urban + theta_2_1)]
  data[,utility_home_1 :=  theta_3_1]
  
  data[,utility_school_2 := beta[3]*urban + theta_1_2]
  data[,utility_work_2 := exp(gamma[4]*urban + theta_2_2)]
  data[,utility_home_2 :=  theta_3_2]
  
  # Solution of model--------------------------------------------------------------------
  
  # create SST
  g = 10:24
  expe = 0:16
  age = 16:32
  id = imin:imax
  SST_full = CJ(id=id, g=g, expe=expe, age=age)
  
  # Start at age 32
  SST=SST_full
  setkey(SST,id,g,expe)
  SST[, `:=` (sum = g + expe, id = as.factor(id))] # Add sum and convert id to factor
  SST = SST[sum <= 26 & age==32] #eliminate combinations that can't enter my model
  
  # include the information from the real data_1
  data_SST <- data[age == 32, .(id, utility_work_1, utility_home_1, utility_work_2, utility_home_2)]
  setkey(data_SST, id)
  SST_32 <- SST[data_SST, on = "id"] #keyed join
  
  SST_32[, `:=`(expe_new = expe,
                age_new = age,
                sum_work_1_33 = 0,
                sum_home_1_33 = 0)]
  
  for (t in 1:32) {
    SST_32[, expe_new := expe_new + 1]
    SST_32[, utility_work_new_1 := utility_work_1*integration_eta*exp(gamma[1]*g + gamma[2]*expe_new + gamma[3]*(expe_new*expe_new))]
    SST_32[, sum_work_1_33 := sum_work_1_33 + discount^(t)*utility_work_new_1]
    
    SST_32[, age_new := age_new + 1]
    SST_32[, utility_home_new_1 := utility_home_1 + delta[1]*age_new + delta[2]*age_new*age_new]
    SST_32[, sum_home_1_33 := sum_home_1_33 + discount^(t)*utility_home_new_1]
  }
  
  
  SST_32[, `:=`(expe_new = expe,
                age_new = age,
                sum_work_2_33 = 0,
                sum_home_2_33 = 0)]
  
  for (t in 1:32) {
    SST_32[, expe_new := expe_new + 1]
    SST_32[, utility_work_new_2 := utility_work_2*integration_eta*exp(gamma[1]*g + gamma[2]*expe_new + gamma[3]*(expe_new*expe_new))]
    SST_32[, sum_work_2_33 := sum_work_2_33 + discount^(t)*utility_work_new_2]
    
    SST_32[, age_new := age_new + 1]
    SST_32[, utility_home_new_2 := utility_home_2 + delta[1]*age_new + delta[2]*age_new*age_new]
    SST_32[, sum_home_2_33 := sum_home_2_33 + discount^(t)*utility_home_new_2]
  }
  
  # Create the utility until age 32
  SST_32[, `:=`(value_function_work_1 = utility_work_1*integration_eta*exp(gamma[1]*g + gamma[2]*expe + gamma[3]*(expe*expe)) + discount*sum_work_1_33,
                value_function_home_1 = utility_home_1 + delta[1]*age + delta[2]*age*age  + discount*sum_home_1_33,
                value_function_work_2 = utility_work_2*integration_eta*exp(gamma[1]*g + gamma[2]*expe + gamma[3]*(expe*expe)) + discount*sum_work_2_33,
                value_function_home_2 = utility_home_2 + delta[1]*age + delta[2]*age*age  + discount*sum_home_2_33)]
  
  
  # For period T
  
  # Now, for each state point I can obtain the expected value function to use in 
  # the backward induction
  
  SST_32[, `:=`(value_function_1_32 = as.numeric(lambda*log(brob(value_function_work_1/lambda)+brob(value_function_home_1/lambda)) + lambda*euler),
                value_function_2_32 = as.numeric(lambda*log(brob(value_function_work_2/lambda)+brob(value_function_home_2/lambda)) + lambda*euler))]
  
  SST_33 = SST_32[, c("id", "age", "g", "expe", "sum_work_1_33", "sum_work_2_33","sum_home_1_33", "sum_home_2_33" )]
  
  SST_32 = SST_32[, c("id", "age", "g", "expe", "value_function_2_32", "value_function_1_32")]
  setkey(SST,id,g,expe)
  
  SST_30_32 = SST_32
  
  counter = 0
  
  for (t in (32-1):30){
    counter = counter + 1 #counter keeps eliminating the state spaces that are not feasible at each age
    
    # Set SST
    SST_t=SST_full
    setkey(SST,id,g,expe)
    SST_t[, `:=` (sum = g + expe, id = as.factor(id))] # Add sum and convert id to factor
    SST_t = SST_t[sum <= 26-counter & age==t] #eliminate combinations that can't enter my model
    setkey(SST_t,id)
    
    
    data_t = data[age==t, c("id", "utility_work_1", "utility_home_1", "utility_work_2", "utility_home_2")]
    setkey(data_t,id)
    
    
    # Work
    data_SST_t_work = data_t[SST_t, nomatch=0]
    data_SST_t_work[, `:=` ( utility_work_1 = utility_work_1*integration_eta*exp(gamma[1]*g + gamma[2]*expe + gamma[3]*(expe*expe)),
                             utility_home_1 = utility_home_1 + delta[1]*age + delta[2]*age*age,
                             utility_work_2 = utility_work_2*integration_eta*exp(gamma[1]*g + gamma[2]*expe + gamma[3]*(expe*expe)),
                             utility_home_2 = utility_home_2 + delta[1]*age + delta[2]*age*age,
                             # Here I will merge the decision of working at age t to the state at t+1 given
                             # given he worked at t
                             g = g + 0,
                             expe = expe + 1,
                             age = age + 1)]
    setkey(data_SST_t_work, id, g, expe)
    
    
    data_SST_t_work <- data_SST_t_work[SST_30_32, on = c("id","age", "g","expe")] #keyed join
    data_SST_t_work = data_SST_t_work[ !is.na(sum)] #eliminate combinations that can't enter my model
    
    data_SST_t_work[, `:=` ( expe = expe - 1,
                             age = age - 1)] #extract the experience gained to return to what is seen in the data_1
    next_period = t+1
    
    
    ## Here it gets messy. The value of working is the utility at time t plus the
    ## discount value function at t+1. Because of the merge above, each row has
    ## been matched with it
    data_SST_t_work = data_SST_t_work
    data_SST_t_work[, `:=` (
      value_work_1 = utility_work_1 + discount * get(paste0("value_function_1_", next_period)),
      value_work_2 = utility_work_2 + discount * get(paste0("value_function_2_", next_period))
    )]
    
    variables_t = c("id", "age", "g", "expe", "value_work_1", "value_work_2")
    data_SST_t_work = data_SST_t_work[, ..variables_t]
    setkey(data_SST_t_work, id, g, expe)
    
    # home
    data_SST_t_home = data_t[SST_t, nomatch=0]
    data_SST_t_home[, `:=` ( utility_work_1 = utility_work_1*integration_eta*exp(gamma[1]*g + gamma[2]*expe + gamma[3]*(expe*expe)),
                             utility_home_1 = utility_home_1 + delta[1]*age + delta[2]*age*age,
                             utility_work_2 = utility_work_2*integration_eta*exp(gamma[1]*g + gamma[2]*expe + gamma[3]*(expe*expe)),
                             utility_home_2 = utility_home_2 + delta[1]*age + delta[2]*age*age,
                             # Here I will merge the decision of homeing at age t to the state at t+1 given
                             # given he homeed at t
                             g = g + 0,
                             expe = expe + 0,
                             age = age + 1)]
    setkey(data_SST_t_home, id, age, g, expe)
    
    
    data_SST_t_home <- data_SST_t_home[SST_30_32, on = c("id","age", "g","expe")] #keyed join
    data_SST_t_home = data_SST_t_home[ !is.na(sum)] #eliminate combinations that can't enter my model
    
    data_SST_t_home[, `:=` ( expe = expe - 0,
                             age = age - 1)] #extract the experience gained to return to what is seen in the data_1
    next_period = t+1
    
    
    ## Here it gets messy. The value of homeing is the utility at time t plus the
    ## discount value function at t+1. Because of the merge above, each row has
    ## been matched with it
    data_SST_t_home = data_SST_t_home
    data_SST_t_home[, `:=` (
      value_home_1 = utility_home_1 + discount * get(paste0("value_function_1_", next_period)),
      value_home_2 = utility_home_2 + discount * get(paste0("value_function_2_", next_period))
    )]
    
    variables_t = c("id", "age", "g", "expe", "value_home_1", "value_home_2")
    data_SST_t_home = data_SST_t_home[, ..variables_t]
    setkey(data_SST_t_home, id, age, g, expe)
    
    
    # Merge
    ## At time t the state space of both data_1 tables above has to be the same, 
    ## hence there shouldn't be any problem in the merge.
    data_SST_t_all = data_SST_t_work[data_SST_t_home, nomatch=0, on=c("id", "age", "g","expe")]
    #data_1_SST_t_all = merge(data_1_SST_t_work, data_1_SST_t_home, by=c("id","g","expe"))
    
    data_SST_t_all[, c(paste0("value_function_1_", t), paste0("value_function_2_", t)) := .(
      as.numeric(lambda * log(brob(value_work_1 / lambda) + brob(value_home_1 / lambda)) + lambda * euler),
      as.numeric(lambda * log(brob(value_work_2 / lambda) + brob(value_home_2 / lambda)) + lambda * euler)
    )]
    
    variables_t = c("id", "age", "g", "expe", paste0("value_function_1_",t), paste0("value_function_2_",t))
    data_SST_t_all = data_SST_t_all[, ..variables_t]
    
    
    
    SST_30_32 = rbindlist(list(SST_30_32, data_SST_t_all), fill=TRUE)
  }
  
  
  
  SST_17_29 = SST_30_32[age==30]
  variables_t = c("id", "age", "g", "expe", "value_function_1_30", "value_function_2_30")
  SST_17_29 = SST_17_29[, ..variables_t]
  
  for (t in 29:17){
    counter = counter + 1 #counter keeps eliminating the state spaces that are not feasible at each age
    
    # Set SST
    SST_t=SST_full
    setkey(SST,id,g,expe)
    SST_t[, `:=` (sum = g + expe, 
                  id = as.factor(id),
                  more_eighteen = ifelse(age>=18,1,0),
                  more_twenty = ifelse(age>=20,1,0))] # Add sum and convert id to factor
    SST_t = SST_t[sum <= 26-counter & age==t] #eliminate combinations that can't enter my model
    setkey(SST_t,id)
    
    
    data_t = data[age==t, c("id", "utility_school_1", "utility_work_1", "utility_home_1", "utility_school_2", "utility_work_2", "utility_home_2")]
    setkey(data_t,id)
    
    # educ
    data_SST_t_educ = data_t[SST_t, nomatch=0]
    data_SST_t_educ[, `:=` ( utility_school_1 = utility_school_1 + beta[1]*more_eighteen + beta[2]*more_twenty,
                             utility_work_1 = utility_work_1*integration_eta*exp(gamma[1]*g + gamma[2]*expe + gamma[3]*(expe*expe)),
                             utility_home_1 = utility_home_1 + delta[1]*age + delta[2]*age*age,
                             utility_school_2 = utility_school_2 + beta[1]*more_eighteen + beta[2]*more_twenty,
                             utility_work_2 = utility_work_2*integration_eta*exp(gamma[1]*g + gamma[2]*expe + gamma[3]*(expe*expe)),
                             utility_home_2 = utility_home_2 + delta[1]*age + delta[2]*age*age,
                             # Here I will merge the decision of educing at age t to the state at t+1 given
                             # given he educed at t
                             g = g + 1,
                             expe = expe + 0,
                             age = age + 1)]
    setkey(data_SST_t_educ, id, g, expe)
    
    
    data_SST_t_educ <- data_SST_t_educ[SST_17_29, on = c("id","age", "g","expe")] #keyed join
    data_SST_t_educ = data_SST_t_educ[ !is.na(sum)] #eliminate combinations that can't enter my model
    
    data_SST_t_educ[, `:=` ( g = g - 1,
                             expe = expe - 0,
                             age = age - 1)] #extract the experience gained to return to what is seen in the data_1
    next_period = t+1
    
    
    ## Here it gets messy. The value of educing is the utility at time t plus the
    ## discount value function at t+1. Because of the merge above, each row has
    ## been matched with it
    data_SST_t_educ = data_SST_t_educ
    data_SST_t_educ[, `:=` (
      value_educ_1 = utility_school_1 + discount * get(paste0("value_function_1_", next_period)),
      value_educ_2 = utility_school_2 + discount * get(paste0("value_function_2_", next_period))
    )]
    
    variables_t = c("id", "age", "g", "expe", "value_educ_1", "value_educ_2")
    data_SST_t_educ = data_SST_t_educ[, ..variables_t]
    setkey(data_SST_t_educ, id, g, expe)
    
    
    # Work
    data_SST_t_work = data_t[SST_t, nomatch=0]
    data_SST_t_work[, `:=` ( utility_school_1 = utility_school_1 + beta[1]*more_eighteen + beta[2]*more_twenty,
                             utility_work_1 = utility_work_1*integration_eta*exp(gamma[1]*g + gamma[2]*expe + gamma[3]*(expe*expe)),
                             utility_home_1 = utility_home_1 + delta[1]*age + delta[2]*age*age,
                             utility_school_2 = utility_school_2 + beta[1]*more_eighteen + beta[2]*more_twenty,
                             utility_work_2 = utility_work_2*integration_eta*exp(gamma[1]*g + gamma[2]*expe + gamma[3]*(expe*expe)),
                             utility_home_2 = utility_home_2 + delta[1]*age + delta[2]*age*age,
                             # Here I will merge the decision of working at age t to the state at t+1 given
                             # given he worked at t
                             g = g + 0,
                             expe = expe + 1,
                             age = age + 1)]
    setkey(data_SST_t_work, id, g, expe)
    
    
    data_SST_t_work <- data_SST_t_work[SST_17_29, on = c("id","age", "g","expe")] #keyed join
    data_SST_t_work = data_SST_t_work[ !is.na(sum)] #eliminate combinations that can't enter my model
    
    data_SST_t_work[, `:=` ( g = g - 0,
                             expe = expe - 1,
                             age = age - 1)] #extract the experience gained to return to what is seen in the data_1
    next_period = t+1
    
    
    ## Here it gets messy. The value of working is the utility at time t plus the
    ## discount value function at t+1. Because of the merge above, each row has
    ## been matched with it
    data_SST_t_work = data_SST_t_work
    data_SST_t_work[, `:=` (
      value_work_1 = utility_work_1 + discount * get(paste0("value_function_1_", next_period)),
      value_work_2 = utility_work_2 + discount * get(paste0("value_function_2_", next_period))
    )]
    
    variables_t = c("id", "age", "g", "expe", "value_work_1", "value_work_2")
    data_SST_t_work = data_SST_t_work[, ..variables_t]
    setkey(data_SST_t_work, id, g, expe)
    
    
    # home
    data_SST_t_home = data_t[SST_t, nomatch=0]
    data_SST_t_home[, `:=` ( utility_school_1 = utility_school_1 + beta[1]*more_eighteen + beta[2]*more_twenty,
                             utility_work_1 = utility_work_1*integration_eta*exp(gamma[1]*g + gamma[2]*expe + gamma[3]*(expe*expe)),
                             utility_home_1 = utility_home_1 + delta[1]*age + delta[2]*age*age,
                             utility_school_2 = utility_school_2 + beta[1]*more_eighteen + beta[2]*more_twenty,
                             utility_work_2 = utility_work_2*integration_eta*exp(gamma[1]*g + gamma[2]*expe + gamma[3]*(expe*expe)),
                             utility_home_2 = utility_home_2 + delta[1]*age + delta[2]*age*age,
                             # Here I will merge the decision of homeing at age t to the state at t+1 given
                             # given he homeed at t
                             g = g + 0,
                             expe = expe + 0,
                             age = age + 1)]
    setkey(data_SST_t_home, id, age, g, expe)
    
    
    data_SST_t_home <- data_SST_t_home[SST_17_29, on = c("id","age", "g","expe")] #keyed join
    data_SST_t_home = data_SST_t_home[ !is.na(sum)] #eliminate combinations that can't enter my model
    
    data_SST_t_home[, `:=` ( g = g - 0,
                             expe = expe - 0,
                             age = age - 1)] #extract the experience gained to return to what is seen in the data_1
    next_period = t+1
    
    
    ## Here it gets messy. The value of homeing is the utility at time t plus the
    ## discount value function at t+1. Because of the merge above, each row has
    ## been matched with it
    data_SST_t_home = data_SST_t_home
    data_SST_t_home[, `:=` (
      value_home_1 = utility_home_1 + discount * get(paste0("value_function_1_", next_period)),
      value_home_2 = utility_home_2 + discount * get(paste0("value_function_2_", next_period))
    )]
    
    variables_t = c("id", "age", "g", "expe", "value_home_1", "value_home_2")
    data_SST_t_home = data_SST_t_home[, ..variables_t]
    setkey(data_SST_t_home, id, age, g, expe)
    
    # Merge
    ## At time t the state space of both data_1 tables above has to be the same, 
    ## hence there shouldn't be any problem in the merge.
    data_SST_t_all = data_SST_t_educ[data_SST_t_work, nomatch=0, on=c("id", "age", "g","expe")]
    data_SST_t_all = data_SST_t_all[data_SST_t_home, nomatch=0, on=c("id", "age", "g","expe")]
    
    data_SST_t_all[, c(paste0("value_function_1_", t), paste0("value_function_2_", t)) := .(
      as.numeric(lambda*log(brob(value_educ_1/lambda)+brob(value_work_1/lambda)+brob(value_home_1/lambda)) + lambda*euler ),
      as.numeric(lambda*log(brob(value_educ_2/lambda)+brob(value_work_2/lambda)+brob(value_home_2/lambda)) + lambda*euler )
    )]
    
    variables_t = c("id", "age", "g", "expe", paste0("value_function_1_",t), paste0("value_function_2_",t))
    data_SST_t_all = data_SST_t_all[, ..variables_t]
    
    SST_17_29 = rbindlist(list(SST_17_29, data_SST_t_all), fill=TRUE)
  }
  
  SST_17_29 = SST_17_29[age<=29]
  SST_all = rbindlist(list(SST_30_32, SST_17_29), fill=TRUE)
  setkey(SST_all, id, g, expe)
  
  SST_all[, identifier := row_number(id)] #create an identifier for each row
  setkey(SST_all, identifier)
  
  # Define the variable names for the value functions
  var_names <- paste0("value_function_1_", 17:32)
  
  # Melt the data.table and subset to the relevant columns
  SST_melted_1 <- melt(SST_all, id.vars = "identifier", measure.vars = var_names,
                       variable.name = "age", value.name = "value_function_1")[!is.na(value_function_1),
                                                                               .(identifier, value_function_1)]
  # Define the variable names for the value functions
  var_names <- paste0("value_function_2_", 17:32)
  
  # Melt the data.table and subset to the relevant columns
  SST_melted_2 <- melt(SST_all, id.vars = "identifier", measure.vars = var_names,
                       variable.name = "age", value.name = "value_function_2")[!is.na(value_function_2),
                                                                               .(identifier, value_function_2)]
  SST_melted = SST_melted_1[SST_melted_2, on=c("identifier")]
  setkey(SST_melted, identifier)
  
  #SST_all = SST_all[SST_2, nomatch=0] #match this new variable to the previous one
  SST_all = SST_all[SST_melted, on=c("identifier")]
  SST_all = SST_all[, c("id","age", "g","expe","value_function_1", "value_function_2")] #eliminate all the variables of value function, whose information is grouped into one variable now
  
  
  # Define functions for renaming columns and setting keys
  rename_cols <- function(df, old_names, new_names) {
    setnames(df, old_names, new_names, skip_absent = TRUE)
  }
  
  set_key <- function(df, key_cols) {
    setkeyv(df, key_cols)
  }
  
  # Rename columns and set key for SST_all
  rename_cols(SST_all, c("age", "g", "expe"), c("age_next", "g_next", "expe_next"))
  set_key(SST_all, c("id", "age_next", "g_next", "expe_next"))
  
  # Rename columns and set key for SST_33
  rename_cols(SST_33, c("age", "g", "expe"), c("age_next", "g_next", "expe_next"))
  set_key(SST_33, c("id", "age_next", "g_next", "expe_next"))
  
  # Chain data manipulation operations for SST_33
  SST_33[, age_next := age_next + 1][, expe_next := expe_next + 1]
  SST_33_work <- SST_33[, .(id, age_next, g_next, expe_next, sum_work_1_33, sum_work_2_33)]
  SST_33_work = SST_33_work[, expe_next := expe_next+1]
  
  SST_33_home <- SST_33[, .(id, age_next, g_next, expe_next, sum_home_1_33, sum_home_2_33)]
  
  ### Value function if he studies
  data_educ = data[age<=29]
  data_educ[, age_next := age + 1]
  data_educ[, g_next := g_+1]
  data_educ[, expe_next := expe_+0]
  setkey(data_educ, id, age_next, g_next, expe_next)
  # Merge. Keep rows in data_1 because otherwise the attrition rows will be left out.
  data_educ = merge(data_educ, SST_all, by=c("id","age_next", "g_next","expe_next"), all.x = TRUE)
  setnames(data_educ, "value_function_1", "value_function_educ_1")
  setnames(data_educ, "value_function_2", "value_function_educ_2")
  data_educ = data_educ[,c("id", "age", "g_", "expe_", "value_function_educ_1", "value_function_educ_2")]
  setkey(data_educ, id, age, g_, expe_)
  
  ### Value function if he works
  data_work = data
  data_work[, age_next := age + 1]
  data_work[, g_next := g_+0]
  data_work[, expe_next := expe_+1]
  setkey(data_work, id, age_next, g_next, expe_next)
  # Merge. Keep rows in data_1 because otherwise the attrition rows will be left out.
  data_work = merge(data_work, SST_all, by=c("id","age_next", "g_next","expe_next"), all.x = TRUE)
  data_work = merge(data_work, SST_33_work, by=c("id","age_next", "g_next","expe_next"), all.x = TRUE)
  data_work = data_work[,`:=`( value_function_1 = ifelse(!is.na(value_function_1), value_function_1, sum_work_1_33),
                               value_function_2 = ifelse(!is.na(value_function_2), value_function_2, sum_work_2_33)
  )]
  setnames(data_work, "value_function_1", "value_function_work_1")
  setnames(data_work, "value_function_2", "value_function_work_2")
  data_work = data_work[,-c("age_next", "g_next", "expe_next", "sum_work_1_33", "sum_work_2_33")]
  setkey(data_work, id, age, g_, expe_)
  
  ### Value function if he stays at home
  data_home = data
  data_home[, age_next := age + 1]
  data_home[, g_next := g_+0]
  data_home[, expe_next := expe_+0]
  setkey(data_home, id, age_next, g_next, expe_next)
  # Merge. Keep rows in data_1 because otherwise the attrition rows will be left out.
  data_home = merge(data_home, SST_all, by=c("id","age_next", "g_next","expe_next"), all.x = TRUE)
  data_home = merge(data_home, SST_33_home, by=c("id","age_next", "g_next","expe_next"), all.x = TRUE)
  data_home = data_home[,`:=`( value_function_1 = ifelse(!is.na(value_function_1), value_function_1, sum_home_1_33),
                               value_function_2 = ifelse(!is.na(value_function_2), value_function_2, sum_home_2_33)
  )]
  setnames(data_home, "value_function_1", "value_function_home_1")
  setnames(data_home, "value_function_2", "value_function_home_2")
  data_home = data_home[,c("id", "age", "g_", "expe_", "value_function_home_1", "value_function_home_2")]
  setkey(data_home, id, age, g_, expe_)
  
  ### Merge all three data
  data_1_1 = data_educ[data_work, nomatch=NA]
  setkey(data_1_1, id, age, g_, expe_)
  data = data_1_1[data_home, nomatch=0]
  
  
  ## Likelihood -------------------------------------------------------------------
  # Realized utilities
  
  # Include the observable g (years of schooling) and experience (expe_)
  
  data[,`:=`( utility_school_1 = beta[1]*more_eighteen + beta[2]*more_twenty + beta[3]*urban + theta_1_1,
              utility_work_1 = exp(gamma[1]*g_ + gamma[2]*expe_ + gamma[3]*(expe_*expe_) + gamma[4]*urban + theta_2_1),
              utility_home_1 = delta[1]*age + delta[2]*age*age + theta_3_1,
              utility_school_2 = beta[1]*more_eighteen + beta[2]*more_twenty + beta[3]*urban + theta_1_2,
              utility_work_2 = exp(gamma[1]*g_ + gamma[2]*expe_ + gamma[3]*(expe_*expe_) + gamma[4]*urban + theta_2_2),
              utility_home_2 = delta[1]*age + delta[2]*age*age + theta_3_2
  )]
  
  
  # Type 1 ----------------------------------------------------------------------------------------------
  
  # group data by choice
  twentynine_younger_1 <- data[age <= 29 & !is.na(choice_)][, `:=` (
    v1 = utility_school_1 + discount * value_function_educ_1,
    v2 = utility_work_1 * integration_eta + discount * value_function_work_1,
    v3 = utility_home_1 + discount * value_function_home_1
  )]
  
  # compute probabilities for each choice
  twentynine_younger_1[, prob_choice_t := as.numeric(log((choice_ == 1) * brob(v1 / lambda) + 
                                                           (choice_ == 2) * brob(v2 / lambda) + 
                                                           (choice_ == 3) * brob(v3 / lambda)) -
                                                       log(sum(brob(v1 / lambda), brob(v2 / lambda), brob(v3 / lambda))))]
  
  # compute residual wages and probability adjustments for work choice
  twentynine_younger_1[choice_ == 2,residual_wage := lwage - log(utility_work_1)]
  twentynine_younger_1[choice_ == 2,denw := dnorm(residual_wage, 0, sd_wages, log = T)]
  twentynine_younger_1[choice_ == 2,prob_choice_t := ifelse(is.na(lwage), prob_choice_t, prob_choice_t + denw)]
  
  #NA 
  twentynine_younger_na_1 = data[age <= 29 & is.na(choice_)]
  
  
  # Joining all
  twentynine_younger_choices_1 = rbindlist(list(twentynine_younger_1,twentynine_younger_na_1), fill = TRUE)
  
  # Likelihood ages 30 to 32
  
  # group data by choice
  thirty_older_1 <- data[age > 29 & !is.na(choice_)][, `:=` (
    v2 = utility_work_1 * integration_eta + discount * value_function_work_1,
    v3 = utility_home_1 + discount * value_function_home_1
  )]
  
  # compute probabilities for each choice
  thirty_older_1[, prob_choice_t := as.numeric(log((choice_ == 2) * brob(v2 / lambda) + 
                                                     (choice_ == 3) * brob(v3 / lambda)) -
                                                 log(sum(brob(v2 / lambda), brob(v3 / lambda))))]
  
  # compute residual wages and probability adjustments for work choice
  thirty_older_1[choice_ == 2,residual_wage := lwage - log(utility_work_1)]
  thirty_older_1[choice_ == 2,denw := dnorm(residual_wage, 0, sd_wages, log = T)]
  thirty_older_1[choice_ == 2,prob_choice_t := ifelse(is.na(lwage), prob_choice_t, prob_choice_t + denw)]
  
  #NA 
  thirty_older_na_1 = data[age > 29 & is.na(choice_)]
  
  
  # Joining all
  thirty_older_choices_1 = rbindlist(list(thirty_older_1,thirty_older_na_1), fill = TRUE)
  
  
  final_likelihood_1 = rbindlist(list(twentynine_younger_choices_1,thirty_older_choices_1), fill = TRUE)
  final_likelihood_1[, likelihood_t_1 := prob_choice_t]
  
  # Type 2 ----------------------------------------------------------------------------------------------
  # group data by choice
  twentynine_younger_2 <- data[age <= 29 & !is.na(choice_)][, `:=` (
    v1 = utility_school_2 + discount * value_function_educ_2,
    v2 = utility_work_2 * integration_eta + discount * value_function_work_2,
    v3 = utility_home_2 + discount * value_function_home_2
  )]
  
  # compute probabilities for each choice
  twentynine_younger_2[, prob_choice_t := as.numeric(log((choice_ == 1) * brob(v1 / lambda) + 
                                                           (choice_ == 2) * brob(v2 / lambda) + 
                                                           (choice_ == 3) * brob(v3 / lambda)) -
                                                       log(sum(brob(v1 / lambda), brob(v2 / lambda), brob(v3 / lambda))))]
  
  # compute residual wages and probability adjustments for work choice
  twentynine_younger_2[choice_ == 2,residual_wage := lwage - log(utility_work_2)]
  twentynine_younger_2[choice_ == 2,denw := dnorm(residual_wage, 0, sd_wages, log = T)]
  twentynine_younger_2[choice_ == 2,prob_choice_t := ifelse(is.na(lwage), prob_choice_t, prob_choice_t + denw)]
  
  #NA 
  twentynine_younger_na_2 = data[age <= 29 & is.na(choice_)]
  
  
  # Joining all
  twentynine_younger_choices_2 = rbindlist(list(twentynine_younger_2,twentynine_younger_na_2), fill = TRUE)
  
  # Likelihood ages 30 to 32
  
  # group data by choice
  thirty_older_2 <- data[age > 29 & !is.na(choice_)][, `:=` (
    v2 = utility_work_2 * integration_eta + discount * value_function_work_2,
    v3 = utility_home_2 + discount * value_function_home_2
  )]
  
  # compute probabilities for each choice
  thirty_older_2[, prob_choice_t := as.numeric(log((choice_ == 2) * brob(v2 / lambda) + 
                                                     (choice_ == 3) * brob(v3 / lambda)) -
                                                 log(sum(brob(v2 / lambda), brob(v3 / lambda))))]
  
  # compute residual wages and probability adjustments for work choice
  thirty_older_2[choice_ == 2,residual_wage := lwage - log(utility_work_2)]
  thirty_older_2[choice_ == 2,denw := dnorm(residual_wage, 0, sd_wages, log = T)]
  thirty_older_2[choice_ == 2,prob_choice_t := ifelse(is.na(lwage), prob_choice_t, prob_choice_t + denw)]
  
  #NA 
  thirty_older_na_2 = data[age > 29 & is.na(choice_)]
  
  
  # Joining all
  thirty_older_choices_2 = rbindlist(list(thirty_older_2,thirty_older_na_2), fill = TRUE)
  
  
  final_likelihood_2 = rbindlist(list(twentynine_younger_choices_2,thirty_older_choices_2), fill = TRUE)
  final_likelihood_2[, likelihood_t_2 := prob_choice_t]
  
  
  # Log likelihood in the sample------------------------------------------------------------------------------------
  
  # Create copies of the data for type 1 and type 2
  type_1 <- copy(final_likelihood_1[, .(id, age, likelihood_t_1, both_parents_french, one_parent_french, 
                                        late_school, occup_father_1, occup_father_2, occup_mother_1, 
                                        occup_mother_2)])
  type_2 <- copy(final_likelihood_2[, .(id, age, likelihood_t_2)])
  
  # Calculate the total likelihood for each individual for type 1 and type 2
  type_1 <- type_1[, .(likelihood_1 = sum(likelihood_t_1, na.rm = TRUE)), 
                   keyby = .(id, both_parents_french, one_parent_french, late_school, 
                             occup_father_1, occup_father_2, occup_mother_1, occup_mother_2)]
  type_2 <- type_2[, .(likelihood_2 = sum(likelihood_t_2, na.rm = TRUE)), keyby = id]
  
  # Merge the two datasets and calculate the final likelihood for each individual
  l_contri_all <- merge(type_1, type_2, by = "id", all = FALSE, allow.cartesian = FALSE)
  # Calculate the probability of belonging to each type
  l_contri_all[, q_1 := omega[1] * both_parents_french + omega[2] * one_parent_french +  
                 omega[3] * late_school + omega[4] * occup_father_1 + omega[5] * occup_father_2 + 
                 omega[6] * occup_mother_1 + omega[7] * occup_mother_2 + omega[8] * 1]
  l_contri_all[, p_1 := as.numeric(brob(q_1) / (1 + brob(q_1)))]
  l_contri_all[, p_2 := 1 - p_1]
  l_contri_all[, log_sum_types := as.numeric(log(p_1 * brob(likelihood_1) + p_2 * brob(likelihood_2))), 
               by = .(id)]
  
  # Calculate the total log-likelihood for the sample
  p <- sum(l_contri_all$log_sum_types)
  
  
  return(-p/n)
  
  
} # close function

# Likelihood function chunks -----------------------------------------------------

likelihood_chunks =  function(param){
  
  
  
  likelihood_sample = foreach (j = 1:no_chunk, .combine="+")  %dopar%{ ### start loop for each chunk
    
    imin=1+(j-1)*no_i_by_chunk
    imax=min(n,j*no_i_by_chunk)
    data_i_chunk = data[id %between% c(as.numeric(imin), as.numeric(imax))]
    data_i_chunk[,id := as.factor(id)]
    contribution_chunk =  complete_likelihood(param, imin, imax, data_i_chunk)
    
  }
  
  
  
  p = likelihood_sample
  
  return(p)
  
}




#Optimization of the code
data[,id := as.numeric(id)]
#anubis
source("/softs/R/createCluster.R")
cl<-createCluster() 
# CPU
# cl<-createCluster(7)
no_chunk=14
no_i_by_chunk=ceiling((n+1)/no_chunk)
registerDoParallel(cl)
clusterExport(cl, list("data", "discount", "euler", "distribution",
                       "n_params", "nb_schooling", "nb_work", "nb_home", "nb_initial", "nb_vip",
                       "n", "no_i_by_chunk",  
                       "param", "no_chunk", "complete_likelihood"))
clusterEvalQ(cl,library("numDeriv"))
clusterEvalQ(cl,library("data.table"))
clusterEvalQ(cl,library("tidyverse"))
clusterEvalQ(cl,library("ucminf"))
clusterEvalQ(cl,library("optimx"))
clusterEvalQ(cl,library("spatstat"))
clusterEvalQ(cl,library("Brobdingnag"))


a = Sys.time()
results = ucminf(param, likelihood_chunks, hessian = 3, control = list(trace=5, grtol = 1e-10))
inv_hessian = solve(results$hessian*n)
standard_errors <- sqrt(diag(inv_hessian))
results$par
results$convergence
warnings()
b = Sys.time()
b
est <- data.table(beta = results$par, stder = standard_errors, p_value = 2*pnorm(-abs(results$par/standard_errors)))
est
stopCluster(cl)

