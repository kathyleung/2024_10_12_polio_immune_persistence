
rm(list=ls());

# Load the given functions
source("mcmcProposal.R");
source("MCMC.R");

# 1. Prior function
# non-informative (flat) priors
logPrior <- function(parms,vecAlpha,vecBeta)
{
  # Flat priors
  logP = -log(vecBeta-vecAlpha);
  return(sum(logP));
}

# 2. Likelihood function
logLikelihood <- function(parms,data,fixedVariables)
{
  # Sanofi data
  
  logL_d4_pv1 = log(dnorm(data$pv1_log[data$dose=='Dose 4'],
                          parms[1]-parms[7]*data$yrs_since_last_dose[data$dose=='Dose 4'],
                          parms[4]*(parms[1]-parms[7]*data$yrs_since_last_dose[data$dose=='Dose 4'])));
  logL_d4_pv2 = log(dnorm(data$pv2_log[data$dose=='Dose 4'],
                          parms[2]-parms[8]*data$yrs_since_last_dose[data$dose=='Dose 4'],
                          parms[5]*(parms[2]-parms[8]*data$yrs_since_last_dose[data$dose=='Dose 4'])));
  logL_d4_pv3 = log(dnorm(data$pv3_log[data$dose=='Dose 4'],
                          parms[3]-parms[9]*data$yrs_since_last_dose[data$dose=='Dose 4'],
                          parms[6]*(parms[3]-parms[9]*data$yrs_since_last_dose[data$dose=='Dose 4'])));
  logL_d5_pv1 = log(dnorm(data$pv1_log[data$dose=='Dose 5'],
                          parms[10]-parms[16]*data$yrs_since_last_dose[data$dose=='Dose 5'],
                          parms[13]*(parms[10]-parms[16]*data$yrs_since_last_dose[data$dose=='Dose 5'])));
  logL_d5_pv2 = log(dnorm(data$pv2_log[data$dose=='Dose 5'],
                          parms[11]-parms[17]*data$yrs_since_last_dose[data$dose=='Dose 5'],
                          parms[14]*(parms[11]-parms[17]*data$yrs_since_last_dose[data$dose=='Dose 5'])));
  logL_d5_pv3 = log(dnorm(data$pv3_log[data$dose=='Dose 5'],
                          parms[12]-parms[18]*data$yrs_since_last_dose[data$dose=='Dose 5'],
                          parms[15]*(parms[12]-parms[18]*data$yrs_since_last_dose[data$dose=='Dose 5'])));
  
  # Normal likelihood 
  logL =  
    sum(logL_d4_pv1,na.rm = TRUE)+sum(logL_d4_pv2,na.rm = TRUE)+sum(logL_d4_pv3,na.rm = TRUE)+
    sum(logL_d5_pv1,na.rm = TRUE)+sum(logL_d5_pv2,na.rm = TRUE)+sum(logL_d5_pv3,na.rm = TRUE);
  return(sum(logL));
}

mcIter = 10000;
burnIn = 0.2;
# Sanofi seroprevalence data
polioData = read.csv('data/Dataset_csv_1.csv',stringsAsFactors = FALSE);
polioData$yrs_since_last_dose = as.numeric((as.Date(polioData$visit.date,'%d-%b-%Y')-as.Date(polioData$last.dose.date,'%d-%b-%Y')))/365.25;
polioData$pv1_log = log2(polioData$pv1.titre);
polioData$pv2_log = log2(polioData$pv2.titre);
polioData$pv3_log = log2(polioData$pv3.titre);
write.csv(polioData,'polio_data.csv',row.names = FALSE);
data = polioData;
# No fixed variables
fixedVariables = NULL;
# After Dose 4
# Mean titre at time 0
d4_pv1_log_zero = mean(data$pv1_log[data$dose == 'Dose 4' & data$yrs_since_last_dose<180/365.25]);
d4_pv2_log_zero = mean(data$pv2_log[data$dose == 'Dose 4' & data$yrs_since_last_dose<180/365.25]);
d4_pv3_log_zero = mean(data$pv3_log[data$dose == 'Dose 4' & data$yrs_since_last_dose<180/365.25]);
# Coefficient of variation of titres
d4_pv1_coeff_var = d4_pv1_log_zero/sd(data$pv1_log[data$dose == 'Dose 4' & data$yrs_since_last_dose<180/365.25]);
d4_pv2_coeff_var = d4_pv2_log_zero/sd(data$pv2_log[data$dose == 'Dose 4' & data$yrs_since_last_dose<180/365.25]);
d4_pv3_coeff_var = d4_pv3_log_zero/sd(data$pv3_log[data$dose == 'Dose 4' & data$yrs_since_last_dose<180/365.25]);
# Waning
d4_pv1_waning = (d4_pv1_log_zero - mean(data$pv1_log[data$dose == 'Dose 4' & data$yrs_since_last_dose>3]))/3;
d4_pv2_waning = (d4_pv2_log_zero - mean(data$pv2_log[data$dose == 'Dose 4' & data$yrs_since_last_dose>3]))/3;
d4_pv3_waning = (d4_pv3_log_zero - mean(data$pv3_log[data$dose == 'Dose 4' & data$yrs_since_last_dose>3]))/3;
# After Dose 5
# Mean titre at time 0
d5_pv1_log_zero = mean(data$pv1_log[data$dose == 'Dose 5' & data$yrs_since_last_dose<180/365.25]);
d5_pv2_log_zero = mean(data$pv2_log[data$dose == 'Dose 5' & data$yrs_since_last_dose<180/365.25]);
d5_pv3_log_zero = mean(data$pv3_log[data$dose == 'Dose 5' & data$yrs_since_last_dose<180/365.25]);
# Coefficient of variation of titres
d5_pv1_coeff_var = d5_pv1_log_zero/sd(data$pv1_log[data$dose == 'Dose 5' & data$yrs_since_last_dose<180/365.25]);
d5_pv2_coeff_var = d5_pv2_log_zero/sd(data$pv2_log[data$dose == 'Dose 5' & data$yrs_since_last_dose<180/365.25]);
d5_pv3_coeff_var = d5_pv3_log_zero/sd(data$pv3_log[data$dose == 'Dose 5' & data$yrs_since_last_dose<180/365.25]);
# Waning
d5_pv1_waning = (d5_pv1_log_zero - mean(data$pv1_log[data$dose == 'Dose 5' & data$yrs_since_last_dose>3]))/3;
d5_pv2_waning = (d5_pv2_log_zero - mean(data$pv2_log[data$dose == 'Dose 5' & data$yrs_since_last_dose>3]))/3;
d5_pv3_waning = (d5_pv3_log_zero - mean(data$pv3_log[data$dose == 'Dose 5' & data$yrs_since_last_dose>3]))/3;

# Group all parameters
parms = c(
  d4_pv1_log_zero,d4_pv2_log_zero,d4_pv3_log_zero,
  d4_pv1_coeff_var,d4_pv2_coeff_var,d4_pv3_coeff_var,
  d4_pv1_waning,d4_pv2_waning,d4_pv3_waning,
  d5_pv1_log_zero,d5_pv2_log_zero,d5_pv3_log_zero,
  d5_pv1_coeff_var,d5_pv2_coeff_var,d5_pv3_coeff_var,
  d5_pv1_waning,d5_pv2_waning,d5_pv3_waning);
# Define the step size for random walk
if (file.exists('mcmc_stepsize.csv')){
  parmStepSize = read.csv('mcmc_stepsize.csv');
  parmStepSize = t(parmStepSize);
}else{
  parmStepSize = 0.02*parms;
}
# Set constraints for parameters
# Lower bound
parmLB = 0*parms;
# Upper bound
parmUB = c(
  100,100,100,
  100,100,100,
  20,20,20,
  100,100,100,
  100,100,100,
  20,20,20);
# Define the constraints for acceptance probability to ensure random walk and convergence
minProbAccept = 0.3;
maxProbAccept = 0.7;
adjStepSize = 0.2;
# Parameters for priors: uniform
priorA = parmLB;
priorB = parmUB;
# Run MCMC
recParm = MCMC(mcIter,burnIn,data,fixedVariables,parms,parmStepSize,parmLB,parmUB,minProbAccept,maxProbAccept,adjStepSize,priorA,priorB)
# Posterior Mean
posteriorMean = colMeans(recParm[(mcIter*burnIn+1):mcIter,])
posteriorMean

# Diagnostics
recParm = read.csv('mcmc_result.csv',header = TRUE);
recParmPosterior = recParm[(burnIn*length(recParm[,1])):length(recParm[,1]),];
hist(recParmPosterior[,1])


