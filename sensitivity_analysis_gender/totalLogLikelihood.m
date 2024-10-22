function totalLogL = totalLogLikelihood(log_titre_male, log_titre_female, ...
    time_from_vaccination_male,time_from_vaccination_female,...
    log_zero_male,coeff_var,waning,log_zero_sex_diff)

% 1. Log likelihood of male
LogL_male = nansum(log(normpdf(table2array(log_titre_male),table2array(log_zero_male-waning.*time_from_vaccination_male), ...
    table2array(coeff_var.*(log_zero_male-waning.*time_from_vaccination_male)))));


% 2. Log Likelihood of female
LogL_female = nansum(log(normpdf(table2array(log_titre_female),table2array((log_zero_male+log_zero_sex_diff)-waning.*time_from_vaccination_female), ...
    table2array(coeff_var.*((log_zero_male+log_zero_sex_diff)-waning.*time_from_vaccination_female)))));


% total loglikelihood
totalLogL = LogL_male+LogL_female;

end

