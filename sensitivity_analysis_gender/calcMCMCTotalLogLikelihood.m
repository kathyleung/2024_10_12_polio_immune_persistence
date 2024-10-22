function out = calcMCMCTotalLogLikelihood(x0,log_titre_male, ...
    log_titre_female,time_from_vaccination_male,time_from_vaccination_female)


log_zero_male = x0(1);
coeff_var = x0(2);
waning = x0(3);
log_zero_sex_diff = x0(4);


totalLogL = totalLogLikelihood(log_titre_male, log_titre_female, ...
    time_from_vaccination_male,time_from_vaccination_female,...
    log_zero_male,coeff_var,waning,log_zero_sex_diff);

out = totalLogL;

end

