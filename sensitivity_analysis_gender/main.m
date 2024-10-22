
clear all;
rng('shuffle');

% 1. Load raw data and titre information
rawData = readtable('data/polioData.csv');
rawData_male = rawData(ismember(rawData.sex, 'M'), :);
rawData_female = rawData(ismember(rawData.sex, 'F'), :);
rawData_d4 = rawData(ismember(rawData.dose, 'Dose 4'), :);
rawData_d5 = rawData(ismember(rawData.dose, 'Dose 5'), :);
rawData_male_d4 = rawData_male(ismember(rawData_male.dose, 'Dose 4'), :);
rawData_male_d5 = rawData_male(ismember(rawData_male.dose, 'Dose 5'), :);
rawData_female_d4 = rawData_female(ismember(rawData_female.dose, 'Dose 4'), :);
rawData_female_d5 = rawData_female(ismember(rawData_female.dose, 'Dose 5'), :);

log_titre_male_d4 = rawData_male_d4(:,{'pv1_log', 'pv2_log', 'pv3_log'});
log_titre_male_d5 = rawData_male_d5(:,{'pv1_log', 'pv2_log', 'pv3_log'});
log_titre_female_d4 = rawData_female_d4(:,{'pv1_log', 'pv2_log', 'pv3_log'});
log_titre_female_d5 = rawData_female_d5(:,{'pv1_log', 'pv2_log', 'pv3_log'});
time_from_vaccination_male_d4 = rawData_male_d4(:,{'yrs_since_last_dose'});
time_from_vaccination_male_d5 = rawData_male_d5(:,{'yrs_since_last_dose'});
time_from_vaccination_female_d4 = rawData_female_d4(:,{'yrs_since_last_dose'});
time_from_vaccination_female_d5 = rawData_female_d5(:,{'yrs_since_last_dose'});

% 3. Estimate Model parameters
log_zero_male_d4 = mean(rawData_male_d4(rawData_male_d4.yrs_since_last_dose < 180/365.25, {'pv1_log','pv2_log','pv3_log'}));
log_zero_male_d5 = mean(rawData_male_d5(rawData_male_d5.yrs_since_last_dose < 180/365.25, {'pv1_log','pv2_log','pv3_log'}));
log_zero_sex_diff_d4 = mean(rawData_female_d4(rawData_female_d4.yrs_since_last_dose < 180/365.25, {'pv1_log','pv2_log','pv3_log'}))-log_zero_male_d4;
log_zero_sex_diff_d5 = mean(rawData_female_d5(rawData_female_d5.yrs_since_last_dose < 180/365.25, {'pv1_log','pv2_log','pv3_log'}))-log_zero_male_d5;
waning_d4 = (mean(rawData_d4(rawData_d4.yrs_since_last_dose < 180/365.25, {'pv1_log','pv2_log','pv3_log'}))- ...
    mean(rawData_d4(rawData_d4.yrs_since_last_dose > 3, {'pv1_log','pv2_log','pv3_log'})))./3;
waning_d5 = (mean(rawData_d5(rawData_d5.yrs_since_last_dose < 180/365.25, {'pv1_log','pv2_log','pv3_log'}))- ...
    mean(rawData_d5(rawData_d5.yrs_since_last_dose > 3, {'pv1_log','pv2_log','pv3_log'})))./3;
coeff_var_d4 = mean(rawData_d4(rawData_d4.yrs_since_last_dose < 180/365.25, {'pv1_log','pv2_log','pv3_log'}))./std( ...
    rawData_d4(rawData_d4.yrs_since_last_dose < 180/365.25, {'pv1_log','pv2_log','pv3_log'}));
coeff_var_d5 = mean(rawData_d5(rawData_d5.yrs_since_last_dose < 180/365.25, {'pv1_log','pv2_log','pv3_log'}))./std( ...
    rawData_d5(rawData_d5.yrs_since_last_dose < 180/365.25, {'pv1_log','pv2_log','pv3_log'}));

serotype = {'PV1' 'PV2' 'PV3'};
dose = {'d4' 'd5'};

for ii = 1:length(serotype)
    for jj = 1:length(dose)
        serotype_choice = string(serotype(ii));
        dose_choice = string(dose(jj));
        if exist(strcat('mcmc_result/',dose_choice,'_',serotype_choice),'dir') == 0
            mkdir(strcat('mcmc_result/',dose_choice,'_',serotype_choice));
        end

        if jj == 1
            log_zero_male_temp = log_zero_male_d4;
            coeff_var_temp = coeff_var_d4;
            waning_temp = waning_d4;
            log_zero_sex_diff_temp = log_zero_sex_diff_d4;
            log_titre_male_temp = log_titre_male_d4;
            log_titre_female_temp = log_titre_female_d4;
            time_from_vaccination_male = time_from_vaccination_male_d4;
            time_from_vaccination_female = time_from_vaccination_female_d4;
        else
            log_zero_male_temp = log_zero_male_d5;
            coeff_var_temp = coeff_var_d5;
            waning_temp = waning_d5;
            log_zero_sex_diff_temp = log_zero_sex_diff_d5;
            log_titre_male_temp = log_titre_male_d5;
            log_titre_female_temp = log_titre_female_d5;
            time_from_vaccination_male = time_from_vaccination_male_d5;
            time_from_vaccination_female = time_from_vaccination_female_d5;
        end

        log_zero_male = table2array(log_zero_male_temp(1,ii));
        coeff_var = table2array(coeff_var_temp(1,ii));
        waning = table2array(waning_temp(1,ii));
        log_zero_sex_diff = table2array(log_zero_sex_diff_temp(1,ii));
        log_titre_male = log_titre_male_temp(:,ii);
        log_titre_female = log_titre_female_temp(:,ii);
    
        totalLogL = totalLogLikelihood(log_titre_male, log_titre_female, ...
        time_from_vaccination_male,time_from_vaccination_female,...
        log_zero_male,coeff_var,waning,log_zero_sex_diff);
        disp(['Starting log likelihood: ',num2str(totalLogL)]);
    
         % Test Likelihood function
        x0 = [log_zero_male,coeff_var,waning,log_zero_sex_diff];
        ori = [log_zero_male,coeff_var,waning,log_zero_sex_diff];
        x0LowerBound = [0,0,0,-5];
        x0UpperBound = [100,100,20,5];
        disp(['Starting neg log likelihood: ',num2str(negTotalLogLikelihood( ...
            x0,log_titre_male,log_titre_female,time_from_vaccination_male,time_from_vaccination_female))]); 
    
        % MCMC
        mcSteps = 20000;
        if exist(strcat('mcmc_result/',dose_choice,'_',serotype_choice,'/','parameter_step.csv'),'file') == 2
            stepSize = csvread(strcat('mcmc_result/',dose_choice,'_',serotype_choice,'/','parameter_step.csv'));
        else
            stepSize = abs(0.02*ori);
        end
        out = mcmcParallel(serotype_choice,dose_choice,mcSteps,log_titre_male,log_titre_female, ...
            time_from_vaccination_male,time_from_vaccination_female,...
            ori,stepSize,x0LowerBound,x0UpperBound);
    end
end 