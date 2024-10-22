
clear all
clf;
% 0. Data
rawData = readtable('data/polio_data.csv');

% % 1. MCMC results
mcmcRes_PV1_d4 = readtable('mcmc_result/d4_PV1/mcmc_res.csv');
mcmcRes_PV1_d4 = mcmcRes_PV1_d4{0.2*length(mcmcRes_PV1_d4{:,1}):end,:};
mcmcRes_PV2_d4 = readtable('mcmc_result/d4_PV2/mcmc_res.csv');
mcmcRes_PV2_d4 = mcmcRes_PV2_d4{0.2*length(mcmcRes_PV2_d4{:,1}):end,:};
mcmcRes_PV3_d4 = readtable('mcmc_result/d4_PV3/mcmc_res.csv');
mcmcRes_PV3_d4 = mcmcRes_PV3_d4{0.2*length(mcmcRes_PV3_d4{:,1}):end,:};
mcmcRes_PV1_d5 = readtable('mcmc_result/d5_PV1/mcmc_res.csv');
mcmcRes_PV1_d5 = mcmcRes_PV1_d5{0.2*length(mcmcRes_PV1_d5{:,1}):end,:};
mcmcRes_PV2_d5 = readtable('mcmc_result/d5_PV2/mcmc_res.csv');
mcmcRes_PV2_d5 = mcmcRes_PV2_d5{0.2*length(mcmcRes_PV2_d5{:,1}):end,:};
mcmcRes_PV3_d5 = readtable('mcmc_result/d5_PV3/mcmc_res.csv');
mcmcRes_PV3_d5 = mcmcRes_PV3_d5{0.2*length(mcmcRes_PV3_d5{:,1}):end,:};

mcmcRes = [mcmcRes_PV1_d4,mcmcRes_PV2_d4,mcmcRes_PV3_d4,mcmcRes_PV1_d5,mcmcRes_PV2_d5,mcmcRes_PV3_d5]
% 
% % 2. Check convergence
% figure(1)
% for ii = 1:length(mcmcRes(1,:))
%     subplot(3,6,ii)
%     histogram(mcmcRes(:,ii))
% end

% % 3. Calculate titre distribution
res_log_mu = mcmcRes(:,[1,5,9,13,17,21]);
res_log_mu_sex = mcmcRes(:,[4,8,12,16,20,24]);
res_CV = mcmcRes(:,[2,6,10,14,18,22]);
res_waning = mcmcRes(:,[3,7,11,15,19,23]);
% 0-15 years after vaccination
yearsAfterVax = 0:0.05:25;
% numRandom = 500;  
log_thres = log2(8);

for ii = 1:length(res_log_mu(1,:))
    disp(ii)
    for jj = 1:length(res_log_mu(:,1))
        disp(jj)
        ij_log_titre_mu_male = res_log_mu(jj,ii)-res_waning(jj,ii)*yearsAfterVax;
        ij_log_titre_mu_female = res_log_mu(jj,ii)+res_log_mu_sex(jj,ii)-res_waning(jj,ii)*yearsAfterVax;
        ij_log_titre_sigma_male = res_CV(jj,ii)*ij_log_titre_mu_male;
        ij_log_titre_sigma_female = res_CV(jj,ii)*ij_log_titre_mu_female;
        res_log_titre_rnd_male(jj,:,ii) = normrnd(ij_log_titre_mu_male,ij_log_titre_sigma_male);
        res_log_titre_rnd_female(jj,:,ii) = normrnd(ij_log_titre_mu_female,ij_log_titre_sigma_female);
        log_titre_mean_male(jj,:,ii) = ij_log_titre_mu_male;
        log_titre_mean_female(jj,:,ii) = ij_log_titre_mu_female;
    end
    prctile_log_titre_mean_male(:,:,ii) = cat(1,nanmean(log_titre_mean_male(:,:,ii)),prctile(log_titre_mean_male(:,:,ii),[2.5,97.5]));
    prctile_log_titre_mean_female(:,:,ii) = cat(1,nanmean(log_titre_mean_female(:,:,ii)),prctile(log_titre_mean_female(:,:,ii),[2.5,97.5]));
end

gender_ratio = 2.^(log_titre_mean_female)./2.^(log_titre_mean_male);
for ii = 1:6
    prctile_log_titre_gender_ratio(:,:,ii) = cat(1,nanmean(gender_ratio(:,:,ii)),prctile(gender_ratio(:,:,ii),[2.5,97.5]));
end


anss=prctile_log_titre_gender_ratio(3,:,6);
anss(:,[1,2,3]);

results_male = []
for ii = 1:6
    for jj = 1:3
        anss=prctile_log_titre_mean_male(jj,:,ii);
        temp = anss(:,[1,81])
        results_male = cat(1,results_male,temp)
    end
end

results_male = 2.^(results_male)

writematrix(results_male, 'results_male.csv')

results_female = []
for ii = 1:6
    for jj = 1:3
        anss=prctile_log_titre_mean_female(jj,:,ii);
        temp = anss(:,[1,81])
        results_female = cat(1,results_female,temp)
    end
end

results_female = 2.^(results_female)

writematrix(results_female, 'results_female.csv')

results_sex_diff = []
for ii = 1:6
    for jj = 1:3
        anss=prctile_log_titre_gender_ratio(jj,:,ii);
        temp = anss(:,[1,2,3])
        results_sex_diff = cat(1,results_sex_diff,temp)
    end
end
writematrix(results_sex_diff, 'results_female.csv')



