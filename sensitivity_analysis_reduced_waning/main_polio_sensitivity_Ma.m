
clear all

clf;
% 0. Data
rawData = readtable('data/polio_data.csv');

% 1. MCMC results
mcmcRes = readtable('../main_polio/mcmc_result.csv');
mcmcRes = mcmcRes{0.2*length(mcmcRes{:,1}):end,:};
prctile(mcmcRes,[50,2.5,97.5])

% 2. Check convergence
figure(1)
for ii = 1:length(mcmcRes(1,:))
    subplot(3,6,ii)
    histogram(mcmcRes(:,ii))
end

% 3. Sensitivity analysis + Ma et al data
titration = 2.^(3:15);
log_titre = log2(titration);
yrs_sample = [10,8,6,4,30/365.25];
idxColumn = 13:17;
pv1_Ma = readtable('../other_chinese_IPV/PV1_wild_IPV_Ma_2023.xlsx');
pv2_Ma = readtable('../other_chinese_IPV/PV2_wild_IPV_Ma_2023.xlsx');
pv3_Ma = readtable('../other_chinese_IPV/PV3_wild_IPV_Ma_2023.xlsx');

% Ma data
for ii = 1:length(idxColumn)
    for jj = 1:length(yrs_sample)
        iiNumPV1 = pv1_Ma{:,idxColumn(jj)};
        iiNumPV2 = pv2_Ma{:,idxColumn(jj)};
        iiNumPV3 = pv3_Ma{:,idxColumn(jj)};
        jjYear = yrs_sample(jj);
        pv1_Ma_temp = [jjYear*ones(sum(iiNumPV1),1),repelem(titration,iiNumPV1)'];
        % 2^(mean(log2(pv1_Ma_temp(:,2))))
        pv2_Ma_temp = [jjYear*ones(sum(iiNumPV2),1),repelem(titration,iiNumPV2)'];
        % 2^(mean(log2(pv2_Ma_temp(:,2))))
        pv3_Ma_temp = [jjYear*ones(sum(iiNumPV3),1),repelem(titration,iiNumPV3)'];
        % 2^(mean(log2(pv3_Ma_temp(:,2))))
        if ii == 1 && jj == 1
            titreMa{1} = pv1_Ma_temp;
            titreMa{2} = pv2_Ma_temp;
            titreMa{3} = pv3_Ma_temp;
        else
            titreMa{1} = [titreMa{1};pv1_Ma_temp];
            titreMa{2} = [titreMa{2};pv2_Ma_temp];
            titreMa{3} = [titreMa{3};pv3_Ma_temp];
        end
    end
end

titreMa{1} = sortrows(titreMa{1},1);
titreMa{2} = sortrows(titreMa{2},1);
titreMa{3} = sortrows(titreMa{3},1);

% non constant waning
log_mu_zero = [12, 12, 11];
waning_zero = [0.5, 0.55, 0.6];
waning_deriv = [0.03, 0.05, 0.05];
sigma_CV = [0.2, 0.15, 0.18];

totalLogLikelihood(log_mu_zero,waning_zero,waning_deriv,sigma_CV,titreMa)

% Test Likelihood function
x0 = [log_mu_zero,waning_zero,waning_deriv,sigma_CV];
x0LowerBound = [0,  0, 0,  0,  0,  0,-1,-1,-1, 0, 0, 0];
x0UpperBound = [50,50,50, 50, 50, 50, 1, 1, 1,10,10,10];
disp(['Starting neg log likelihood: ',num2str(negTotalLogLikelihood(x0,titreMa))]); 

% 5. MLE
% Point estimates from fmincon
redeffun = @(x)negTotalLogLikelihood(x,titreMa);
options = optimoptions(@fmincon,'Display','iter','MaxFunEvals',30000);
if exist(strcat('mcmc_result/','mle.csv'),'file') == 2
    x0 = readmatrix(strcat('mcmc_result/','mle.csv'));
    xfmin = fmincon(redeffun,x0,[],[],[],[],x0LowerBound,x0UpperBound,[],options);
    write_matrix_new(xfmin,strcat('mcmc_result/','mle.csv'),'w',',','dec');
else
    xfmin = fmincon(redeffun,x0,[],[],[],[],x0LowerBound,x0UpperBound,[],options);
    write_matrix_new(xfmin,strcat('mcmc_result/','mle.csv'),'w',',','dec');
end
write_matrix_new(xfmin,strcat('mcmc_result/','mle.csv'),'w',',','dec');
disp(['MLE neg log likelihood: ',num2str(negTotalLogLikelihood(xfmin,titreMa))]); 

% 6. MCMC
mcSteps = 20000;
if exist(strcat('mcmc_result/','parameter_step.csv'),'file') == 2
    stepSize = readmatrix(strcat('mcmc_result/','parameter_step.csv'));
else
    stepSize = abs(0.05*(x0UpperBound-x0LowerBound));
end
out = mcmcParallel(mcSteps,titreMa,x0,stepSize,x0LowerBound,x0UpperBound);
write_matrix_new(out,strcat('mcmc_result/','mcmc_res.csv'),'w',',','dec');

% 7. Check convergence
mcmcRes = readtable('mcmc_result/mcmc_res.csv');
out = mcmcRes{0.2*length(mcmcRes{:,1}):end,:};

figure(2)
for ii = 1:length(out(1,:))
    subplot(3,4,ii)
    histogram(out(:,ii))
end

prctile(out,[50,2.5,97.5])

dt = 0.01;
numYearsDt = 20/dt;
for iiMC = 1:length(out(:,1))
    for iiPolioType = 1:length(titreMa)
        waning_zero = out(iiMC,3+iiPolioType);
        waning_deriv = out(iiMC,6+iiPolioType);
        waningRate(:,iiPolioType,iiMC) = waning_zero+waning_deriv*(1:numYearsDt)*dt;
        waningRateCum(:,iiPolioType,iiMC) = cumsum(waningRate(:,iiPolioType,iiMC)*dt);
    end
end


prctile(waningRateCum(4/dt,:,:),[50,2.5,97.5],3)

% prctile(2.^out(:,1),[50,2.5,97.5],1)

prctile(2.^out(:,1:3),[50,2.5,97.5],1)
prctile(2.^(out(:,1:3)-reshape(waningRateCum(2.5/dt,:,:),[3,length(out(:,1))])'),[50,2.5,97.5],1)
prctile(2.^(out(:,1:3)-reshape(waningRateCum(4.5/dt,:,:),[3,length(out(:,1))])'),[50,2.5,97.5],1)
prctile(2.^(out(:,1:3)-reshape(waningRateCum(6.5/dt,:,:),[3,length(out(:,1))])'),[50,2.5,97.5],1)
prctile(2.^(out(:,1:3)-reshape(waningRateCum(8.5/dt,:,:),[3,length(out(:,1))])'),[50,2.5,97.5],1)

prctile((2.^out(:,1:3))./(2.^(out(:,1:3)-reshape(waningRateCum(4/dt,:,:),[3,length(out(:,1))])')),[50,2.5,97.5],1)
prctile((2.^out(:,1:3))./(2.^(out(:,1:3)-reshape(waningRateCum(8/dt,:,:),[3,length(out(:,1))])')),[50,2.5,97.5],1)

prctile(out(:,7:9)./out(:,4:6)*100,[50,2.5,97.5],1)