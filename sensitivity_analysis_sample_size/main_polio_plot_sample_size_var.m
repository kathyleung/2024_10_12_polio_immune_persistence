
clear all
clf;
% 0. Data
rawData = readtable('data/polioData.csv');

% % 1. MCMC results
mcmcRes = readtable('600data/mcmc_result.csv');
mcmcRes = mcmcRes{0.2*length(mcmcRes{:,1}):end,:};

% 2. Check convergence
figure(1)
for ii = 1:length(mcmcRes(1,:))
    subplot(3,6,ii)
    histogram(mcmcRes(:,ii))
end

% % 3. Calculate titre distribution
res_log_mu = mcmcRes(:,[1:3,10:12]);
res_CV = mcmcRes(:,[4:6,13:15]);
res_waning = mcmcRes(:,[7:9,16:18]);
% 0-15 years after vaccination
yearsAfterVax = 0:0.05:25;
numRandom = 500;
log_thres = log2(8);

adjustment = [[ones(51, 1); 0.6933*ones(450, 1)],[ones(91, 1); 0.1579*ones(410, 1)],...
    [ones(91, 1); 0.1463*ones(410, 1)]];
adjustment = [adjustment,adjustment];

% parfor
for ii = 1:length(res_log_mu(1,:))
    %disp(ii)
    for jj = 1:length(res_log_mu(:,1))
        disp(jj)
        for mm = 1:length(yearsAfterVax)
            ij_log_titre_mu(mm) = res_log_mu(jj,ii)-res_waning(jj,ii)*yearsAfterVax(mm)*adjustment(mm,ii);
        end 
        % ij_log_titre_mu = res_log_mu(jj,ii)-res_waning(jj,ii)*yearsAfterVax;
        ij_log_titre_sigma = res_CV(jj,ii)*ij_log_titre_mu;
        res_log_titre_rnd(jj,:,ii) = normrnd(ij_log_titre_mu,ij_log_titre_sigma);
        % seroprotection
        parfor kk = 1:numRandom
            % disp(kk)
            iijj_log_titre_rnd(:,:,kk) = normrnd(ij_log_titre_mu,ij_log_titre_sigma);
        end
        res_seroprotect(jj,:,ii) = sum(iijj_log_titre_rnd>log_thres,3)/numRandom;
    end
    prctile_log_titre(:,:,ii) = cat(1,nanmean(res_log_titre_rnd(:,:,ii)),prctile(res_log_titre_rnd(:,:,ii),[2.5,97.5]));
    prctile_seroprotect(:,:,ii) = cat(1,nanmean(res_seroprotect(:,:,ii)),prctile(res_seroprotect(:,:,ii),[2.5,97.5]));
end
save('600data/prctile_log_titre.mat','prctile_log_titre');
save('600data/prctile_seroprotect.mat','prctile_seroprotect');


% 4. Plot figures
% Colors
colorMp = [
    204, 0,   0;
    0, 153,   76;
    255, 178, 102;
    51, 153, 255]/255;
arrTitle = {'PV1','PV2','PV3'};
arrData = {'pv1_log','pv2_log','pv3_log'};

load('600data/prctile_log_titre.mat','prctile_log_titre');
load('600data/prctile_seroprotect.mat','prctile_seroprotect');

figure(2)
% left column
for ii = 1:3
    subplot(3,2,2*(ii-1)+1)
    % Dose 4
    hDose4Patch = patch('XData',[yearsAfterVax,fliplr(yearsAfterVax)],...
        'YData',[smooth(prctile_log_titre(2,:,ii))',smooth(fliplr(prctile_log_titre(3,:,ii)))'],...
        'FaceColor',colorMp(2,:),...
        'EdgeColor',colorMp(2,:),...
        'FaceAlpha',0.1, ...
        'LineStyle','-.');
    hold on
    hDose4 = plot(yearsAfterVax,prctile_log_titre(1,:,ii),...
        'LineWidth',1,'Color',colorMp(2,:));
    hold on
    hDataDose4 = scatter(rawData{strcmp(rawData.dose,{'Dose 4'}),'yrs_since_last_dose'},...
        rawData{strcmp(rawData.dose,{'Dose 4'}),arrData{ii}},8,...
        'filled',...
        'MarkerFaceColor',colorMp(2,:));
    % Dose 5
    hDose5Patch = patch('XData',[yearsAfterVax,fliplr(yearsAfterVax)]+5,...
        'YData',[smooth(prctile_log_titre(2,:,ii+3))',smooth(fliplr(prctile_log_titre(3,:,ii+3)))'],...
        'FaceColor',colorMp(3,:),...
        'EdgeColor',colorMp(3,:),...
        'FaceAlpha',0.1, ...
        'LineStyle','-.');
    hold on
    hDose5 = plot(yearsAfterVax+5,prctile_log_titre(1,:,ii+3),...
        'LineWidth',1,'Color',colorMp(3,:));
    hold on
    hDataDose5 = scatter(rawData{strcmp(rawData.dose,{'Dose 5'}),'yrs_since_last_dose'},...
        rawData{strcmp(rawData.dose,{'Dose 5'}),arrData{ii}},8,...
        'filled',...
        'MarkerFaceColor',colorMp(3,:));

    xlim([0,9])
    xlabel('Years after receiving the 4th dose')
    ylim([0,17.2])
    ylabel('Neutralising antibody titre')
    set(gca,...
        'XTick',0:9,...
        'YTick',log2(2.^(3:4:20)),...
        'YTickLabels',{'8','128','2048','32768','524288'})

    % Threshold
    hLine = line([0,15],[log2(8),log2(8)],'Color','black','LineWidth',1,'LineStyle','--');
    
    % Legend
    hLegend = legend([hDose4,hDose5],{'4th dose','5th dose'},'Location','southwest');
    set(hLegend,'box','off')
    % Title
    title(arrTitle{ii})
end

% right column
for ii = 1:3
    subplot(3,2,2*ii)
    % Dose 4
    hDose4Patch = patch('XData',[yearsAfterVax,fliplr(yearsAfterVax)],...
        'YData',[prctile_seroprotect(2,:,ii),fliplr(prctile_seroprotect(3,:,ii))],...
        'FaceColor',colorMp(2,:),...
        'EdgeColor',colorMp(2,:),...
        'FaceAlpha',0.1, ...
        'LineStyle','-.');
    hold on
    hDose4 = plot(yearsAfterVax,prctile_seroprotect(1,:,ii),...
        'LineWidth',1,'Color',colorMp(2,:));
    hold on

    % Dose 5
    hDose5Patch = patch('XData',[yearsAfterVax,fliplr(yearsAfterVax)],...
        'YData',[prctile_seroprotect(2,:,ii+3),fliplr(prctile_seroprotect(3,:,ii+3))],...
        'FaceColor',colorMp(3,:),...
        'EdgeColor',colorMp(3,:),...
        'FaceAlpha',0.1, ...
        'LineStyle','-.');
    hold on
    hDose5 = plot(yearsAfterVax,prctile_seroprotect(1,:,ii+3),...
        'LineWidth',1,'Color',colorMp(3,:));
    hold on

    xlim([0,20])
    xlabel('Years after receiving the vaccination')
    ylim([0,1])
    ylabel('Seroprotection rate')
    
    % Legend
    hLegend = legend([hDose4,hDose5],{'4th dose','5th dose'},'Location','southwest');
    set(hLegend,'box','off')
    % Title
    title(arrTitle{ii})
end