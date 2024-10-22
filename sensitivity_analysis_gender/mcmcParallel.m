function [iPar,iiLogL] = mcmcParallel(serotype_choice,dose_choice,mcSteps,log_titre_male, ...
    log_titre_female,time_from_vaccination_male,time_from_vaccination_female, ...
    parameters,parameterSteps,lowerLimit,upperLimit)

goodMC = false;
pAccept = 0.5*ones(1,length(parameters));
numWhileLoop = 0;
numPars = length(parameters);

while goodMC == false
    iPar=zeros(mcSteps,numPars); 
    iiLogL=zeros(mcSteps,1);
    nAccept = zeros(1,numPars);
    parameterSteps(pAccept > 0.7) = parameterSteps(pAccept > 0.7)*1.2;
    parameterSteps(pAccept < 0.3) = parameterSteps(pAccept < 0.3)/1.25;
    parameterSteps = min([parameterSteps;upperLimit-lowerLimit]);
    % Current loglikelihood
    parameters_c = parameters;
    logLikelihood = calcMCMCTotalLogLikelihood(parameters_c,log_titre_male, ...
        log_titre_female,time_from_vaccination_male,time_from_vaccination_female);
    for tt = 1:mcSteps
        for ii = 1:numPars
            tempNew = mcmcProposal(parameters_c(ii),parameterSteps(ii),lowerLimit(ii),upperLimit(ii));
            parameters_mc = parameters_c;
            parameters_mc(ii) = tempNew;
            logLikelihood_new = calcMCMCTotalLogLikelihood(parameters_mc,log_titre_male, ...
                log_titre_female,time_from_vaccination_male,time_from_vaccination_female);
            alpha = min(1,exp(logLikelihood_new-logLikelihood));
            uu = rand;
            if uu <= alpha
                parameters_c = parameters_mc;
                logLikelihood = logLikelihood_new;
                nAccept(ii) = nAccept(ii)+1;
            end
        end
        iiLogL(tt,:) = logLikelihood;
        iPar(tt,:) = parameters_c;
        % Check MCMC
        if mod(tt,mcSteps/10) == 0
            pAccept = nAccept/tt;
            write_matrix_new(parameterSteps,strcat('mcmc_result/',dose_choice,'_',serotype_choice,'/','parameter_step.csv'),'w',',','dec');
            write_matrix_new(iPar,strcat('mcmc_result/',dose_choice,'_',serotype_choice,'/','mcmc_res.csv'),'w',',','dec');
            write_matrix_new(iiLogL,strcat('mcmc_result/',dose_choice,'_',serotype_choice,'/','log_likelihood.csv'),'w',',','dec');
            disp('MCMC')
            display(tt/mcSteps);
            disp('Acceptance probability');
            display(pAccept);
            if any(pAccept > 0.7) || any(pAccept < 0.3)
                goodMC = false;
                if numWhileLoop <= 30
                    if tt/mcSteps >= 0.2
                        numWhileLoop = numWhileLoop+1;
                       break;
                    else
                       goodMC = true;
                    end
                else
                    goodMC = true;
                end
            else
                goodMC = true;
            end
        end     
    end
end

end



