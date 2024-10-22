function totalLogL = totalLogLikelihood(log_mu_zero,waning_zero,waning_deriv,sigma_CV,titreMa)

dt = 0.01;
numYearsDt = 20/dt;
for iiPolioType = 1:length(titreMa)
    waningRate(:,iiPolioType) = waning_zero(iiPolioType)+waning_deriv(iiPolioType)*(1:numYearsDt)*dt;
    waningRateCum(:,iiPolioType) = cumsum(waningRate(:,iiPolioType)*dt);
end

for iiPolioType = 1:length(titreMa)
    iiTitreMa = titreMa{iiPolioType};
    iiTitreMa = iiTitreMa(~isnan(iiTitreMa(:,2)),:);
    iiWaning = waningRateCum(ceil(iiTitreMa(:,1)/dt),iiPolioType);
    iiMu = log_mu_zero(iiPolioType)-iiWaning;
    iiSigma = iiMu*sigma_CV(iiPolioType);
    iiLogL = log(normpdf(log2(iiTitreMa(:,2)),iiMu,iiSigma.*ones(size(iiMu))));
    iiLogLSum(iiPolioType) = sum(iiLogL);
end
iiLogLSum(isnan(iiLogLSum)) = -1e9;

totalLogL = sum(iiLogLSum);

end

