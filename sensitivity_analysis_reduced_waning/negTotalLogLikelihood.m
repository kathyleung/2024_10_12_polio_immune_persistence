function out = negTotalLogLikelihood(x0,titreMa)

numType = length(titreMa);

log_mu_zero = x0(1:numType);
waning_zero = x0(numType+(1:numType));
waning_deriv = x0(2*numType+(1:numType));
sigma_CV = x0(3*numType+(1:numType));

totalLogL = totalLogLikelihood(log_mu_zero,waning_zero,waning_deriv,sigma_CV,titreMa);

out = -totalLogL;

end

