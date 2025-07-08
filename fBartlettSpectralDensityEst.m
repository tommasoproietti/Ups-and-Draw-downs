function sdf = fBartlettSpectralDensityEst(y, m, omega)
% Time domain estimator of spectral density using a Bartlett window
%   
gamma = fAutocovarianceFunction(y, m);
sdf = gamma(1);
for k=1:(m-1)
    sdf = sdf + 2 * (m-k) * gamma(k+1) .* cos(omega * k) / m;
end
sdf = sdf/(2*pi);
end

