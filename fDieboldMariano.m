function [ dDMtest, dpval ] = fDieboldMariano(ve1, ve2, ctrunc )
    vd = ve1.^2-ve2.^2; % loss differential
    cn = length(vd);
    dLRv = 2* pi * fBartlettSpectralDensityEst(vd, ctrunc, 0);
    dDMtest  = mean(vd)/sqrt(dLRv/cn); 
    dpval  = normcdf(dDMtest);
end

