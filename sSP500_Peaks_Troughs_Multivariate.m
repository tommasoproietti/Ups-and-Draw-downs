%% Dating Bull and Bear Markets for SP500 stocks (starting not later than 2016) 
clear vars
load('sp500.mat', 'PT') 
%%
mPrice = table2array(PT);
mY     = price2ret(mPrice); % 
vdates = PT.Index
%[cn,cN]= size(mPrice); %  
%mp     = rmmissing(log(mPrice),2); 
vsel = sum(ismissing(mPrice)) < (cn - 1928); % series must at least start in 2016
mp =  log(mPrice(:,vsel));
[cn,cN]= size(mp); % 
%% 
ctau    = 65; cK = 22;          
mpeak = nan(cn,cN); mtrough = nan(cn,cN);
mbull = nan(cn,cN); mbear   = nan(cn,cN);

add = NaN(cn, ctau, cN); adu = NaN(cn, ctau, cN);
for k = 1:cN
    vy = mp(:,k) + 10e-15*randn(cn,1); 
    disp(['Series =', num2str(k)]);
    [vpeak, vtrough,vbull_ind, vbear_ind, vbull_dur, vbear_dur] = fBearBullDating(ctau, cK, vy); 
    mpeak(ctau+1:end,k) = vpeak; mtrough(ctau+1:end,k) = vtrough; 
    mbull(ctau+1:end,k) = vbull_ind; mbear(ctau+1:end,k) = vbear_ind; 
end
mpeak = mpeak(:,1:cN);mtrough = mtrough(:,1:cN); 
mbull = mbull(:,1:cN);mbear = mbear(:,1:cN); 

vbear_diffusion = mean(mbear,2)
plot(vdates , vbear_diffusion)
vbull_diffusion = mean(mbull,2)
plot(vdates , vbull_diffusion)
vtrough_diffusion = mean(mtrough,2)
plot(vdates , vtrough_diffusion)
save sSP500_Peaks_Troughs_Multivariate.mat


