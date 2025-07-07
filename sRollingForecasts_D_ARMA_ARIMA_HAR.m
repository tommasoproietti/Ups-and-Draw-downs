%% Rolling forecasts D-ARMA, D-ARIMA(0,1,1) and D-HAR 
clear vars;
%% Graphic parameters
set(0,'DefaultAxesFontName','Times New Roman');
set(0,'DefaultTextFontName','Times New Roman');
set(0,'DefaultAxesFontSize',15);
set(0,'DefaultTextFontSize',15);
%% Read data
TT0     = readtimetable(append('^GSPC.csv'));
PT      = TT0(:,4);  
vPrice  = table2array(PT);
vdates  = PT.Index;
[cn,cN] = size(vPrice); %  
vp      = log(vPrice); 
rng     = 5211314;  % random number seed
vpe     = vp  + 10e-15*randn(cn,1); % add tiny random error to solve ties
vy      = diff(vpe);     % returns
vdatesy = vdates(2:end); 
cny     = cn -1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% drawdown
ctau    = 22;   % max filter 
[vMax, vGap_max, mS_max, mT_max, vpi_max] = fMaxFilter(vpe, ctau);
[vMin, vGap_min, mS_min ,mT_min, vpi_min] = fMaxFilter(-vpe, ctau);
vdd     = NaN(cn,1);    % drawdown 
vdd(ctau+1:end) = -vGap_max;  % observed value
vdu     = NaN(cn,1);    % drawdown 
vdu(ctau+1:end) = -vGap_min;  % observed value
vc = vdu-vdd;
%% Specifications 
% ARIMA(0,1,1)
Mdl_arima = arima(0,1,1); Mdl_arima.Constant = 0;
cL    = 22;         %   forecast horizon
cM    = 2;         %   number of MC replications
cT = 4000;      % length of window
%% select model for dd using the first cT obs
vdd0 = vdd(ctau+1:cT+ctau); 
cpmax = 5;
cqmax = 2;
maic = NaN(cpmax+1,cqmax+1);  mbic =  NaN(cpmax+1,cqmax+1);
for  p = 0:cpmax
    cp = p;
    for q =0:2
        cq = q;
        Mdl_arma = arima(cp,0,cq); Mdl.Constant = NaN;
        [EstMdl_arma, ~, dlogL,info] = estimate(Mdl_arma, vdd0);
        [maic(p+1,q+1),mbic(p+1,q+1)] = aicbic(dlogL,cp+cq+1,cT);
    end
end
minAIC = min(maic,[],'all'); [min1,min2] = find(minAIC == maic);
Mdl_arma = arima(min1-1,0,min2-1); Mdl.Constant = NaN;
minBIC = min(mbic,[],'all'); [min1,min2] = find(minBIC == mbic);
% ARMA selected by aic 
% initialize
cr = cny-cT-cL; % number of rolling window
mPar = NaN(5,cr);   
mfError_arma = NaN(cr, cL); mfError_arma = NaN(cr, cL);
mf_arma = NaN(cr,cL);     mf_armac = NaN(cr,cL);
mf_har  = NaN(cr,cL);     mfError_har=NaN(cr,cL);
%%
for i = 0:cr-1
    clc;
    disp(['Processing rolling forecast ', num2str(i), ' out of ', num2str(cr)])
    vdd_obs = vdd(cT+i+2:cT+i+1+cL);   % observed values
    vys  = vy(i+1:cT+i);        % rolling sample of returns
    if (i<=ctau) 
        vdds = vdd(i+1+ctau:cT+i);  % rolling sample of dd
    else 
        vdds = vdd(i+1:cT+i); 
    end 
    %% estimation  
    % HAR specification 
    vma5   = movavg(vdds,'simple',5); vma22 = movavg(vdds,'simple',22);
    mx     = [vdds(22:end-1) vma5(22:end-1) vma22(22:end-1)]; 
    mdl    = fitlm(mx, vdds(23:end)); 
    dc     = mdl.Coefficients.Estimate(1);
    dphi1  = mdl.Coefficients.Estimate(2); 
    dphi5  = mdl.Coefficients.Estimate(3); 
    dphi22 = mdl.Coefficients.Estimate(4); 
    vphi   = dphi1*[1; zeros(21,1)] + dphi5 * [ones(5,1); zeros(17,1)]/5 + dphi22/22;
    Mdlhar = arima('Constant', dc, 'ARLags', 1:22, 'AR', num2cell(vphi,22));  
    Mdlhar.Variance = mdl.MSE; 
    if any(i==0:10:(cr-1))
        [EstMdl_arma, ~, ~,info] = estimate(Mdl_arma, vdds,'Display','off');
        [EstMdl_arima, ~, ~,info] = estimate(Mdl_arima, vdds,'Display','off');
        % HAR estimation  and AR representation
        vma5 = movavg(vdds,'simple',5); vma22 = movavg(vdds,'simple',22);
        mx  = [vdds(22:end-1) vma5(22:end-1) vma22(22:end-1)];
        mdl = fitlm(mx, vdds(23:end)); 
        dc    = mdl.Coefficients.Estimate(1);  dphi1  = mdl.Coefficients.Estimate(2); 
        dphi5 = mdl.Coefficients.Estimate(3);  dphi22 = mdl.Coefficients.Estimate(4);
        vphi = dphi1*[1; zeros(21,1)] + dphi5 * [ones(5,1); zeros(17,1)]/5 + dphi22/22;
        EstMdl_har = arima('Constant', dc, 'ARLags', 1:22, 'AR', num2cell(vphi,22));
        EstMdl_har.Variance = mdl.MSE;
    end
    %% forecasting
    [vdd_f_arma,~] = forecast(EstMdl_arma,cL,vdds);
    [vdd_f_arima,~] = forecast(EstMdl_arima,cL,vdds);
    [vdd_f_har, ~] = forecast(EstMdl_har, cL,vdds);

    mf_arma(i+1,:) = vdd_f_arma';   mfError_arma(i+1,:) = (vdd_obs-vdd_f_arma)';
    mf_arima(i+1,:) = vdd_f_arima'; mfError_arima(i+1,:) = (vdd_obs-vdd_f_arima)';
    mf_har(i+1,:) = vdd_f_har';     mfError_har(i+1,:) = (vdd_obs-vdd_f_har)';
 
end
save sRollingForecasts_arima_har.mat
%load sRollingForecasts_arima_har.mat
