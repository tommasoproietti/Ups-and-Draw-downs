%% Rolling forecasts C-HAR and C-factor HAR
clear vars;  %  
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
%% multiperiod returns
mr = NaN(cn-1,ctau);
for i = 1:ctau
    mr(i:end,i) =  vp(i+1:end)- vp(1:end-i);
end
mar = abs(mr)
%% Specifications 
cL = 22;         %   forecast horizon
cT = 4000;      % length of window
cr = cny-cT-cL; % number of rolling window
mPar = NaN(4,cr);   
mfError_cond = NaN(cr, cL); mf_cond = NaN(cr,cL); 
%%
for i = 0:cr-1
    clc;
    disp(['Processing rolling forecast ', num2str(i), ' out of ', num2str(cr)])
    vdd_obs = vdd(cT+i+2:cT+i+1+cL);   % observed values
    vps  = vp(i+1:cT+i);        % rolling sample of returns
    if (i<=ctau) 
         mars = mar(i+1+ctau:cT+i,:);        % rolling sample of returns
    else 
         mars = mar(i+1:cT+i,:);        % rolling sample of returns
    end 
    %% estimation and forecasting
    mvar_f_har = zeros(ctau+1,cL); vdd_f_cond = NaN(cL,1); 
    % factor model
    cK = 3; % number of factors
    mf_f_har = zeros(cL,cK);
    mS = cov(mars);
    [mV, mLambda] = eigs(mS, cK);
    vLambda = diag(mLambda)
    mFactors =  mars*mV;
    for t = 1:cK
        vars = mFactors(:,t);
        if any(i==0:10:(cr-1))
            % HAR specification
            vma5 = movavg(vars,'simple',5); vma22 = movavg(vars,'simple',22);
            mx  = [vars(22:end-1) vma5(22:end-1) vma22(22:end-1)];
            mdl = fitlm(mx, vars(23:end));
            dc = mdl.Coefficients.Estimate(1);  dphi1 = mdl.Coefficients.Estimate(2);
            dphi5 = mdl.Coefficients.Estimate(3);     dphi22 = mdl.Coefficients.Estimate(4);
            vphi = dphi1*[1; zeros(21,1)] + dphi5 * [ones(5,1); zeros(17,1)]/5 + dphi22/22;
            EstMdl_har = arima('Constant', dc, 'ARLags', 1:22, 'AR', num2cell(vphi,22));
            EstMdl_har.Variance = mdl.MSE;
        end
        [vf_f_har, ~] = forecast(EstMdl_har, cL,vars);
        mf_f_har(:,t) =  vf_f_har;
    end
    mvar_f_har(2:end,:) = mV*mf_f_har';
    [~, ~, mS_max, mT_max, vpi_max] = fMaxFilter(vps+10e-15*randn(cT,1), ctau);
    mTs = eye(ctau+1); vs = mS_max(:,end);
    for l = 1:cL
        mTs = mT_max * mTs;
        vdd_f_cond(l) = vs'*mTs*mvar_f_har(:,l);
    end
    %% forecasting
    
    mf_cond(i+1,:) = vdd_f_cond';   mfError_cond(i+1,:) = (vdd_obs-vdd_f_cond)';
 
 
end
 
save sRollingForecasts_C_Factor_HAR.m
 