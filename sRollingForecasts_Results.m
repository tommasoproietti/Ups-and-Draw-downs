%% Rolling forecasting results Tables 5-6
load('sRollingForecasts_PMC.mat', 'mfError', 'mdd_sim', 'vdates');
load('sRollingForecasts_D_ARMA_ARIMA_HAR.mat', 'mfError_arima','mfError_arma','mfError_har');
load('sRollingForecasts_C_HAR.mat', 'mfError_cond')
mfError_cond_har = mfError_cond;
load('sRollingForecasts_C_Factor_HAR.mat', 'mfError_cond')
%% Accuracy metrics (Table 4)
vnumeraire = sqrt(mean(mfError.^2));
mRMSE       = [vnumeraire; sqrt(mean(mfError_arma.^2)); sqrt(mean(mfError_arima.^2)); sqrt(mean(mfError_har.^2)); ...
              sqrt(mean(mfError_cond_har.^2)); sqrt(mean(mfError_cond.^2))]; 
mRMSEratios = [sqrt(mean(mfError_arma.^2)); sqrt(mean(mfError_arima.^2)); sqrt(mean(mfError_har.^2)); ...
               sqrt(mean(mfError_cond_har.^2)); sqrt(mean(mfError_cond.^2))] ./ vnumeraire  ;
mRes        = [vnumeraire; mRMSEratios];
mA = [mRes(:,1:11)'; mRes(:,12:end)'];
writematrix(mA, 'rRollingForecast_results.xlsx')
%% Diebold Mariano Test, Table 5
mDMtests  = NaN(22,5);
for h =1:22
    vdm = NaN(1,5);
    ctrunc = h + 1;
    ve2 = mfError(:,h);
    ve1 = mfError_arma(:,h); 
    [dDMtest, ~] = fDieboldMariano(ve1, ve2, ctrunc); vdm(1) = dDMtest;
    ve1 = mfError_arima(:,h); 
    [dDMtest, ~] = fDieboldMariano(ve1, ve2, ctrunc); vdm(2) = dDMtest;
        ve1 = mfError_har(:,h); 
    [dDMtest, ~] = fDieboldMariano(ve1, ve2, ctrunc); vdm(3) = dDMtest;
    ve1 = mfError_cond_har(:,h); 
    [dDMtest, ~] = fDieboldMariano(ve1, ve2, ctrunc); vdm(4) = dDMtest;
    ve1 = mfError_cond(:,h); 
    [dDMtest, dpval] = fDieboldMariano(ve1, ve2, ctrunc); vdm(5) = dDMtest;
    mDMtests(h,:) = vdm;
end
writematrix(mDMtests, 'rRollingForecast_results_DM.xlsx')
 