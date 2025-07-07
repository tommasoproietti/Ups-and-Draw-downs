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
vp      = log(vPrice); 
vdates  = PT.Index;
[cn,cN] = size(vPrice); %  
EstMdl = estimate(arima(0,1,1), vp); dtheta = EstMdl.MA{1};
%% Model based filtering using cutoff1=2*pi/5 and cutoff2 = 2*pi/22
dcutoff1 =  2*pi/5; dcutoff2 =  2*pi/22; cs = 6; 
[mZ, mGG, mT, mHH, mW, va, mP, mW0] = fSsfMBfiltering([], dtheta, cs, dcutoff1);
[mAs, ~, ~, ~] = fAugmentedStateSmoother(vp, mZ, mGG, mT, mHH, mW, va, mP, mW0);
%plot(vdates , vp, '.'); hold on; plot(vdates,  mAs(1,:)'); hold off
%plot(vdates,  vp-mAs(1,:)');  
vps1 = mAs(1,:)';
[mZ, mGG, mT, mHH, mW, va, mP, mW0] = fSsfMBfiltering([], dtheta, cs, dcutoff2);
[mAs, mPs, vBeta, mVarBeta] = fAugmentedStateSmoother(vp, mZ, mGG, mT, mHH, mW, va, mP, mW0);
 
vps2 = mAs(1,:)';
%% Calculate drawdown
ctau  = 22;
[vMaxs1, vGap_maxs1, mS_maxs1, mT_maxs1, vpi_maxs1] = fMaxFilter(vps1, ctau);
[vMaxs2, vGap_maxs2, mS_maxs2, mT_maxs2, vpi_maxs2] = fMaxFilter(vps2, ctau);

[vMax, vGap_max, mS_max, mT_max, vpi_max] = fMaxFilter(vp, ctau);
plot( vdates(23:end), [vps2(23:end), vMaxs2 ])
plot( vdates(end-165:end), [vp(end-165:end), vps1(end-165:end), vps2(end-165:end)]); hold on;
plot( vdates(end-165:end), [vMax(end-165:end), vMaxs1(end-165:end), vMaxs2(end-165:end)]); hold off;

%% Figure 9 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
g= figure();
plot( vdates(end-165:end), -vGap_max(end-165:end), LineStyle= ":", ...
    LineWidth = 2, Color="green" );  hold on;
plot( vdates(end-165:end), -vGap_maxs1(end-165:end), LineStyle= "--", ...
    LineWidth = 2, Color="red" );  hold on; 
plot( vdates(end-165:end), -vGap_maxs2(end-165:end), LineStyle= "-", ...
    LineWidth = 2, Color="blue" );  hold off; 
legend("$d_t$", "Lowpass (cutoff $2\pi/5$)", "Lowpass  (cutoff $2\pi/22$)", Interpreter = "Latex") 
TT = timetable(vp, vps1, vps2, [NaN(ctau,1); -vGap_max],[NaN(ctau,1); -vGap_maxs1],[NaN(22,1); -vGap_maxs2] , 'RowTimes', vdates);
writetimetable(TT, "sp500_robust.xlsx" );
