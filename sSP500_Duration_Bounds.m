%% script used for producing table 3 and figure 8 
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
vp_high = log(table2array(TT0(:,2)));   
vp_low  = log(table2array(TT0(:,3)));
vp_open = log(table2array(TT0(:,1))); 
vy      = price2ret(vPrice); % 
vdates  = PT.Index;
[cn,cN] = size(vPrice); %  
plot(vdates, [vp_open, vp vp_high, vp_low])
%% max filter, drawdown and drawup
ctau    = 22;
rng     = 5211314;
vpe = vp  + 10e-15*randn(cn,1);
[vMax, vGap_max, mS_max, mT_max, vpi_max] = fMaxFilter(vpe, ctau);
[vMin, vGap_min, mS_min ,mT_min, vpi_min] = fMaxFilter(-vpe, ctau);
vdd     = NaN(cn,1);    % drawdown
vdd(ctau+1:end)     = -vGap_max;
vdu     = NaN(cn,1);    % drawup
vdu(ctau+1:end)     = -vGap_min;
vMax  = [NaN(ctau,1);  vMax];
vMin  = [NaN(ctau,1); -vMin];
vrange  = vdd-vdu;
%% descriptive stats
vsplus = (0:ctau)*mS_max; vsminus = (0:ctau)*mS_min;
mY = [vdd(ctau+1:end), vdu(ctau+1:end), vsplus', vsminus'];
mDescr = [min(mY); quantile(mY, [0.25, 0.50, 0.75]); max(mY); mean(mY); std(mY); skewness(mY); kurtosis(mY)];
%% Duration analysis (table 3)
vpe = vp  + 10e-15*randn(cn,1);
[vMax, vGap_max, mS_max, mT_max, vpi_max] = fMaxFilter(vpe, ctau);
[vMin, vGap_min, mS_min ,mT_min, vpi_min] = fMaxFilter(-vpe, ctau);
dp = 1; dpe = 1;
vpdd = NaN(ctau+1,1);  %   duration probs
vpdu = NaN(ctau+1,1);  %  
vpdd(1) = dp * mT_max(1,1); vpdu(1) = dpe * mT_min(1,1);
for r = 2:(ctau+1)
    dp = dp * mT_max(r-1,r);   dpe = dpe * mT_min(r-1,r);
    vpdd(r) = dp*mT_max(r,1); vpdu(r) = dpe*mT_min(r,1);
end
vpdd(ctau+2) = 1-sum(vpdd(1:ctau+1)); vpdu(ctau+2) = 1-sum(vpdu(1:ctau+1));
vpdds = 1-cumsum(vpdd(1:ctau+1));
vpdus = 1-cumsum(vpdu(1:ctau+1));
writematrix([vpdd(1:ctau+1) vpdds vpdu(1:ctau+1) vpdus], "sp500_duration.xlsx" );
[vpdd(1:ctau+1) vpdds vpdu(1:ctau+1) vpdus]
%% Evaluate the dd on the 4 series of daily prices: open close low high
mp = [vp_open, vp vp_high, vp_low];
mdd = NaN(cn,4);    % drawdown
for i = 1:4
    vpe = mp(:,i)  + 10e-15*randn(cn,1);
    [vMax, vGap_max, mS_max, mT_max, vpi_max] = fMaxFilter(vpe, ctau);
    mdd(ctau+1:end,i)     = -vGap_max;
end
% plot(vdates, mdd)
% plot(vdates, range(mdd'))
% plot(vdates, mdd(:,2)); hold on;
% plot(vdates, mdd(:,2)+.5*range(mdd')');  hold off;
%% upper and lower bound for dd 
vpe = vp_high+ 10e-15*randn(cn,1);
[vMax, vGap_max, mS_max, mT_max, vpi_max] = fMaxFilter(vpe, ctau);
vMax  = [NaN(ctau,1);  vMax];
vdd_ub = vMax-vp_low;  % upper bound
% lower bound for dd 
vpe = vp_low+ 10e-15*randn(cn,1);
[vMax, vGap_max, mS_max, mT_max, vpi_max] = fMaxFilter(vpe, ctau);
vMax  = [NaN(ctau,1);  vMax];
vdd_lb = max(vMax-vp_high,0);  % lower bound
mY = [vdd, vdd_ub, vdd_lb];
mDescr = [min(mY); quantile(mY, [0.25, 0.50, 0.75]); max(mY); nanmean(mY); nanstd(mY); skewness(mY); kurtosis(mY)];
figure()
cm = 165
plot(vdates(end-cm:end),  vdd(end-cm:end), LineStyle="--", Marker=".", Color='red' ); hold on;
plot(vdates(end-cm:end),  vdd_ub(end-cm:end), LineStyle="-", Color='blue' ); hold on;
plot(vdates(end-cm:end),  vdd_lb(end-cm:end), LineStyle="-", Color='green' ); hold off;
legend("$d_t$","Upper bound", "Lower bound", 'Location', 'northeast', Interpreter = 'latex' )
% save results
TT = timetable(vdd, vdd_ub, vdd_lb, 'RowTimes', vdates);
writetimetable(TT,  "sp500_dd_bounds.xlsx" );
 
 
