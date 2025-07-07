%% Figure 1 and Table 1 SP500
clear all;
%% Read data
TT0     = readtimetable(append('^GSPC.csv'), VariableNamingRule='preserve');
PT      = TT0(:,4);  
vPrice  = table2array(PT);
vp      = log(vPrice); 
vp_high = log(table2array(TT0(:,2)));   
vp_low  = log(table2array(TT0(:,3)));
vp_open = log(table2array(TT0(:,1))); 
vy      = price2ret(vPrice); % 
vdates  = PT.Index;
[cn,cN] = size(vPrice); %  
%% max filter, drawdown and drawup
ctau    = 22;
rng     = 5211314;
vpe = vp  + 10e-15*randn(cn,1);  % add tiny random error to break ties
[vMax, vGap_max, mS_max, ~, ~] = fMaxFilter(vpe, ctau);
[vMin, vGap_min, mS_min ,~, ~] = fMaxFilter(-vpe, ctau);
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
%% Summary statistics printed in Table 1 Columns 1-5 (tau = 22)
mDescr = [min(mY); quantile(mY, [0.25, 0.50, 0.75]); max(mY); mean(mY); std(mY); skewness(mY); kurtosis(mY)]
writematrix(mDescr, "sp500_Table1_1.xlsx" );
%% save data used for figure 1 (The published graph is done in OxMetrics)
TT = timetable(vp, vMax, vMin, vdd, vdu, vrange, [NaN(ctau, 1); vsplus'], [NaN(ctau, 1); vsminus'],'RowTimes', vdates);
writetimetable(TT, "sp500_Data4Figure1.xlsx" );
%% Select period
t1      = '01-Jan-2023'; 
t2      = '30-Aug-2023';
S       = timerange(t1, t2);
TTs     = TT(S,:);
%% Graphic parameters
set(0,'DefaultAxesFontName','Times New Roman');
set(0,'DefaultTextFontName','Times New Roman');
set(0,'DefaultAxesFontSize',15);
set(0,'DefaultTextFontSize',15);
g = figure("Name","Figure 1"); 
subplot(3,1,1);
    plot(TTs.Time , [TTs.vp , TTs.vMax, TTs.vMin])
    legend("$p_t$", "$m_t^+$", "$m_t^-$", 'Location', 'southeast',   Interpreter="latex") 
    title("i. Log-prices, $\tau$-maxima and minima", Interpreter="latex");
subplot(3,1,2);
    plot(TTs.Time , [TTs.vdd, TTs.vdu])
    legend("$d_t$", "$u_t$", Interpreter="latex")     
    title("ii. Drawdown and Drawup", Interpreter="latex");
subplot(3,1,3);
    plot(TTs.Time , [TTs.Var7, TTs.Var8], '*')
    legend("$s_t^+$", "$s_t^-$", Interpreter="latex")     
    title("iii. Current lead time from maximum and minimum", Interpreter="latex");    
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% max filter, drawdown and drawup for tau = 65
ctau    = 65;
rng     = 5211314;
[vMax, vGap_max, mS_max, mT_max, vpi_max] = fMaxFilter(vpe, ctau);
[vMin, vGap_min, mS_min ,mT_min, vpi_min] = fMaxFilter(-vpe, ctau);
vdd65     = NaN(cn,1);    % drawdown
vdd65(ctau+1:end)     = -vGap_max;
vdu65     = NaN(cn,1);    % drawup
vdu65(ctau+1:end)     = -vGap_min;
vMax65  = [NaN(ctau,1);  vMax];
vMin65  = [NaN(ctau,1); -vMin];
vrange65  = vdd-vdu;

%%  Summary statistics printed in Table 1 Columns 1-5 (tau = 22)
vsplus65 = (0:ctau)*mS_max; vsminus65 = (0:ctau)*mS_min;
mY = [vdd65(ctau+1:end), vdu65(ctau+1:end), vsplus65', vsminus65'];
mDescr65 = [min(mY); quantile(mY, [0.25, 0.50, 0.75]); max(mY); mean(mY); std(mY); skewness(mY); kurtosis(mY)]
writematrix(mDescr65, "sp500_Table1_2.xlsx" );