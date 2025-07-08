%% SP500 index illustrations: Max, min, drawdown (DD) etc (Figure 1) and 
%  descriptive statistics (Table 1) 
%  DD and DU duration probabilities (Table 3)
%  Upper and lower bounds for DD (Figure 8)
clear vars;
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
[vMax, vGap_max, mS_max, mT_max, ~] = fMaxFilter(vpe, ctau);
[vMin, vGap_min, mS_min ,mT_min, ~] = fMaxFilter(-vpe, ctau);
vdd     = NaN(cn,1);    % drawdown
vdd(ctau+1:end)     = -vGap_max;
vdu     = NaN(cn,1);    % drawup
vdu(ctau+1:end)     = -vGap_min;
vMax  = [NaN(ctau,1);  vMax];
vMin  = [NaN(ctau,1); -vMin];
vrange  = vdd-vdu;
%% Summary statistics printed in Table 1 Columns 1-5 (tau = 22)
vsplus = (0:ctau)*mS_max; vsminus = (0:ctau)*mS_min;
mY = [vdd(ctau+1:end), vdu(ctau+1:end), vsplus', vsminus'];
% 
mDescr = [min(mY); quantile(mY, [0.25, 0.50, 0.75]); max(mY); mean(mY); std(mY); skewness(mY); kurtosis(mY)];
writematrix(mDescr, "sp500_Table1_1.xlsx" );
disp("------------------ Table 1 part 1 (tau=22) ---------------")
column_names = {'drawdown', 'drawup', 's_t_plus', 's_t_minus'};
my_table = array2table(mDescr, 'VariableNames', column_names);

% Display the table
disp(my_table);
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
%% Duration analysis (table 3)
dp = 1; dpe = 1;
vpdd = NaN(ctau+1,1);  %   duration probs
vpdu = NaN(ctau+1,1);  %  
vpdd(1) = dp * mT_max(1,1);   vpdu(1) = dpe * mT_min(1,1);
for r = 2:(ctau+1)
    dp = dp * mT_max(r-1,r);      dpe = dpe * mT_min(r-1,r);
    vpdd(r) = dp*mT_max(r,1); vpdu(r) = dpe*mT_min(r,1);
end
vpdd(ctau+2) = 1-sum(vpdd(1:ctau+1)); vpdu(ctau+2) = 1-sum(vpdu(1:ctau+1));
vpdds = 1-cumsum(vpdd(1:ctau+1));
vpdus = 1-cumsum(vpdu(1:ctau+1));
mProbs =[vpdd(1:ctau+1) vpdds vpdu(1:ctau+1) vpdus];
writematrix(mProbs, "sp500_duration.xlsx" );
disp("------------------------- Table 3 ------------------------------")
column_names = {'P(D_d=k)', 'P(D_d>k)', 'P(D_u=k)', 'P(D_u>k)'};
my_table = array2table(mProbs, 'VariableNames', column_names);
disp(my_table);
%% Upper and lower bound for dd (Figure 8)
vpe_h = vp_high+ 10e-15*randn(cn,1);
[vMax, ~, ~, ~, ~] = fMaxFilter(vpe_h, ctau);
vMax  = [NaN(ctau,1);  vMax];
vdd_ub = vMax-vp_low;  % upper bound
% lower bound for dd 
vpe_l = vp_low+ 10e-15*randn(cn,1);
[vMax, ~, ~, ~, ~] = fMaxFilter(vpe_l, ctau);
vMax  = [NaN(ctau,1);  vMax];
vdd_lb = max(vMax-vp_high,0);  % lower bound
mY = [vdd, vdd_ub, vdd_lb];
mDescr = [min(mY); quantile(mY, [0.25, 0.50, 0.75]); max(mY); nanmean(mY); nanstd(mY); skewness(mY); kurtosis(mY)];
figure("Name", "Figure 8")
cm = 165;
plot(vdates(end-cm:end),  vdd(end-cm:end), LineStyle="--", Marker=".", Color='red' ); hold on;
plot(vdates(end-cm:end),  vdd_ub(end-cm:end), LineStyle="-", Color='blue' ); hold on;
plot(vdates(end-cm:end),  vdd_lb(end-cm:end), LineStyle="-", Color='green' ); hold off;
legend("$d_t$","Upper bound", "Lower bound", 'Location', 'northeast', Interpreter = 'latex' )
% save results
TT = timetable(vdd, vdd_ub, vdd_lb, 'RowTimes', vdates);
writetimetable(TT,  "sp500_dd_bounds.xlsx" );
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% max filter, drawdown and drawup for tau = 65
ctau    = 65;
[vMax, vGap_max, mS_max, mT_max, vpi_max] = fMaxFilter(vpe, ctau);
[vMin, vGap_min, mS_min ,mT_min, vpi_min] = fMaxFilter(-vpe, ctau);
vdd65     = NaN(cn,1);    % drawdown
vdd65(ctau+1:end)     = -vGap_max;
vdu65     = NaN(cn,1);    % drawup
vdu65(ctau+1:end)     = -vGap_min;
vMax65  = [NaN(ctau,1);  vMax];
vMin65  = [NaN(ctau,1); -vMin];
%%  Summary statistics printed in Table 1 Columns 1-5 (tau = 22)
vsplus65 = (0:ctau)*mS_max; vsminus65 = (0:ctau)*mS_min;
mY = [vdd65(ctau+1:end), vdu65(ctau+1:end), vsplus65', vsminus65'];
mDescr65 = [min(mY); quantile(mY, [0.25, 0.50, 0.75]); max(mY); mean(mY); std(mY); skewness(mY); kurtosis(mY)];
disp("--------------- Table 1 part 2 (tau=65) ------------------")
column_names = {'drawdown', 'drawup', 's_t_plus', 's_t_minus'};
my_table = array2table(mDescr65, 'VariableNames', column_names);

% Display the table
disp(my_table);
writematrix(mDescr65, "sp500_Table1_2.xlsx" );