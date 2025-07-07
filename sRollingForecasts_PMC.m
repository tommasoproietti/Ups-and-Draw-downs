%% Rolling forecasts Particle Monte Carlo (GARCH model) 
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
figure('Name', 'SP log returns');
plot(vdatesy, vy), title('Time series') 
%% drawdown
ctau    = 22;   % max filter 
[vMax, vGap_max, mS_max, mT_max, vpi_max] = fMaxFilter(vpe, ctau);
[vMin, vGap_min, mS_min ,mT_min, vpi_min] = fMaxFilter(-vpe, ctau);
vdd     = NaN(cn,1);    % drawdown 
vdd(ctau+1:end) = -vGap_max;  % observed value
vdu     = NaN(cn,1);    % drawdown 
vdu(ctau+1:end) = vGap_min;  % observed value
vc = vdu+vdd;
%% Specifications 
Mdl = garch('GARCHLags',1,'ARCHLags',1,'Offset',NaN);  % garch model
Mdl.Distribution = "t";
cL    = 22;         %   forecast horizon
cM    = 2000;         %   number of MC replications
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cT = 4000;      % length of window
[EstMdl_garch, ~, ~,info] = estimate(Mdl, vy(1:cT),'Display','off');
cr = cny-cT-cL; % number of rolling window
mPar = NaN(5,cr);   
mfError = NaN(cr, cL); mfRelError = NaN(cr, cL); mfDens = NaN(cr,cL);
mPIT = NaN(cr,cL); mPITe = NaN(cr,cL); mf = NaN(cr,cL)
for i = 0:cr-1
    clc;
    disp(['Processing rolling forecast ', num2str(i), ' out of ', num2str(cr)])
    vys = vy(i+1:cT+i);        % rolling sample
    %% Garch model
    if any(i==0:10:(cr-1))
         [EstMdl_garch0, ~, ~,info] = estimate(Mdl, vys,'Display','off');
         if (info.exitflag == 2) 
             EstMdl_garch = EstMdl_garch0;
         end
    end
    dalpha0 = EstMdl_garch.Constant;    dalpha1 = EstMdl_garch.ARCH{1};
    dbeta1  = EstMdl_garch.GARCH{1};  dnu     = EstMdl_garch.Distribution.DoF;
    dcmean  = EstMdl_garch.Offset;
    mPar(:,i+1) = [dalpha0; dalpha1; dbeta1; dnu; dcmean];
    vCond_var    = infer(EstMdl_garch,vys); % conditional variances
    vCond_stdev  = sqrt(vCond_var);
    %% Monte Carlo Predictions 
    myfor = NaN(cM, cL);   % returns MC predictions
    dh    = vCond_var(end);
    mh    = NaN(cM, cL);
    dy_s  = vys(end);
    %
    vdd_obs = vdd(cT+i+2:cT+i+1+cL);   % observed values
    mdd_sim = NaN(cL,cM);   % dd predictions
    for m = 1:cM
        for l=1:cL
            dh = dalpha0 + dalpha1 * (dy_s-dcmean)^2 + dbeta1 * dh;
            dy_s = dcmean + sqrt(dh) * trnd(dnu,1)*sqrt((dnu-2)/dnu);
            myfor(m,l) = dy_s;
            mh(m,l)= dh;
        end
        [~, vgs, ~, ~, ~] = fMaxFilter([vpe(cT+i+1-ctau+1:cT+i+1);  cumsum(myfor(m,:))'+vpe(cT+i+1)], ctau);  % filter extended series
        mdd_sim(:,m)  = -vgs;
    end

    % mean forecasts
    vdd_f = mean(mdd_sim,2);
 
    mf(i+1,:) = vdd_f';
    mfError(i+1,:) = (vdd_obs-vdd_f)';
end 
save sRollingForecasts_PMC.mat
%load sRollingForecasts.mat
