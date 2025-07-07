%% Script for dating phases of SP500 and for producing Figure 2 
clear all;
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
vy      = price2ret(vPrice); % 
vdates  = PT.Index;
[cn,cN] = size(vPrice); %  
%% max filter, drawdown and drawup
ctau    = 65;
rng     = 5211314;
vpe = vp  + 10e-15*randn(cn,1);
[vMax, vGap_max, mS_max, mT_max, vpi_max] = fMaxFilter(vpe, ctau);
[vMin, vGap_min, mS_min ,mT_min, vpi_min] = fMaxFilter(-vpe, ctau);
vSplus = (0:ctau) * mS_max
vSminus = (0:ctau) * mS_min
[vSminus; vSplus]
cK = 22;
cns = length(vSminus);
%% Preliminary identification of turning points
vpeak = zeros(cns,1); vtrough = zeros(cns,1);
for t = cK+1:cns-cK 
    if     and(all(vSplus(t:t+cK) == (0:cK)),   all(vSminus(t-1:t-cK) > 0) )  
 
        vpeak(t) = 1;
    elseif  and(all(vSminus(t:t+cK) == (0:cK)),   all(vSplus(t-1:t-cK) > 0) ) 
     vtrough(t) = 1;
    end
end
figure("Name","Preliminary turning points")
plot(vdates(ctau+1:end), zscore(vp(ctau+1:end ))); hold on;
b1 = bar(vdates(ctau+1:end), -2+5*vpeak, -2); hold on
b2 = bar(vdates(ctau+1:end), -2+5*vtrough, -2); hold off
b1(1).BaseValue = -2;b2(1).BaseValue = -2;
%% Enforce alternation of turning points
vloc_t = find(vtrough==1); % indices of troughs
vloc_p = find(vpeak == 1); % indices of peaks
vp_loc_t = vpe(vloc_t+ctau); vp_loc_p = vpe(vloc_p+ctau);
mA_t = [vloc_t, vp_loc_t, zeros(length(vloc_t),1)];
mA_p = [vloc_p, vp_loc_p,  ones(length(vloc_p),1)];
[a cselpt] = min([vloc_t(1),vloc_p(1)]);
if cselpt == 2
    mA  = sortrows([mA_p; mA_t]);
else mA = sortrows([mA_t; mA_p]);
end
cruns = (diff(find([1;diff(mA(:,3));1]))); % number of runs
ccruns = [0;cumsum(cruns)] ;
[cruns diff(ccruns), ccruns(2:end)]
cm = length(cruns);
mAc = nan(cm,3);
for i = 1:cm
    if (cruns(i)== 1)
        mAc(i,:) = mA(ccruns(i+1),:);
    else
        ma = mA(ccruns(i)+1:ccruns(i+1),:);
        if (ma(1,3) == 0)
            [~,imin] = min(ma(:,2));
            mAc(i,:) = ma(imin,:);
        else
            [~,imax] = max(ma(:,2));
            mAc(i,:) = ma(imax,:);
        end
    end
end
mA_peak = mAc(2:2:cm, :); mA_trough = mAc(1:2:cm, :);
vpeak = zeros(cns,1); vtrough = zeros(cns,1);
vpeak(mA_peak(:,1)) = 1; vtrough(mA_trough(:,1)) = 1;
%% bull and bear indicators
vbull_ind = zeros(cns,1); vbear_ind = zeros(cns,1);
vtp_ind =  [0; mAc(:,1); cns];
for i = 1:2:cm
    if cselpt == 2
        vbull_ind(vtp_ind(i)+1:vtp_ind(i+1)) = 1;
        vbear_ind(vtp_ind(i+1)+1:vtp_ind(i+2)) = 1;
    else
        vbear_ind(vtp_ind(i)+1:vtp_ind(i+1)) = 1;
        vbull_ind(vtp_ind(i+1)+1:vtp_ind(i+2)) = 1;
    end
end
figure("Name","Bull and Bear Phases")
plot(vdates(ctau+1:end), zscore(vp(ctau+1:end )), LineWidth=2, Color=[.9, .6, 0]); hold on;
plot(vdates(ctau+1:end), zscore(vp(ctau+1:end )), LineWidth=2, Color=[.9, .6, 1]); hold on;
b1 = bar(vdates(ctau+1:end), -2+5*vpeak, "BarWidth",2); hold on
b2 = bar(vdates(ctau+1:end), -2+5*vtrough,   "BarWidth",2); hold on
b1(1).BaseValue = -2;b2(1).BaseValue = -2; 
b1.FaceColor = [0.1 0.1 0.9]; b2.FaceColor = [0.1 0.1 0.9];
b3 = bar(vdates(ctau+1:end), -2+5*vbear_ind,  "BarWidth",1 ); hold off
b3(1).BaseValue = -2 ;  b3.FaceColor = [0.9 0.9 0.9];
%% Figure 2
gsd = figure("Name","Bull and bull phases and diffusion index");
b1 = bar(vdates(ctau+1:end),   vpeak, "BarWidth",2); hold on;
b2 = bar(vdates(ctau+1:end),  vtrough,   "BarWidth",2); hold on;
%b1.FaceColor = [0.1 0.1 0.9]; b2.FaceColor = [0.1 0.1 0.9];
b1.FaceColor = 'b'; b2.FaceColor = 'k';
b3 = bar(vdates(ctau+1:end),  vbear_ind,  "BarWidth",1 ); hold on;
plot(vdates(ctau+1:end),  (vp(ctau+1:end)-min(vp))./range(vp), LineWidth=.5, Color='r'); hold on;
b3.FaceColor = [0.9 0.9 0.9];
load("sSP500_Peaks_Troughs_Multivariate.mat", "vbear_diffusion")
plot(vdates , vbear_diffusion, 'b', LineWidth=1);
orient(gsd,'landscape')
gsd.Position = [-50 -50 1350 850];  
print(gsd,'gSP500_BullBear','-dpdf',  '-r600')
%% Comparison bull and bear phases and drawdowns
vdrawup   = -vGap_min; % 
vdrawdown = - vGap_max; 
figure("Name","Bull and bear phases and drawdown");
b1 = bar(vdates(ctau+1:end),   vpeak, "BarWidth",2); hold on;
b2 = bar(vdates(ctau+1:end),  vtrough,   "BarWidth",2); hold on;
%b1.FaceColor = [0.1 0.1 0.9]; b2.FaceColor = [0.1 0.1 0.9];
b1.FaceColor = 'b'; b2.FaceColor = 'k';
b3 = bar(vdates(ctau+1:end),  vbear_ind,  "BarWidth",1 ); hold on;
plot(vdates(ctau+1:end),  (vp(ctau+1:end)-min(vp))./range(vp), LineWidth=.5, Color='r'); hold on;
b3.FaceColor = [0.9 0.9 0.9];
plot(vdates(ctau+1:end) ,  (vdrawup /max(vdrawup)), 'g', LineWidth=1); hold on;
plot(vdates(ctau+1:end) ,  (vdrawdown/max(vdrawdown)), 'b', LineWidth=1); hold off;

crosstab(vdrawup>0,vdrawdown>0,  vbear_ind,  vbull_ind) 
crosstab(vdrawup>0,  vbull_ind) 
crosstab(vdrawup>0,  vbear_ind) 
crosstab(vdrawdown>0,  vbear_ind) 
crosstab(vbear_ind,  vbull_ind) 
tabulate(vbear_ind)
tabulate(vbull_ind)
tabulate(vdrawdown>0)
tabulate(vdrawup>0)