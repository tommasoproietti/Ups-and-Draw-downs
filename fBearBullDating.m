function [vpeak, vtrough,vbull_ind, vbear_ind, vbull_dur, vbear_dur] = fBearBullDating(ctau, cK, vy)
%% 
cn = length(vy); %  
%% max filter, drawdown and drawup
vpe = vy + 10e-15*randn(cn,1);
[~, ~, mS_max, ~, ~] = fMaxFilter(vpe, ctau);
[~, ~, mS_min ,~, ~] = fMaxFilter(-vpe, ctau);
vSplus = (0:ctau) * mS_max;
vSminus = (0:ctau) * mS_min;
%[vSminus; vSplus];
cns = length(vSminus);
%% step 1: identification of candidate peaks and troughs
vpeak = zeros(cns,1); vtrough = zeros(cns,1);
for t = cK+1:cns-cK
    if     and(all(vSplus(t:t+cK) == (0:cK)),   all(vSminus(t-1:t-cK) > 0) )

        vpeak(t) = 1;
    elseif  and(all(vSminus(t:t+cK) == (0:cK)),   all(vSplus(t-1:t-cK) > 0) )
        vtrough(t) = 1;
    end
end
%% Step 2: Enforce alternation of turning points
vloc_t = find(vtrough==1); % indices of troughs
vloc_p = find(vpeak == 1); % indices of peaks
vp_loc_t = vpe(vloc_t+ctau); vp_loc_p = vpe(vloc_p+ctau);
mA_t = [vloc_t, vp_loc_t, zeros(length(vloc_t),1)];
mA_p = [vloc_p, vp_loc_p,  ones(length(vloc_p),1)];
[~, cselpt] = min([vloc_t(1),vloc_p(1)]);
if cselpt == 2
    mA  = sortrows([mA_p; mA_t]);
else
    mA = sortrows([mA_t; mA_p]);
end
%
cruns = (diff(find([1;diff(mA(:,3));1]))); % number of runs
ccruns = [0;cumsum(cruns)] ;
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

%% duration analysis
vdur = diff(find(vtrough+vpeak ==1));
cd = length(vdur);
if cselpt == 2
    vbull_dur = vdur(2:2:cd);
    vbear_dur = vdur(1:2:cd);
else
    vbull_dur = vdur(1:2:cd);
    vbear_dur = vdur(2:2:cd);
end

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

end 