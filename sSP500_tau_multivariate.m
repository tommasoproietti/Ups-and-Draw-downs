%% Computes tau-drawdowns and tau-drawups of component series for tau ranging from 1 to 65
clear vars
load('sp500.mat', 'PT') 
%%
mPrice = table2array(PT);
mY     = price2ret(mPrice); % 
vdates = PT.Index;
[cn,cN]= size(mPrice); %  
mp     = log(mPrice); 
%% 
ctau    = 65;                       % 3* (365.25-2*(365.25/7))/12; % n days in a quarter
rng     = 5211314;
add = NaN(cn, ctau, cN); adu = NaN(cn, ctau, cN);
parfor k = 1:cN
    disp(['Series =', num2str(k)]) 
    vpe = mp(:,k)+10e-15*randn(cn,1);
    mdd = NaN(cn, ctau); mdu = NaN(cn, ctau);
    for tau = 1:ctau
        disp(['Series =', num2str(k), ' ** tau value =', num2str(tau)]);
        [vMax, vGap_max, mS_max, mT_max, vpi_max] = fMaxFilter(vpe, tau);
        [vMin, vGap_min, mS_min ,mT_min, vpi_min] = fMaxFilter(-vpe, tau);
        vdd     = NaN(cn,1);    % drawdown
        vdd(tau+1:end)     = -vGap_max;
        vdu     = NaN(cn,1);    % drawup
        vdu(tau+1:end)     = -vGap_min;
        mdd(:,tau) = vdd;  mdu(:,tau) = vdu;
    end
    add(:,:,k) = mdd; adu(:,:,k) = mdu;
end
%%
save sSP500_tau_multivariate.mat
%%

% mmeandd = squeeze(nanmean(add(ctau+1:end,:,:),1));
% mmeandu = squeeze(nanmean(adu(ctau+1:end,:,:),1));
% mvardd = squeeze(nanvar(add(ctau+1:end,:,:),1));
% mvardu = squeeze(nanvar(adu(ctau+1:end,:,:),1));
% mmeandd ./ sqrt(mvardd)
% plot(1:ctau, mmeandd ./ sqrt(mvardd)) 
% plot(1:ctau, sqrt(mvardd)./mmeandd) 
% plot(1:ctau,diff([-mmeandd+mmeandu])) 
% plot(1:ctau, [-mmeandd+mmeandu]) 
% mesh(1:ctau, 1:501, -mmeandd'); hold on;
% mesh(1:ctau, 1:501, mmeandu'); hold off;
% mesh(1:ctau, 1:501, mmeandu'-mmeandd');  
% mmaxdd = squeeze(max(add(ctau+1:end,:,:)));
% plot(log(1:ctau),  [log(mmeandd)] ); hold on;
% mycolors = [1 0 0; 0 1 0; 0 0 1];
% ax = gca; 
% ax.ColorOrder = mycolors;
% plot(log(1:ctau),  [log(mmeandu)]); hold off;
% ax = gca; 
% ax.ColorOrder = mycolors;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%
% load sSP500_tau_multivariate.mat
% ctau = 60
% vx  = log(1:ctau)'- mean(log(1:ctau)'); 
% mx1 = vx; 
% mxo = vx;
% cp =5;
% for k=2:cp
%     vxadd = vx.^k;
%     mod  = fitlm(mx1, vxadd); ve = mod.Residuals.Raw; 
%     mx1 = [mx1, vxadd];     mxo = [mxo, ve];
% end
% mmeandd = squeeze(nanmean(add(ctau+1:end,:,:),1));
% mmeandu = squeeze(nanmean(adu(ctau+1:end,:,:),1));
% 
% mmeandd = rmmissing(mmeandd,2);
% cN = size(mmeandd,2)
% mbeta_d = nan(cN, cp+1);
% mbeta_u = nan(cN, cp+1);
% for i = 1:cN
%     mod = fitlm(mxo, log(mmeandd(:,i)));
%     % mod = fitlm(mxo, log(mmaxdd(:,i)));
%     mbeta_d(i,:) = mod.Coefficients.Estimate';
%     mod = fitlm(mxo, log(mmeandu(:,i)));
%     mbeta_u(i,:) = mod.Coefficients.Estimate';
% end
% 
% histogram(mbeta_u(:,2), Normalization="pdf"); hold on;
% histogram(mbeta_d(:,2), Normalization="pdf"); hold off;
% % biplot(prc)
% 
% histogram(mbeta_u(:,2)-mbeta_d(:,2)  )
% 
% histogram2(mbeta_d(:,2), mbeta_u(:,2)  )
% plot(mbeta_d(:,2), mbeta_u(:,2),'.')
% exp(mbeta_u(:,1)-mbeta_d(:,1))
% mZ = [mbeta_d mbeta_u]
% [coefs,score, v] = pca(zscore(mZ));
% biplot(coefs(:,1:2),'Scores',score(:,1:2) );
% 
% histogram(mbeta_u(:,2)-mbeta_d(:,2), Normalization="pdf")