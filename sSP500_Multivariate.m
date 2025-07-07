%% Some stylized facts: The time series of drawdowns, dt, and drawup, ut, with τ = 22 have been computed for the
%  N = 501 constituent stocks of S&P500, along with the states s_t_plus and s_t_minus, 
clear vars
load('sp500.mat', 'PT') 
mPrice = table2array(PT);
mY     = price2ret(mPrice); % 
vdates = PT.Index;
[cn,cN]= size(mPrice); %  
mp     = log(mPrice); 
%% compute drawdown, and other variables for constituent series with tau = 22
ctau     = 22;   % horizon tau = 22
mGap     = NaN(cn,cN); 
mGap_max = NaN(cn,cN); mGap_min = NaN(cn,cN);
mMax     = NaN(cn,cN); mMin     = NaN(cn,cN);
mPi_max  = NaN(ctau+1,cN); mPi_min = NaN(ctau+1,cN);
vGapMean = NaN(1,cN);
mC       = NaN(cn,cN);
mRange   = NaN(cn,cN); mprec = NaN(ctau+2,cN); mpexp = NaN(ctau+2,cN);
vTmax00  = NaN(1,cN); vTmaxqq = NaN(1,cN);
vTmin00  = NaN(1,cN); vTminqq = NaN(1,cN);
aS       = NaN(ctau+1,cn-ctau,cN);
[cn,cN]= size(mPrice);
for k = 1:cN
    disp(['Series =', num2str(k)]) 
    vy = log(mPrice(:,k)) + 10e-15*randn(cn,1);  % note: ties must be broken at random
    [vMax, vGap_max, mS_max, mT_max, vpi_max] = fMaxFilter(vy, ctau);
    aS(:,:,k)       = mS_max;
    [vMin, vGap_min, mS_min ,mT_min, vpi_min] = fMaxFilter(-vy, ctau);
    vTmax00(k) = mT_max(1,1); vTmaxqq(k) = mT_max(end,end); 
    vTmin00(k) = mT_min(1,1); vTminqq(k) = mT_min(end,end);
    mPi_max(:,k) = vpi_max; mPi_min(:,k) = vpi_min;
    mMax(ctau+1:end,k) = vMax; mMin(ctau+1:end,k) = vMin;
    % note vG^- = -vGap_min (change sign)
    vm  = zeros(ctau+1,1);
    for i=1:ctau
        vdy_i   = vy(i+1:end) - vy(1:end-i); % 
        vm(i+1) = mean(abs(vdy_i));
    end
    dGapMean = dot(vm,  vpi_min- vpi_max);
    vGapMean(k) = dGapMean;
    [vpi_max vpi_min];
    vT    =  0.5*(vMax-vMin);
    vC    =  0.5*(vGap_max-vGap_min) ;
    vGap = 0.5*(vGap_max-vGap_min  ) ;
    mGap(ctau+1:end,k) = vGap; mGap_max(ctau+1:end,k) = vGap_max; mGap_min(ctau+1:end,k) = vGap_min;
    mC(ctau+1:end,k) = vC;
    mRange(ctau+1:end,k) = vMax+vMin;
    %
    dp = 1; dpe = 1;
    vprec = NaN(ctau+1,1);  % recession duration probs
    vpexp = NaN(ctau+1,1);  % expansion duration probs
    vprec(1) = dp * mT_max(1,1); vpexp(1) = dpe * mT_min(1,1);
    for r = 2:(ctau+1)
        dp = dp * mT_max(r-1,r);   dpe = dpe * mT_min(r-1,r);
        vprec(r) = dp*mT_max(r,1); vpexp(r) = dpe*mT_min(r,1);
    end
    vprec(ctau+2) = 1-sum(vprec(1:ctau+1)); vpexp(ctau+2) = 1-sum(vpexp(1:ctau+1));
    mprec(:,k)=vprec; mpexp(:,k)=vpexp;
end
%% drawdown dup and maximum dd
mddown = -mGap_max;
mdup   = -mGap_min;

%% cross-sectional correlation analysis of features
vsel = sum(ismissing(mPrice)) == 0; % select complete series
vmaxdd = max(mddown);
vsr    = mean(mY(:,vsel),1)./sqrt(var(mY(:,vsel),1));  % Sharpe ratios
%  
vpi0   = mPi_max(1,vsel); vpi0m = mPi_min(1,vsel); 
%
vmaxdd = max(mddown(:,vsel)); vmeandd = mean(mddown(ctau+1:end,vsel)); vmediandd = median(mddown(ctau+1:end,vsel));
vmaxdu = max(mdup(:,vsel));   vmeandu = mean(mdup(ctau+1:end,vsel)); vmediandu = median(mdup(ctau+1:end,vsel)); 
vCalmar = mean(mY(:,vsel),1)./vmaxdd;
vmean   = mean(mY(:,vsel),1);
vstdev  = sqrt(var(mY(:,vsel),1));
%
vdurdd = (0:ctau)*(mprec(1:ctau+1,vsel)./(1-mprec(end,vsel)));
vdurdu = (0:ctau)*(mpexp(1:ctau+1,vsel)./(1-mpexp(end,vsel)));
msurvdd =1-cumsum(mprec,1);
msurvdu =1-cumsum(mpexp,1);
%% Concordance and correlation analysis
%  Cramer's V  and concordance Goodman-Kruskal
vs  = 0:ctau;
mVCramer = eye(cN); mcorr_s = eye(cN); mGK= eye(cN);
for i = 1:cN
    disp(['Cramer and GK for series =', num2str(i)]) 
    mSi = aS(:,:,i);
    vsi = vs * mSi;
    for j = i+1:cN
        mSj   = aS(:,:,j);
        vsj   = vs * mSj;
        mCTab = mSi*mSj'; dtot = sum(sum(mCTab));
        vctot = sum(mCTab,1); vrtot = sum(mCTab,2);
        mE    = (mCTab-vrtot*vctot/dtot)./sqrt(vrtot*vctot/dtot);
        dv    = sqrt((sum(sum(mE.^2))/dtot)/ctau);
        mVCramer(i,j) = dv ;
        mVCramer(j,i) = dv ;
        dc            = corr(vsi',vsj');
        mcorr_s(i,j)  =  dc;
        mcorr_s(j,i)  = dc;
        dgamma  = fGoodmanKruskal(mCTab);
        mGK(i,j) = dgamma ;
        mGK(j,i) = dgamma ;
    end
end
mCgk = tril(mGK(vsel,vsel),-1);  vcgk = mCgk(:);  vcgk(find(vcgk==0))=[]; 
%  Pearson, Spearman  
mcorr = corr(mddown(23:end,vsel), "type", 'Pearson','rows','pairwise')
mcorr_sp = corr(mddown(23:end,vsel), "type", 'Spearman' , 'rows','pairwise')
mcorr_s = mcorr_s(vsel, vsel) % Pearson correlation of S_t^+ 
mC = tril(mcorr,-1);  vc = mC(:);  vc(find(vc==0))=[]; 
mCs = tril(mcorr_s,-1);  vcs = mCs(:);  vcs(find(vcs==0))=[]; 
mCsp = tril(mcorr_sp,-1);  vcsp = mCsp(:);  vcsp(find(vcsp==0))=[]; 
mCv = tril(mVCramer(vsel, vsel),-1);  vcv = mCv(:);  vcv(find(vcv==0))=[]; 
%% Figure 3
mX =[vc vcsp vcs vcv vcgk];
sc = 'Component '
for k = 1:5
    x = mX(:,k);
    for j=1:5
        y = mX(:,j);
        subplot(5, 5, 5*(k-1)+j)
        if ne(k,j)
             plot(x,y,'b.')
         else 
            histogram(mX(:,k), Normalization="pdf")
        end
          
    end
end
gsd = figure() 
[~,ax]=plotmatrix(mX); 
ax(1,1).YLabel.String='Pearson d_t';    ax(2,1).YLabel.String='Spearman'; 
ax(3,1).YLabel.String='Pearson s_t^+';  ax(4,1).YLabel.String='Cramer'; 
ax(5,1).YLabel.String='Goodman-Kruskal';
ax(5,1).XLabel.String='Pearson d_t';    ax(5,2).XLabel.String='Spearman'; 
ax(5,3).XLabel.String='Pearson s_t^+';  ax(5,4).XLabel.String='Cramer'; 
ax(5,5).XLabel.String='Goodman-Kruskal'; 
orient(gsd,'landscape')
gsd.Position = [0 -0 1275 875];  
print(gsd,'gSP500_Association','-dpdf',   '-r250')     
%% Characteristics and heatmap Fig 4
mC = [vmeandd', vmediandd', vmaxdd', vmeandu', vmediandu', vmaxdu', ...
      vmean', vstdev', vsr', vCalmar', ...
      vpi0',  vTmax00(vsel)', vpi0m', vTmin00(vsel)', msurvdd(22,vsel)', msurvdu(22,vsel)'];
slab ={'CDaR(22,0)', 'CDaR(22,0.5)', 'CDaR(22,1)', 'CUaR(22,0)', 'CUaR(22,0.5)', 'CUaR(22,1)', '$\hat{\mu}$', '$\hat{\sigma}$', 'Sharpe', 'Calmar', '$\hat{\pi}_0^+$', ...
    '$\hat{p}_{00}+$', '$\hat{\pi}_0^-$', '$\hat{p}_{00}^-$', '$P(\mathcal{D}_d>22)$', '$P(\mathcal{D}_u>22)$' }
% figure 
gsd = figure("Name","Heat of the moment")
h = heatmap(slab, slab, corr(mC, 'type','Spearman'))
h.Interpreter = 'latex'; 
h.FontSize = 10
h.CellLabelFormat = '%.2f'
h.Colormap = sky 
orient(gsd,'landscape')
gsd.Position = [-10 -10 1275 875]; 
print(gsd,'gSP500_heatmap','-dpdf',    '-r600')     