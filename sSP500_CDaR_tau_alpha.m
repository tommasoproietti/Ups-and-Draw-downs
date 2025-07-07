%% Analisis of conditional drawdown at risk - Figures 5 and 6
clear all
load("sSP500_tau_multivariate.mat", "add", "adu", "vdates");
[cn ctau cN] = size(add);
load('sp500.mat', 'PT') 
mPrice = table2array(PT);
vsel = sum(ismissing(mPrice)) < (cn - 1928); % series must at least start in 2016
%% conditional dd at risk
valpha = 0:0.025:.99;
cq = length(valpha);
aCDaR = NaN(ctau, cq, cN); aCUaR = NaN(ctau, cq, cN);
for i = 1:cN
    disp(['Series =', num2str(i)]);
    mdd = add(:,:,i); mdu = adu(:,:,i);
    mqd  = quantile(mdd,valpha,1);  mqu  = quantile(mdu,valpha,1);
    mCDaR = NaN(ctau, cq);  mCUaR = NaN(ctau, cq);
    for tau = 1:ctau
        vdd = mdd(ctau+1:end, tau);vdu = mdu(ctau+1:end, tau);
        for q = 1:cq
            mCDaR(tau,q) = nanmean(vdd(vdd>=mqd(q,tau)));
            mCUaR(tau,q) = nanmean(vdu(vdu>=mqu(q,tau)));
        end
    end
    aCDaR(:,:,i) = mCDaR; aCUaR(:,:,i) = mCUaR;
end

mUd = NaN(ctau, cN); mVd = NaN(cq,cN); vlambda_d = NaN(1, cN);
mUu = NaN(ctau, cN); mVu = NaN(cq,cN); vlambda_u = NaN(1, cN);
vGoF = nan(1,cN);
for i = 1:cN
    disp(['Series =', num2str(i)]);
    mCDaR = aCDaR(:,:,i); mCUaR = aCUaR(:,:,i);
    % 
    vs = svd(mCDaR);
    vGoF(i) = vs(1)/sum(vs); vGoF2(i) = vs(1)^2/sum(vs.^2);
    [vu, dlambda, vv] = svds((mCDaR),1);
    mUd(:,i) = vu; mVd(:,i)=vv; vlambda_d(i)=dlambda;
    %
    [vu, dlambda, vv] = svds((mCUaR),1);
    mUu(:,i) = vu; mVu(:,i)=vv; vlambda_u(i)=dlambda;
end
%% CDaR comparison - Contour plot - Figure 5
itic   = 269; % Stock A
itic2  = 363; % Stock B 
gsd = figure("Name","CDaR comparison")
subplot(1,2,1)
    plot(vdates, log(mPrice(:,itic)), "-b", LineWidth=0.8); hold on;
      plot(vdates, log(mPrice(:,itic2)), "-.r", LineWidth=0.8); hold off;
    title("Logarithmic prices",Interpreter="latex")
    legend("Stock A","Stock B", Location="southeast",Orientation="vertical")
    legend("boxoff")
subplot(1,2,2)
    contour(1:ctau, valpha,  aCDaR(:,:,itic)',  "-b",   "LevelList",[0.05, 0.10, 0.15, 0.20, 0.25], ...
        ShowText="on", LineWidth=1); hold on 
    contour(1:ctau, valpha, aCDaR(:,:,itic2)', "-.r", "LevelList",[0.05, 0.10, 0.15, 0.20, 0.25],...
        ShowText="on", LineWidth=1); hold off
    title("CDaR($\tau,\alpha)$ contours", Interpreter="latex");
    xlabel("$\tau$",Interpreter="latex");
    ylabel("$\alpha$", Interpreter="latex")
    %legend("Stock A","Stock B", Location="southwest",Orientation="vertical")
    %legend("boxoff")
orient(gsd,'landscape')
gsd.Position = [-80 -100 1310 850];  
print(gsd,'gSP500_CDaRcomparison','-dpdf',   '-r600') 

%% CDaR Singular value decomposition 
line_color = [vlambda_d(vsel)'/max(vlambda_d(vsel)) .6*ones(sum(vsel),2)];
[vdens,vl] = ksdensity(vlambda_d);
gsd = figure("Name","CDaR comparison")
subplot(1,3,1)
    h = histogram(vlambda_d(vsel), 30, 'Normalization','pdf',   'FaceAlpha',.6); hold on
    h.FaceColor = [1 0.6 0.6];
    h.EdgeColor = [.1 0.6 0.6];
    plot(vl,vdens, Color=[.2 0.6 0.6], LineWidth=2);
    title("Risk factor $\lambda_i$", Interpreter="latex")
    xlabel("$\lambda$", Interpreter="latex")
    ylabel("Density", Interpreter="latex") 
    hold off
    c = colormap(line_color)
subplot(1,3,2)
    plot((1:ctau), mUd(:,vsel)  )
    xlabel("$\tau$",Interpreter="latex");
    ylabel("$\theta_i(\tau)$", Interpreter="latex") 
    title("$\tau$-profile: $\theta_i(\tau)$", Interpreter="latex")
    xlim([0 66])
    colororder(line_color);
subplot(1,3,3)
    b = plot(valpha, mVd(:,vsel) )
    xlabel("$\alpha$", Interpreter="latex")
     ylabel("$\nu_i(\tau)$", Interpreter="latex") 
    title("Level $\alpha$-profile: $\nu_i(\alpha)$", Interpreter="latex")
    colororder(line_color);
orient(gsd,'landscape')
gsd.Position = [-80 -100 1210 750];  
print(gsd,'gSP500_CDaRdecomposition','-dpdf',   '-r600')  