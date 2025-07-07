function [mZ, mGG, mT, mHH, mW, va, mP, mW0] = fSsfMBfiltering(vphi, vtheta, cs, dcutoff)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MB filtering for ARIMA MODEL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cp = length(vphi); cq= length(vtheta);
dlambda =  power((1+cos(dcutoff))/(1-cos(dcutoff)),cs);
vGamma  = fACGF(cs, dlambda); 
[vVarPhi, ds2] = fSpectralFactorisation(vGamma);
vSum = zeros(cs+1,1); vDif = zeros(cs+1,1); 
for j=0:cs
    vSum(j+1) = nchoosek(cs, j);
    vDif(j+1) = (-1)^j * nchoosek(cs, j);
end
vARpoly = conv([1;(-vphi)], [1; (vVarPhi)]);   
vMApoly_lp = conv([1; vtheta], vSum); vMApoly_hp = conv([1; vtheta], vDif);
vAR = -vARpoly(2:end);
cms = max(cp+cs,cq+cs+1);
if cms > cp+cs
    vARe = [vAR; zeros(cms-cp-cs, 1)];
end
cm  = 2*cms+1;
vMA_lp = zeros(cms-1,1);  vMA_hp = zeros(cms-1,1); 
vMA_lp(1:cq+cs) = vMApoly_lp(2:end); vMA_hp(1:cq+cs) = vMApoly_hp(2:end);
mZarma = [1, zeros(1, cms-1)];
% Measurement equation: % y_{t} = Z_t \alpha_{t} +  G_t \epsilon_t
mZ = [1, 0, zeros(1, cms-1), mZarma];     mGG = 0; % G_t * G_t'
% Transition equation:  % \alpha_{t+1} = T_t \alpha_{t} + H_t \epsilon_t
mTarma =  [vARe,  [eye(cms-1); zeros(1, cms-1)]];
mT =  blkdiag(1,kron(eye(2), mTarma));   
mT(1, 2:cms+1) = mZarma*mTarma;
mHarma = [[1;vMA_lp], zeros(cms,1);
      zeros(cms,1), [1;vMA_hp]] * diag([1/sqrt(ds2), sqrt(dlambda/ds2)]);
mH = [mZarma * mHarma(1:cms,:); mHarma];
mHH = mH*mH'; % H_t * H_t' 
mW  = zeros(cm,2); mW(1,2) = 1;   % drift
%%%%%%%%%%%%%%% Initial conditions               %%%%%%%%%%%%%%%%%%%%%%%%%
va = zeros(cm,1);
mParma = reshape(inv(eye((cm-1)^2)-kron(mT(2:end,2:end),mT(2:end,2:end))) * reshape(mHarma*mHarma', (cm-1)^2,1) , cm-1, cm-1);    %        
mP = blkdiag(0,mParma);
mW0 =  zeros(cm,2); mW0(1,1) = 1;   % initial trend value
end


