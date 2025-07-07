function [mAs, mPs, vBeta, mVarBeta] = fAugmentedStateSmoother(vy, mZ, mGG, mT, mHH, mW, va, mP, mW0)
%--------------------------------------------------------------------------
cT    = length(vy);  % number of time points 
cm    = size(mT,1);
[cm, ck] = size(mW0);   % ck is number of diffuse state elements
%--------------------------------------------------------------------------
% Initialisation of matrices and scalars
av  = NaN([1 1 cT]);     aV = NaN([1 ck cT]);    af  = NaN([1 1 cT]);     aK = NaN([cm 1 cT]);
aaf = NaN([cm 1 cT]);    aAf = NaN([cm ck cT]);  aP  = NaN([cm cm cT]);   vs = 0; mS = 0;    
cn = 0;  % observations counter
% -------------------------------------------------------------------------
mA = -mW0;
for i = 1:cT
 
	aaf(:,:,i) = va; 	aAf(:,:,i) = mA; 	 aP(:,:,i) = mP;	  
    if  (isnan(vy(i)) == 0)  % y(i) is not missing
 		dv = vy(i) - mZ * va;           vV =  - mZ * mA;			
       % disp(dv);
		df= mZ * mP * mZ' + mGG;
		%vK =( mT * mP * mZ' +  mHGt)/ df;
        vK =( mT * mP * mZ' )/ df;
        av(:,:,i) = dv; 		aV(:,:,i) = vV;	  af(:,:,i) = df; 		aK(:,:,i) = vK;
		va = mT * va + vK * dv; 		mA = mW + mT * mA + vK * vV;
		mP = mT * mP * mT' + mHH - vK * vK' * df ;
	 	vs = vs + dv * vV' / df ; 		mS = mS + vV' * vV/df;
 		cn = cn + 1;  % observation counter
	else % y(i) is   missing
		va = mT * va; 		mA = mW + mT * mA; 		mP = mT * mP * mT' + mHH; 
    end 
end
[mQ,mR]  = qr(mS); opts.UT  = true;
mRinv    = linsolve(mR,eye(size(mR)),opts);
mVarBeta = mRinv*mQ';  	 vBeta    = mVarBeta * vs;
%%	
vr = zeros(cm,1);  mR = zeros(cm, ck); mN = zeros(cm,cm); 
mAs = NaN(cm, cT); 	mPs = NaN(cm,cm, cT); 
for i = cT:-1.0:1
	%disp(['Smoothing observation: ', num2str(i)]);
    if (isnan(vy(i)) == 0)
	 	dfinv =  1 / af(:,:,i);
		mL = mT - aK(:,:,i) * mZ;			
		vr = mZ' * dfinv * av(:,:,i) + mL' * vr;	
		mR = mZ' * dfinv * aV(:,:,i) + mL' * mR;
	 	mN = mZ' * dfinv * mZ + mL' * mN * mL;
    else 
        vr = mT' * vr; 		mR = mT' * mR;  		mN = mT' * mN * mT; 
    end
	mAs(:,i) = aaf(:,:,i) - aAf(:,:,i)*vBeta + aP(:,:,i) * (vr - mR * vBeta );
	mAstar = aAf(:,:,i) + aP(:,:,i) * mR; 
    mPs(:,:,i) =  (aP(:,:,i) - aP(:,:,i) * mN * aP(:,:,i) + mAstar * mVarBeta * mAstar' );
end
end