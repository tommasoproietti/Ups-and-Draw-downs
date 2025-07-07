function [vTheta, ds2] = fSpectralFactorisation(vGamma)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
%disp("Spectral factorisation, based on Sayed and Kailath, 2001. ---");
cr = length(vGamma)-1;
mZ = [1 zeros(1, cr-1)]; mT = zeros(cr,cr); mT(1:end-1, 2:end) =eye(cr-1);
ds2 = 1;  
dconv = 1; 
dg0 = vGamma(1); vG=vGamma(2:end); 
vTheta = [1; zeros(cr-1,1)]; mP = eye(cr); % initial values
for i = 0:1000
	vOld = [vTheta; ds2];
	mP   = mT*mP*mT' + ds2 * vTheta * vTheta';
	ds2  = dg0 - mZ * mP * mZ';
	vTheta = (vG - mT*mP*mZ')/ds2;
	vNew   = [vTheta; ds2];
 	dconv = max(abs((vOld - vNew) ./ vNew));
	if (dconv < 10^-9)
	%    disp(["Convergence after ", num2str(i+1), " iterations. "]);
        break;
    end
   
end
end