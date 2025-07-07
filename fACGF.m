function  vGamma = fACGF(cs, dlambda)
% |\varphi(L)|^2; lambda = variance 
vGamma = zeros(cs+1, 1);
for k=0:cs
    da = nchoosek(2 * cs, cs + k);
    db = (-1)^k * nchoosek(2 * cs, cs + k);
    vGamma(k+1) = da + dlambda * db;
end
end