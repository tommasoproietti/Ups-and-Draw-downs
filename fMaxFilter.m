function [vMax, vGap, mS, mT, vpi] = fMaxFilter(vy, cq)
%%
cn = length(vy);
vdy = diff(vy);
cn = length(vdy);
mO = ones(cq,cq);
mC0 = tril(ones(cq,cq));
aC = NaN(cq, cq, cq+1);
aC(:,:,1) = mC0;
for j = 2:cq 
    mCi = blkdiag(-triu(ones(j-1, j-1)), tril(ones(cq-j+1,cq-j+1)));
    aC(:,:,j) = mCi;
end
aC(:,:,cq+1) = -triu(ones(cq, cq));
%%
mS = zeros(cq+1,cn-cq+1);
mvGap = zeros(cq+1,cn-cq+1);
for t = 1:(cn-cq+1)
    vdyt = flip(vdy(t:cq+t-1));
    for j =1:cq+1
        mC = aC(:,:,j);
        cInd = (sum(mC * vdyt > 0) == cq);
        mS(j,t) = cInd;
        if (cInd == 1) & (j==1) mvGap(j,t) = 0;
        else mvGap(j,t) =  abs(sum(vdyt(1:j-1))) * cInd;
        end
    end
end
vGap = -sum(mvGap)';
% Transition matrix
mT = mS(:,1:cn-cq)*mS(:,2:cn-cq+1)' ; % transition matrix
mT = diag(1./sum(mT,2)) *    mT;
% ergodic probabilities
vpi = sum(mS,2)/(cn-cq+1);
%bar( 0:cp, vpi)
%figure(); plot(  vGap, '-');
vMax = vy(cq+1:end)-vGap;
% plot([vy; vMax]');
end
 
 