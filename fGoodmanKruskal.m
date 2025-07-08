function dgamma  = fGoodmanKruskal(mN)
[I,J] = size(mN);
% Inizializzo C e D 
cC=0;
cD=0;
for i=1:I-1
    for j=1:J
        % coppie concordanti solo quando j<J
        if j<J
            xsel=mN(i+1:I,j+1:J);
            cC=cC+mN(i,j)*sum(xsel(:));
        end
        % coppie discordanti solo quando j>1
        if j>1
            xsel=mN(i+1:I,1:j-1);
            cD=cD+mN(i,j)*sum(xsel(:));
        end
    end
end

% Indice gamma
dgamma = (cC - cD)/(cC + cD);
%disp(['Indice gamma= ' num2str(dgamma)])
end