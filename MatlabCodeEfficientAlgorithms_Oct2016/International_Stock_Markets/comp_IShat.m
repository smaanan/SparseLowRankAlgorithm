function [YY,flag] = comp_IShat(n,phat,Sigmahat,Bhat)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%input:
% n                 = no. of variables
% phat              = estimated order for AR model
% Sigmahat          = estimated covariance matrix noise
% Bhat              = estimated matrix coefficients B
%output
% YY                = causal part of estimated matrix polynomial $S^{-1}$ (IShat)
% err_spectral_fact = estimation error 
% flag              = is 0 in case of an error
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[V,D] = eig(Sigmahat);
if min(diag(D))>0,
    flag = 1;
    temp = zeros(n,n);
    for ind=1:n,
        temp(ind,:) = V(:,ind)'/sqrt(D(ind,ind));
    end
    AA{1} = temp;
    
    for ind=1:phat,
        AA{ind+1} = -AA{1}*Bhat(:,(ind-1)*n+1:ind*n);
    end
    
    for k = 1 : phat+1  
        YY{k} = zeros(n);
        for i = 1 : phat+2-k,
            YY{k} = YY{k} + AA{i}' * AA{i+k-1};
        end
    end   
else
    YY = Inf;
    flag = 0;
end

end %function

