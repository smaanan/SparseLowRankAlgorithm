function [B, Sigma] = gen_IS_model(n, p, mask)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%AR model as in in Songsiri PhD thesis, Experiment 2, page 70
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%input:
% n           = no. of variables
% p           = order of the true AR model
% mask        = mask of 0's
%output
% B           = coefficients of AR model
% Sigma       = cov. matrix driven noise
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

flag_stab_check=0;
while flag_stab_check==0,
        % Generate IS.
        Y = gen_IS(n, p, mask);
        % Perform spectral factorization.
        [A, ~] = spectral_fact(n,p,Y);
        % Pre-multiply the coefficients of the model
        B = zeros(n,n*p);
        for ind=1:p,
            B(:,(ind-1)*n+1:ind*n) = -A{1}\A{ind+1};
        end
        Sigma = A{1}\eye(size(A{1}));
        Sigma = Sigma*Sigma';
        % Check the stability.
        comp = [B; eye(n*(p-1)) zeros(n*(p-1),n)];
        flag_stab_check = (max(abs(eig(comp)))<1);
end

end% function gen_IS_model

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Y = gen_IS(n, p, mask)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%input:
% n           = no. of variables
% p           = order of the true AR model
% mask        = mask of 0's
%output
% Y           = causal part of "true" matrix polynomial $S^{-1}$ (IS),
%               generated as in Songsiri PhD thesis, Experiment 2, page 70 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Following settings can be adjusted by the user:
%approximate magnitude of nonzero entries
mu = 0.2;
%control for how close are to mu the nonzero entries
factor = 0.01;
%threshold for the eigenvalues of IS (test IS if it is pos. def.)
th = 10^(-1);

%Init:
j = sqrt(-1);
vec_om = 0:(pi/1000):pi;
    
Y{1} = (mu + factor*randn(n,n)).*mask;
Y{1} = (Y{1}+Y{1}')/2;

for k=2:(p+1),
    Y{k} = (mu + factor*randn(n,n)).*mask;
end
    
vec_ind = zeros(size(vec_om));
while prod(vec_ind)==0,
    Y{1} = Y{1}+eye(n,n);
    vec_ind = zeros(size(vec_om));
    
     for ind=1:length(vec_ind),
        om = vec_om(ind);
        IS = Y{1};
        for k=1:p,
            IS = IS+Y{k+1}*exp(-j*om*k) + Y{k+1}'*exp(j*om*k);
        end
        if min(real(eig(IS)))>th,
            vec_ind(ind)=1;
        end 
     end %for 
end %while prod

end %function gen_IS
