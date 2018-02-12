function val_crit = eval_crit(crit, data, n, phat, Yout, noz)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input:
%crit   = (1) SBC; (2) FPE; (3) RNML; (4) AICC
%data   = measurements
%n      = no. of components of time series
%phat   = estimated VAR order
%Yout   = causal part of estimated matrix polynomial $S^{-1}$ (IShat)
%noz    = number of zeros in the mask (below the main diagonal)
%Output:
%val_crit   = criterion value
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Perform spectral factorization
[A,~] = spectral_fact(n, phat, Yout);

%Pre-multiply the coefficients of the model
B = zeros(n, (phat+1)*n);
for ind = 1:phat+1
    B(:,(ind-1)*n+1:ind*n) = A{1}\A{ind};
end

%Compute log(det(Sigma))
N = size(data, 1);
T = N - phat;
Ut= zeros(T, n*(phat+1));
for j = 1:phat+1
    Ut(:,n*(j-1)+1:j*n) = data(phat+2-j:N-j+1,:);
end
temp  = Ut*B';
Sigma = temp'*temp/T;
logdetSigma = log(det(Sigma));

%Compute R2 (modified R)
Y = data(phat+1:end,:);
R2= trace(Y'*Y - T*Sigma)/T;

%No. of param. [Songsiri et al., book chapter, p. 106]
nop = n*(n+1)/2 - noz + phat*(n^2-2*noz);
nop1= nop/n;
 
%Here we calculate the values of the criteria
%SBC
if crit == 1,
    val_crit = T*logdetSigma + nop*log(T);
%FPE    
elseif crit == 2,
    val_crit = n*log((T+nop1)/(T-nop1))+logdetSigma;    
%RNML    
elseif crit == 3,
    if R2 <= 0
        val_crit = Inf;
    else
        ell     = nop1;
        term    = zeros(4,1);
        term(1) = ((T-ell-n+1)/2)*logdetSigma;
        term(2) = -sum(gammaln((T-ell:-1:T-ell+1-n)/2));
        term(3) = -gammaln(ell*n/2);
        term(4) = ((ell*n)/2)*log(R2);
        val_crit = sum(term);
    end
%AICC    
elseif crit == 4,
val_crit = T*logdetSigma + (2*nop*T)/(T-nop-1);
else
    fprintf('Error: Incorrect ITC!\n');
    val_crit = Inf;
end

end %function


