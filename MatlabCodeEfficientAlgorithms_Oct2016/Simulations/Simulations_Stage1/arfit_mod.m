function  vec_phat = arfit_mod(data, pmin, pmax)
% function  [vec_phat, A_SBC, Sigma_SBC, A_FPE, Sigma_FPE, A_RNML, Sigma_RNML]=arfit_mod(data, pmin, pmax)
%ARFIT	Stepwise least squares estimation of multivariate AR model.

% n:   number of time steps (per realization)
% m:   number of variables (dimension of state vectors) 

[n,m] = size(data); 
mcor = 0;                       %no intercept vector
ne  	= n-pmax;               % number of block equations of size m
% npmax	= m*pmax+mcor;          % maximum number of parameter vectors of length m
crit    = {'sbc','fpe','rnml','aic','aicc','kic','kicc'};
vec_phat = zeros(1,7);

% compute QR factorization for model of order pmax
[R, scale]   = arqr(data, pmax, mcor);

% compute approximate order selection criteria for models 
% of order pmin:pmax
[sbc, fpe, rnml, aic, aicc, kic, kicc]   = arord_mod(data, R, m, mcor, ne, pmin, pmax);

for ind_crit=1:7,
    [~, iopt]  = min(eval(crit{1,ind_crit})); 
    vec_phat(ind_crit) = pmin + iopt-1; % estimated order 
end
  
end %function 

  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   function [w, A, C] = comp_A(popt, R, scale, m, mcor, npmax, ne)
%   
%   np           = m*popt + mcor; % number of parameter vectors of length m
% 
%   % decompose R for the optimal model order popt according to 
%   %
%   %   | R11  R12 |
%   % R=|          |
%   %   | 0    R22 |
%   %
%   R11   = R(1:np, 1:np);
%   R12   = R(1:np, npmax+1:npmax+m);    
%   R22   = R(np+1:npmax+m, npmax+1:npmax+m);
% 
%   % get augmented parameter matrix Aaug=[w A] if mcor=1 and Aaug=A if mcor=0
%   if (np > 0)   
%     if (mcor == 1)
%       % improve condition of R11 by re-scaling first column
%       con 	= max(scale(2:npmax+m)) / scale(1); 
%       R11(:,1)	= R11(:,1)*con; 
%     end;
%     Aaug = (R11\R12)';
%     
%     %  return coefficient matrix A and intercept vector w separately
%     if (mcor == 1)
%       % intercept vector w is first column of Aaug, rest of Aaug is 
%       % coefficient matrix A
%       w = Aaug(:,1)*con;        % undo condition-improving scaling
%       A = Aaug(:,2:np);
%     else
%       % return an intercept vector of zeros 
%       w = zeros(m,1);
%       A = Aaug;
%     end
%   else
%     % no parameters have been estimated 
%     % => return only covariance matrix estimate and order selection 
%     % criteria for ``zeroth order model''  
%     w   = zeros(m,1);
%     A   = [];
%   end
%   
%    % return covariance matrix
%    dof   = ne-np;                % number of block degrees of freedom
%    C     = R22'*R22./dof;        % bias-corrected estimate of covariance matrix
%  
%     end %function comp_coeff  
% 
% 
% 
%   % for later computation of confidence intervals return in th: 
%   % (i)  the inverse of U=R11'*R11, which appears in the asymptotic 
%   %      covariance matrix of the least squares estimator
%   % (ii) the number of degrees of freedom of the residual covariance matrix 
% %   invR11 = inv(R11);
% %   if (mcor == 1)
% %     % undo condition improving scaling
% %     invR11(1, :) = invR11(1, :) * con;
% %   end
% %   Uinv   = invR11*invR11';
% %   th     = [dof zeros(1,size(Uinv,2)-1); Uinv];
% 
% 
% 
