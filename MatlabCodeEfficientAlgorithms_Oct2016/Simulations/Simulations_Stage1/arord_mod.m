function [sbc, fpe, rnml, aic, aicc, kic, kicc] = arord_mod(data, R, m, ~, ne, pmin, pmax)

% [sbc, fpe, logdp, np] = arord(R, m, mcor, ne, pmin, pmax)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%MOST OF THE CODE IS FROM ARFIT; NEW CODE IS MARKED AS "My code" 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % n:   number of time steps (per realization)
% % m:   number of variables (dimension of state vectors) 
% [n,m] = size(data); 
% ne = ntr*(n-pmax);         % number of block equations of size m
% mcor = 0;                  %no intercept vector
% 
% % compute QR factorization for model of order pmax
% [R, ~]   = arqr(data, pmax, mcor);

 imax 	  = pmax-pmin+1;        % maximum index of output vectors
  
  % initialize output vectors
  sbc     = zeros(1, imax);     % Schwarz's Bayesian Criterion
  fpe     = zeros(1, imax);     % log of Akaike's Final Prediction Error
  rnml    = zeros(1, imax);     % 
  aic     = zeros(1, imax);     % 
  aicc    = zeros(1, imax);     % 
  kic     = zeros(1, imax);     % 
  kicc    = zeros(1, imax);     % 
  logdp   = zeros(1, imax);     % determinant of (scaled) covariance matrix
  np      = zeros(1, imax);     % number of parameter vectors of length m
  np(imax)= m*pmax;

  % Get lower right triangle R22 of R: 
  %
  %   | R11  R12 |
  % R=|          |
  %   | 0    R22 |
  %
  R22     = R(np(imax)+1 : np(imax)+m, np(imax)+1 : np(imax)+m);
  Delta = R22'*R22; %My code

  % From R22, get inverse of residual cross-product matrix for model
  % of order pmax
  invR22  = inv(R22);
  Mp      = invR22*invR22';
  
  % For order selection, get determinant of residual cross-product matrix
  %       logdp = log det(residual cross-product matrix)
  logdp(imax) = 2.*log(abs(prod(diag(R22))));

  % Compute approximate order selection criteria for models of 
  % order pmin:pmax
  i = imax;
  
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  %END OF CODE FROM ARFIT
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
  Y = data(pmax+1:end,:);   %My code
  trYY = trace(Y'*Y);       %My code
  
  for p = pmax:-1:pmin,
      
      np(i) = m*p;	% number of parameter vectors of length m (ARFIT)
      
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      %CODE FROM ARFIT 
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
   if p < pmax
      % Downdate determinant of residual cross-product matrix
      % Rp: Part of R to be added to Cholesky factor of covariance matrix
      Rp       = R(np(i)+1:np(i)+m, np(imax)+1:np(imax)+m);
      Delta_prime = Delta + Rp'*Rp; %My code

      % Get Mp, the downdated inverse of the residual cross-product
      % matrix, using the Woodbury formula
      L        = chol(eye(m) + Rp*Mp*Rp')';
      N        = L \ Rp*Mp;
      Mp       = Mp - N'*N;

      % Get downdated logarithm of determinant
      logdp(i) = logdp(i+1) + 2.* log(abs(prod(diag(L))));
   end

   % Schwarz's Bayesian Criterion
   sbc(i) = logdp(i)/m - log(ne) * (ne-np(i))/ne;

   % logarithm of Akaike's Final Prediction Error
   fpe(i) = logdp(i)/m - log(ne*(ne-np(i))/(ne+np(i)));
   
    %Compute RNML  - My code
    T = ne;
    ell = np(i);
    term = zeros(4,1);
    term(1) = ((T-ell-m+1)/2)*( logdp(i) - m*log(T) );
    term(2) = -sum(gammaln( (T-ell:-1:T-ell+1-m)/2 ));
    term(3) = -gammaln(ell*m/2);
    F = (trYY - trace(Delta))/T;
    if p < pmax
        Delta = Delta_prime;
    end
    term(4) = ((ell*m)/2)*log(F); 
    rnml(i) = sum(term);
    
    %Compute AIC
    aic(i) = T*logdp(i) + 2*(p*m^2+m*(m+1)/2);
    
    %Compute AICc
    b = T/(T-(p*m+m+1));
    aicc(i) = T*logdp(i) + 2*b*(p*m^2+m*(m+1)/2);
    
    %Compute KIC
    kic(i) = T*logdp(i) + 3*(p*m^2+m*(m+1)/2);
    
    %Compute KICc
    kicc(i) = T*logdp(i) + T*m*(2*p*m+m+1)/(T-p*m-m-1) + ...
        (T*m)/(T-m*p-(m-1)/2) + (2*m^2*p+m^2-m)/2;
    
    i      = i-1;                % go to next lower order
   
  end %for p

end %function

