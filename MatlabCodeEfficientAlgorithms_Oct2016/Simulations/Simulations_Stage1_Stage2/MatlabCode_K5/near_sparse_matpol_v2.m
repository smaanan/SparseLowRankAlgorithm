function Yout = near_sparse_matpol_v2(n, p, Yestim, mask)
% Find the nearest positive matrix polynomial with sparsity pattern
% Second version, using simple positivity constraints

% BD 6.08.2015

% data sizes
% n = 3;      % no. of variables
% p = 4;      % order of the (true) AR model
% N = 1000;   % sample size
% 
% % generate target matrix
% [Ytrue,Yestim] = est_IShat(n, p, N);  % use Yestim as target
% mask = zeros(n);        % sparsity mask
% mask(2,1) = 1;          % indicator for the zero positions
% mask(1,2) = 1;
% 
% % convert target matrix (Yestim is its causal part) from cell to vector form
% % (as required by Pos3Poly)
ncoef_sim = n*(n+1)/2 + p*n*n;  % number of coefficients for the symmetric part
ncoef = (2*p+1)*n*n;            % total number of coefficients

Yevecs = zeros(ncoef_sim,1);    % vectorized symmetric part (not needed!)
i = n*(n+1)/2;
Yevecs(1:i) = vecs(Yestim{1});
for k = 1:p
  Yevecs(i+1:i+n*n) = vec(Yestim{k+1});
  i = i + n*n;
end

% formulate Pos3Poly problem
pars.fid = 0;                   % don't show any message
Xp = [p n];                     % descriptor for the symmetric matrix polynomial
Xptype.trigonometric = 1;
imask = find(vec(mask));        % zero positions in vectorized coefficient
ismask = find(vecs(mask));
cvx_begin
  cvx_quiet('true');
  variable X(ncoef_sim);        % desired matrix polynomial, symmetric by default
  variable err_norm;            % infinity norm of the error
  minimize err_norm             % to be minimized
  subject to
    % positivity constraint
    X == sos_pol(Xp, Xptype, pars);   % X is positive matrix polynomial
    
    % zero mask constraints
    X(ismask) == 0;             % for the symmetric coefficient
    i = n*(n+1)/2;
    for k = 1:p                 % for the other coefficients
      X(i+imask) == 0;
      i = i + n*n;
    end
    
    % error minimization constraints
    -Yevecs + X + err_norm*unitpol(p,n) == sos_pol(Xp, Xptype, pars);
    Yevecs - X + err_norm*unitpol(p,n) == sos_pol(Xp, Xptype, pars);
cvx_end

% convert result back to cell form
Yout{1} = zeros(n);
j = 0;
len = n;
for i = 1:n
  Yout{1}(i:n,i) = X(j+1:j+len);
  Yout{1}(i,i:n) = X(j+1:j+len)';
  j = j+len;
  len = len-1;
end
for k = 1:p
  Yout{k+1} = mat(X(j+1:j+n*n));
  j = j + n*n;
end

% the error is Yout - Yestim
for k = 1 : p+1
  E{k} = Yout{k} - Yestim{k};
end

% its max norm on the unit circle is
% max_error = err_norm
