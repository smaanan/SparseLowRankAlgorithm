function v = min_pol_eigenvalue(Xc)
% Compute smallest eigenvalue over [0,pi] of a symmetric matrix polynomial Xc
% given in cell form

% BD 24.08.2015

% data sizes
p = length(Xc)-1;    % degree of the matrix polynomial
n = size(Xc{1},1);   % coefficient size (assume it's square, as it should)

% convert from cell to vector form (as required by Pos3Poly)
ncoef_sim = n*(n+1)/2 + p*n*n;  % number of coefficients for the symmetric part
ncoef = (2*p+1)*n*n;            % total number of coefficients

Xv = zeros(ncoef_sim,1);        % vectorized symmetric part
i = n*(n+1)/2;
Xv(1:i) = vecs(Xc{1});
for k = 1:p
  Xv(i+1:i+n*n) = vec(Xc{k+1});
  i = i + n*n;
end

% formulate Pos3Poly problem
pars.fid = 0;                   % don't show any message
Xp = [p n];                     % descriptor for the symmetric matrix polynomial
Xptype.trigonometric = 1;
cvx_begin
  cvx_quiet('true');
  variable v;                   % minimum eigenvalue
  maximize v                    % to be maximized
  subject to
    % positivity constraint
    Xv - v*unitpol(p,n) == sos_pol(Xp, Xptype, pars);
cvx_end