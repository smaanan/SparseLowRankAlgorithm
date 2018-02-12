% Regularized maximum entropy solution for covariance matrix polynomial...
% Problem from Songsiri & Vandenberghe 2010

% BD 1.09.2016

% data sizes
n = 3;      % no. of variables
p = 4;      % order of the (true) AR model
N = 10000;   % sample size
lam = 0.1;  % regularization parameter

% generate target matrix
[Ytrue,Yestim,Rest] = est_IShat_R(n, p, N);  % Rest are covariance matrix estimates

% solve problem
Q = max_entropy_reg(n, p, Rest, lam);