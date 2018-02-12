function Q = max_entropy(n, p, mask, Rest)
% Maximum entropy solution for covariance matrix polynomial...
% Problem from Avventi et al 2013

% BD 15.10.2015

% data sizes
% n = 3;      % no. of variables
% p = 4;      % order of the (true) AR model
% N = 10000;   % sample size
% 
% % generate target matrix
% [Ytrue,Yestim,Rest] = est_IShat_R(n, p, N);  % Rest are covariance matrix estimates
% mask = zeros(n);        % sparsity mask
% mask(2,1) = 1;          % indicator for the zero positions
% mask(1,2) = 1;


% build block Toeplitz matrix
d = n*(p+1);            % size of all interesting matrices
T = zeros(d);
for i = 0 : p           % lazy solution...
  for j = i : p
    T(i*n+1:(i+1)*n, j*n+1:(j+1)*n) = Rest{j-i+1};
    T(j*n+1:(j+1)*n, i*n+1:(i+1)*n) = Rest{j-i+1}';
  end
end

% formulate CVX problem
cvx_begin
  cvx_expert('true');
  cvx_quiet('true');
  variable X(d,d);              % the Gram matrix
  minimize( trace(T * X) - log_det(X(1:n,1:n)) )
  subject to
    X == semidefinite(d);
    
    % zero mask constraints
    for i = 1 : n
      for j = 1 : n
        if mask(i,j)    % if the coefficient must be zero
          mm = zeros(n);
          mm(i,j) = 1;
          for k = 0 : p
            v = zeros(p+1,1);
            v(k+1) = 1;      % diagonal
            T = triu(toeplitz(v));
            trace(kron(T,mm)'*X) == 0;
          end
        end
      end
    end
    
cvx_end

% extract estimated inverse spectrum
for i = 0 : p
  Q{i+1} = block_trace(X,n,i);
end
