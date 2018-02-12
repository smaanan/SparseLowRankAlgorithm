function Q = max_entropy_reg(n, p, Rest, lam)
% Regularized maximum entropy solution for covariance matrix polynomial...
% Problem from Songsiri & Vandenberghe 2010

% BD 1.09.2016

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
  variable b(n*(n-1));          % bounds for the infinity norm of the coefficients vector on a given position
  minimize( trace(T * X) - log_det(X(1:n,1:n)) + lam * sum(b) )
  subject to
    X == semidefinite(d);
    
    % enforce bound for the infinity norms
    ii = 1;
    for i = 1 : n
      for j = 1 : n
        if i ~= j   % only on off-diagonal entries
          mm = zeros(n);
          mm(i,j) = 1;
          for k = 0 : p
            v = zeros(p+1,1);
            v(k+1) = 1;      % diagonal number
            T = triu(toeplitz(v));
            trace(kron(T,mm)'*X) <= b(ii);
            trace(kron(T,mm)'*X) >= -b(ii);
          end
          ii = ii+1;
        end
      end
    end

cvx_end

% extract estimated inverse spectrum
for i = 0 : p
  Q{i+1} = block_trace(X,n,i);
end
