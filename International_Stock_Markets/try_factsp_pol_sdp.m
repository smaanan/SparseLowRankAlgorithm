% positive matrix polynomial coefficients, only the causal part
n = 2;  % size of the matrix coefficient
p = 2;  % degree
R{1} = [2 0; 0 2];
R{2} = [0.1 0.2; -0.1 0.1];
R{3} = [-0.2 0.1; 0.1 -0.1];

% solve CVX problem (B.7)
s = n*(p+1);    % size of matrix Q
cvx_begin
  variable Q(s,s) symmetric;
  maximize trace(Q(1:n,1:n));
  subject to
    Q == semidefinite(s);
    % equality constraints
    for k = 0 : p     % for each block diagonal
      ind_diag = diag(ones(p+1-k,1),k);
      for i = 1 : n     % for each element of the block
        for j = 1 : n
          ind_bl = zeros(n);
          ind_bl(i,j) = 1;
          %kron(ind_diag,ind_bl)
          trace( kron(ind_diag,ind_bl) * Q ) == R{k+1}(i,j);
        end
      end
    end
cvx_end

% rank should be n (numerically)
[V,D] = eig(Q);
%diag(D)   % visual check

% build spectral factor as a matrix
[xx,i] = sort(diag(D), 'descend');
i = i(1:n);     % assume rank is n !!!
H = sqrt(D(i,i)) * V(:,i)';

%Q - H'*H  % should be zero

% build coefficients of spectral factor polynomial
for k = 1 : p+1
  A{k} = zeros(n);
  A{k} = H(:,(k-1)*n+1:k*n);
end

% check if indeed spectral factorization (i.e. compute error)
e = 0;
for k = 1 : p+1  % use convolution, to be sure
  RR{k} = zeros(n);
  for i = 1 : p+2-k
    RR{k} = RR{k} + A{i}' * A{i+k-1};
  end
  e = e + norm( RR{k} - R{k} )^2;
end
e = sqrt(e)  % should be zero
