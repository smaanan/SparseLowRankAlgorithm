function B = block_trace(A, n, id)

% Block trace of matrix A on block diagonal id (positive), with block size n

% BD 15.10.2005

B = zeros(n);

[m,m] = size(A);  % must be square
p = m / n;        % the number of blocks must be integer :)
v = zeros(p,1);
v(id+1) = 1;      % diagonal
T = triu(toeplitz(v));
for i = 1 : n
  for j = 1 : n
    mask = zeros(n);
    mask(i,j) = 1;
    B(i,j) = trace(kron(T,mask)'*A);
  end
end
