%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [A,Sigma] = solve_yw(Rest, n, p)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Solves YW equations
%Input:
% Rest = estimated covariances
% n    = no. of variables
% p    = AR order
%Output:
% A{1},...,A{p}, Sigma
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% build block Toeplitz matrix
d = n*p;            
T = zeros(d);
t = zeros(d,n);
for i = 0 : p-1           % lazy solution...
  t(i*n+1:(i+1)*n,:) = Rest{i+2}';  
  for j = i : p-1
    T(i*n+1:(i+1)*n, j*n+1:(j+1)*n) = Rest{j-i+1};
    T(j*n+1:(j+1)*n, i*n+1:(i+1)*n) = Rest{j-i+1}';
  end
end

MA = T\t;

Sigma = Rest{1};
for i=0:p-1,
    A{i+1} = MA(i*n+1:(i+1)*n,:)';
    Sigma = Sigma - Rest{i+2}*A{i+1}';
end

end %function

