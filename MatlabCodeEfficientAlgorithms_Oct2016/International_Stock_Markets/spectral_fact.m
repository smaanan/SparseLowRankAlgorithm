function [A, err_spectral_fact] = spectral_fact(n,p,Y)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%input:
% n                 = no. of variables
% p                 = order of the true AR model
% Y                 = causal part of "true" matrix polynomial $S^{-1}$ (IS)
%output
% A                 = spectral factor
% err_spectral_fact = error spectral factorization             
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

R = Y;

s = n*(p+1);
cvx_begin quiet
    variable Q(s,s) symmetric;
    maximize trace(Q(1:n,1:n));
    subject to
    Q == semidefinite(s);
    for k = 0 : p
        ind_diag = diag(ones(p+1-k,1),k);
        for i = 1 : n
            for j = 1 : n
                ind_bl = zeros(n);
                ind_bl(i,j) = 1;
                trace(kron(ind_diag,ind_bl)*Q) == R{k+1}(i,j);
            end
        end
    end
    cvx_end
    %-------------------------------------%
    %             'A' FACTORS             %
    %-------------------------------------%
    opts.tol = 1e-2;
    kk = n;
    [V,D] = eigs(full(Q), kk, 'la', opts);
    [~,i] = sort(diag(D),'descend');
    i = i(1:n);
    H = sqrt(D(i,i)) * V(:,i)';
    for k = 1 : p + 1
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
err_spectral_fact = sqrt(e);  % should be zero
    
end %function

