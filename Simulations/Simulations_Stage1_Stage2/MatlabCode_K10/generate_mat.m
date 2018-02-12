function NewMat = generate_mat()
K = 10;
K_bar = K*(K-1)/2;
N = [0 1 5:5:40];
NewMat = zeros(K, K*K);
n = size(NewMat, 2)/K;
for ind_mask = 1:n
    A = triu(ones(K),1);
    B = ones(1, K_bar);
    rp= randperm(K_bar);
    B(rp(1:N(ind_mask))) = 0;
    A(A==1) = B;
    NewMat(:,(ind_mask-1)*K+1:ind_mask*K) = A + A' + eye(K);
end
end