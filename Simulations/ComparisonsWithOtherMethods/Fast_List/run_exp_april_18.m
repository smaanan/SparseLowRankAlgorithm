function run_exp_april_18
warning ('off','all');
NoCrit = 4;
vec_N  = [600:100:1000 3000 9000 27000];
Nmax   = max(vec_N);
p      = 5;
pmin   = 1;
pmax   = 10;
K      = 10;
Kbar   = K*(K-1)/2;
vec_Ntot = Kbar:-1:1;
Ntr    = 10;
mask_set = generate_mat();
n      = size(mask_set, 2)/K;
Ngrid  = 1024;
flag_max = 0;
% List
vec_N0 = ones(size(vec_Ntot));

parfor ind_mask = 1:n
    mask_ns = zeros(Ntr, length(vec_N), NoCrit);
    mask_me = zeros(Ntr, length(vec_N), NoCrit);
    phatvec = cell(Ntr, length(vec_N));
    mask    = mask_set(:,(ind_mask-1)*K+1:ind_mask*K);
    true_mask = mask;
    [~, A, Sigma] = gen_IS_model(K, p, mask);
    for ind_tr = 1:Ntr
        all_data = arsim(zeros(K,1), A, Sigma, Nmax);
        for ind_N = 1:length(vec_N)
            N = vec_N(ind_N);
            data = all_data(1:N,:);
            [vec_phat, ~, Ah, Sh] = arfit_mod(data, pmin, pmax);
            phatvec{ind_tr, ind_N} = vec_phat;
            Rest = comp_R(data, max(vec_phat));
            for crit = 1:NoCrit
                phat = vec_phat(crit);
                [Yest,~] = comp_IShat(K, phat, Sh{crit}, Ah{crit});
                Yout = Yest;
                cur_mask = ones(K,K);
                noz = 0;
                best_crit = eval_crit(crit, data, K, phat, Yout, noz);
                M = -Inf*ones(K,K);
                for ii=1:K
                    for jj=ii:K
                        M(ii,jj) = Inf;
                    end
                end
                best_mask_ns = estimate_ns_me('near_sparse_matpol_v2', M, K, Yest, Ngrid, cur_mask, phat, data, crit, Rest, best_crit, vec_N0, flag_max);
                mask_ns(ind_tr, ind_N, crit) = sum(sum(xor(true_mask, best_mask_ns)))/2;
                best_mask_me = estimate_ns_me('max_entropy', M, K, Yest, Ngrid, cur_mask, phat, data, crit, Rest, best_crit, vec_N0, flag_max);
                mask_me(ind_tr, ind_N, crit) = sum(sum(xor(true_mask, best_mask_me)))/2;
                parsave(sprintf('june16_4_%d.mat',ind_mask), mask_ns, mask_me, phatvec);
            end
        end
    end
end
end

function best_mask = estimate_ns_me(fname, M, k, Yest, Ngrid, cur_mask, phat, data, crit, Rest, best_crit, vec_N0, flag_max)
Yout = Yest;
best_mask = cur_mask;
c_mask = cur_mask;
ind_vec_N0 = 1;
kbar = k*(k-1)/2;
while min(min(M))<Inf
    if( vec_N0(ind_vec_N0)+ind_vec_N0 == kbar+1),
        mat_pos = get_pos_zeros(k, M, vec_N0(ind_vec_N0));
    else
        mat_pos = get_pos_N0_zeros(k, Yout, Ngrid, M, vec_N0(ind_vec_N0), flag_max);
    end
    crit2 = Inf;
    for counter=1:size(mat_pos,1),
        row_zero = mat_pos(counter,1);
        col_zero = mat_pos(counter,2);
        cur_mask_copy = c_mask;
        cur_mask_copy(row_zero, col_zero) = 0;
        cur_mask_copy(col_zero, row_zero) = 0;
        mask_estim = xor(cur_mask_copy, ones(k,k));
        if strcmp(fname,'near_sparse_matpol_v2'),
            Yout = feval(fname, k, phat, Yest, mask_estim);
        elseif strcmp(fname,'max_entropy'),
            Yout = feval(fname, k, phat, mask_estim, Rest);
        else
            best_mask = Inf*ones(k,k);
            fprintf('Wrong function name!\n');
            return;
        end
        noz = sum(sum(mask_estim))/2;
        val_crit = eval_crit(crit, data, k, phat, Yout, noz);
        if val_crit<=crit2, 
            crit2 = val_crit;
            mask2 = cur_mask_copy;
            rz = row_zero;
            cz = col_zero;
        end
    end 
    M(rz, cz) = Inf;
    c_mask = mask2;    
    if crit2 < best_crit
        best_crit = crit2;
        best_mask = mask2;
    end
    ind_vec_N0 = ind_vec_N0+1;
end 
end 

function parsave(fname, a, b, c)
save(fname, 'a', 'b', 'c')
end