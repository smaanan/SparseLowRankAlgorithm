function run_exp_SelectSP_NS_ME
warning('off','all');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
NoCrit = 4;                                 %no. of ITC used in simulations
%ITC: (1) SBC; (2) FPE; (3) RNML; (4) AICC
vec_N = [600:100:1000 3000 9000 27000];     %sample sizes
Nmax = max(vec_N);
p = 10;                                     %true order VAR
pmin = 1;                                   %min VAR order in selection
pmax = 20;                                  %max VAR order in selection
k = 5;                                      %no. components time series
kbar = k*(k-1)/2;
vec_Ntot = kbar:-1:1;                       %max. no of tested 0's at each step
Ntr = 10;                                   %no. trials
mask_set = generate_mask_set(k, 9);         %true SP's
n = size(mask_set, 2)/k;
Ngrid = 1024;                               %no. points unif. grid for \omega
% flag_max = 1; %use max when searching for the entries of ISDM which should be turned to zero
flag_max = 0; %use norm-2 when searching for the entries of ISDM which should be turned to zero

%No. of 0's tested at each step of the algorithm
vec_N0 = 4*ones(size(vec_Ntot));
vec_N0 = min([vec_N0; vec_Ntot],[],1)
%Greedy:
%vec_N0 = vec_Ntot
%List:
%vec_N0 = ones(size(vec_Ntot));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for ind_mask = 1:n,
    mask_ns = zeros(Ntr, length(vec_N), NoCrit);
    mask_me = zeros(Ntr, length(vec_N), NoCrit);
    phatvec = cell(Ntr, length(vec_N));
    mask = mask_set(:,(ind_mask-1)*k+1:ind_mask*k);
    true_mask = mask;    
    [~, A, Sigma] = gen_IS_model(k, p, mask);
    for ind_tr = 1:Ntr,
        all_data = arsim(zeros(k,1), A, Sigma, Nmax);
        for ind_N = 1:length(vec_N)
            N = vec_N(ind_N);
            data = all_data(1:N,:);
            [vec_phat, ~, Ah, Sh] = arfit_mod(data, pmin, pmax);
            phatvec{ind_tr, ind_N} = vec_phat;
            Rest = comp_R(data, max(vec_phat));
            for crit = 1:NoCrit,
                phat = vec_phat(crit);
                [Yest,~] = comp_IShat(k, phat, Sh{crit}, Ah{crit});
                Yout = Yest;
                cur_mask = ones(k,k);
                noz = 0;
                best_crit = eval_crit(crit, data, k, phat, Yout, noz);
                M = -Inf*ones(k,k);
                for ii=1:k
                    for jj=ii:k
                        M(ii,jj) = Inf;
                    end
                end
                best_mask_ns = estimate_ns_me('near_sparse_matpol_v2', M, k, Yest, Ngrid, cur_mask, phat, data, crit, Rest, best_crit, vec_N0, flag_max);
                mask_ns(ind_tr, ind_N, crit) = sum(sum(xor(true_mask, best_mask_ns)))/2;
                % fprintf('NS: mask=%2i res=%2i\n',ind_mask, mask_ns(ind_tr, ind_N, crit));
                best_mask_me = estimate_ns_me('max_entropy', M, k, Yest, Ngrid, cur_mask, phat, data, crit, Rest, best_crit, vec_N0, flag_max);
                mask_me(ind_tr, ind_N, crit) = sum(sum(xor(true_mask, best_mask_me)))/2;
                % fprintf('ME: mask=%2i res=%2i\n',ind_mask ,mask_me(ind_tr, ind_N, crit));
                % parsave(sprintf('new%d.mat',ind_mask), mask_ns, mask_me, phatvec);
            end %crit
        end %ind_N
    end %ind_tr
end % %ind_mask
end % function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function best_mask = estimate_ns_me(fname, M, k, Yest, Ngrid, cur_mask, phat, data, crit, Rest, best_crit, vec_N0, flag_max)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input:     Obvious from the function run_exp_SelectSP_NS_ME
%Output:    Mask selected by the ctiterion crit
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Yout = Yest;
best_mask = cur_mask;
c_mask = cur_mask;
ind_vec_N0 = 1;
kbar = k*(k-1)/2;
while min(min(M))<Inf
    if( vec_N0(ind_vec_N0)+ind_vec_N0 == kbar+1),
        %get all possible positions
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
        %Compute Yout
        if strcmp(fname,'near_sparse_matpol_v2'),
            Yout = feval(fname, k, phat, Yest, mask_estim);
        elseif strcmp(fname,'max_entropy'),
            Yout = feval(fname, k, phat, mask_estim, Rest);
        else
            best_mask = Inf*ones(k,k);
            fprintf('Wrong function name!\n');
            return;
        end
        %Evaluate criterion
        noz = sum(sum(mask_estim))/2;
        val_crit = eval_crit(crit, data, k, phat, Yout, noz);
        %Comparison
        if val_crit<=crit2, %The condition is satisfied if val_crit=crit2=Inf
            crit2 = val_crit;
            mask2 = cur_mask_copy;
            rz = row_zero;
            cz = col_zero;
        end
    end %for counter
    %New SP
    M(rz, cz) = Inf;
    c_mask = mask2;
    
    %Comparison
    if crit2 < best_crit
        best_crit = crit2;
        best_mask = mask2;
    end
    ind_vec_N0 = ind_vec_N0+1;
end %while
end %function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function parsave(fname, a, b, c)
% save(fname, 'a', 'b', 'c')
% end