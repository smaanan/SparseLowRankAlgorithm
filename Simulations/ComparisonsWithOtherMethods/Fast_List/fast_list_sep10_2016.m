function fast_list_sep10_2016

%criteria:
criteria = {'SBC', 'FPE', 'RNML', 'AICC'};
NoCrit = 4;

%sample sizes:
% vec_N  = [600:100:1000 3000 9000 27000];
vec_N = 1000;
Nmax   = max(vec_N);

%model orders; in estimation, we use the true order
vec_p  = [1 5];
%no. of variables:
K      = 10;
Kbar   = K*(K-1)/2;
% Fast List
vec_N0 = 3*ones(1,Kbar);
%no. of trials:
Ntr    = 10;
%set of true masks:
load('.\Data\mask_set.mat');
n      = size(mask_set, 2)/K;
%param. used in estimation:
Ngrid  = 1024;
flag_max = 0;

for ind_mask=[4,n],
    for ind_p = 1:length(vec_p),
        p = vec_p(ind_p);
    for ind_tr = 1:Ntr,
        fdata = strcat('.\Data\Data_mask',num2str(ind_mask),'_','p',num2str(p),'_','tr',num2str(ind_tr));
        load(fdata);
        for ind_N = 1:length(vec_N),
            N = vec_N(ind_N);
            
            data = all_data(1:N,:);
            Rest = comp_R(data, p);
            %Mask without zeros
            [~, ~, Ah, Sh] = arfit_mod(data, p, p);
            [Yest,~] = comp_IShat(K, p, Sh{1}, Ah{1});
            cur_mask = ones(K,K);
            noz = 0;
            %crit   = (1) SBC; (2) FPE; (3) RNML; (4) AICC;
            for ind_crit=1:4,
                vec_crit(ind_crit,1) = eval_crit(ind_crit, data, K, p, Yest, noz);
            end
            
            M = -Inf*ones(K,K);
            for ii=1:K
                for jj=ii:K
                    M(ii,jj) = Inf;
                end
            end
            %NS
            [vec_crit_ns, vec_mask_ns] = estimate_ns_me('near_sparse_matpol_v2', M, K, Yest, Ngrid, cur_mask, p, data, vec_crit, Rest, vec_N0, flag_max);
            fres = strcat('.\Results_NS3\Res_mask',num2str(ind_mask),'_','p',num2str(p),'_','N',num2str(N),'_','tr',num2str(ind_tr));
            save(fres, 'criteria', 'mask_true', 'vec_crit_ns', 'vec_mask_ns');    
            %ME
%             [vec_crit_me, vec_mask_me] = estimate_ns_me('max_entropy', M, K, Yest, Ngrid, cur_mask, p, data, vec_crit, Rest, vec_N0, flag_max);
%             fres = strcat('.\Results_ME3\Res_mask',num2str(ind_mask),'_','p',num2str(p),'_','N',num2str(N),'_','tr',num2str(ind_tr));
%             save(fres, 'criteria', 'mask_true', 'vec_crit_me', 'vec_mask_me');               
        end % N
    fprintf('\n');
    end %tr
    end %p
end %mask


end %function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [vec_crit_out, vec_mask_out] = estimate_ns_me(fname, M, k, Yest, Ngrid, cur_mask, p, data, vec_crit_in, Rest, vec_N0, flag_max)

Yout = Yest;
vec_crit_out = vec_crit_in;
vec_mask_out{1} = cur_mask;
ind_vec_N0 = 1;
kbar = k*(k-1)/2;
noz = 0;

while min(min(M))<Inf,
    if( vec_N0(ind_vec_N0)+noz > kbar),
        mat_pos = get_pos_zeros(k, M, (kbar-noz));
    else
        mat_pos = get_pos_N0_zeros(k, Yout, Ngrid, M, vec_N0(ind_vec_N0), flag_max);
    end
    for counter=1:size(mat_pos,1),
        row_zero = mat_pos(counter,1);
        col_zero = mat_pos(counter,2);
        cur_mask(row_zero, col_zero) = 0;
        cur_mask(col_zero, row_zero) = 0;
        M(row_zero, col_zero) = Inf;
    end %counter
    
    vec_mask_out{ind_vec_N0+1} = cur_mask;
    mask_estim = xor(cur_mask, ones(k,k));
    if strcmp(fname,'near_sparse_matpol_v2'),
        Yout = feval(fname, k, p, Yest, mask_estim);
        fprintf('<');
    elseif strcmp(fname,'max_entropy'),
        Yout = feval(fname, k, p, mask_estim, Rest);
        fprintf('>');
    else
         fprintf('Wrong function name!\n');
         return;
    end
    noz = sum(sum(mask_estim))/2;
    for crit=1:4,
        vec_crit_out(crit,ind_vec_N0+1) = eval_crit(crit, data, k, p, Yout, noz);
    end
        
    ind_vec_N0 = ind_vec_N0+1;
end 

end %function

