function run_exp_SelectSP_NS_ME
% warning('off','all');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
NoCrit = 4;                                 %no. of ITC used in simulations
ITC = {'SBC', 'FPE', 'RNML', 'AICC'};       %criteria applied
pmin = 1;                                   %min VAR order in selection
pmax = 9;                                   %max VAR order in selection
k = 10;                                      %no. components time series
kbar = k*(k-1)/2;
vec_Ntot = kbar:-1:1;                       %max. no of tested 0's at each step
Ngrid = 1024;                               %no. points unif. grid for \omega
% flag_max = 1; %use max when searching for the entries of ISDM which should be turned to zero
flag_max = 0; %use norm-2 when searching for the entries of ISDM which should be turned to zero

TH_EIG_NS = 10^(-6); %min eig value should be larger than TH_EIG_NS (for near sparse)
TH_EIG_ME = []; %min eig value should be larger than TH_EIG_ME (for max. entropy)

results_file = strcat('Results_FinData','.mat');   %file to save the results
%Best SP's selected by NS-method
SP_NS = cell(3,4); %where
%Cols: (1) SBC; (2) FPE; (3) RNML; (4) AICC
%Rows: (1) List; (2) Mixed List-Greedy; (3) GREEDY
%Best SP's selected by ME-method
SP_ME = cell(3,4); %significance of cols and rows is the same as above

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load FinDataConv.mat
% data = newsrm2;
% clear newsrm2

load FinData.mat
data = newsrm;
clear newsrm

%Zero-mean
data = data-kron(mean(data,1),ones(size(data,1),1));

%Test for valid data
N = size(data,1);
if N==0,
    return;
end

%Estimate the order of the model
[vec_phat, ~, Ah, Sh] = arfit_mod(data, pmin, pmax);
save(results_file,'vec_phat'); %save the estimated orders

%Estimate the autocovariances
Rest = comp_R(data, max(vec_phat));

%Estimate the SP
for crit = 1:NoCrit,
    phat = vec_phat(crit);
    [Yest,~] = comp_IShat(k, phat, Sh{crit}, Ah{crit});
    Yout = Yest;
    cur_mask = ones(k,k);
    noz = 0;
    %Evaluate criterion
    temp_crit = eval_crit(crit, data, k, phat, Yout, noz); 
    if min_pol_eigenvalue(Yout) > TH_EIG_NS,
        best_crit_NS = temp_crit;
    else
        best_crit_NS = Inf;
    end
    best_crit_ME = temp_crit;

    %Initialize M-values
    M = -Inf*ones(k,k);
    for ii=1:k
        for jj=ii:k
            M(ii,jj) = Inf;
        end
    end
   for method=1:2:3,
        if method==1,
            vec_N0 = ones(size(vec_Ntot));  %LIST
            fprintf('%s, %s %s \n', ITC{crit}, 'List', 'NS');
            SP_NS{method,crit} = estimate_ns_me('near_sparse_matpol_v2', M, k, Yest, Ngrid, cur_mask, phat, data, crit, Rest, best_crit_NS, vec_N0, flag_max, TH_EIG_NS);            
            fprintf('%s, %s %s \n', ITC{crit}, 'List', 'ME');
            SP_ME{method,crit} = estimate_ns_me('max_entropy', M, k, Yest, Ngrid, cur_mask, phat, data, crit, Rest, best_crit_ME, vec_N0, flag_max, TH_EIG_ME);
        elseif method==2, 
            %This part is not used
            fprintf('%s, %s \n', ITC{crit}, 'Inter');
            vec_N0 = 4*ones(size(vec_Ntot)); 
            vec_N0 = min([vec_N0; vec_Ntot],[],1); %Mixed List-Greedy
            fprintf('%s, %s \n', ITC{crit}, 'Mixed List-Greedy');
            SP_NS{method,crit} = estimate_ns_me('near_sparse_matpol_v2', M, k, Yest, Ngrid, cur_mask, phat, data, crit, Rest, best_crit_NS, vec_N0, flag_max, TH_EIG_NS);            
            SP_ME{method,crit} = estimate_ns_me('max_entropy', M, k, Yest, Ngrid, cur_mask, phat, data, crit, Rest, best_crit_ME, vec_N0, flag_max, TH_EIG_ME);
        elseif method==3,
            vec_N0 = vec_Ntot;  %Greedy
            fprintf('%s, %s %s \n', ITC{crit}, 'Greedy', 'NS');
            SP_NS{method,crit} = estimate_ns_me('near_sparse_matpol_v2', M, k, Yest, Ngrid, cur_mask, phat, data, crit, Rest, best_crit_NS, vec_N0, flag_max, TH_EIG_NS);            
            fprintf('%s, %s %s \n', ITC{crit}, 'Greedy', 'ME');
            SP_ME{method,crit} = estimate_ns_me('max_entropy', M, k, Yest, Ngrid, cur_mask, phat, data, crit, Rest, best_crit_ME, vec_N0, flag_max, TH_EIG_ME);
        else
        end
    end %method
    save(results_file,'SP_NS','-append');
    save(results_file,'SP_ME','-append');
end %crit
      
end % function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function best_mask = estimate_ns_me(fname, M, k, Yest, Ngrid, cur_mask, phat, data, crit, Rest, best_crit, vec_N0, flag_max, TH_EIG)

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
        noz = sum(sum(mask_estim))/2;
        %Compute Yout and crit
        if strcmp(fname,'near_sparse_matpol_v2'),
            Yout = feval(fname, k, phat, Yest, mask_estim);
            if min_pol_eigenvalue(Yout) > TH_EIG,
                val_crit = eval_crit(crit, data, k, phat, Yout, noz);
            else
                val_crit = Inf;
            end
        elseif strcmp(fname,'max_entropy'),
            Yout = feval(fname, k, phat, mask_estim, Rest);
            val_crit = eval_crit(crit, data, k, phat, Yout, noz);
        else
            best_mask = Inf*ones(k,k);
            fprintf('Wrong function name!\n');
            return;
        end
        
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

