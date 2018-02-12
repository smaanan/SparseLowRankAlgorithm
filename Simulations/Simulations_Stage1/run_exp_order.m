function run_exp_order

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function produces the results for Example 1, reported in Figure 1.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Initialization
vec_N = [600:100:1000 3000 9000 27000]; %sample sizes
vec_p = [1 5 10 15];                    %orders AR model
pmin  = 1;                              %min. AR order in estimation
pmax  = 20;                             %max. AR order in estimation 
n     = 5;                              %no. of time series
k     = 9;                              %max. no. of 0's in the upper triangular part
Ntr   = 100;                            %no. of realizations for each model

Nmax  = vec_N(end);                     %max. sample size
fname = strcat('Result_Order_Estimation'); %file results 

% Crit1 = SBC
% Crit2 = FPE
% Crit3 = RNML
% Crit4 = AIC
% Crit5 = AICc
% Crit6 = KIC
% Crit7 = KICc
%
ZE = zeros(length(vec_p),k+1);
%
ord_eq = cell(7,length(vec_N));
for crit=1:7,
    for ind_N=1:length(vec_N),
        ord_eq{crit,ind_N}=ZE;
    end
end
%
ord_under = cell(7,length(vec_N));
for crit=1:7,
    for ind_N=1:length(vec_N),
        ord_under{crit,ind_N}=ZE;
    end
end
clear crit;

%Mask_set contains all the masks we use
[mask_set, flag_err] = generate_mask_set(n,k);
if flag_err,
    return;
end

fprintf('\n');

for ind_p=1:length(vec_p),
    %True order
    p = vec_p(ind_p);
    for ind_mask=1:k+1,
        ZEc = ZE;
        ZEc(ind_p,ind_mask)=1;
        %True mask
        no_true = ind_mask;
        mask = mask_set(:,(no_true-1)*n+1:no_true*n);
        [A, Sigma] = gen_IS_model(n, p, mask);
        for ind_tr=1:Ntr,     
            %Generate data
            all_data = arsim(zeros(n,1),A,Sigma,Nmax);
            for ind_N = 1:length(vec_N),
                %Sample size
                N = vec_N(ind_N);
                data = all_data(1:N,:);
                %Estimate order
                vec_phat = arfit_mod(data, pmin, pmax);
                %Results for order estimation
                for crit=1:7,
                   if vec_phat(crit)==p,
                        ord_eq{crit,ind_N} = ord_eq{crit,ind_N} + ZEc;
                    elseif vec_phat(crit)<p,
                        ord_under{crit,ind_N} = ord_under{crit,ind_N} + ZEc;  
                    else
                    end 
                end %for crit 
            end %sample size
        end %trial
    end % mask
 
save(fname,'ord_eq','ord_under','vec_N','Ntr','vec_p');
fprintf('*');
end %order

fprintf('\n');
end %function