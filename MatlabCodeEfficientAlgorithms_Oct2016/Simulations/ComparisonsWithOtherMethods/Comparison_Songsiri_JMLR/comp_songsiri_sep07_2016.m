function comp_songsiri_sep07_2016

%criteria:
criteria = {'SBC', 'FPE', 'RNML', 'AICC'};
NoCrit = 4;

%sample sizes:
% vec_N  = [600:100:1000 3000 9000 27000];
vec_N = 1000;
Nmax   = max(vec_N);

%model orders; in estimation, we use the true order
% vec_p  = [1 5];
vec_p = 5;
%no. of variables:
K      = 10;
%no. of trials:
Ntr    = 10;
%set of true masks:
load('.\Data\mask_set.mat');
n      = size(mask_set, 2)/K;
%param. used in estimation:
Ngrid  = 1024;
flag_max = 1;
Th = 0.01;


for ind_mask=[4,n],
    for ind_p = 1:length(vec_p),
        p = vec_p(ind_p);
    for ind_tr = 1:Ntr,
        fdata = strcat('.\Data\Data_mask',num2str(ind_mask),'_','p',num2str(p),'_','tr',num2str(ind_tr));
        load(fdata);
        for ind_N = 1:length(vec_N),
            N = vec_N(ind_N);
            fres = strcat('.\Results_Sep13\Res_mask',num2str(ind_mask),'_','p',num2str(p),'_','N',num2str(N),'_','tr',num2str(ind_tr));
            
            %Values lambda
            lamin = sqrt(log(K)/N);
            vec_lam1 = get_vec_lam(lamin);
            vec_lam2 = vec_lam1+max(vec_lam1);
            vec_lam = [vec_lam1 vec_lam2];
   
            data = all_data(1:N,:);
            Rest = comp_R(data, p);
            vec_crit = zeros(4,length(vec_lam)+1);
            %Mask without zeros
            [~, ~, Ah, Sh] = arfit_mod(data, p, p);
            [ISDM,~] = comp_IShat(K, p, Sh{1}, Ah{1});
            vec_mask{1} = ones(K,K);
            mask_old = ones(K,K);
            noz = 0;
            %crit   = (1) SBC; (2) FPE; (3) RNML; (4) AICC;
            for ind_crit=1:4,
                vec_crit(ind_crit,1) = eval_crit(ind_crit, data, K, p, ISDM, noz);
            end
             
            %Loop for lambda
            for ind_lam = 1:length(vec_lam),
                lam = vec_lam(ind_lam);
                %Regularized ME
                ISDM = max_entropy_reg(K, p, Rest, lam);
                % fprintf('1');
                %Get mask
                mask = get_mask(K, ISDM, Ngrid, flag_max, Th);
                mask_estim = xor(mask, ones(K,K));
                %Evaluate criteria
                if ( sum(sum(mask-mask_old)) ~=0 ),
                    vec_mask{ind_lam+1} = mask;
                    mask_old = mask;
                    ISDM = max_entropy(K, p, mask_estim, Rest);
                    fprintf('2');
                    noz = sum(sum(mask_estim))/2;
                    %crit   = (1) SBC; (2) FPE; (3) RNML; (4) AICC;
                    for ind_crit=1:4,
                        vec_crit(ind_crit,ind_lam+1) = eval_crit(ind_crit, data, K, p, ISDM, noz);
                    end
                else
                    for ind_crit=1:4,
                        vec_crit(ind_crit,ind_lam+1) = vec_crit(ind_crit,ind_lam);
                    end
                    vec_mask{ind_lam+1} = vec_mask{ind_lam};
                end %if
                if mask(1:K,1:K)==eye(K,K),
                    for ii=(ind_lam+2):(length(vec_lam)+1),
                        for ind_crit=1:4,
                            vec_crit(ind_crit,ii) = vec_crit(ind_crit,ind_lam+1);
                        end
                        vec_mask{ii} = mask;    
                    end
                    break;
                end %mask 
                % save(fres, 'criteria', 'mask_true', 'vec_crit', 'vec_mask', 'vec_lam');
            end %ind_lam
            save(fres, 'criteria', 'mask_true', 'vec_crit', 'vec_mask', 'vec_lam');
            fprintf('\n');
        end % N
    end %tr
    end %p
end %mask



end %function


