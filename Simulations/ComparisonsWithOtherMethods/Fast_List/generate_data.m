function generate_data

vec_N  = [600:100:1000 3000 9000 27000];
Nmax   = max(vec_N);
vec_p  = [1 5]; 
K      = 10;
Ntr    = 10;

temp = dir('.\Data\mask_set.mat');
if numel(temp)==0,
    mask_set = generate_mat();
    save('.\Data\mask_set.mat','mask_set');
else
    load('.\Data\mask_set.mat');
end
clear temp;
n      = size(mask_set, 2)/K;

for ind_mask=1:n,
    mask = mask_set(:,(ind_mask-1)*K+1:ind_mask*K);
    mask_true = mask;
    for ind_p = 1:length(vec_p),
        p = vec_p(ind_p);
        [~, A, Sigma] = gen_IS_model(K, p, mask);
        for ind_tr = 1:Ntr,
            fdata = strcat('.\Data\Data_mask',num2str(ind_mask),'_','p',num2str(p),'_','tr',num2str(ind_tr));
            all_data = arsim(zeros(K,1), A, Sigma, Nmax);
            save(fdata, 'all_data', 'A', 'Sigma', 'mask_true');
        end
    end
end
   
end %function

