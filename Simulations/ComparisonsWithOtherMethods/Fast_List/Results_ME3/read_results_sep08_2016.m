function read_results_sep08_2016


NoCrit = 4;

vec_N = 1000;
Nmax   = max(vec_N);

%model orders; in estimation, we use the true order
vec_p  = [1 5];
%no. of variables:
K      = 10;
kbar = K*(K-1)/2;
%no. of trials:
Ntr    = 10;
%set of true masks:
load('..\Data\mask_set.mat');
n      = size(mask_set, 2)/K;

ii = 0;
for ind_mask=[4,n],
% for ind_mask=4, 
    ii = ii+1;
    Mat{ii} = zeros(2,5);
    for ind_p = 1:length(vec_p),
        p = vec_p(ind_p);
        for ind_tr = 1:Ntr,
            for ind_N = 1:length(vec_N), %one single value for N
                N = vec_N(ind_N);
                fres = strcat('Res_mask',num2str(ind_mask),'_','p',num2str(p),'_','N',num2str(N),'_','tr',num2str(ind_tr));
                load(fres);
                for cc=1:4,
                    [~,w] = min(vec_crit_me(cc,:));
                    Mat{ii}(ind_p,cc) = Mat{ii}(ind_p,cc)+sum(sum(xor(mask_true,vec_mask_me{w})))/2;
                end %cc
                oracle = Inf;
                for jj=1:size(vec_mask_me,2);
                    mask = vec_mask_me{jj};
                    temp = sum(sum(xor(mask,mask_true)))/2;
                    oracle = min(oracle,temp);
                end
                Mat{ii}(ind_p,5) = oracle;
            end %N
        end %tr
    end %p
end %mask

ii = 0;
for ind_mask=[4,n],
% for ind_mask=4,
    mask_true = mask_set(:,(ind_mask-1)*K+1:ind_mask*K);
    no_ones = (sum(sum(mask_true))-K)/2;
    noz = kbar-no_ones;
    ii = ii+1;
    Mat{ii} = ones(size(Mat{ii})) - (Mat{ii}/kbar/Ntr);
    figure(ii);
    bar(vec_p, Mat{ii});
    xlabel('Model order (p^\circ)');
    ylabel('Similarity Index SIM');
    tit = sprintf('SP for which #0=%i', noz);
    title(tit);
    legend('SBC','FPE','RNML','AICC','Oracle','Location','NorthWest');
end

        
end %function

