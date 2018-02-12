function mask = get_mask(K, Yestim, Ngrid, flag_max, Th)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input:
%K          = no. of components for time series (manifest+latent)
%Yestim     = causal part of estimated matrix polynomial $S^{-1}$ (IShat)
%Ngrid      = no. of points for the uniform grid on [0,\pi]
%flag_max = 1 : use max when searching for the entries of ISDM which should be turned to zero
%flag_max = 0 :use norm-2 when searching for the entries of ISDM which should be turned to zero
%Th         = threshold for deciding which entries are zero
%Output:
%mask       = sparsity pattern
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

vec_om = 0:(pi/Ngrid):pi;
i = sqrt(-1);

%Compute partial coherence for \omega on the grid
Sinv0 = Yestim{1};
Sinv = Sinv0;
Pcoh_tot = zeros(K,K,length(vec_om));
for ind=1:length(vec_om),
    om = vec_om(ind);
    for ell=2:size(Yestim,2),
        Sinv = Sinv + exp(-i*om*(ell-1))*Yestim{ell} +exp(i*om*(ell-1))*transpose(Yestim{ell});
    end
    Pcoh_tot(:,:,ind) = abs(diag(1./sqrt(diag(Sinv)))*Sinv*diag(1./sqrt(diag(Sinv))));
    Sinv = Sinv0;
end

%Compute M-values 
M = zeros(K,K);
if flag_max==1, %max
    for ii=1:K,
        for jj=1:ii-1,
            M(ii,jj)=max(squeeze(Pcoh_tot(ii,jj,:)));
        end
    end
elseif flag_max==0, %norm-2
    for ii=1:K,
        for jj=1:ii-1,
            M(ii,jj)=norm(squeeze(Pcoh_tot(ii,jj,:)));
        end
    end
else
    mask = Inf;
    fprintf('Wrong value for flag_max!\n');
    return;
end

mask = M>Th;
mask = mask + mask' + eye(size(mask));
% mask

end %function




