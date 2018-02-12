function mat_pos = get_pos_N0_zeros(K, Yestim, Ngrid, M, N0, flag_max)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input:
%K          = no. of components for time series
%Yestim     = causal part of estimated matrix polynomial $S^{-1}$ (IShat)
%Ngrid      = no. of points for the uniform grid on [0,\pi]
%M          = matrix which contains sparsity pattern
%N0         = no. of zeros 
%flag_max = 1 : use max when searching for the entries of ISDM which should be turned to zero
%flag_max = 0 :use norm-2 when searching for the entries of ISDM which should be turned to zero
%Output:
%{(row_zero, col_zero)}  = positions of the next N0 zeros
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
Mcopy = M;
if flag_max==1, %max
    for ii=1:K,
        for jj=1:ii-1,
            Mcopy(ii,jj)=max( Mcopy(ii,jj), max(squeeze(Pcoh_tot(ii,jj,:))) );
        end
    end
elseif flag_max==0, %norm-2
    for ii=1:K,
        for jj=1:ii-1,
            Mcopy(ii,jj)=max( Mcopy(ii,jj), norm(squeeze(Pcoh_tot(ii,jj,:))) );
        end
    end
else
    mat_pos = Inf;
    fprintf('Wrong value for flag_max!\n');
    return;
end

%Get the positions
mat_pos = zeros(N0,2);
for counter=1:N0,
    [mc,imin] = min(Mcopy);
    [~,jmin] = min(mc);
    row_zero = imin(jmin);
    col_zero = jmin;
    mat_pos(counter,:) = [row_zero col_zero];
    Mcopy(row_zero,col_zero) = Inf;
end %for

end %function




