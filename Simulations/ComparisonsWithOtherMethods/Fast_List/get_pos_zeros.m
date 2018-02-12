function mat_pos = get_pos_zeros(K, M, N0)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input:
%K          = no. of components for time series
%M          = matrix which contains sparsity pattern
%N0         = no. of zeros 
%Output:
%{(row_zero, col_zero)}  = positions of the next N0 zeros
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mat_pos = zeros(N0,2);
counter = 1;
for ii=1:K,
    for jj=1:ii-1,
        if M(ii,jj)~=Inf,
            mat_pos(counter,:) = [ii jj];
            counter = counter+1; 
        end         
    end
end

end %function




