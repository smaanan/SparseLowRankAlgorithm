%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [mask_set flag_err] = generate_mask_set(n, k)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%input:
% n = no. of variables
% k = max. no. of 0's in the upper triangular part
%output:
% mask_set - Contains k+1 masks
% flag_err = Takes value 1 when an error is found
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
flag_err = 0;
mask_set = ones(n,n*(k+1)); 
no = 1;
i = 1;
j = 1;
mask = ones(n,n);
while no<=k+1,
    % mask
    mask_set(:,(no-1)*n+1:no*n) = mask;
    no = no+1;
    if j<n,
        j = j+1;
    elseif i<n-1
        i=i+1;
        j = i+1;
    else
        flag_err = 1;
        fprintf('k is too large\n');
        return;
    end
    mask(i,j) = 0;
    mask(j,i) = 0;
end

end %function generate_mask_set 

