function vec_lam = get_vec_lam( lamin )

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Input:
% lamin     = value of lambda calculated on theoretical grounds
%Output:
% vec_lam   = grid for lambda
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if lamin>=1,
    fprintf('Init val. for lamda is larger than one!\n');
    vec_lam = lamin;
else
    q = 0;
    lam = lamin;
    while lam<1,
        lam = 10*lam;
        q = q+1;
    end
    vec_lam = (10^(-q):10^(-q):10^(-q+1));
end

end %function

