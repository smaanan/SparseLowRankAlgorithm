function comp_sim_index

[true_mat_weak, ~] = build_true_mat;
k = size(true_mat_weak,1);
kbar = k*(k-1)/2;

load Results_FinData

SIM = zeros(3,4);

for i=1:2:3,
    for j=1:4,
        SIM(i,j) = 1-sum(sum(xor(true_mat_weak,SP_ME{i,j})))/2/kbar;
    end 
end

fprintf('Results for ME:\n');

SIM

SIM = zeros(3,4);

for i=1:2:3,
    for j=1:4,
        SIM(i,j) = 1-sum(sum(xor(true_mat_weak,SP_NS{i,j})))/2/kbar;
    end 
end

fprintf('Results for  NS:\n');

SIM

end %function

