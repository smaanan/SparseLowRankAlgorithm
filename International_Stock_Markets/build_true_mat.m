function [true_mat_weak, true_mat_strong] = build_true_mat


true_mat_weak = [
%US 
zeros(1,10)
%UK
zeros(1,10)
%JA
zeros(1,10)
%AU
1 0 1 zeros(1,7)
%CA
1 zeros(1,9)
%SW
0 1 0 0 zeros(1,6)
%FR
0 1 0 0 0 1 zeros(1,4)
%GE
1 1 1 0 0 1 1 zeros(1,3)
%HK
0 0 1 1 0 1 0 0 zeros(1,2)
%IT
0 1 0 0 0 0 1 1 0 zeros(1,1)
];

true_mat_weak = true_mat_weak + true_mat_weak' + eye(10,10);

true_mat_strong = [
%US 
zeros(1,10)
%UK
zeros(1,10)
%JA
zeros(1,10)
%AU
1 0 0 zeros(1,7)
%CA
1 zeros(1,9)
%SW
0 1 0 0 zeros(1,6)
%FR
0 1 0 0 0 1 zeros(1,4)
%GE
1 0 0 0 0 1 1 zeros(1,3)
%HK
0 0 1 0 0 0 0 0 zeros(1,2)
%IT
0 1 0 0 0 0 1 0 0 zeros(1,1)
];

true_mat_strong = true_mat_strong + true_mat_strong' + eye(10,10);

end

