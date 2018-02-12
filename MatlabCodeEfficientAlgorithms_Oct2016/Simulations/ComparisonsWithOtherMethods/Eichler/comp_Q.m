
function YY = comp_Q(A, n, p)

AA{1} = eye(size(A{1}));
for k = 2 : p+1,
    AA{k} = -A{k-1};
end
    
for k = 1 : p+1  
        YY{k} = zeros(n);
        for i = 1 : p+2-k,
            YY{k} = YY{k} + AA{i}' * AA{i+k-1};
        end
end
    
end %function
