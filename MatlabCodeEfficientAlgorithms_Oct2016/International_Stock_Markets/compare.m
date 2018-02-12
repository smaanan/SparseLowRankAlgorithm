function compare(true_SP,est_SP )

mark = {'US','UK','JA','AU','CA','SW','FR','GE','HK','IT'};

k = size(true_SP,1);

for i=1:k,
    for j=1:i-1,
        if xor(true_SP(i,j),est_SP(i,j))
            if (true_SP(i,j)==1) && (est_SP(i,j)==0),
                fprintf('%s %s %s %s \n','FN:',mark{i},',',mark{j});
            else
                fprintf('%s %s %s %s \n','FP:',mark{i},',',mark{j});
            end
        end
    end
end


end

