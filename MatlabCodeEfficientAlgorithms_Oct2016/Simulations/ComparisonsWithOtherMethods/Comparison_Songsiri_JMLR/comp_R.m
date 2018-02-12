function Rest = comp_R(data,p)

% Estimate covariance matrices from the data

N = size(data,1);
for k = 1 : p+1,
  Rest{k} = data(k:end,:)' * data(1:end-k+1,:) / N;
end


end

