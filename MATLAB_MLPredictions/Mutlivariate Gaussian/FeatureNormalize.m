function matrix = FeatureNormalize(data_matrix)

%matrix = zeros(size(data_matrix));
count = 0;
for i = 1:size(data_matrix,2);
    
    feature = data_matrix(:,i); 
    % Check NaN
    if sum(isnan(feature))/numel(feature) <= 0;    % threshold
        count = count + 1;
        new_feature_array = (feature - min(feature))/(max(feature)-min(feature));
        matrix(:,count) = new_feature_array;
end

end
    

