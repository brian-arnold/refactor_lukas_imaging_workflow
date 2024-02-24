function F_mean_responses_norm = normalize_mean_responses2(F_mean_responses)

num_animals = length(F_mean_responses);

F_mean_responses_norm = cell(num_animals,1);
for(i=1:num_animals)
       F_mean_responses_norm{i} = zscore(F_mean_responses{i},[],1); 
       %F_mean_responses_norm{i} = min_max_normalize(F_mean_responses{i});
end