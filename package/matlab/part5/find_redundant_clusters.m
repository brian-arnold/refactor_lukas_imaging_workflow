%find_redundant_clusters(F_mean_responses{1}, 0.001)
%
function [unique_clusters, redundant, non_redundant, representative, to_merge] = find_redundant_clusters(F_mean_responses, threshold)

F_mean_responses_norm = zscore(F_mean_responses,[],2);
rng(42);

%num_glomeruli = 34; %number of glomeruli in the Orco atlas
%clustering = kmeans(F_mean_responses_norm, num_glomeruli,'Replicates',10);

%determine number of clusters automatically (the number for which the improvement in clustering quality is lower than the given thrshold)
score=[]; clusterings=[];
for(i=2:40)
    [idx, C, sumD] = kmeans(F_mean_responses_norm, i,'Replicates',10);
    score = cat(1,score,mean(sumD));
    clusterings = cat(2,clusterings,idx);
end
score=score/sum(score);
clustering = clusterings(:,min(find(score<=threshold)));

%get cluster members:
clusters = unique(clustering);
members  = zeros(length(clusters),1);
for(i=1:length(clusters))
    members(i) = length(find(clustering==clusters(i)));
end
more_than_one = clusters(find(members>1));

redundant = [];
for(i=1:length(more_than_one))
    redundant = cat(1, redundant, find(clustering==more_than_one(i)));
end
non_redundant = setdiff(1:size(F_mean_responses,1), redundant);

representative = [];
to_merge = cell(length(more_than_one),1);
for(i=1:length(more_than_one))
    in_cluster = find(clustering==more_than_one(i));
    to_merge{i} = uint8(in_cluster);
    representative = cat(1, representative, in_cluster(1));
end

unique_clusters = [non_redundant, representative'];

%show_heatmap(F_time_series.odor_names, F_mean_responses_norm(non_redundant,:))
%show_heatmap(F_time_series.odor_names, F_mean_responses_norm(redundant,:))
%show_heatmap(F_time_series.odor_names, F_mean_responses_norm(representative,:))
%show_heatmap(F_time_series.odor_names, F_mean_responses_norm(unique_clusters,:))


    