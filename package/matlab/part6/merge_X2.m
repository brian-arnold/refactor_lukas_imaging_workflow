function [redundant, non_redundant, to_merge] = merge_X2(X_vis, F_mean_responses, correlation_threshold, dice_threshold)

%detect relevant overlap between two clusters
num_clusters = size(X_vis,2);
C = zeros(num_clusters, num_clusters);

for(i=1:size(X_vis,2))
    for(j=1:size(X_vis,2))
        in_clus1 = find(X_vis(:,i)==1);
        in_clus2 = find(X_vis(:,j)==1);
        overlap  = length(intersect(in_clus1, in_clus2)); 
        dice     = (2*overlap) / (length(in_clus1)+length(in_clus2));
        
        corr = corrcoef( F_mean_responses(i,:), F_mean_responses(j,:) ); 
        pair_corr = corr(1,2);
        if(dice>=dice_threshold && pair_corr>=correlation_threshold)   
             C(i,j)=1;
        end
    end
end
%figure,imagesc(C)

%add sets of overlapping clusters to 'to_merge'
to_merge  = cell(1,1);
redundant = [];
non_redundant = [];

already_found = [];
counter = 1;
for(i=1:num_clusters)
    [current, already_found] = get_neighbours(C,i,already_found);
    if(length(current)>0)
        to_merge{counter} = current;
        counter = counter + 1;
        if(length(current)==1)
            non_redundant = [non_redundant, current];
        else
            redundant = [redundant, current];
        end
    end
end