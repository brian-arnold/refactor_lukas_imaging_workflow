function output = save_cluster(X, w, h, num_slices)

output = zeros(w, h, num_slices);
num_clusters = size(X, 2);
% remove very large cluster, which is likely to be background noise
size_threshold = 10000;
% count the number of clusters that pass the filter
cluster_count = 1;
for idx=1:num_clusters
    temp = reshape(X(:,idx),w,h,num_slices);
    loc_idx = find(temp==1);
%     if(size(loc_idx,1)>size_threshold)
%         continue
%     end
    [x,y,z] = ind2sub(size(temp), loc_idx);
    % remove clusters whose mode is in the first or last stack
    % likely to be artifacts
%     if(mode(z)==1 || mode(z)==num_slices)
%         continue
%     end
    output(loc_idx) = cluster_count;
    cluster_count = cluster_count + 1;
end
