%1) Compute glomerulus maps from saved .mat files:
%-------------------------------------------------
out_path  = 'D:\zhilei_data\human_nonhuman_reg\';
meta_data = load_meta_data('humanNonhumanOrco_reg_meta_data.mat');

num_clusters = 40; 
selected_odors = '[A-Z][0-9]_m1'; %perform clustering only on the first response of the reference odors, i.e. U1_m1, T1_m1, ...
F_cc = compute_glomerulus_maps2(out_path,num_clusters, meta_data, selected_odors);

%Plot the 3D glomerulus maps:
for(i=1:length(F_cc))
    figure, volume_plot(F_cc{i}.X_vis, F_cc{1}.parameters);
end

%Compute time series and mean odor responses:
F_time_series = apply_odor_sorting(F_cc, meta_data);
regime1 = 1:28;                         
regime1_intervals.reference = [10:20];  
regime1_intervals.signal    =  [30:45];    
regime2_intervals.reference = [10:100];
regime2_intervals.signal    = [190:225];


F_mean_responses = compute_mean_responses(F_time_series, regime1, regime1_intervals, regime2_intervals);

%Plot heatmaps:
for(i=1:length(F_mean_responses))
    show_heatmap(F_time_series{i}.odor_names, F_mean_responses{i});
end


%2) Merge clusters:
%------------------

%clusters with an overlap >= threshold (overlap: between 0 (no overlap) and 1) and a response correlation >= threshold will be merged
%correlation_threshold = 0.7; overlap_threshold=0.05;
correlation_threshold = 0.8; overlap_threshold=0.1;
for(i=1:length(F_cc))
    [F_time_series{i}, F_mean_responses{i}, F_cc{i}] = merge_redundant_clusters2(F_time_series{i}, F_mean_responses{i}, F_cc{i}, correlation_threshold, overlap_threshold);
end


for(i=1:length(F_mean_responses))
    show_heatmap(F_time_series{i}.odor_names, F_mean_responses{i});
end



%3) Glomerulus matching:
%-----------------------

%3a) variant 1: match all animals to one target:
target_index = 4;
constraints.functional_threshold = "off";
constraints.spatial_threshold = "off";

%3b) variant 2: find the best target to match as the one with the minimum overall matching cost:
constraints.functional_threshold = "off";
constraints.spatial_threshold = "off";
cost=zeros(5,1);
for(i=1:5)
    [mapping_all, mapped_indices_all, this_cost] = register_all_to_target(normalize_mean_responses2(extract_odors(F_mean_responses,1:28)), F_cc, constraints, i);
    cost(i)=sum(this_cost);
end
cost
[minimum, min_pos] = min(cost)
target_index = min_pos;

%3c) the actual matching:
responses_train = normalize_mean_responses2(extract_odors(F_mean_responses,1:28)); %match using only the reference odors
responses_all = normalize_mean_responses2(F_mean_responses); %then apply to all odors (could use different normalization here as well)
[mapping_all, mapped_indices_all, cost] = register_all_to_target(responses_train, F_cc, constraints, target_index);
matched = apply_glomerulus_matching_to_all(F_time_series, responses_all, F_cc, mapping_all, mapped_indices_all, target_index);

%normalization variant: we can either normalize glomeruli or odors
%responses_train = normalize_mean_responses(extract_odors(F_mean_responses,1:28));
%responses_all = normalize_mean_responses2(F_mean_responses);


%4) Visualize matching results:
%------------------------------

[responses_2D, responses_3D] = combine_responses(matched);
%Heatmap for all animals concatenated (with NAs for glomeruli that are not matched)
show_heatmap(responses_2D.odor_names, responses_2D.responses)
%Heatmap for the average odor responses over all animals (NAs ignored):
average_heatmap = mean(responses_3D.responses,3,'omitnan'); 
show_heatmap(responses_3D.odor_names, average_heatmap);
%PCA plots:
pca_plot(average_heatmap, responses_3D.odor_names);
pca_plot(responses_2D.responses, responses_2D.odor_names);

%AL maps with matched colours:
for(i=1:length(matched))
   figure, volume_plot(matched{i}.F_cc.X_vis, matched{i}.F_cc.parameters); 
   hold on; volume_plot(matched{i}.F_cc.X_vis_unmatched, matched{i}.F_cc.parameters,{-37.5, 30, 'gray',1});
end 


