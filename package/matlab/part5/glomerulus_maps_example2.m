%Compute glomerulus maps from saved .mat files:
out_path  = 'D:\zhilei_data\humanNonhumanOrco_new_preprocessing\odorEvoked\'
meta_data = load_meta_data('humanNonhumanOrco_meta_data.mat');
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
regime1_intervals.signal =  [30:45];    
%regime2_intervals.reference = [10:150]; 
regime2_intervals.reference = [10:100]; 
regime2_intervals.signal = [190:225];

F_mean_responses = compute_mean_responses(F_time_series, regime1, regime1_intervals, regime2_intervals);
%Optional: remove an odor (here the odor labeled as 'mistake' in animal 1).
%Sorting and removing surplus odors is relevant for glomerulus matching, where we need the same odor sequence for all animals.
[F_time_series{1}, F_mean_responses{1}] = remove_odor(F_time_series{1}, F_mean_responses{1}, 37); 

for(i=1:length(F_mean_responses))
    show_heatmap(F_time_series{i}.odor_names, F_mean_responses{i});
end


%optional: merge clusters
correlation_threshold = 0.7; %adjacent clusters with a response correlation >= threshold will be merged
for(i=1:length(F_cc))
    [F_time_series{i}, F_mean_responses{i}, F_cc{i}] = merge_redundant_clusters(F_time_series{i}, F_mean_responses{i}, F_cc{i}, correlation_threshold);
end

%merging creates bigger clusters that might now be adjacent --> it can make sense to merge again:
correlation_threshold = 0.9; %more strict threshold criterion to avoid arbitrarly large clusters
for(i=1:length(F_cc))
    [F_time_series{i}, F_mean_responses{i}, F_cc{i}] = merge_redundant_clusters(F_time_series{i}, F_mean_responses{i}, F_cc{i}, correlation_threshold);
end


%Match all animals to the first:
target_index = 1;
constraints.functional_threshold = "off";
constraints.spatial_threshold = "off";

target_index = 1;
constraints.functional_threshold = 0;
constraints.spatial_threshold = 20;

[mapping_all, mapped_indices_all] = register_all_to_target(normalize_mean_responses(F_mean_responses), F_cc, constraints, target_index);
matched = apply_glomerulus_matching_to_all(F_time_series, normalize_mean_responses2(F_mean_responses), F_cc, mapping_all, mapped_indices_all, target_index);

%Visualize results:
[responses_2D, responses_3D] = combine_responses(matched);
%Heatmap for all animals concatenated (with NAs for glomeruli that are not matched)
show_heatmap(responses_2D.odor_names, responses_2D.responses)
%Heatmap for the average odor responses over all animals (NAs ignored):
average_heatmap = mean(responses_3D.responses,3,'omitnan');
show_heatmap(responses_3D.odor_names, average_heatmap);
pca_plot(average_heatmap, responses_3D.odor_names);
pca_plot(responses_2D.responses, responses_2D.odor_names);

for(i=1:length(matched))
   figure, volume_plot(matched{i}.F_cc.X_vis, matched{i}.F_cc.parameters); 
   hold on; volume_plot(matched{i}.F_cc.X_vis_unmatched, matched{i}.F_cc.parameters,{-37.5, 30, 'gray',1});
end


%Mean animal odor and mean human odor mapped to the target AL:
animal_mean = mean(average_heatmap(:,29:33),2); 
human_mean = mean(average_heatmap(:,37:44),2);
show_3D_heatmap(matched{target_index}.F_cc, [animal_mean, human_mean], 1, {'transparent', 60, 30});
show_3D_heatmap(matched{target_index}.F_cc, [animal_mean, human_mean], 2, {'transparent', 60, 30});

matched = apply_glomerulus_matching_to_all(F_time_series, normalize_mean_responses(F_mean_responses), F_cc, mapping_all, mapped_indices_all, target_index);
