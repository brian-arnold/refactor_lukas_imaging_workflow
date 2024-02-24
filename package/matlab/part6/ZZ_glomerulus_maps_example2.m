% run Martin's glomerulus matching algorithm
% this should be run after the convex cone algorithm on segmentation

%ZZ: experiment specific settings
% experiment_name = 'humanNonhumanOrco';
% brains = {'humanNonhumanOrcoF1', 'humanNonhumanOrcoF2', 'humanNonhumanOrcoF3', 'humanNonhumanOrcoF4', 'humanNonhumanOrcoF5'};
% datasets = {'20191019_humanNonhumanOrco', '20191022_humanNonhumanOrco', '20191023_humanNonhumanOrco', '20191025_humanNonhumanOrco', '20191027_humanNonhumanOrco'};
% meta_file = 'humanNonhumanOrco_meta_data.mat';
experiment_name = 'humanNonhumanOrcoNew2CMTK';
brains = {'humanNonhumanOrcoF1' ,'humanNonhumanOrcoF2' ,'humanNonhumanOrcoF3' ,'humanNonhumanOrcoF4' ,'humanNonhumanOrcoF5'};
datasets = {'20191019_humanNonhumanOrcoNew2' ,'20191022_humanNonhumanOrcoNew2' ,'20191023_humanNonhumanOrcoNew2' ,'20191025_humanNonhumanOrcoNew2' ,'20191027_humanNonhumanOrcoNew2'};
meta_file = 'humanNonhumanOrcoNew2CMTK_meta_data.mat';
% all movies are used for segmentation; some movies may need to be removed
% for matching 
% remove_odor_brain = {1};
% remove_odor_odor = {37};
remove_odor_brain = {}; 
remove_odor_odor = {};

% more universal setting
dim_x = 128; dim_y = 128; dim_z = 24;

% folder settings
% base folder in scotty
folder_home = '/mnt/bucket/labs/mcbride/zhileizhao/Scotty/Analysis/CalciumImaging/AnalysisPipelineOrco/';
% if running script on lab desktop
%folder_home = '/Volumes/mcbride/zhileizhao/Scotty/Analysis/CalciumImaging/AnalysisPipelineOrco/';
folder_working = fullfile(folder_home, 'MartinStrauch', 'glomerulus_matching_scotty');
cd(folder_working);
addpath(genpath(folder_working));
% save results to folders
out_path_svd = fullfile(folder_home, 'MartinStrauch', 'SVDResults', experiment_name);
out_path_segment = fullfile(folder_home, 'MartinStrauch', 'SegmentationResults', experiment_name);
if exist(out_path_segment, 'dir')
    rmdir(out_path_segment, 's');
end
mkdir(out_path_segment);

meta_data_fn = fullfile(folder_home, 'MartinStrauch', 'MetaData', meta_file);
meta_data = load_meta_data(meta_data_fn);

%linux path and windows path name handled differently
out_path_svd_temp = strcat(out_path_svd, '/');

%Set the desired number of clusters (glomeruli):
num_clusters = 40; 
selected_odors = '[A-Z][0-9]_m1'; %perform clustering only on the first response of the reference odors, i.e. U1_m1, T1_m1, ...
F_cc = compute_glomerulus_maps2(out_path_svd_temp,num_clusters, meta_data, selected_odors);

%ZZ: save the convex cone results, so can be used later
fn_cc = fullfile(out_path_segment, strcat(experiment_name,'_convex_cone_res.mat'));
save(fn_cc, 'F_cc', '-v7.3');

%Plot and save the 3D glomerulus maps:
for(i=1:length(F_cc))
    figure, volume_plot(F_cc{i}.X_vis, F_cc{1}.parameters);
    %save figure
    fig_fn = fullfile(out_path_segment, strcat(brains(i),'_segment.fig'));
    savefig(gcf, string(fig_fn));
end

% save segmentation as tif; not flipped, original orientation
for brain_idx=1:size(F_cc,1)
    cc_res = F_cc{brain_idx};
    cc_X = cc_res.X;
    % binarize the result
    threshold = 0; 
    cc_bi = binarize_X(cc_X, threshold);
    cc_tif = save_cluster(cc_bi, dim_x, dim_y, dim_z);
    % save as tif
    options.overwrite = true;
    fname = fullfile(out_path_segment, strcat(brains{brain_idx}, '_segment.tif'));
    saveastiff(single(cc_tif), fname, options);
end

% save segmentation as tif; flipped orientation
for brain_idx=1:size(F_cc,1)
    cc_res = F_cc{brain_idx};
    cc_bi = cc_res.X_vis;
    cc_tif = save_cluster(cc_bi, dim_x, dim_y, dim_z);
    % save as tif
    options.overwrite = true;
    fname = fullfile(out_path_segment, strcat(brains{brain_idx}, '_segmentFlipped.tif'));
    saveastiff(single(cc_tif), fname, options);
end



% folder where segmentation results are stored
out_path_segment = fullfile(folder_home, 'MartinStrauch', 'SegmentationResults', experiment_name);
% folder to store matching results
path_matching = fullfile(folder_home, 'MartinStrauch', 'MatchingResults', experiment_name);
if exist(path_matching, 'dir')
    rmdir(path_matching, 's');
end
mkdir(path_matching);


%Compute time series and mean odor responses:
F_time_series = apply_odor_sorting(F_cc, meta_data);
regime1 = 1:28;                         
regime1_intervals.reference = [10:20];  
regime1_intervals.signal =  [30:45];    
regime2_intervals.reference = [10:150]; 
regime2_intervals.signal = [190:225];
F_mean_responses = compute_mean_responses(F_time_series, regime1, regime1_intervals, regime2_intervals);
%Optional: remove an odor (here the odor labeled as 'mistake' in animal 1).
%Sorting and removing surplus odors is relevant for glomerulus matching, where we need the same odor sequence for all animals.
%[F_time_series{1}, F_mean_responses{1}] = remove_odor(F_time_series{1}, F_mean_responses{1}, 37); 

for(i=1:length(F_mean_responses))
    show_heatmap(F_time_series{i}.odor_names, F_mean_responses{i}(sort_by_max_response(F_mean_responses{i}),:));
end


%Match all animals to the first:
target_index = 3;
%constraints.functional_threshold = "off";
%constraints.spatial_threshold = "off";

constraints.functional_threshold = 0.3;
constraints.spatial_threshold = 0.5;

[mapping_all, mapped_indices_all] = register_all_to_target(normalize_mean_responses(extract_odors(F_mean_responses,1:14)), F_cc, constraints, target_index);
matched = apply_glomerulus_matching_to_all(F_time_series, normalize_mean_responses(F_mean_responses), F_cc, mapping_all, mapped_indices_all, target_index);

%ZZ: save the matching results, so can be used later
%for target_index = 1:size(brains,2)
    %matching_index = all_mapped{target_index};
    matching_index = cell2mat(mapping_all);
    fn_save = fullfile(path_matching, strcat(experiment_name,'_', brains(target_index), '_asRefMatchingIndex.mat'));
    save(fn_save{1}, 'matching_index');
%end



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
   hold on; volume_plot(matched{i}.F_cc.X_vis_unmatched, matched{i}.F_cc.parameters,{-37.5, 30, 'gray'});
end


%Mean animal odor and mean human odor mapped to the target AL:
animal_mean = mean(average_heatmap(:,29:33),2); 
human_mean = mean(average_heatmap(:,37:44),2);
show_3D_heatmap(matched{target_index}.F_cc, [animal_mean, human_mean], 1, {'transparent', 60, 30});
show_3D_heatmap(matched{target_index}.F_cc, [animal_mean, human_mean], 2, {'transparent', 60, 30});

