%code by: Martin.Strauch@lfb.rwth-aachen.de
%
%Computing 3D glomerulus maps, time series signals and mean responses


%1) Process (reading, filtering, SVD) all animals in a folder and save results as .mat files. 
%--------------------------------------------------------------------------------------------
%This step has high computation time and memory requirements. You may want to run this on a cluster computer. 
%Assumed path structure: <path>\animal_name\<subfolder>\, where <subfolder> is the experiment name, e.g. 'odorEvoked'

%Adapt paths/file name endings:
in_path       = 'D:\humanNonhumanOrco\';
subfolder     = 'odorEvoked';
files_to_read = '*-mc.tif.between.tif'; %point to the files that should be read, e.g. only motion-corrected tif-files.
out_path      = 'D:\zhilei_data\humanNonhumanOrco\odorEvoked\'
dim.x=128; dim.y=128; dim.z=24;
process_folder_orco(in_path, subfolder, files_to_read, out_path, dim.x, dim.y, dim.z, 1:20, 13, 5); 



%2) Load the saved .mat files from <out_path> and compute the 3D glomerulus maps:
%--------------------------------------------------------------------------------

%Load meta data (file names matched to odor names, information on whether the AL should be flipped)
meta_data = load_meta_data('humanNonhumanOrco_meta_data.mat');
%Set the desired number of clusters (glomeruli):
num_clusters = 40; 
%border_width.x/.y voxels at the borders of the volume will be cut out (to remove motion artifacts):
border_width.x = 9; border_width.y = 12;

F_cc = compute_glomerulus_maps(out_path, num_clusters, meta_data, border_width);
%if flipping ALs/cutting the border is not desired: F_cc = compute_glomerulus_maps(out_path, num_clusters, [], []);

%Plot the 3D glomerulus maps:
for(i=1:length(F_cc))
    figure, volume_plot(F_cc{i}.X_vis, F_cc{1}.parameters);
end



%3) Compute mean odor responses
%------------------------------

%Glomerulus time series reordered s.t. we have the same odor sequence in each animal. The sorted name list is in .odor_names
F_time_series = apply_odor_sorting(F_cc, meta_data);

%Compute mean odor responses from the time series:
regime1 = 1:28;                         %in the example, the first 28 odors are single odorants
regime1_intervals.reference = [10:20];  %baseline before odor stimulation
regime1_intervals.signal =  [30:45];    %time points at which the odor response is extracted
%Use different intervals for regime2 (comprises all odor numbers not in regime1): In the example, these are odor blends
regime2_intervals.reference = [10:150]; 
regime2_intervals.signal = [190:225];

F_mean_responses = compute_mean_responses(F_time_series, regime1, regime1_intervals, regime2_intervals);

%Optional: remove an odor (here the odor labeled as 'mistake' in animal 1).
%Sorting and removing surplus odors is relevant for glomerulus matching, where we need the same odor sequence for all animals.
[F_time_series{1}, F_mean_responses{1}] = remove_odor(F_time_series{1}, F_mean_responses{1}, 37); 

%Visualize odor reponses with a heatmap:
animal_nr=5
show_heatmap(F_time_series{animal_nr}.odor_names, F_mean_responses{animal_nr});

%3D heatmap of the mean odor response to odor <odor_nr> in animal <animal_nr>  
animal_nr=1; odor_nr=1;
show_3D_heatmap(F_cc{animal_nr}, F_mean_responses{animal_nr}, odor_nr);
show_3D_heatmap(F_cc{animal_nr}, F_mean_responses{animal_nr}, odor_nr, {'transparent', 60, 30});



%4) Glomerulus matching 
%-----------------------

%4a) Finding an assignment of glomeruli between animals:

%Animal F_cc{subject_index} will be mapped to animal F_cc{target_index}:
subject_index = 3; 
target_index  = 1;

%Optimize functional and spatial distance with these constraints. 
%The stricter the constraints, the fewer glomeruli can be matched.
constraints.functional_threshold = 0.3; %Functional distances (odor response correlations) < functional_threshold will be penalized
constraints.spatial_threshold = 10;     %Spatial distances (Euclidean distances in the 3D volume) > spatial_threshold will be penalized

%If we do not use constraints, all glomeruli will be matched:
%constraints.functional_threshold = "off";
%constraints.spatial_threshold = "off";

%Match glomeruli using mean odor responses and the 3D glomerulus maps in F_cc:
[subject, target, mapping, mapped_indices] = register_animals(F_mean_responses, F_cc, constraints, target_index, subject_index);
%Alternative: use only the responses to the reference odors (here: nr. 1 to 28):
[subject, target, mapping, mapped_indices] = register_animals(extract_odors(F_mean_responses,1:28), F_cc, constraints, target_index, subject_index);

%4b) Visualizations:

%Apply glomerulus matching consistently to time series, mean responses and 3D maps:
matched_subject = apply_glomerulus_matching(F_time_series, F_mean_responses, F_cc, mapping, mapped_indices, 'subject', subject_index);
matched_target  = apply_glomerulus_matching(F_time_series, F_mean_responses, F_cc, mapping, mapped_indices, 'target', target_index);

%Plot the 3D glomerulus maps for subject and target animal with matched colors:
figure, volume_plot(matched_target.F_cc.X_vis, matched_target.F_cc.parameters);
figure, volume_plot(matched_subject.F_cc.X_vis, matched_subject.F_cc.parameters); 

%Also plot the unmatched glomeruli (transparent gray):
figure, volume_plot(matched_subject.F_cc.X_vis, matched_subject.F_cc.parameters);
hold on; volume_plot(matched_subject.F_cc.X_vis_unmatched, matched_subject.F_cc.parameters,{-37.5, 30, 'gray',15});

%Compute and visualize distances between the spatial centroids/the odor response profiles of the matched glomeruli:
spatial_distances    = evaluate_spatial_match(matched_target.F_cc, matched_subject.F_cc);
functional_distances = evaluate_functional_match(zscore(matched_target.F_mean_responses,[],2), zscore(matched_subject.F_mean_responses,[],2))

%Heatmaps with matched glomeruli 
show_heatmap(matched_target.F_time_series.odor_names, matched_target.F_mean_responses);
show_heatmap(matched_subject.F_time_series.odor_names, matched_subject.F_mean_responses);

odor_nr1 = 30;
show_3D_heatmap(matched_subject.F_cc, matched_subject.F_mean_responses, odor_nr1, {'transparent', 60, 30});
show_3D_heatmap(matched_target.F_cc, matched_target.F_mean_responses, odor_nr1, {'transparent', 60, 30});
odor_nr2 = 40;
show_3D_heatmap(matched_subject.F_cc, matched_subject.F_mean_responses, odor_nr2, {'transparent', 60, 30});
show_3D_heatmap(matched_target.F_cc, matched_target.F_mean_responses, odor_nr2, {'transparent', 60, 30});


%4c) Match all animals to the same target animal:

target_index = 3;

%optional: merge clusters
%I. glomerulus merging based on odor response similarity:
clustering_threshold = 0.001; % the higher, the more glomeruli will be merged
for(i=1:length(F_cc))
    [F_time_series{i}, F_mean_responses{i}, F_cc{i}] = remove_redundant_clusters(F_time_series{i}, F_mean_responses{i}, F_cc{i}, clustering_threshold);
end

%II. alternative: merge only adjacent clusters:
correlation_threshold = 0.7; %adjacent clusters with a response correlation >= threshold will be merged
for(i=1:length(F_cc))
    [F_time_series{i}, F_mean_responses{i}, F_cc{i}] = merge_redundant_clusters(F_time_series{i}, F_mean_responses{i}, F_cc{i}, correlation_threshold);
end

constraints.functional_threshold = 0; 
constraints.spatial_threshold = 10;

%Match using all odors:
[mapping_all, mapped_indices_all] = register_all_to_target(normalize_mean_responses(F_mean_responses), F_cc, constraints, target_index);
%Match using only reference odors (1 to 28):
[mapping_all, mapped_indices_all] = register_all_to_target(normalize_mean_responses(extract_odors(F_mean_responses,1:28)), F_cc, constraints, target_index);
matched = apply_glomerulus_matching_to_all(F_time_series, normalize_mean_responses(F_mean_responses), F_cc, mapping_all, mapped_indices_all, target_index);

%Some visualizations:

%3D glomerulus maps with matched colors:
for(i=1:length(matched))
   figure, volume_plot(matched{i}.F_cc.X_vis, matched{i}.F_cc.parameters); 
   hold on; volume_plot(matched{i}.F_cc.X_vis_unmatched, matched{i}.F_cc.parameters,{-37.5, 30, 'gray',15});
end

for(i=1:length(matched))
   figure, volume_plot(matched{i}.F_cc.X_vis, matched{i}.F_cc.parameters); 
   spatial_distances  = evaluate_spatial_match(matched{target_index}.F_cc, matched{i}.F_cc);
end

[responses_2D, responses_3D] = combine_responses(matched);
%Heatmap for all animals concatenated (with NAs for glomeruli that are not matched)
show_heatmap(responses_2D.odor_names, responses_2D.responses)
%Heatmap for the average odor responses over all animals (NAs ignored):
average_heatmap = mean(responses_3D.responses,3,'omitnan');
show_heatmap(responses_3D.odor_names, average_heatmap);
pca_plot(average_heatmap, responses_3D.odor_names);

%Mean animal odor and mean human odor mapped to the target AL:
animal_mean = mean(average_heatmap(:,29:33),2); 
human_mean = mean(average_heatmap(:,37:44),2);
show_3D_heatmap(matched{target_index}.F_cc, [animal_mean, human_mean], 1, {'transparent', 60, 30});
show_3D_heatmap(matched{target_index}.F_cc, [animal_mean, human_mean], 2, {'transparent', 60, 30});



%Orco atlas:
load('atlas_orco_rotated.mat')
figure, volume_plot(atlas_orco_rotated, atlas_parameters)
figure, volume_plot(atlas_orco_rotated, atlas_parameters, {90, 30})
figure, volume_plot(F_cc{1}.X_vis, F_cc{1}.parameters)

