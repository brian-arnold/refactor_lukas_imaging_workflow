%code by: Martin.Strauch@lfb.rwth-aachen.de
%
%Computing 3D glomerulus maps, time series signals and mean responses


%1) Process (reading, filtering, SVD) all animals in a folder and save results as .mat files. 
%--------------------------------------------------------------------------------------------
%This step has high computation time and memory requirements. You may want to run this on a cluster computer. 
%Assumed path structure: <path>\animal_name\<subfolder>\, where <subfolder> is the experiment name, e.g. 'odorEvoked'

%Image stack dimensions:
dim_x = 128; dim_y = 128; dim_z = 24;

%Adapt paths/file name endings:
in_path       = '/mnt/bucket/labs/mcbride/zhileizhao/Scotty/Analysis/CalciumImaging/AnalysisPipelineOrco/MartinStrauch/testData/humanNonhumanOrco/';
subfolder     = 'odorEvoked';
files_to_read = '*-mc.tif.between.tif'; %point to the files that should be read, e.g. only motion-corrected tif-files.
out_path      = '/mnt/bucket/labs/mcbride/zhileizhao/Scotty/Analysis/CalciumImaging/AnalysisPipelineOrco/MartinStrauch/testData/humanNonhumanOrco_res/';
process_folder_orco(in_path, subfolder, files_to_read, out_path, dim_x, dim_y, dim_z, 1:20, 13, 5); 
%ZZ: save other results in a different folder
out_path_other = '/mnt/bucket/labs/mcbride/zhileizhao/Scotty/Analysis/CalciumImaging/AnalysisPipelineOrco/MartinStrauch/testData/humanNonhumanOrco_resOther/';
mkdir(out_path_other);

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

%ZZ: save the convex cone results, so can be used later
fn_cc = fullfile(out_path_other, 'convex_cone_res.mat');
save(fn_cc, 'F_cc', '-v7.3');


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

%Plot heatmaps:
for(i = 1:length(F_time_series))
    names = F_time_series{i}.odor_names;
    figure, imagesc(global_scaling(F_mean_responses{i}, 4096)); 
    xticks(1:length(names)); xticklabels(names); xtickangle(90); set(gca,'TickLabelInterpreter','none'); set(gcf,'color','w');
end


%4) Glomerulus matching 
%-----------------------

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

%Match glomeruli using normalized mean odor responses and the 3D glomerulus maps in F_cc:
F_mean_responses_norm = normalize_mean_responses(F_mean_responses);
[subject, target, mapping, mapped_indices] = register_animals(F_mean_responses_norm, F_cc, constraints, subject_index, target_index);

%Plot the 3D glomerulus maps for subject and target animal with matched colors:
figure, volume_plot(subject.movie.X_vis(:,mapping(mapped_indices)), subject.movie.parameters); 
figure, volume_plot(target.movie.X_vis(:,mapped_indices), subject.movie.parameters);

