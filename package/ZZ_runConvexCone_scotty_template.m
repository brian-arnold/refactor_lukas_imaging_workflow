%code by: Martin.Strauch@lfb.rwth-aachen.de
%Modified by Zhilei Zhao to run on PNI scotty
%run the segmentation (convex cone) for odor evoked movies


%1) Input parameters
%Assumed path structure: <path>\animal_name\<subfolder>\, where <subfolder> is the experiment name, e.g. 'odorEvoked'

%Image stack dimensions:
dim_x = 128; dim_y = 128; dim_z = 24;

%ZZ: experiment specific settings
% experiment_name = 'humanNonhumanOrco';
% brains = {'humanNonhumanOrcoF1', 'humanNonhumanOrcoF2', 'humanNonhumanOrcoF3', 'humanNonhumanOrcoF4', 'humanNonhumanOrcoF5'};
% datasets = {'20191019_humanNonhumanOrco', '20191022_humanNonhumanOrco', '20191023_humanNonhumanOrco', '20191025_humanNonhumanOrco', '20191027_humanNonhumanOrco'};
% meta_file = 'humanNonhumanOrco_meta_data.mat';
experiment_name = 'EXPERIMENTNAME';
brains = {BRAINS};
datasets = {DATASETS};
meta_file = 'METAFILE';

folder_home = 'FOLDER_HOME';
cd('FOLDERWORKING');
addpath(genpath('FOLDERWORKING'));
addpath(genpath('FOLDER_CODE_MATLAB'));

files_to_read = '*-mc.tif.between.tif'; %point to the files that should be read, e.g. only motion-corrected tif-files.
subfolder     = 'odorEvoked';  %what type of data to analyze

out_path_svd = fullfile(folder_home, 'MartinStrauch_hybridNewMerging', 'SVDResults', experiment_name);
out_path_segment = fullfile(folder_home, 'MartinStrauch_hybridNewMerging', 'SegmentationResults', experiment_name);
if exist(out_path_segment, 'dir')
    rmdir(out_path_segment, 's');
end
mkdir(out_path_segment);

%2) Load the saved .mat files from <out_path> and compute the 3D glomerulus maps:
%--------------------------------------------------------------------------------
%Save results in the SegmentationResults folder
%Load meta data (file names matched to odor names, information on whether the AL should be flipped)
meta_data_fn = fullfile(folder_home, 'MartinStrauch_hybridNewMerging', 'MetaData', meta_file);
meta_data = load_meta_data(meta_data_fn);
%Set the desired number of clusters (glomeruli):
num_clusters = 15; 
%num_clusters = 35; 
%border_width.x/.y voxels at the borders of the volume will be cut out (to remove motion artifacts):
border_width.x = 9; border_width.y = 12;
%linux path and windows path name handled differently
out_path_svd_temp = strcat(out_path_svd, '/');
%F_cc = compute_glomerulus_maps(out_path_svd_temp, num_clusters, meta_data, border_width);
selected_odors = '[A-Z][0-9]_m1'; %perform clustering only on the first response of the reference odors, i.e. U1_m1, T1_m1, ...
%F_cc = compute_glomerulus_maps2(out_path_svd_temp,num_clusters, meta_data, selected_odors);
[F_cc, C_index] = compute_glomerulus_maps2_returnCIndex(out_path_svd_temp,num_clusters, meta_data, selected_odors);
%if flipping ALs/cutting the border is not desired: F_cc = compute_glomerulus_maps(out_path, num_clusters, [], []);

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

%ZZ: save the index
for(i=1:length(F_cc))
    [x,y,z]=ind2sub([dim_x dim_y dim_z], C_index{i,1});
    pure_idx = [x;y;z]';
    fn_idx = fullfile(out_path_segment, strcat(brains{i},'_pureIndex.mat'));
    save(fn_idx, 'pure_idx');
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

exit(0);