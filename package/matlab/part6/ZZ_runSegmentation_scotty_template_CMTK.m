%code by: Martin.Strauch@lfb.rwth-aachen.de
%Modified by Zhilei Zhao to run on scotty
%run the preprocessing and segmentation (convex cone) for odor evoked movies


%1) Process (reading, filtering, SVD) all animals in a folder and save results as .mat files. 
%--------------------------------------------------------------------------------------------
%This step has high computation time and memory requirements. You may want to run this on a cluster computer. 
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

folder_home = '/mnt/bucket/labs/mcbride/zhileizhao/Scotty/Analysis/CalciumImaging/AnalysisPipelineOrco/';
cd('FOLDERWORKING');
addpath(genpath('FOLDERWORKING'));
files_to_read = '*-mc.tif.registered.tif'; %point to the files that should be read, e.g. only motion-corrected tif-files.
subfolder     = 'odorEvoked';  %what type of data to analyze
in_path = fullfile(folder_home, 'MartinStrauch', 'MotionCorrectedData', experiment_name);
if exist(in_path, 'dir')
    rmdir(in_path, 's');
end
mkdir(in_path);
out_path_svd = fullfile(folder_home, 'MartinStrauch', 'SVDResults', experiment_name);
if exist(out_path_svd, 'dir')
    rmdir(out_path_svd, 's');
end
mkdir(out_path_svd);
out_path_segment = fullfile(folder_home, 'MartinStrauch', 'SegmentationResults', experiment_name);
if exist(out_path_segment, 'dir')
    rmdir(out_path_segment, 's');
end
mkdir(out_path_segment);

%copy files from MotionCorrection folder to the working folder
for brain_idx = 1:length(brains)
    in_path_source = fullfile(folder_home, 'Registration', 'RegisteredMovies', datasets(brain_idx), brains(brain_idx));
    in_path_target = fullfile(in_path, brains(brain_idx), subfolder);
    in_path_source = string(in_path_source);
    in_path_target = string(in_path_target);
    mkdir(in_path_target);
    copyfile(in_path_source, in_path_target);
end

%run SVD, save results in the 'SVDResults' folder
in_path = strcat(in_path, '/');
ZZ_process_folder_orco(in_path, subfolder, files_to_read, out_path_svd, dim_x, dim_y, dim_z, 1:20, 13, 5); 



%2) Load the saved .mat files from <out_path> and compute the 3D glomerulus maps:
%--------------------------------------------------------------------------------
%Save results in the SegmentationResults folder
%Load meta data (file names matched to odor names, information on whether the AL should be flipped)
meta_data_fn = fullfile(folder_home, 'MartinStrauch', 'MetaData', meta_file);
meta_data = load_meta_data(meta_data_fn);
%Set the desired number of clusters (glomeruli):
num_clusters = 40; 
%border_width.x/.y voxels at the borders of the volume will be cut out (to remove motion artifacts):
border_width.x = 9; border_width.y = 12;
%linux path and windows path name handled differently
out_path_svd_temp = strcat(out_path_svd, '/');
F_cc = compute_glomerulus_maps(out_path_svd_temp, num_clusters, meta_data, border_width);
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
