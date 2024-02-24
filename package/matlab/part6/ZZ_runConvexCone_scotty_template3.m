% modified from template2, by supplementing a baseline mask

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

folder_home = '/mnt/cup/labs/mcbride/lukas/CalciumImaging/FCV/Analysis/';
cd('FOLDERWORKING');
addpath(genpath('FOLDERWORKING'));

out_path_svd = fullfile(folder_home, 'MartinStrauch', 'SVDResults', experiment_name);
out_path_segment = fullfile(folder_home, 'MartinStrauch', 'SegmentationResults', experiment_name);
if exist(out_path_segment, 'dir')
    rmdir(out_path_segment, 's');
end
mkdir(out_path_segment);

%2) Load the saved .mat files from <out_path> and compute the 3D glomerulus maps:
%--------------------------------------------------------------------------------
%Save results in the SegmentationResults folder
%Load meta data (file names matched to odor names, information on whether the AL should be flipped)
meta_data_fn = fullfile(folder_home, 'MartinStrauch', 'MetaData', meta_file);
meta_data = load_meta_data(meta_data_fn);

% load the mask files
path_baseline = fullfile(folder_home, 'MartinStrauch', 'BaselineFluoData', experiment_name);
masks = {};
threshold = 7; 
for ii=1:length(brains)
    fn = fullfile(path_baseline, strcat(brains{ii}, '_baseline.tif'));
    m = read_tiff_stack(fn);
    p = single(m>=threshold);
    p_filter = zeros(size(p));
    for z=1:size(p,3)
        temp = imgaussfilt(p(:,:,z),2);
        p_filter(:,:,z) = temp;
    end
    masks{ii} = p_filter;
end

% 
% options.overwrite = true;
% fname = fullfile(folder_home, 'brainF3_maskMatlab.tif');
% saveastiff(single(masks{3}), fname, options);

%Set the desired number of clusters (glomeruli):
num_clusters = 40; 
selected_odors = '[A-Z][0-9]_m1'; %perform clustering only on the first response of the reference odors, i.e. U1_m1, T1_m1, ...
%linux path and windows path name handled differently
out_path_svd_temp = strcat(out_path_svd, '/');
F_cc = ZZ_compute_glomerulus_maps2(out_path_svd_temp,num_clusters, meta_data, selected_odors, masks);

%ZZ: save the convex cone results, so can be used later
% fn_cc = fullfile(out_path_segment, strcat(experiment_name,'_convex_cone_res.mat'));
% save(fn_cc, 'F_cc', '-v7.3');

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

exit(0);