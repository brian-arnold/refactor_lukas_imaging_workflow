% Zhilei Zhao
% Merge redundant glomeruli
% this should be run after the convex cone algorithm on segmentation

% dataset specific information, change this part to your dataset information
experiment_name = 'FCV_GCaMP6f';
brains = {'220505_FCV_orco_GCaMP6f_F1', '220506_FCV_orco_GCaMP6f_F4', '220506_FCV_orco_GCaMP6f_F5', '220506_FCV_orco_GCaMP6f_F6'};
datasets = {'220505', '220506', '220506', '220506'};
meta_file = 'FCV_GCaMP6f_meta_data.mat';
% end of dataset specific information

% sometimes we may need to remove an odor puff from the data, e.g. a
% mistake or low quality puff
% remove_odor_brain = {1};
% remove_odor_odor = {37};
remove_odor_brain = {}; 
remove_odor_odor = {};

% more universal setting
dim_x = 128; dim_y = 128; dim_z = 24;
%Set the desired number of clusters (glomeruli):
num_clusters = 30;

% folder settings: may need to change to your PATH
% base folder in Scotty server
folder_home = '/mnt/cup/labs/mcbride/lukas/CalciumImaging/FCV/Analysis/';
% if running script on lab desktop
%folder_home = '/Volumes/mcbride/zhileizhao/Scotty/Analysis/CalciumImaging/AnalysisPipelineOrco/';

%cd('Volumes/mcbride/lukas/CalciumImaging/FCV/Analysis/MartinStrauch_hybridNewMerging/glomerulus_matching_scotty');
addpath(genpath('/mnt/cup/labs/mcbride/lukas/CalciumImaging/FCV/Analysis/MartinStrauch_hybridNewMerging/glomerulus_matching_scotty/'));

% folder where segmentation results are stored
out_path_segment = fullfile(folder_home, 'MartinStrauch_hybridNewMerging', 'SegmentationResults', experiment_name);
%Load meta data (file names matched to odor names, information on whether the AL should be flipped)
meta_data_fn = fullfile(folder_home, 'MartinStrauch_hybridNewMerging', 'MetaData', meta_file);
meta_data = load_meta_data(meta_data_fn);
fn_cc = fullfile(out_path_segment, strcat(experiment_name,'_convex_cone_res.mat'));
%load the convex cone results
load(fn_cc);


%3) Compute mean odor responses
%-----------------------------
%Glomerulus time series reordered s.t. we have the same odor sequence in each animal. The sorted name list is in .odor_names
F_time_series = apply_odor_sorting(F_cc, meta_data);

% Optional: remove an odor (here the odor labeled as 'mistake' in animal 1).
% Sorting and removing surplus odors is relevant for glomerulus matching, where we need the same odor sequence for all animals.
if length(remove_odor_brain)>0
    % note that for the same brain, removing one odor change the overall
    % index, need to count how many has been removed
    count_removed = zeros(length(brains),1);
    for ii=1:length(remove_odor_brain)
        brain_idx = remove_odor_brain{ii};
        odor_idx = remove_odor_odor{ii} - count_removed(brain_idx);
        [F_time_series{brain_idx}] = ZZ_remove_odor(F_time_series{brain_idx}, odor_idx);
        count_removed(brain_idx) = count_removed(brain_idx) + 1;
    end
end

%Compute mean odor responses from the time series:
regime1 = 1:20;                         %in the example, the first 28 odors are single odorants
% regime1_intervals.reference = [10:20];  %baseline before odor stimulation
% regime1_intervals.signal =  [30:45];    %time points at which the odor response is extracted
regime1_intervals.reference = [10:20];  %baseline before odor stimulation
regime1_intervals.signal =  [27:40];   %time points at which the odor response is extracted
%Use different intervals for regime2 (comprises all odor numbers not in regime1): In the example, these are odor blends
% regime2_intervals.reference = [10:150]; 
% regime2_intervals.signal = [190:225];
%ZZ extended interval
regime2_intervals.reference = []; 
regime2_intervals.signal = [];

F_mean_responses = compute_mean_responses(F_time_series, regime1, regime1_intervals, regime2_intervals);

%Remove redundant (w.r.t. the reference odor responses) clusters in the target AL: (note: this function is still preliminary)
clustering_threshold = 0.001; % the higher, the more glomeruli will be merged
for(i=1:length(F_cc))
    [F_time_series{i}, F_mean_responses{i}, F_cc{i}] = remove_redundant_clusters(F_time_series{i}, F_mean_responses{i}, F_cc{i}, clustering_threshold);
    % save new segmentation
    figure, volume_plot(F_cc{i}.X_vis, F_cc{1}.parameters);
    %save figure
    fig_fn = fullfile(out_path_segment, strcat(brains(i),'_segmentRemoved.fig'));
    savefig(gcf, string(fig_fn));
    % save the new segmentation for the target brain as tif
    cc_res = F_cc{i};
    cc_bi = cc_res.X_vis;
    cc_tif = save_cluster(cc_bi, dim_x, dim_y, dim_z);
    % save as tif
    options.overwrite = true;
    fname = fullfile(out_path_segment, strcat(brains{i}, '_segmentRemoved.tif'));
    saveastiff(single(cc_tif), fname, options);
end

exit(0);