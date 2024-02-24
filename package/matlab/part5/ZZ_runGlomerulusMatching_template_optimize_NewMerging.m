% run Martin's glomerulus matching algorithm
% this should be run after the convex cone algorithm on segmentation

%ZZ: experiment specific settings
% experiment_name = 'humanNonhumanOrco';
% brains = {'humanNonhumanOrcoF1', 'humanNonhumanOrcoF2', 'humanNonhumanOrcoF3', 'humanNonhumanOrcoF4', 'humanNonhumanOrcoF5'};
% datasets = {'20191019_humanNonhumanOrco', '20191022_humanNonhumanOrco', '20191023_humanNonhumanOrco', '20191025_humanNonhumanOrco', '20191027_humanNonhumanOrco'};
% meta_file = 'humanNonhumanOrco_meta_data.mat';
experiment_name = 'EXPERIMENTNAME';
brains = {BRAINS};
datasets = {DATASETS};
meta_file = 'METAFILE';
% all movies are used for segmentation; some movies may need to be removed
% for matching 
% remove_odor_brain = {1};
% remove_odor_odor = {37};
remove_odor_brain = {REMOVE_ODOR_BRAIN}; 
remove_odor_odor = {REMOVE_ODOR_ODOR};

% more universal setting
dim_x = 128; dim_y = 128; dim_z = 24;
%Set the desired number of clusters (glomeruli):
num_clusters = 35;

% folder settings
% base folder in scotty
folder_home = '/mnt/bucket/labs/mcbride/zhileizhao/Scotty/Analysis/CalciumImaging/AnalysisPipelineOrco/';
% if running script on lab desktop
%folder_home = '/Volumes/mcbride/zhileizhao/Scotty/Analysis/CalciumImaging/AnalysisPipelineOrco/';
cd('FOLDERWORKING');
addpath(genpath('FOLDERWORKING'));

% folder where segmentation results are stored
out_path_segment = fullfile(folder_home, 'MartinStrauch_hybridNewMerging', 'SegmentationResults', experiment_name);
% folder to store matching results
path_matching = fullfile(folder_home, 'MartinStrauch_hybridNewMerging', 'MatchingResults', experiment_name);
if exist(path_matching, 'dir')
    rmdir(path_matching, 's');
end
mkdir(path_matching);
%Load meta data (file names matched to odor names, information on whether the AL should be flipped)
meta_data_fn = fullfile(folder_home, 'MartinStrauch_hybridNewMerging', 'MetaData', meta_file);
meta_data = load_meta_data(meta_data_fn);
fn_cc = fullfile(out_path_segment, strcat(experiment_name,'_convex_cone_res.mat'));
%load the convex cone results
load(fn_cc);


%3) Compute mean odor responses
%------------------------------

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
regime1 = 1:28;                         %in the example, the first 28 odors are single odorants
% regime1_intervals.reference = [10:20];  %baseline before odor stimulation
% regime1_intervals.signal =  [30:45];    %time points at which the odor response is extracted
regime1_intervals.reference = [SINGLE_BASELINE];  %baseline before odor stimulation
regime1_intervals.signal =  [SINGLE_PEAK];    %time points at which the odor response is extracted
%Use different intervals for regime2 (comprises all odor numbers not in regime1): In the example, these are odor blends
% regime2_intervals.reference = [10:150]; 
% regime2_intervals.signal = [190:225];
%ZZ extended interval
regime2_intervals.reference = [MARKES_BASELINE]; 
regime2_intervals.signal = [MARKES_PEAK];

F_mean_responses = compute_mean_responses(F_time_series, regime1, regime1_intervals, regime2_intervals);

% save the original time-series and mean response
fn_save = fullfile(path_matching, strcat(experiment_name,'_F_time_series.mat'));
save(fn_save, 'F_time_series');
fn_save = fullfile(path_matching, strcat(experiment_name,'_F_mean_response.mat'));
save(fn_save, 'F_mean_responses');


%4) Glomerulus matching 
%-----------------------

%4a) Finding an assignment of glomeruli between animals:

%Optimize functional and spatial distance with these constraints. 
%The stricter the constraints, the fewer glomeruli can be matched.
constraints.functional_threshold = FUNCTIONAL; %Functional distances (odor response correlations) < functional_threshold will be penalized
constraints.spatial_threshold = SPATIAL;     %Spatial distances (Euclidean distances in the 3D volume) > spatial_threshold will be penalized

% make a copy of the original results, because it will change later
F_cc_original = F_cc; 
F_time_series_original = F_time_series;
F_mean_responses_original = F_mean_responses;

%Remove redundant (w.r.t. the reference odor responses) clusters in the target AL: (note: this function is still preliminary)
% clustering_threshold = MERGING; % the higher, the more glomeruli will be merged
% for(i=1:length(F_cc))
%     %[F_time_series{i}, F_mean_responses{i}, F_cc{i}] = remove_redundant_clusters(F_time_series{i}, F_mean_responses{i}, F_cc{i}, clustering_threshold);
%     [F_time_series{i}, F_mean_responses{i}, F_cc{i}] = ZZ_remove_redundant_clusters(F_time_series{i}, F_mean_responses{i}, F_cc{i}, clustering_threshold, num_clusters);
%     % save new segmentation
%     figure, volume_plot(F_cc{i}.X_vis, F_cc{1}.parameters);
%     %save figure
%     fig_fn = fullfile(out_path_segment, strcat(brains(i),'_segmentRemoved.fig'));
%     savefig(gcf, string(fig_fn));
%     % save the new segmentation for the target brain as tif
%     cc_res = F_cc{i};
%     cc_bi = cc_res.X_vis;
%     cc_tif = save_cluster(cc_bi, dim_x, dim_y, dim_z);
%     % save as tif
%     options.overwrite = true;
%     fname = fullfile(out_path_segment, strcat(brains{i}, '_segmentRemoved.tif'));
%     saveastiff(single(cc_tif), fname, options);
% end

%new merging approach
%optional: merge clusters
correlation_threshold = MERGING; %adjacent clusters with a response correlation >= threshold will be merged
overlap_threshold = OVERLAP;
for(i=1:length(F_cc))
    [F_time_series{i}, F_mean_responses{i}, F_cc{i}] = merge_redundant_clusters2(F_time_series{i}, F_mean_responses{i}, F_cc{i}, correlation_threshold, overlap_threshold);
end

% %merging creates bigger clusters that might now be adjacent --> it can make sense to merge again:
% correlation_threshold = 0.9; %more strict threshold criterion to avoid arbitrarly large clusters
% for(i=1:length(F_cc))
%     [F_time_series{i}, F_mean_responses{i}, F_cc{i}] = merge_redundant_clusters(F_time_series{i}, F_mean_responses{i}, F_cc{i}, correlation_threshold);
% end

% save merged glomeruli segmentation
for(i=1:length(F_cc))
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

%ZZ: map all brains to one reference brain 
%ZZ: choose different brain as reference, see how the numbers change 
% save matching relationship for each reference brain
all_mapped = {};
for target_index = 1:size(brains,2)
%for target_index = 3:3
% run the matching algorithm
    % save all mapped index into one matrix
    if REFONLY==0
        [mapping_all, mapped_indices_all] = register_all_to_target(normalize_mean_responses(F_mean_responses), F_cc, constraints, target_index);
    else
        [mapping_all, mapped_indices_all] = register_all_to_target(normalize_mean_responses(extract_odors(F_mean_responses,1:28)), F_cc, constraints, target_index);
    end
    % apply matching to time series and response data
    matched = apply_glomerulus_matching_to_all(F_time_series, normalize_mean_responses(F_mean_responses), F_cc, mapping_all, mapped_indices_all, target_index);
    % save matching relationship
    temp = cell2mat(mapping_all);
    all_mapped{target_index} = temp;
    % save 3D glomerulus maps with matched colors:
    figure;
    ha = tight_subplot(1, length(matched), [.01 .03], [.1 .01], [.01 .01]);
    for(i=1:length(matched))
        axes(ha(i));
        volume_plot(matched{i}.F_cc.X_vis, matched{i}.F_cc.parameters); 
    end
    fig_fn = fullfile(path_matching, strcat(brains(target_index),'_asRef'));
    fig_fn = fig_fn{1};
    fig = gcf;
    fig.PaperUnits = 'inches';
    fig.PaperPosition = [0 0 10*length(matched) 10];
    print(fig_fn,'-dpng','-r300')
end

%ZZ: save the matching results, so can be used later
for target_index = 1:size(brains,2)
    matching_index = all_mapped{target_index};
    fn_save = fullfile(path_matching, strcat(experiment_name,'_', brains(target_index), '_asRefMatchingIndex.mat'));
    save(fn_save{1}, 'matching_index');
end

%count how many are mapped, plot out the results
figure;
for target_index = 1:size(brains,2)
    count_mapped = [];
    mapped_temp = all_mapped{target_index};
    count_this  = sum(mapped_temp~=0, 1);
    for count_thre = 2:size(brains,2)
        temp =  sum(count_this>=count_thre);
        count_mapped = [count_mapped, temp];
    end
    plot(2:size(brains,2), count_mapped, '-o', 'LineWidth', 2, 'MarkerSize',8);
    xlabel("Matched in at least No. brains", 'FontSize', 18);
    ylabel("No. glomeruli", 'FontSize', 18);
    ylim([0,40]);
    hold on;
end
legend(brains);
%ZZ: save the figure
fig_fn = fullfile(path_matching, strcat(experiment_name,'_matching.fig'));
savefig(gcf, string(fig_fn));

exit(0);