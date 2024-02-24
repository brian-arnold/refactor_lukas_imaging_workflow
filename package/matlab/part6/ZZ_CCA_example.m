% experiment specific parameters
experiment_name = 'humanNonhumanOrco';
brains = {'humanNonhumanOrcoF1' ,'humanNonhumanOrcoF2' ,'humanNonhumanOrcoF3' ,'humanNonhumanOrcoF4' ,'humanNonhumanOrcoF5'};
%datasets = {'20191019_humanNonhumanOrcoNew2' ,'20191022_humanNonhumanOrcoNew2' ,'20191023_humanNonhumanOrcoNew2' ,'20191025_humanNonhumanOrcoNew2' ,'20191027_humanNonhumanOrcoNew2'};
meta_file = 'humanNonhumanOrco_meta_data.mat';

folder_home = '/mnt/bucket/labs/mcbride/zhileizhao/Scotty/Analysis/CalciumImaging/AnalysisPipelineOrco/';
cd('/mnt/bucket/labs/mcbride/zhileizhao/Scotty/Analysis/CalciumImaging/AnalysisPipelineOrco/MartinStrauch_hybridNewMerging/glomerulus_matching_scotty');
addpath(genpath('/mnt/bucket/labs/mcbride/zhileizhao/Scotty/Analysis/CalciumImaging/AnalysisPipelineOrco/MartinStrauch_hybridNewMerging/glomerulus_matching_scotty'));
files_to_read = '*-mc.tif.between.tif'; %point to the files that should be read, e.g. only motion-corrected tif-files.
subfolder     = 'odorEvoked';  %what type of data to analyze

out_path_cca = fullfile(folder_home, 'MartinStrauch_hybridNewMerging', 'CCAResults', experiment_name);
if ~exist(out_path_cca, 'dir')
    mkdir(out_path_cca);
end
out_file_cca = fullfile(out_path_cca, 'figure_data.csv');

data_path = fullfile(folder_home, 'MartinStrauch_hybridNewMerging', 'SVDResults', experiment_name);
data_path = strcat(data_path, '/');
meta_data_file = fullfile(folder_home, 'MartinStrauch_hybridNewMerging', 'MetaData', meta_file);

%%%
%data_path      = 'D:\zhilei_data\humanNonhumanOrco_new_preprocessing\odorEvoked\';
%meta_data_file = 'humanNonhumanOrco_meta_data.mat';

%regime1: single odor measurements named e.g. "S1_m1", regime 2: blends: all other measurements                    
intervals.regime1_odors = '[A-Z][0-9]_m'; 
%signal and reference (baseline) intervals for computing mean odor responses:
intervals.regime1_reference = [10:20];  
intervals.regime1_signal    = [30:45];    
intervals.regime2_reference = [10:100]; 
intervals.regime2_signal    = [190:225];
%%%

%compute mean odor responses for all voxels:
[F_time_series, F_mean_responses, masks] = compute_voxel_responses(data_path, meta_data_file, intervals);

%this only works if we have the same odors in all data sets
%--> here, we need to remove an additional odor (nr. 37) from data set 1:
[F_time_series{1}, F_mean_responses{1}] = remove_odor(F_time_series{1}, F_mean_responses{1}, 37);

%restrict the analysis to these odor measurements:
odor_subset=29:44; %here: only blends
F_time_series{1}.odor_names(odor_subset)

%plot figure and write data into the specified csv file:
%[coordinates, names] = cca_figure(F_time_series, F_mean_responses, masks, odor_subset, 'D:\figure_data.csv');
[coordinates, names] = cca_figure(F_time_series, F_mean_responses, masks, odor_subset, out_file_cca);