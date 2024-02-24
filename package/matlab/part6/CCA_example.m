%%%
data_path      = 'D:\zhilei_data\humanNonhumanOrco_new_preprocessing\odorEvoked\';
meta_data_file = 'humanNonhumanOrco_meta_data.mat';

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
[coordinates, names] = cca_figure(F_time_series, F_mean_responses, masks, odor_subset, 'D:\figure_data.csv');