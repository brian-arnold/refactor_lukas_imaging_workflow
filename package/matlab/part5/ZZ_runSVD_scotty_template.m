%code by: Martin.Strauch@lfb.rwth-aachen.de
%Modified by Zhilei Zhao to run on PNI scotty
%run the preprocessing (SVD) for odor evoked movies


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
experiment_name = 'EXPERIMENTNAME';
brains = {BRAINS};
datasets = {DATASETS};

folder_home = '/mnt/bucket/labs/mcbride/zhileizhao/Scotty/Analysis/CalciumImaging/AnalysisPipelineOrco/';
cd('FOLDERWORKING');
addpath(genpath('FOLDERWORKING'));
files_to_read = '*-mc.tif.between.tif'; %point to the files that should be read, e.g. only motion-corrected tif-files.
subfolder     = 'odorEvoked';  %what type of data to analyze
in_path = fullfile(folder_home, 'MartinStrauch', 'MotionCorrectedData', experiment_name);
out_path_svd = fullfile(folder_home, 'MartinStrauch', 'SVDResults', experiment_name);
if exist(out_path_svd, 'dir')
    rmdir(out_path_svd, 's');
end
mkdir(out_path_svd);

%run SVD, save results in the 'SVDResults' folder
in_path = strcat(in_path, '/');
ZZ_process_folder_orco(in_path, subfolder, files_to_read, out_path_svd, dim_x, dim_y, dim_z, 1:20, 13, 5); 

exit(0);