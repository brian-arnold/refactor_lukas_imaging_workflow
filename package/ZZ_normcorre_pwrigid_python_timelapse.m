%% rigid motion correction with NormCorre for time-lapse movies (spontaneous & odor-evoked)
cd('FOLDER_HOME/MotionCorrection');
addpath(genpath('FOLDER_HOME/MotionCorrection/NoRMCorre'));
addpath(genpath('FOLDER_CODE'));
parpool('local',22);
dirData = 'RAWDATAFOLDER/DATASETNAME/BRAINNAME/MOVIETYPE';
dirRes = 'SAVEFOLDER/DATASETNAME/BRAINNAME/MOVIETYPE';
if(~exist(dirRes,'dir'))
    mkdir(dirRes);
end

%% info for the tiff volume
% timePoint doesn't matter here, since it's infered based on matrix size
volume_info = struct('xRes', 128, 'yRes', 128, 'nChannel', 1, ...
    'zStacks', NUMSTACKS, 'timePoint', 118, 'dtype', 'int16', 'offset',[0]);

%% parameters for normcorre
mc_option.max_shift = [30 30 3];
% mc_option.max_shift = [15 15 3];
mc_option.us_fac = 20;  % up-sampling factor
mc_option.init_batch = 20;  % num of volumes use to estimate the initial template
mc_option.iter = 5;  % num of iterations
mc_option.n_pad = 1; % how many zero stacks to pad at the start and end of volume
% to make grid, N*grid_size + overlap_pre = dim
% the patch will be 56X56X7, num of patches will be 3X3X4
mc_option.grid_size = [36 36 5]; % size of non-overlapping regions
mc_option.overlap_pre = [20 20 2]; % size of overlapping region
mc_option.min_patch_size = [24 24 4];  % minimum size of patch (default: [32,32,16])
mc_option.mot_uf = [4 4 1]; % degree of patches upsampling
mc_option.max_dev = [3 3 1]; % maximum deviation of patch shift from rigid shift
mc_option.use_parallel = true;

%% get file names for the data dir
fileNamesTiff = dir(fullfile(dirData,'*.tif'));
fileNames = fullfile(dirData, {fileNamesTiff.name}');

[shifts] = ZZ_pwrigid_normcorre_batch(fileNames, volume_info, 1, mc_option, dirRes, true);
