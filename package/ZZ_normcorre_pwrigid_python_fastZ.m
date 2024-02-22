%% rigid motion correction with NormCorre for fastZ
cd('/mnt/cup/labs/mcbride/lukas/CalciumImaging/FCV/Analysis/MotionCorrection');
addpath(genpath('/mnt/cup/labs/mcbride/lukas/CalciumImaging/FCV/Analysis/MotionCorrection/NoRMCorre'));
parpool('local',22);
dirData = 'BeforeCorrection/DATASETNAME/BRAINNAME/fastZ';
dirRes = 'SAVEFOLDER/DATASETNAME/BRAINNAME/fastZ';

if(~exist(dirRes,'dir'))
    mkdir(dirRes);
end

%% info for the tiff volume
volume_info = struct('xRes', 256, 'yRes', 256, 'nChannel', 1, ...
    'zStacks', NUMSTACKS, 'timePoint', 20, 'dtype', 'int16', 'offset',[0]);

%% parameters for normcorre for fastZ high-res
mc_option.max_shift = [30 30 4];
mc_option.us_fac = 20;  % up-sampling factor
mc_option.init_batch = 3;  % num of volumes use to estimate the initial template
mc_option.iter =10;  % num of iterations
mc_option.n_pad = 1; % how many zero stacks to pad at the start and end of volume 
% to make grid, N*grid_size + overlap_pre = dim
% the patch will be 56X56X7, num of patches will be 3X3X4
mc_option.grid_size = [72 72 23]; % size of non-overlapping regions
mc_option.overlap_pre = [40 40 9]; % size of overlapping region
mc_option.min_patch_size = [56 56 12];  % minimum size of patch (default: [32,32,16])
mc_option.mot_uf = [4 4 1]; % degree of patches upsampling
mc_option.max_dev = [3 3 1]; % maximum deviation of patch shift from rigid shift
mc_option.use_parallel = true; 

%% get file names for the data dir
fileNamesTiff = dir(fullfile(dirData,'*.tif'));
fileNames = fullfile(dirData, {fileNamesTiff.name}');

[shifts] = ZZ_pwrigid_normcorre_batch(fileNames, volume_info, 1, mc_option, dirRes, true);

%% make a mean volume from the motion corrected volume
% read tiff and reshape into 4D (xyczt)
num_to_ave = 5;
for i=1:size(fileNames,1)
    tempNames = strsplit(fileNames{i},'/');
    mcFileName = strcat(fullfile(dirRes,tempNames{end}), '-mc.tif');
    mc_volume_info = volume_info;
    mc_volume_info.zStacks = mc_volume_info.zStacks + 2*mc_option.n_pad;
    v = ZZ_read_reshape_tiff(mcFileName, mc_volume_info);
    cor = csvread(strcat(fullfile(dirRes,tempNames{end}), '-correlation.csv'));
    % get the volumes with the top 3 correlation
    [val ind] = sort(cor,'descend');
    I = ind(1:num_to_ave);
    v_top = squeeze(v(:,:,1,:,I));
    % make a mean volume
    v_mean = mean(v_top, 4);
    % save mean volume
    options.overwrite = true;
    saveastiff(single(v_mean), strcat(fullfile(dirRes,tempNames{end}), '-mean.tif'), options);
    % adjust contrast and convert to uint8
    v_mean_int8 = scanImageToUint8(v_mean, 0.001, 0.05);
    saveastiff(v_mean_int8, strcat(fullfile(dirRes,tempNames{end}), '-mean-int8.tif'), options);
end
