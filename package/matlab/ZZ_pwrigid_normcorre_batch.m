function [dummy] = ZZ_pwrigid_normcorre_batch(fileNames, volume_info, channel, mc_option, save_dir, plot_metrics)
% do motion correction on time-lapse data
% input a list of file names, info of the volume and which channel to use
% save motion correct data in tiff if specified
% returns nothing
gcp;

options.overwrite = true;
numFiles = size(fileNames, 1);

%% get parameters for normcorre rigid
max_shift = mc_option.max_shift;
us_fac = mc_option.us_fac;  %up-sampling factor
init_batch = mc_option.init_batch;  %num of volumes use to estimate the initial template
iter = mc_option.iter;  %num of iterations
n_pad = mc_option.n_pad;  %how many zero stacks to pad at the start and end of volume
grid_size = mc_option.grid_size; % size of non-overlapping regions
overlap_pre = mc_option.overlap_pre; % size of overlapping region
mot_uf = mc_option.mot_uf; % degree of patches upsampling
max_dev = mc_option.max_dev; % maximum deviation of patch shift from rigid shift
use_parallel = mc_option.use_parallel;

%% prelocate the arrays
% mc = zeros(volume_info.xRes, volume_info.yRes, volume_info.zStacks+2*n_pad, volume_info.timePoint, numFiles, 'single');
% template = zeros(volume_info.xRes, volume_info.yRes, volume_info.zStacks+2*n_pad, numFiles, 'single');
% shifts = struct('shifts',{},'shifts_up',{}, 'diff',{});

%% go through each file
for i=1:numFiles
    display(fileNames{i});
    %% read tiff and reshape into 4D (xyczt)
    v = ZZ_read_reshape_tiff(fileNames{i}, volume_info);
    % pad zero stacks to volume
    v_padded = zeros(volume_info.xRes, volume_info.yRes, volume_info.nChannel, ...
        volume_info.zStacks+2*n_pad, size(v,5), volume_info.dtype);
    v_padded(:,:,:,(1+n_pad):(volume_info.zStacks+n_pad),:) = v;
    mc = v_padded;

    %% calculate correction only on one channel
    Y_padded = squeeze(v_padded(:,:,channel,:,:));
    % convert to single precision
    Y_padded = single(Y_padded);

    % substract the minimal pixel value to make it positive to make the
    % algorithm work, but apply the shifts to the original data
    substracted_min = single(min(v(:)));
    Y_padded(:,:,(1+n_pad):(volume_info.zStacks+n_pad),:) = Y_padded(:,:,(1+n_pad):(volume_info.zStacks+n_pad),:) - substracted_min;
    disp('size of Y_padded');
    disp(size(Y_padded));
    %% set up parameter for normcorre rigid
    options_pwrigid = NoRMCorreSetParms('d1',size(Y_padded,1),'d2',size(Y_padded,2),...
        'd3',size(Y_padded,3), 'bin_width',size(v,5),'max_shift',max_shift,...
        'us_fac',us_fac,'init_batch',init_batch,'shifts_method','cubic', 'iter', iter, ...
        'grid_size',grid_size, 'overlap_pre', overlap_pre, 'mot_uf', mot_uf', 'max_dev', max_dev,...
        'use_parallel', use_parallel);

    %% perform the motion correction
    [mc_temp,shifts_temp,template_temp,options_pwrigid] = normcorre(Y_padded,options_pwrigid);

    %% apply shifts to each channel, add offset if specified
    for c=1:size(v_padded, 3)
        mc(:,:,c,:,:) = apply_shifts(squeeze(v_padded(:,:,c,:,:)), shifts_temp, options_pwrigid);
        mc(:,:,c,:,:) = mc(:,:,c,:,:) + volume_info.offset(c);
    end

%     shifts(:,i) = shifts_temp;

    %% save mc tiff and shift values if required
    if ~isempty(save_dir)
        % save tiff
        tempNames = strsplit(fileNames{i},'/');
        % save the template volume
        % add the substracted value back
        template_temp_added = template_temp + substracted_min;
        saveastiff(single(template_temp_added), strcat(fullfile(save_dir,tempNames{end}), '-template.tif'), options);
        temp = reshape(mc, volume_info.xRes, volume_info.yRes, []);
        saveastiff(temp, strcat(fullfile(save_dir,tempNames{end}), '-mc.tif'), options);
        % save shifts in N*3 table, corresponding to xyz
        shift_xyz = [];
        for j=1:size(shifts_temp,1)
            this_shifts = reshape(shifts_temp(j).shifts, [], 3);
            nrow = size(this_shifts,1);
            this_shifts = [this_shifts (1:nrow)' j*ones(nrow,1)];
            shift_xyz = [shift_xyz; this_shifts];
        end
        csvwrite(strcat(fullfile(save_dir,tempNames{end}), '-shifts.csv'),shift_xyz);
    end

    % calculate correlation and shifts
    %% plot metrics if required
    if plot_metrics
        tempNames = strsplit(fileNames{i},'/');
        fig_name =  strcat(fullfile(save_dir,tempNames{end}), '-correlation.fig');
        [correlation, cripness] = ZZ_plotpw_normcore(Y_padded, mc_temp, 10, fig_name);
        csvwrite(strcat(fullfile(save_dir,tempNames{end}), '-correlation.csv'),correlation);
        csvwrite(strcat(fullfile(save_dir,tempNames{end}), '-cripness.csv'),cripness);
    end
    %% choose the volumes with top 50% correlation to make a new template
    correlation_after = correlation(:,3);
    [sorted, sorted_idx] = sort(correlation_after, 'descend');
    top50_idx = sorted_idx(1:(int16(length(sorted_idx)/2)));
    mc_top50 = mc(:,:,:,:,top50_idx);
    mean_template = mean(mc_top50, 5);
    mean_template = reshape(mean_template, volume_info.xRes, volume_info.yRes, []);
    saveastiff(single(mean_template), strcat(fullfile(save_dir,tempNames{end}), '-template2.tif'), options);
end
dummy = 1;
end
