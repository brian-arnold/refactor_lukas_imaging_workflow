# a set of functions that are used to calculate calcium response
# Zhilei Zhao

import tifffile as tif
import numpy as np
import cv2
import scipy
import scipy.ndimage
from collections import OrderedDict
from ZZ_tiff_file import *
import glob
import os
import nrrd
from scipy import ndimage
from skimage import morphology


def ZZ_calculate_dff(volume, background_range=[10, 20], time_range=[], offset=0, smooth='', smooth_args={}):
    """calculate df/f for xyzt movie
    flatten negative values to 0
    smooth the results if specified based on smooth_args
    support smooth method: Gaussian, median, mean
    return both dff and background volume"""
    # don't use numpy broadcasting for all 4 dim, since memory may explode
    # substract offset, default offset is zero, since offset of PMT has been tuned
    # to have negative value when laser shutter is off
    v_o = volume - offset
    # add a small value to every pixel to avoid divide by 0
    # since the raw pixel values are int16, adding 0.5 will remove all zeros
    v_o = v_o.astype('float64')
    v_o += 0.5
    # select the volumes over time_range
    # if time_range is not specified, calculate for all volumes
    if not time_range:
        time_range = [0, volume.shape[3]]
    v_t = v_o[:, :, :, range(time_range[0], time_range[1])]
    # calculate a background volume by averaging
    v_b = v_o[:, :, :, range(background_range[0], background_range[1])]
    v_b = np.mean(v_b, axis=3)
    v_dff = np.zeros(v_t.shape)
    v_b_abs = np.abs(v_b)
    # loop through each time point
    for t in range(v_t.shape[3]):
        # substract the background volume from each volume
        v_dff[:, :, :, t] = v_t[:, :, :, t] - v_b
        # devide by the background volume, use absolute value, since some pixels
        # may have small negative values due to PMT fluctuation
        v_dff[:, :, :, t] /= v_b_abs
    # do 2D xy image smoothing if specificied
    if smooth:
        v_dff_s = np.zeros(v_dff.shape)
        if smooth == 'gaussian':
            func_smooth = ZZ_gaussian_smooth
        elif smooth == 'median':
            func_smooth = ZZ_median_smooth
        elif smooth == 'mean':
            func_smooth = ZZ_mean_smooth
        else:
            print('specified smooth method is not supported! return the original dff')
            return v_dff
        # loop through each time point
        for t in range(v_dff.shape[3]):
            vol = v_dff[:, :, :, t]
            vol_list = [vol[:, :, i] for i in range(vol.shape[2])]
            vol_s = map(lambda img: func_smooth(img, **smooth_args), vol_list)
            v_dff_s[:, :, :, t] = np.dstack(vol_s)
        return v_dff_s, v_b

    else:
        return v_dff, v_b


def ZZ_gaussian_smooth(img, size=(3, 3), sigmaX=3, sigmaY=3):
    """warp the cv2 gaussian smooth method for 2D image"""
    img_s = cv2.GaussianBlur(img, size, sigmaX, sigmaY)
    return img_s


def ZZ_median_smooth(img, size=3):
    """warp the cv2 median smooth method for 2D image"""
    img_s = cv2.medianBlur(img.astype('single'), size)
    return img_s.astype('float32')


def ZZ_mean_smooth(img, size=(3, 3)):
    """warp the cv2 mean smooth method for 2D image"""
    img_s = cv2.blur(img, size)
    return img_s


def ZZ_movie_blur(movie, gaussian_sigma=0, median_size=3):
    """do smoothing on the time-lapse movie
    use median smooth on XY, use Gaussian smooth on time if specified
    meddian smooth will remove shot noise"""
    # do median filtering on XY for each time point
    after_median = np.zeros(movie.shape)
    for t in range(movie.shape[3]):
        vol = movie[:, :, :, t]
        vol_list = [vol[:, :, i] for i in range(vol.shape[2])]
        vol_s = map(lambda img: ZZ_median_smooth(img, median_size), vol_list)
        after_median[:, :, :, t] = np.dstack(vol_s)
    # do gaussian filtering on time if specified
    if gaussian_sigma:
        after_gaussian = scipy.ndimage.filters.gaussian_filter1d(
            after_median, gaussian_sigma, axis=3)
        return after_gaussian
    else:
        return after_median


def ZZ_movie_threshold(movie, threshold):
    """hard thresholding to remove artifacts of the movie
    assign all elements that <= threshold to zero, while keep
    other elements unchanged"""
    movie_thre = movie
    movie_thre[movie_thre <= threshold] = 0
    return movie_thre


def ZZ_parse_odor_string(odor_pattern_string, file_range, channels=[]):
    """given odor pattern string and file num range, return file num that corresponds
    to channels"""
    odor_channels = [chr(c) for c in list(range(ord('A'), ord('Y') + 1)) +
                     list(range(ord('a'), ord('y') + 1))]
    string_info = [a[0] for a in odor_pattern_string.split('_') if a[0] in odor_channels]
    odor_file = OrderedDict()
    for c in channels:
        idx = [file_range[i] for i, x in enumerate(string_info) if x == c]
#         idx = [file_range[i] for i, x in enumerate(string_info) if c in x]
        odor_file[c] = idx
    return odor_file


def ZZ_calculate_peak_dff(file_movie, movie_info, channel, range_baseline, range_peak, threshold, filter_size_median, filter_size_gaussian):
    """calculate dff for peak volumes, ignore dynamics"""
    # what's most informative in calcium imaging studies is the peak df/f
    # this could be positive or negative, i.e. excitation or inhibition
    # this ignores the dynamics, also ignores off-response and tonic response
    movie = ZZ_read_reshape_tiff(file_movie, movie_info)
    movie_channel = movie[:, :, channel, :, :]
    movie_channel = np.squeeze(movie_channel)
    # calculate the baseline and peak volume, take average over the range
    # this already acts as a mean filter on the time axis
    volume_baseline = np.mean(movie_channel[:, :, :, range_baseline], axis=3)
    volume_peak = np.mean(movie_channel[:, :, :, range_peak], axis=3)
    # add a 3D median filter to remove salt-pepper noise
    volume_baseline = scipy.ndimage.median_filter(volume_baseline, filter_size_median)
    volume_peak = scipy.ndimage.median_filter(volume_peak, filter_size_median)
    # calculate df: change in fluorescence
    volume_df = volume_peak - volume_baseline
    # do thresholding based on the baseline volume, only care voxels that pass the threshold
    idx_good = volume_baseline > threshold
    # pre-locate the dff volume, set voxels that failed the threshold as dff=0
    threshold_dff = np.zeros(volume_baseline.shape)
    # calculate dff for elements that pass the threshold
    threshold_dff[idx_good] = volume_df[idx_good] / volume_baseline[idx_good]
    # smooth the final dff
    # final_dff = scipy.ndimage.median_filter(threshold_dff, filter_size)
    final_dff = scipy.ndimage.gaussian_filter(threshold_dff, filter_size_gaussian)
    return final_dff


def ZZ_plot_peak_dff_rgb(volume, value_range, color_range, filename=''):
    """map peak dff values to a diverging color map
    value_range = [min, middle, max]
    color_range = [color_min, color_middle, color_max], where color_X is a tuple
    with values from 0 to 255, corresponds to RGB"""
    # calculate the slopes and intercepts, separately for lower and upper half
    slope_upper = np.array([(color_range[1][i] - color_range[2][i]) /
                            (value_range[1] - value_range[2]) for i in range(3)])
    intercept_upper = np.array([(value_range[1] * color_range[2][i] - value_range[2]
                                 * color_range[1][i]) / (value_range[1] - value_range[2]) for i in range(3)])
    slope_lower = np.array([(color_range[0][i] - color_range[1][i]) /
                            (value_range[0] - value_range[1]) for i in range(3)])
    intercept_lower = np.array([(value_range[0] * color_range[1][i] - value_range[1]
                                 * color_range[0][i]) / (value_range[0] - value_range[1]) for i in range(3)])
    # define a function that maps value to a RGB color tuple

    def ZZ_numpy_to_rgb_diverge_map(x):
        if x >= value_range[1]:
            # belongs to the higher half
            if x >= value_range[2]:
                return color_range[2]
            else:
                res = (slope_upper * x + intercept_upper).astype('uint8')
        else:
            if x <= value_range[0]:
                return color_range[0]
            else:
                res = (slope_lower * x + intercept_lower).astype('uint8')
        return tuple(res[:])
    # vectorize this function to speed up
    func = np.vectorize(ZZ_numpy_to_rgb_diverge_map)
    volume_rgb = func(volume)
    volume_rgb_int8 = np.asarray(volume_rgb)
    # save as RGB tiff if specified
    if filename:
        volume_rgb_int8 = volume_rgb_int8.swapaxes(0, -1)
        volume_rgb_int8 = np.fliplr(np.rot90(volume_rgb_int8, axes=(1, 2)))
        tif.imsave(filename, volume_rgb_int8, photometric='rgb')
    return volume_rgb_int8


def ZZ_combine_dff_rgb(odor_channel, odor_pattern_strings, file_ranges, brain_names, folder_save_dff, tiff_type, folder_save_combined=''):
    """combine dff from different brains into one volume tiff
    so response could be compared, odor_pattern_string, file_ranges and brain_names are lists
    with the same number of elements"""
    # get the odor session file names
    odor_files = [ZZ_parse_odor_string(odor_pattern_strings[i], file_ranges[i], odor_channel)[
        odor_channel] for i in range(len(file_ranges))]
    # go through each session file for each brain, read the tiff
    dff_tiffs = []
    for i in range(len(odor_files)):
        temp = []
        for j in range(len(odor_files[i])):
            file_base_glob = f'{brain_names[i]}f{odor_files[i][j]:05d}*{tiff_type}'
            file_name = glob.glob(folder_save_dff + file_base_glob)
            try:
                fn = file_name[0]
            except:
                temp.append('NA')
            else:
                dff_volume = tif.imread(fn)
                temp.append(dff_volume)
        # combine the dff from one brain
        temp2 = np.concatenate(temp, axis=1)
        dff_tiffs.append(temp2)
    # combine dff from different brains
    dff_combined = np.concatenate(dff_tiffs, axis=2)
    # save if specified
    if folder_save_combined:
        tif.imsave(folder_save_combined + 'dff_combined_channel' +
                   odor_channel + '-' + tiff_type, dff_combined, photometric='rgb')
    return dff_combined


def ZZ_combine_dff_rgb2(odor_channel, odor_pattern_string, channel_odorant, file_range, brain_name, folder_save_dff, tiff_type, folder_save_combined=''):
    """combine data for one channel from one brain"""
    odor_files = ZZ_parse_odor_string(odor_pattern_string, file_range, odor_channel)[odor_channel]
    # go through each session file for each brain, read the tiff
    dff_tiffs = []
    for i in range(len(odor_files)):
        file_base_glob = f'{brain_name}_{odor_files[i]:05d}*{tiff_type}'
        file_name = glob.glob(os.path.join(folder_save_dff, file_base_glob))
        try:
            fn = file_name[0]
        except:
            dff_tiffs.append('NA')
        else:
            dff_volume = tif.imread(fn)
            dff_tiffs.append(dff_volume)
        # combine the dff from one brain
    dff_combined = np.concatenate(dff_tiffs, axis=2)
    # save if specified
    if folder_save_combined:
        tif.imsave(os.path.join(folder_save_combined,
                                f'{brain_name}_{channel_odorant[odor_channel]}_dff_rgb_combined.tif'), dff_combined, photometric='rgb')
    return dff_combined


def ZZ_combine_dff_rgb3(odor_channel, odor_pattern_string, channel_odorant, file_range, brain_name, folder_save_dff, tiff_type, folder_save_combined=''):
    """combine data for one channel from one brain, to make it work with channel name like 'A1' """
    odor_files = ZZ_parse_odor_stringNew(odor_pattern_string, file_range, [odor_channel])
    odor_files = odor_files[odor_channel]
    # go through each session file for each brain, read the tiff
    dff_tiffs = []
    for i in range(len(odor_files)):
        file_base_glob = f'{brain_name}_{odor_files[i]:05d}*{tiff_type}'
        file_name = glob.glob(os.path.join(folder_save_dff, file_base_glob))
        try:
            fn = file_name[0]
        except:
            dff_tiffs.append('NA')
        else:
            dff_volume = tif.imread(fn)
            dff_tiffs.append(dff_volume)
        # combine the dff from one brain
    dff_combined = np.concatenate(dff_tiffs, axis=2)
    # save if specified
    if folder_save_combined:
        tif.imsave(os.path.join(folder_save_combined,
                                f'{brain_name}_{channel_odorant[odor_channel]}_dff_rgb_combined.tif'), dff_combined, photometric='rgb')
    return dff_combined


def ZZ_combine_dff_raw(odor_channel, odor_pattern_string, channel_odorant, file_range, brain_name, folder_save_dff, tiff_type, folder_save_combined=''):
    """combine data for one channel from one brain"""
    odor_files = ZZ_parse_odor_string(odor_pattern_string, file_range, odor_channel)[odor_channel]
    # go through each session file for each brain, read the tiff
    dff_tiffs = []
    for i in range(len(odor_files)):
        file_base_glob = f'{brain_name}_{odor_files[i]:05d}*{tiff_type}'
        file_name = glob.glob(os.path.join(folder_save_dff, file_base_glob))
        try:
            fn = file_name[0]
        except:
            dff_tiffs.append('NA')
        else:
            dff_volume = tif.imread(fn)
            dff_tiffs.append(dff_volume)
        # combine the dff from one brain
    dff_combined = np.concatenate(dff_tiffs, axis=2)
    # save if specified
    if folder_save_combined:
        tif.imsave(os.path.join(folder_save_combined,
                                f'{brain_name}_{channel_odorant[odor_channel]}_{tiff_type}_combined.tif'), dff_combined)
    return dff_combined


def ZZ_combine_dff_raw3(odor_channel, odor_pattern_string, channel_odorant, file_range, brain_name, folder_save_dff, tiff_type, folder_save_combined=''):
    """combine data for one channel from one brain, to make it work with channel name like 'A1' """
    odor_files = ZZ_parse_odor_stringNew(odor_pattern_string, file_range, [odor_channel])
    odor_files = odor_files[odor_channel]
    # go through each session file for each brain, read the tiff
    dff_tiffs = []
    for i in range(len(odor_files)):
        file_base_glob = f'{brain_name}_{odor_files[i]:05d}*{tiff_type}'
        file_name = glob.glob(os.path.join(folder_save_dff, file_base_glob))
        try:
            fn = file_name[0]
        except:
            dff_tiffs.append('NA')
        else:
            dff_volume = tif.imread(fn)
            dff_tiffs.append(dff_volume)
        # combine the dff from one brain
    dff_combined = np.concatenate(dff_tiffs, axis=2)
    # save if specified
    if folder_save_combined:
        tif.imsave(os.path.join(folder_save_combined,
                                f'{brain_name}_{channel_odorant[odor_channel]}_{tiff_type}_combined.tif'), dff_combined)
    return dff_combined


def ZZ_merge_across_dataset(folder_home, brains_merge, name_dataset_array, folder_dff, merge_subfolder, tif_types, odorants, folder_save, tif_format):
    """input:
    brains_merge: what brain names to merge, format as as the final merged layout
    name_dataset_array: what dataset are those brains located
    folder_dff: the subfolder name under 'Response' that store dff for different datasets
    merge_subfolder: what subfolder in 'RawPeakDff/dataset' to merge
    tif_types: what tif type to merge (use for glob function)
    odorants: what odorants to merge
    folder_save: what folder to save merged tif"""
    # determine the number of rows and columns for merged tif
    num_row = len(brains_merge)
    num_column = max([len(a) for a in brains_merge])
    # loop through each tif type to be merged
    for tif_type in tif_types:
        # loop through each odorant
        for odorant in odorants:
            dff_volume_merged_lists = []
            # loop through the brains
            for ii in range(len(brains_merge)):
                for jj in range(len(brains_merge[ii])):
                    # merge the rgb volumes
                    # read the dff volume
                    if tif_format == 'tif':  # this is unregistered volume
                        fname_dff_brain = glob.glob(os.path.join(
                            folder_home, 'Response', folder_dff, name_dataset_array[ii][jj], merge_subfolder, f'{brains_merge[ii][jj]}*{odorant}*{tif_type}.{tif_format}'))[0]
                        dff_volume_brain = tif.imread(fname_dff_brain)
                    if 'nrrd' in tif_format:  # this is registered volume
                        #                         print(os.path.join(folder_home,'Response', folder_dff,name_dataset_array[ii][jj], brains_merge[ii][jj], merge_subfolder, f'{brains_merge[ii][jj]}*{odorant}*{tif_type}.{tif_format}'))
                        fname_dff_brain = glob.glob(os.path.join(folder_home, 'Response', folder_dff, name_dataset_array[ii][
                                                    jj], brains_merge[ii][jj], merge_subfolder, f'{brains_merge[ii][jj]}*{odorant}*{tif_type}.{tif_format}'))[0]
                        dff_volume_brain, header_info = nrrd.read(fname_dff_brain)
                        dff_volume_brain = np.swapaxes(dff_volume_brain, -1, 0)
                    dff_volume_merged_lists.append(dff_volume_brain)
            shape = dff_volume_merged_lists[0].shape
            if len(shape) == 4:  # rgb volume
                dff_volume_merged = np.zeros(
                    (shape[0], shape[1] * num_row, shape[2] * num_column, shape[3]), dtype='int8')
                # if the layout is not symmetrical, the empty block will be white
                dff_volume_merged.fill(255)
            else:
                dff_volume_merged = np.zeros(
                    (shape[0], shape[1] * num_row, shape[2] * num_column), dtype='float32')
            brain_count = 0
            # go through each brain, fill out the merged dff volume
            for ii in range(len(brains_merge)):
                for jj in range(len(brains_merge[ii])):
                    if len(shape) == 4:  # rgb volume
                        dff_volume_merged[:, ii * shape[1]:(ii + 1) * shape[1], jj * shape[2]:(
                            jj + 1) * shape[2], :] = dff_volume_merged_lists[brain_count]
                    else:
                        dff_volume_merged[:, ii * shape[1]:(ii + 1) * shape[1], jj * shape[2]:(
                            jj + 1) * shape[2]] = dff_volume_merged_lists[brain_count]
                    brain_count += 1
            fname_save = os.path.join(folder_save, f'{odorant}_{tif_type}.tif')
            tif.imsave(fname_save, dff_volume_merged)


def ZZ_segmentVOI_iterative(dff, segment_quantile=0.9995, glomeruli_size=10, dilation_iter=15, stop_threshold=1, save_fname=''):
    """Segment VOI out from responding glomeruli by iterative thresholding the dff volume
    Inputs:
        dff: the response volume
        segment_quantile: threshold to segment glomeruli
        glomeruli_size: filter out regions that are smaller than this size (num of voxels)
        dilation_iter: num of iterations to dilate the glomeruli center to mask in the dff volume
        stop_threshold: terminate the iteration when current region has dff smaller than this value
        save_fname: if not empty, save intermediate volume as tif
    Outputs:
        ndarray with the same size as dff, where glomerular centers are labeled"""
    # a volume to record the labels
    d_label = np.zeros(dff.shape, dtype='uint8')
    # a volume that keeps getting updated
    d_now = dff
    # how many regions have been labeled
    label_count = 0
    idx = 0
    while(1):
        thre = np.quantile(d_now, segment_quantile)
        # break the loop if no more regions passed the stop threshold
        if thre < stop_threshold:
            break
        # get regions that pass the quantile threshold
        this_pass = np.zeros(d_now.shape, dtype="float32")
        this_pass[d_now > thre] = d_now[d_now > thre]
        # dilate the labeled region then mask on the dff volume
        this_pass_bool = d_now > thre
        this_pass_bool_dilate = ndimage.binary_dilation(this_pass_bool, iterations=dilation_iter)
        this_pass_mask = d_now.copy()
        this_pass_mask[this_pass_bool_dilate] = 0
        # get labels for the thresholded region
        this_label, num_label = morphology.label(this_pass_bool, return_num=True)
        # sort the labels based on mean dff value
        region_mean = [np.mean(dff[this_label == v]) for v in range(1, num_label + 1)]
        # sort the region mean
        labels_old = [v for v in range(1, num_label + 1)]
        labels_new = np.flip(np.argsort(region_mean)) + 1
        # re-assign the labels
        this_label_sorted = np.zeros(this_label.shape, dtype='uint8')
        for jdx, label_old in enumerate(labels_old):
            label_new = labels_new[jdx]
            this_label_sorted[this_label == label_old] = label_new
        # filter out regions that are too small
        label_count_failed = 0
        for v in range(1, num_label + 1):
            label_size = np.sum(this_label_sorted == v)
            if label_size >= glomeruli_size:
                d_label[this_label_sorted == v] = v + label_count + label_count_failed
            else:
                label_count_failed += -1
        # update current data volume after masking
        d_now = this_pass_mask
        label_count += num_label + label_count_failed
        # save intermediate volumes if specified
        if save_fname:
            ZZ_save_as_tiff(f'{save_fname}.pass{idx}.tif', this_pass)
            ZZ_save_as_tiff(f'{save_fname}.dilate{idx}.tif', this_pass_bool_dilate)
            ZZ_save_as_tiff(f'{save_fname}.mask{idx}.tif', this_pass_mask)
            ZZ_save_as_tiff(f'{save_fname}.label{idx}.tif', this_label_sorted)
    return d_label


def ZZ_calculate_area(file_movie, movie_info, channel, baseline_range, area_range, threshold, filter_size_gaussian):
    """calculate area under df/f curve as a measure of total neural activity"""
    vol = ZZ_read_reshape_tiff(file_movie, movie_info)
    vol = np.squeeze(vol[:,:,channel,:,:])
    # make a running average movie to smooth out the temporal element
    # def runningAve(x):
    #     res = np.convolve(x, np.ones((5,))/5, mode='valid')
    #     return res
    # ave = np.apply_along_axis(runningAve, axis=3, arr=vol)
    # smooth out the temporal element
    def smoothTemporal(x):
        res = ndimage.gaussian_filter1d(x, sigma=3)
        return res
    ave = np.apply_along_axis(smoothTemporal, axis=3, arr=vol)
    # calculate a baseline volume
    baseline = np.mean(ave[:,:,:,baseline_range], axis=3)
    # remove non-AL pixels by thresholding
    idx_good = baseline > threshold
    # calculate df/f for each time point
    dff_all = np.zeros(ave.shape)
    # reshape baseline volume for broadcasting
    baseline = baseline[idx_good].reshape(-1,1)
    dff_all[idx_good,:] = ave[idx_good, :] / baseline - 1
    # calculate area
    area = np.trapz(dff_all[:,:,:,area_range], axis=3)
    # smooth the volume with Gaussian filter to remove noise
    filter_size_gaussian = [4,4,2]
    area = scipy.ndimage.gaussian_filter(area, filter_size_gaussian)
    return area


def ZZ_parse_odor_stringNew(odor_string, file_ranges, channels):
    """a new function to parse the odor string, when channel is letter+number, e.g. (A0)"""
    info = odor_string.split('_')
    odor_file = OrderedDict()
    for channel in channels:
        file_idx = []
        for idx, ch in enumerate(info):
            if ch==f'({channel})':
                file_idx.append(file_ranges[idx])
        odor_file[channel] = file_idx
    return(odor_file)


def ZZ_combine_dff_rawNew(odor_channel, odor_pattern_string, channel_odorant, file_range, brain_name, folder_save_dff, tiff_type, folder_save_combined=''):
    """combine data for one channel from one brain, modified on ZZ_combine_dff_raw, to coupe with new channel representation, e.g. (A0)"""
    odor_files = ZZ_parse_odor_stringNew(odor_pattern_string, file_range, [odor_channel])[odor_channel]
    # go through each session file for each brain, read the tiff
    dff_tiffs = []
    for i in range(len(odor_files)):
        file_base_glob = f'{brain_name}_{odor_files[i]:05d}*{tiff_type}'
        file_name = glob.glob(os.path.join(folder_save_dff, file_base_glob))
        try:
            fn = file_name[0]
        except:
            dff_tiffs.append('NA')
        else:
            dff_volume = tif.imread(fn)
            dff_tiffs.append(dff_volume)
        # combine the dff from one brain
    dff_combined = np.concatenate(dff_tiffs, axis=2)
    # save if specified
    if folder_save_combined:
        tif.imsave(os.path.join(folder_save_combined,
                                f'{brain_name}_{channel_odorant[odor_channel]}_{tiff_type}_combined.tif'), dff_combined)
    return dff_combined


def ZZ_time_trace_glomerulus(fn, center, onset_point = 26, radius=7, dim_x=128, dim_y=128, dim_z=24, subtract_background=True, baseline_range=[5,20]):
    """get the time trace given the file name and the center of a glomerulus
    center should be formated [z, x, y], values as indicated in Fiji"""
    # define the averaging region according to the circle center
    x = np.arange(0, dim_x)
    y = np.arange(0, dim_y)
    xy_idx = (x[np.newaxis,:]-center[1])**2 + (y[:,np.newaxis]-center[2])**2 <= radius**2
    vol = tif.imread(fn)
    vol = np.reshape(vol, (-1, dim_z, dim_x, dim_y))
    timeseries = vol[:, center[0], :, :]
    # loop through each timepoint, calculate the average area of the glomerulus
    ave_all = []
    for t in range(timeseries.shape[0]):
        stack = np.squeeze(timeseries[t, :, :])
        ave = np.mean(stack[xy_idx])
        ave_all.append(ave)
    # subtract background if specified
    if subtract_background:
        background = np.mean(ave_all[baseline_range[0]:baseline_range[1]])
        ave_all = [(v-background)/background for v in ave_all]
    return ave_all


def ZZ_time_trace_glomerulus3D(fn, center, onset_point = 26, radius=7, dim_x=128, dim_y=128, dim_z=24, subtract_background=True, baseline_range=[5,20]):
    """get the time trace given the file name and the center of a glomerulus
    center should be formated [z, x, y], values as indicated in Fiji
    response is averaged over 3 z-stacks: centered around z
    note z should be larger than 0, and smaller than the dim_z"""
    # define the averaging region according to the circle center
    # for the upper and lower z-stack, take a smaller circle
    radius_diff = 2
    x = np.arange(0, dim_x)
    y = np.arange(0, dim_y)
    # center stack
    xy_idx = (x[np.newaxis,:]-center[1])**2 + (y[:,np.newaxis]-center[2])**2 <= radius**2
    # upper and lower stack
    xy_idx_upper = (x[np.newaxis,:]-center[1])**2 + (y[:,np.newaxis]-center[2])**2 <= (radius-radius_diff)**2
    # read the tif movie
    vol = tif.imread(fn)
    vol = np.reshape(vol, (-1, dim_z, dim_x, dim_y))
    # define a local function to calculate averaged value for each time point
    def calculateAverageTimePoint(movie, zstack, xy_range):
        timeseries = movie[:, zstack, :, :]
        # loop through each timepoint, calculate the average area of the glomerulus
        ave_all = []
        for t in range(timeseries.shape[0]):
            stack = np.squeeze(timeseries[t, :, :])
            ave = np.mean(stack[xy_range])
            ave_all.append(ave)
        return(ave_all)
    # calculate the averaged fluroscence value
    ave_all_center = calculateAverageTimePoint(vol, center[0], xy_idx)
    ave_all_upper = calculateAverageTimePoint(vol, center[0]-1, xy_idx_upper)
    ave_all_lower = calculateAverageTimePoint(vol, center[0]+1, xy_idx_upper)
    # average over these three stacks
    res = np.stack([ave_all_upper, ave_all_center, ave_all_lower])
    res = np.mean(res, axis=0)
    # subtract background and divide background if specified; return is actually df/f
    if subtract_background:
        background = np.mean(res[baseline_range[0]:baseline_range[1]])
        res = [(v-background)/background for v in res]
    return res


def ZZ_combine_time_traces(fns, center, onset_point, radius=7, dim_x=128, dim_y=128, dim_z=24, imaging_freq=3.76, baseline_range=[5,20]):
    """combine time traces for different puffs of the same odorant
    center is the center of glomerulus, in [z, x, y] order
    onset_point is at what volume did the valve open
    imaging_freq is the volumetric imaging rate"""
    # get the time trace from each file
    glo_ave_all = []
    for fn in fns:
        glo_ave =  ZZ_time_trace_glomerulus3D(fn, center, onset_point, radius, dim_x, dim_y, dim_z, True, baseline_range)
        glo_ave_all.append(glo_ave)
    # note the timepoints may vary by 1 due to trigger signal difference, so take the smallest
    min_time = np.min([len(a) for a in glo_ave_all])
    glo_ave_all_array = [a[0:min_time] for a in glo_ave_all]
    # smooth the trace of each puff with gaussian before averaging
    glo_ave_all_array_smooth = [ndimage.gaussian_filter1d(a, sigma=2) for a in glo_ave_all_array]
    glo_ave_all_array = np.stack(glo_ave_all_array)
    glo_ave_all_array_smooth = np.stack(glo_ave_all_array_smooth)
    # take the average of the puffs for the same odorant
    glo_ave_all_ave = np.mean(glo_ave_all_array_smooth, axis=0)
    # append the averaged trace to end of the array
    res = np.append(glo_ave_all_array, glo_ave_all_ave.reshape(1,-1), axis=0)
    # append the real time to the end of the array
    real_time = (np.arange(res.shape[1]) - onset_point) / imaging_freq
    res = np.append(res, real_time.reshape(1, -1), axis=0)
    return res


def ZZ_align_traces_by_time(trace_value, time_value, imaging_freq=3.76):
    """combine traces by align them to time=0
    both trace_value and time_value are lists of the same length"""
    # find the index for time=0
    zero_idx = []
    for t in time_value:
        zero_idx.append(np.argmin(np.abs(t)))
    # calculate the max index for left and right side
    left_long = np.max(zero_idx)
    right_side = [time_value[ii].shape[0] - zero_idx[ii] for ii in range(len(trace_value))]
    right_long = np.max(right_side)
    # put the trace values into a new numpy array, according to the new index
    # the last row being the average
    comb = np.zeros((len(trace_value)+2, left_long+right_long))
    # calculate offset of index for each trace
    offsets = [left_long - zero_idx[ii] for ii in range(len(trace_value))]
    # fill the values in
    for ii in range(len(trace_value)):
        comb[ii, offsets[ii]:(offsets[ii]+time_value[ii].shape[0])] = trace_value[ii]
    # calculate the average
    ave = np.mean(comb[0:len(trace_value), :], axis=0)
    comb[-2, :] = ave
    # append the real time
    comb[-1, :] = (np.arange(comb.shape[1]) - left_long)/imaging_freq
    return comb


def ZZ_calculate_area_by_peak(value, area_idx, drop_level = 0.95, plot_out=False):
    """given a time trace (value) and where to search peak (area_idx)
    find the peak, define the peak, and calculate area under curve"""
    # first find the peaks
    search_region = area_idx[0]
    value_search = value[area_idx]
    idx_max = np.argmax(value_search) + search_region[0]
    # then check when the peak height drops to the specified value
    baseline = np.mean(value[0:search_region[0]])
    peak_value = value[idx_max]
    peak_height = peak_value - baseline
    boundary_value = peak_value - drop_level*peak_height
    # find the left boundary
    for left in range(idx_max, 0, -1):
        if value[left] < boundary_value:
            break
    # find the right boundary
    for right in range(idx_max, len(value), 1):
        if value[right] < boundary_value:
            break
    # calculate the area
    area = np.trapz(value[left:right], axis=0)
    # plot to check
    if plot_out:
        fig, ax = plt.subplots()
        ax.plot(value)
        ax.plot(idx_max, value[idx_max], 'ro')
        ax.plot(left, value[left], 'bo')
        ax.plot(right, value[right], 'bo')
    # return area and peak
    return area, value[idx_max]


def ZZ_calculate_area_peakDrop(fluo, real_time, interval=[0,60], drop_level=0.95, plot_out=False):
    """given a calcium trace (fluo and real_time) and where to search for peak (interval)
    find the peak, define the peak by extending from highest point until noise level, and calculate area under curve
    unit for fluo is df/f
    unit for real_time and interval is seconds"""
    # first find the location of the peak
    left = np.where(real_time>=interval[0])
    right = np.where(real_time<=interval[1])
    area_idx = np.intersect1d(left, right)
    fluo_search = fluo[area_idx]
    # the peak could be excitation or inhibition peak
    # find the point with maximal absolute value
    idx_max_abs = np.argmax(np.abs(fluo_search))
    idx_max = idx_max_abs + left[0][0]
    # check if the peak is excitation or inhibition
    peak_value = fluo[idx_max]
    boundary_value = peak_value - drop_level*peak_value
    # define the peaks by extending from the max point until drops to the noise level
    if peak_value > 0:
        peak_type = 'excitation'
        # find the left boundary
        for left_bound in range(idx_max, 0, -1):
            if fluo[left_bound] < boundary_value:
                break
        # find the right boundary
        for right_bound in range(idx_max, len(fluo), 1):
            if fluo[right_bound] < boundary_value:
                break
    else:
        peak_type = 'inhibition'
        # find the left boundary
        for left_bound in range(idx_max, 0, -1):
            if fluo[left_bound] > boundary_value:
                break
        # find the right boundary
        for right_bound in range(idx_max, len(fluo), 1):
            if fluo[right_bound] > boundary_value:
                break
    # calculate area by integrating under curve
    area = np.trapz(fluo[left_bound:right_bound], real_time[left_bound:right_bound], axis=0)
    # plot to check
    if plot_out:
        fig, ax = plt.subplots()
        ax.plot(real_time, fluo)
        ax.plot(real_time[area_idx], fluo[area_idx], color='red')
        ax.plot(real_time[idx_max], fluo[idx_max], marker='o', color='blue')
        ax.plot(real_time[left_bound], fluo[left_bound], marker='o', color='green')
        ax.plot(real_time[right_bound], fluo[right_bound], marker='o', color='green')
    return(area, peak_value)


def ZZ_quantify_dff_area_brain(data, selected_odor_area_single, selected_odor_area_markes, area_range_single, area_range_markes, normalization_factor=1, plot_out=False):
    """quantify area under curve for df/f time trace of selected odors
    data is a dict that has odor as keys and calcium traces as values
    calcium trace is a np array that has averaged trace in the 2nd last row and
    real time in the last row"""
    temp_area = {}
    temp_peak = {}
    for odor in selected_odor_area_single+selected_odor_area_markes:
        trace = data[odor]
        # get the averaged response
        fluo = trace[-2, :]
        # get the real time
        real_time = trace[-1, :]
        # normalize the raw response
        fluo = fluo / normalization_factor
        # decide what interval to find peak; different between single odorant and markes
        if odor in selected_odor_area_single:
            interval = area_range_single
        else:
            interval = area_range_markes
        # calculate area under curve
        area, peak = ZZ_calculate_area_peakDrop(fluo, real_time, interval, plot_out=plot_out)
        temp_area[odor] = area
        temp_peak[odor] = peak
    return(temp_area, temp_peak)
