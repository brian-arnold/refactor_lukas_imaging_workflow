# a set of functions that are used to calculate calcium response
# Zhilei Zhao

import tifffile as tif
import numpy as np
import cv2
import scipy
from collections import OrderedDict
from ZZ_tiff_file import *
import glob
import os
import nrrd
from scipy import ndimage
from skimage import morphology
from scipy import stats
import pandas as pd



def ZZ_parse_odor_stringNew(odor_string, file_ranges, channels):
    """parse the odor string with the new format, where different odorants may under the same channel
    odor_string: String, the order of puffs (single odorants or Markes), e.g. (A0)_(B0)_(B0)_(A1)
    file_ranges: List, what files corresponds to the odor_string, e.g. [2,3,5,6]
    channels: List, what channels to return results, e.g. ['B0', 'A1']
    return: OrderedDict, e.g. {'B0':[3, 5], 'A1':[6]} """
    info = odor_string.split('_')
    odor_file = OrderedDict()
    for channel in channels:
        file_idx = []
        for idx, ch in enumerate(info):
            if ch==f'({channel})':
                file_idx.append(file_ranges[idx])
        odor_file[channel] = file_idx
    return(odor_file)



def ZZ_reverse_parse_odorString(odor_string, file_range, fidx):
    """reverse parse the odor pattern string, i.e. from file index to channel name
    odor_string: String
    file_range: List or range
    fidx: Int, the movie file index"""
    info = odor_string.split('_')
    file_range = list(file_range)
    ii = file_range.index(fidx)
    channel = info[ii].replace('(','').replace(')','')
    return(channel)



def ZZ_smooth_movies(movie, gaussian_spatial=[4,4,2], gaussian_temporal=3):
    """smooth odor-evoked movie with gaussian filter to remove noise
    movie: 4-d array, [x, y, z, t]
    gaussian_spatial: List, filter size on the spatial domain
    gaussian_temporal: Int, filter size on the temporal domain
    return: 4-d array, smoothed movie"""
    # smooth on the spatial domain
    smoothed1 = np.zeros(movie.shape)
    for tt in range(movie.shape[3]):
        v = np.squeeze(movie[:,:,:,tt])
        v = ndimage.gaussian_filter(v, gaussian_spatial)
        smoothed1[:,:,:,tt] = v
    # smooth on the temporal domain
    def smoothTemporal(x):
        res = ndimage.gaussian_filter1d(x, sigma=gaussian_temporal)
        return res
    smoothed2 = np.apply_along_axis(smoothTemporal, axis=3, arr=smoothed1)
    return(smoothed2)



def ZZ_calculate_area_peakDrop_once(fluo, real_time, search_idx, baseline_sd):
    """calculate area for calcium df/f trace using the peak-drop method
    fluo_search: 1-d array, df/f curve
    real_time_search: 1-d array, corresponding time points, unit is seconds
    search_idx: where to expect the peaks
    baseline_sd: float, standard deviation of the baseline, used to define boundaries
    return: List, [area, peak_idx, left_bound, right_bound]"""
    
    # define the region of searching 
    fluo_search = fluo[search_idx]
    real_time_search = real_time[search_idx]
    
    # the peak could be excitation or inhibition peak
    # find the point with maximal absolute value
    idx_max = np.argmax(np.abs(fluo_search))
    # convert to original scale
    peak_idx = idx_max + search_idx[0]
    # check if the peak is excitation or inhibition
    peak_value = fluo[peak_idx]

    ## define the entire peak by extending from the highest/lowest point
    # define the peaks by extending from the max point until drops to the noise level (1.96*sd)
    if peak_value >= 0:
        peak_type = 'excitation'
        # find the left boundary
        boundary_value = 1.96*baseline_sd
        for left_bound in range(peak_idx, -1, -1):
            if fluo[left_bound] < boundary_value:
                break
        # find the right boundary
        for right_bound in range(peak_idx, len(fluo), 1):
            if fluo[right_bound] < boundary_value:
                break
    else:
        peak_type = 'inhibition'
        # find the left boundary
        boundary_value = -1.96*baseline_sd
        for left_bound in range(peak_idx, -1, -1):
            if fluo[left_bound] > boundary_value:
                break
        # find the right boundary
        for right_bound in range(peak_idx, len(fluo), 1):
            if fluo[right_bound] > boundary_value:
                break

    # check if the fluroscence drops or not after the highest/lowest point
    # if doesn't drop at all, it's likely to be a fake peak
    # use linear fit to check if it drops or not
    x = np.arange(0, right_bound-peak_idx)
    yn = fluo[peak_idx:right_bound]
    flag_fake_peak = 0
    if len(x)<=1:
        flag_fake_peak = 1
    else:
        slope, intercept, r_value, p_value, std_err = stats.linregress(x,yn)
        if peak_type == 'excitation':
            if slope>=0 or p_value>=0.1:
                flag_fake_peak = 1
        if peak_type == 'inhibition':
            if slope<=0 or p_value>=0.1:
                flag_fake_peak = 1

    ## calculate area by integrating the defined peak
    # return 0 if consider peak as fake
    if flag_fake_peak==1:
        area = 0
    else:
        area = np.trapz(fluo[left_bound:right_bound], real_time[left_bound:right_bound], axis=0)
   
    return([area, peak_idx, left_bound, right_bound])


def ZZ_calculate_area_peakDrop(fluo, real_time, interval=[0,60], markes=False, plot_out=False, return_points=False):
    """given a calcium trace, calculate the area under curve by the peak drop method
    fluo: 1-d array, fluorescence trace, unit is df/f, baseline should be 0
    real_time: 1-d array, real time corresponds to each point in the fluo trace
    interval: List, where to search for the highest/lowest point in a peak, unit is seconds in real time
    markes: Logical, whether it's markes puff, if so check possibility of multiple peaks
    plot_out: Logical, whether to plot the trace and landmark points
    return: Float, total area of identified peaks, unit for area is df/f * second"""
    # first find the major peak; then mask it, run the same algorithm on its left or right to check for possible peaks
    # define the initial searching region
    left = np.where(real_time>=interval[0])
    right = np.where(real_time<=interval[1])
    search_idx = np.intersect1d(left, right)

    # calculate the sd of baseline 
    baseline_range = real_time<0
    baseline_sd = np.std(fluo[baseline_range])

    # search for the 1st peak
    [area1, peak_idx1, left_bound1, right_bound1] = ZZ_calculate_area_peakDrop_once(fluo, real_time, search_idx, baseline_sd)

    if markes:
        # for Markes, mask the 1st peak; 
        # run the same algorithm again for remaing search region on the left and right
        flag_left = 0
        flag_right = 0
        area2 = 0
        area3 = 0
        # left side
        if left_bound1>search_idx[0]:
            # left side is bounded by time zero and left_bound of 1st peak
            flag_left = 1
            fluo_left = fluo[search_idx[0]:left_bound1]
            real_time_left = real_time[search_idx[0]:left_bound1]
            search_idx_left = np.arange(len(fluo_left))
            [area2, peak_idx2, left_bound2, right_bound2] = ZZ_calculate_area_peakDrop_once(fluo_left, real_time_left, search_idx_left, baseline_sd)
            # convert the index into the original scale
            peak_idx2 += search_idx[0]
            left_bound2 += search_idx[0] 
            right_bound2 += search_idx[0]
        # right side
        if right_bound1<search_idx[-1]:
            # right left is bounded by right_bound of 1st peak and last time point
            flag_right = 1
            fluo_right = fluo[right_bound1:]
            real_time_right = real_time[right_bound1:]
            search_idx_right = np.arange(0, search_idx[-1]-right_bound1, 1)
            [area3, peak_idx3, left_bound3, right_bound3] = ZZ_calculate_area_peakDrop_once(fluo_right, real_time_right, search_idx_right, baseline_sd)
            # convert the index into the original scale
            peak_idx3 += right_bound1
            left_bound3 += right_bound1 
            right_bound3 += right_bound1

        if plot_out:
            fig, ax = plt.subplots()
            ax.plot(real_time, fluo)
            #     ax.plot(real_time_search[area_idx], fluo_search[area_idx], color='red')
            ax.plot(real_time[peak_idx1], fluo[peak_idx1], marker='o', color='blue')
            ax.plot(real_time[left_bound1], fluo[left_bound1], marker='o', color='green')
            ax.plot(real_time[right_bound1], fluo[right_bound1], marker='o', color='green')
            if flag_left:
                ax.plot(real_time[peak_idx2], fluo[peak_idx2], marker='+', color='blue')
                ax.plot(real_time[left_bound2], fluo[left_bound2], marker='+', color='green')
                ax.plot(real_time[right_bound2], fluo[right_bound2], marker='+', color='green')
            if flag_right:
                ax.plot(real_time[peak_idx3], fluo[peak_idx3], marker='x', color='blue')
                ax.plot(real_time[left_bound3], fluo[left_bound3], marker='x', color='green')
                ax.plot(real_time[right_bound3], fluo[right_bound3], marker='x', color='green')
#         print(area1, area2, area3)
        return(area1+area2+area3)
    else:
        # for single odorants, just return 1st peak
        if plot_out:
            fig, ax = plt.subplots()
            ax.plot(real_time, fluo)
            #     ax.plot(real_time_search[area_idx], fluo_search[area_idx], color='red')
            ax.plot(real_time[peak_idx1], fluo[peak_idx1], marker='o', color='blue')
            ax.plot(real_time[left_bound1], fluo[left_bound1], marker='o', color='green')
            ax.plot(real_time[right_bound1], fluo[right_bound1], marker='o', color='green')
        if return_points:
            return(area1, [peak_idx1, left_bound1, right_bound1])
        else:
            return(area1)



def ZZ_calculate_area_volume_peakDrop(fn, movie_info, markes=False, channel=0, imaging_freq=3.76, pre_puff=7, search_interval=[0,15], threshold=10, plot_out=False, remove_ends=5):
    """from calcium imaging movie, calculate an area volume by 'peak drop' method
    fn: String, file name of the movie, full path
    movie_info: Dict, xyz size, number of channels and etc, used to reshape the np array
    imaging_freq: Float, volumetric imaging rate, used to calculate real time
    pre_puff: Int, how many seconds before puffing were recorded, used to identify time zero
    search_interval: List, when to expect the highest point in a peak
    drop_level: Float, extend from the highest point until fluorescence drops to this level
    plot_out: Logical, whether to plot the fluorescence trace and identification of landmark points
    remove_ends: Int, how many volumes to ignore when calculate baseline
    return: 3d-array, area volume"""

    # parameters that are relatively fixed
    # filter size to smooth the data
    gaussian_spatial = [4,4,2]
    gaussian_temporal = 3

    ## read in the movie
    vol = ZZ_read_reshape_tiff(fn, movie_info)
    # only use the green channel
    movie = np.squeeze(vol[:,:,channel,:,:])

    ## smoothing to remove noise
    smoothed2 = ZZ_smooth_movies(movie, gaussian_spatial=gaussian_spatial, gaussian_temporal=gaussian_temporal)

    ## define a local function to call the peak drop method
    def ZZ_calculate_area_local(trace):
        # calculate df/f
        # at what time point valve open
        time_zero = int(imaging_freq * pre_puff)
        # calculate baseline fluorescence, excluding the 1st and last 5 volumes
        baseline = np.mean(trace[remove_ends:(time_zero-remove_ends)])
        # ignore voxels those baseline fluorescence is smaller than threshold
        if baseline>=threshold:
            dff = trace/baseline - 1
            # calculate real time for each volume based on volumetric imaging rate
            real_time = np.arange(trace.shape[0])/imaging_freq - pre_puff
            # calculate area by find the highest point, then drop to noise level
            area = ZZ_calculate_area_peakDrop(dff, real_time, interval=search_interval, markes=markes, plot_out=plot_out)
        else:
            area = 0
        return(area)

    ## calculate the area volume 
    area = np.zeros(list(smoothed2.shape[0:3]))
    # loop through each pixel
    for xx in range(smoothed2.shape[0]):
        for yy in range(smoothed2.shape[1]):
            for zz in range(smoothed2.shape[2]):
                trace = smoothed2[xx, yy, zz, :]
                area_temp = ZZ_calculate_area_local(trace)
                area[xx, yy, zz] = area_temp

    # smooth the final area volume
    area_smoothed = ndimage.gaussian_filter(area, gaussian_spatial)
    return(area_smoothed)



def ZZ_calculate_area_volume_fixedInterval(fn, movie_info, channel=0, imaging_freq=3.76, pre_puff=30, interval_baseline_markes=[-25,-5], interval_peak_markes=[0,120], threshold=10):
    """given an odor-evoked movie, calculate an area volume by integrating over a fixed time interval
    fn: String, file name of the movie, full path
    movie_info: Dict, xyz size, number of channels and etc, used to reshape the np array
    imaging_freq: Float, volumetric imaging rate, used to calculate real time
    pre_puff: Int, how many seconds before puffing were recorded, used to identify time zero
    interval_baseline_markes: List, what time interval to calculate baseline fluorescence, unit is seconds
    interval_peak_markes: List, what time interval to integrate for area, unit is seconds
    threshold: Int, ignore voxels whose baseline is smaller than threshold
    return: 3d-array, area volume"""

    # parameters that are relatively fixed
    # filter size to smooth the data
    gaussian_spatial = [4,4,2]
    gaussian_temporal = 3

    ## read in the movie
    vol = ZZ_read_reshape_tiff(fn, movie_info)
    # only use the green channel
    movie = np.squeeze(vol[:,:,channel,:,:])

    ## smoothing to remove noise
    smoothed2 = ZZ_smooth_movies(movie, gaussian_spatial=gaussian_spatial, gaussian_temporal=gaussian_temporal)

    ## calculate df/f 
    # calculate real time
    real_time = np.arange(smoothed2.shape[3])/imaging_freq - pre_puff
    # calculate a baseline volume
    # convert real time into index
    left = np.where(real_time>=interval_baseline_markes[0])
    right = np.where(real_time<=interval_baseline_markes[1])
    baseline_range = np.intersect1d(left, right)
    left = np.where(real_time>=interval_peak_markes[0])
    right = np.where(real_time<=interval_peak_markes[1])
    area_range = np.intersect1d(left, right)
    baseline = np.mean(smoothed2[:,:,:,baseline_range], axis=3)
    # remove non-AL pixels by thresholding
    idx_good = baseline > threshold
    # calculate df/f for each time point
    dff_all = np.zeros(smoothed2.shape)
    # reshape baseline volume for broadcasting
    baseline = baseline[idx_good].reshape(-1,1)
    dff_all[idx_good,:] = smoothed2[idx_good, :] / baseline - 1

    ## calculate area by integrating
    area = np.trapz(dff_all[:,:,:,area_range], axis=3)
    # divide area with imaging frequency, so unit for area is df/f * second
    area = area / imaging_freq
    # smooth the area volume with Gaussian filter 
    area_smoothed = ndimage.gaussian_filter(area, gaussian_spatial)
    return(area_smoothed)

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

def ZZ_combine_dff_rawNew(channel, odor_pattern_string, channel_odorant, file_range, brain_name, folder_dff, tiff_type, folder_save_merged=''):
    """combine dff area volume of the same odorant into one master volume to check similarity between puffs
    channel: String, what channels to merge, e.g. 'A0'
    odor_pattern_string: String, order of delivery for each channel
    channel_odorant: Dict, what odorant in what channel, e.g. {'A0':'decanal'}
    file_range: List, file numbers that corresponds to the odor_pattern_string
    brain_name: String, brain name
    folder_dff: String, folder that stores the raw dff area volume
    tiff_type: List, what tif types to merge, e.g. ['area.tif', 'areaNormalized.tif']
    folder_save_merged: String, folder to store the merged tif
    return: 3d-array, merged volume of the same odorant"""
    odor_files = ZZ_parse_odor_stringNew(odor_pattern_string, file_range, [channel])[channel]
    # go through each puff, read in the volumes
    dff_tiffs = []
    for i in range(len(odor_files)):
        file_base_glob = f'{brain_name}_{odor_files[i]:05d}*{tiff_type}'
        fn = glob.glob(os.path.join(folder_dff, file_base_glob))[0]
        dff_volume = tif.imread(fn)
        dff_tiffs.append(dff_volume)
    # combine the dff from one brain; default is combine horizontally
    dff_combined = np.concatenate(dff_tiffs, axis=2)
    # save if specified
    if folder_save_merged:
        tif.imsave(os.path.join(folder_save_merged,
                                f'{brain_name}_{channel_odorant[channel]}_{tiff_type}_combined.tif'), dff_combined)
    return(dff_combined)



def ZZ_average_dff_rawNew(channel, odor_pattern_string, channel_odorant, file_range, brain_name, folder_dff, tiff_type, folder_save_average=''):
    """average multiple puffs of the same odorant into one volume
    channel: String, what channels to merge, e.g. 'A0'
    odor_pattern_string: String, order of delivery for each channel
    channel_odorant: Dict, what odorant in what channel, e.g. {'A0':'decanal'}
    file_range: List, file numbers that corresponds to the odor_pattern_string
    brain_name: String, brain name
    folder_dff: String, folder that stores the raw dff area volume
    tiff_type: List, what tif types to merge, e.g. ['area.tif', 'areaNormalized.tif']
    folder_save_average: String, folder to store the averaged tif
    return: 3d-array, average volume of the same odorant"""
    odor_files = ZZ_parse_odor_stringNew(odor_pattern_string, file_range, [channel])[channel]
    # go through each puff, read in the volumes
    dff_tiffs = []
    for i in range(len(odor_files)):
        file_base_glob = f'{brain_name}_{odor_files[i]:05d}*{tiff_type}'
        fn = glob.glob(os.path.join(folder_dff, file_base_glob))[0]
        dff_volume = tif.imread(fn)
        dff_tiffs.append(dff_volume)
    # calculate average
    dff_tiffs = np.stack(dff_tiffs)
    area_ave = np.mean(dff_tiffs, axis=0)
    # save if specified
    if folder_save_average:
        tif.imsave(os.path.join(folder_save_average, f'{brain_name}_{channel_odorant[channel]}_{tiff_type}_average.tif'), area_ave)
    return(area_ave)



def ZZ_time_trace_glomerulus(fn, center, onset_point = 26, radius=7, dim_x=128, dim_y=128, dim_z=24, subtract_background=True, remove_ends=3):
    """get the time trace given the file name and the center of a glomerulus
    fn: String, file name of the movie
    center: List, formated [z, x, y], values as indicated in Fiji
    onset_point: Int, odor got puffed at what volume
    radius: Int, what large is the circle to average over
    dim_x, dim_y, dim_z: Int, dimension of the volume
    subtract_background: Logical, whether to subtract the background, if True returns df/f
    remove_ends: Int, how many volumes to exclude when calculate baseline fluorescence
    return: 1-d array, time trace of the glomerulus"""
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
    baseline_range = range(remove_ends, onset_point-remove_ends)
    if subtract_background:
        background = np.mean(ave_all[baseline_range])
        ave_all = [(v-background)/background for v in ave_all]
    return ave_all



def ZZ_time_trace_glomerulus3D(fn, center, onset_point = 26, radius=7, dim_x=128, dim_y=128, dim_z=24, subtract_background=True, remove_ends=3):
    """get the time trace given the file name and the center of a glomerulus
    fn: String, file name of the movie
    center: List, formated [z, x, y], values as indicated in Fiji
    onset_point: Int, odor got puffed at what volume
    radius: Int, what large is the circle to average over
    dim_x, dim_y, dim_z: Int, dimension of the volume
    subtract_background: Logical, whether to subtract the background, if True returns df/f
    remove_ends: Int, how many volumes to exclude when calculate baseline fluorescence
    return: 1-d array, time trace of the glomerulus"""
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
    baseline_range = range(remove_ends, onset_point-remove_ends)
    if subtract_background:
        background = np.mean(res[baseline_range])
        res = [(v-background)/background for v in res]
    return res



def ZZ_combine_time_traces(fns, center, onset_point, radius=7, dim_x=128, dim_y=128, dim_z=24, imaging_freq=3.76, remove_ends=3):
    """combine df/f time traces for different puffs of the same odorant
    fns: List of Strings, file names of different puffs
    center: List, formated [z, x, y], values as indicated in Fiji
    onset_point: Int, odor got puffed at what volume
    radius: Int, what large is the circle to average over
    dim_x, dim_y, dim_z: Int, dimension of the volume
    imaging_freq: Float, volumetric imaging rate
    remove_ends: Int, how many volumes to exclude when calculate baseline fluorescence
    return: 2-d array, last row: real time; 2nd last row: averaged df/f trace; remaining: df/f trace of each puff"""
    # get the time trace from each file
    glo_ave_all = []
    for fn in fns:
        glo_ave =  ZZ_time_trace_glomerulus3D(fn, center, onset_point, radius, dim_x, dim_y, dim_z, True, remove_ends)
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
    trace_value: Lists of 1-d array, calcium traces need to be aligned
    time_value: Lists of 1-d array, real time of corresponding calcium traces
    imaging_freq: Float, volumetric imaging rate
    return: 2-d array, last row: real time; 2nd last row: averaged df/f trace; remaining: aligned df/f trace of each puff"""
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
        comb[ii, offsets[ii]:(offsets[ii]+trace_value[ii].shape[0])] = trace_value[ii]
    # calculate the average
    ave = np.mean(comb[0:len(trace_value), :], axis=0)
    comb[-2, :] = ave
    # append the real time
    comb[-1, :] = (np.arange(comb.shape[1]) - left_long)/imaging_freq
    return comb



def ZZ_quantify_dff_area_brain(data, odor_single, odor_markes, odor_markesFake, search_interval_single=[0,15], search_interval_markes=[0, 60], normalization_factor=1, plot_out=False):
    """quantify area under curve for df/f time trace of selected odors
    data: Dict, odor name as keys and calcium traces as values
        calcium trace is a 2-d array that has averaged trace in the 2nd last row and real time in the last row
    odor_single, odor_markes, odor_fakeMarkes: Lists of odorant names
    search_interval_single: List, for single odorants, what time interval to expect peak
    search_interval_markes: List, for markes, what time interval to expect peak
    normalization factor: Float, factor to normalize between brains
    plot_out: Logical, whether to plot the calcium trace and landmark points
    return: Dict, odor name as keys, and peak area as values"""
    # store results in dict
    temp_area = {}
    for odor in odor_single + odor_markes + odor_markesFake:
        trace = data[odor]
        # get the averaged df/f trace
        fluo = trace[-2, :]
        # get the real time
        real_time = trace[-1, :]
        # normalize the raw response
        fluo = fluo / normalization_factor
        if odor in odor_single:
            area = ZZ_calculate_area_peakDrop(fluo, real_time, interval=search_interval_single, markes=False, plot_out=plot_out)
        elif odor in odor_markes:
#             left = np.where(real_time>=interval_peak_markes[0])
#             right = np.where(real_time<=interval_peak_markes[1])
#             area_range = np.intersect1d(left, right)
#             area = np.trapz(fluo[area_range], real_time[area_range])
            area = ZZ_calculate_area_peakDrop(fluo, real_time, interval=search_interval_markes, markes=True, plot_out=plot_out)
        else:
            area = ZZ_calculate_area_peakDrop(fluo, real_time, interval=search_interval_single, markes=False, plot_out=plot_out)
        temp_area[odor] = area
    return(temp_area)



def ZZ_quantify_dff_area_brain_fromRaw(data, selected_odor_single, selected_odor_markes, selected_odor_markesFake, search_interval_single=[0,15], search_interval_markes=[0, 60], normalization_factor=1, plot_out=False):
    """quantify area under curve for df/f time trace of selected odors
    data: Dict, odor name as keys and calcium traces as values
        calcium trace is a 2-d array that has averaged trace in the 2nd last row and real time in the last row
    selected_odor_single, selected_odor_markes, selected_odor_fakeMarkes: Lists of odorant names
    search_interval_single: List, for single odorants, what time interval to expect peak
    search_interval_markes: List, for markes, what time interval to expect peak
    normalization factor: Float, factor to normalize between brains
    plot_out: Logical, whether to plot the calcium trace and landmark points
    return: Dict, odor name as keys, and peak area as values"""
    # get the name of Null puffs and solvent puffs
    odor_single_all = selected_odor_single + [a+'Null' for a in selected_odor_single] + ['paraffin', 'paraffinNull'] + [a+'Null' for a in selected_odor_markesFake]
    odor_markes_all = selected_odor_markes + ['conditioned']
    # calculate the area from raw traces
    area_raw = ZZ_quantify_dff_area_brain(data, odor_single=odor_single_all, odor_markes=odor_markes_all, odor_markesFake=selected_odor_markesFake, search_interval_single=search_interval_single, search_interval_markes=search_interval_markes, normalization_factor=normalization_factor, plot_out=plot_out)
    # subtract the area for null puff
    area_nullSubtracted = {}
    for odor in selected_odor_single + selected_odor_markesFake + ['paraffin']:
        area_odor = area_raw[odor]
        area_null = area_raw[odor+'Null']
        area_nullSubtracted[odor] = area_odor - area_null
    for odor in selected_odor_markes + ['conditioned']:
        area_nullSubtracted[odor] = area_raw[odor]
    # subtract the area for solvent
    area_solventSubtracted = {}
    for odor in selected_odor_single + selected_odor_markesFake:
        area_odor = area_nullSubtracted[odor]
        area_solvent = area_nullSubtracted['paraffin']
        area_solventSubtracted[odor] = area_odor - area_solvent
    for odor in selected_odor_markes:
        area_odor = area_nullSubtracted[odor]
        area_solvent = area_nullSubtracted['conditioned']
        area_solventSubtracted[odor] = area_odor - area_solvent
    return(area_solventSubtracted)



def ZZ_martin_calculate_area_matrix(option_ds, refbrain, pre_puff_single=25, pre_puff_markes=30, imaging_freq=3.76, search_interval_single=[0,25], interval_peak_markes=[0,140]):
    """calculate df/f area for segmented glomeruli
    option_ds: Dict, dataset-specific parameters
    refbrain: String, what brain is the target brain
    pre_puff_single: Int, how many seconds recorded before the puff for single odorants
    pre_puff_markes: Int, how many seconds recorded before trap fire for Markes puffs
    imaging_freq: Float, volumetric imaging rate
    search_interval_single: List, for single-odorant puffs, where to expect the highest/lowest point, unit is seconds
    interval_peak_markes: List, what time interval to integrate df/f for Markes puffs
    return: (area_matrix, ave_area, time_traces)
    area_matrix: List of 2-d arrays, response matrix, where rows are file index, columns are glomeruli id
    ave_area: List of pandas dataframes, averaged response of different puffs of the same odorant
    time_traces: List of Lists of 2-d array, time traces of each glomerulus of each brain"""

    # unpack the dataset-specific parameters
    experiment_name = option_ds['experiment_name']
    brains = option_ds['brains']
    datasets = option_ds['datasets']
    odor_pattern_strings_single = option_ds['odor_pattern_strings_single']
    channel_odorant_single = option_ds['channel_odorant_single']
    channel_solvent_single = option_ds['channel_solvent_single']
    file_ranges_odorEvoked_single = option_ds['file_ranges_odorEvoked_single']
    odor_pattern_strings_markes = option_ds['odor_pattern_strings_markes']
    channel_odorant_markes = option_ds['channel_odorant_markes']
    channel_solvent_markes = option_ds['channel_solvent_markes']
    channel_fakeMarkes = option_ds['channel_fakeMarkes']
    file_ranges_odorEvoked_markes = option_ds['file_ranges_odorEvoked_markes']

    # universal parameters
    folder_home = '/mnt/bucket/labs/mcbride/zhileizhao/Scotty/Analysis/CalciumImaging/AnalysisPipelineOrco'
    folder_martin = os.path.join(folder_home, 'MartinStrauch')
    movie_info = {'xRes': 128, 'yRes': 128, 'nChannel': 1, 'zStacks': 24, 'dtype': np.dtype('float32'), 'zoomin': 10}

    # read the segmentation results from Convex cone algorithm
    folder_seg = os.path.join(folder_home, 'MartinStrauch', 'SegmentationResults', experiment_name)
    # put all segmentation into one list
    seg_all = []
    for brain in brains:
        # for refbrain need to read the new segmentation where some glomeruli are merged
        if brain == refbrain:
            fn = os.path.join(folder_seg, f'{brain}_segmentRemoved.tif')
        else:
            fn = os.path.join(folder_seg, f'{brain}_segment.tif')
        vol = ZZ_read_reshape_tiff(fn, movie_info)
        vol = np.squeeze(vol)
        seg_all.append(vol)
    
    # save results in list
    area_matrix = []
    time_traces = []
    # loop through each brain    
    for bi in range(len(brains)):
        seg = seg_all[bi]
        num_glo = int(np.amax(seg))
        traces_brain = []
        all_file_range = sorted([v for v in file_ranges_odorEvoked_single[bi]] + [v for v in file_ranges_odorEvoked_markes[bi]])
        area_brain = np.zeros([len(all_file_range), num_glo])
        # loop through each file idx of single odorants
        for count_fi, fi in enumerate(all_file_range):
            fn = glob.glob(os.path.join(folder_home, 'MartinStrauch', 'MotionCorrectedData', experiment_name, brains[bi], 'odorEvoked', f'{brains[bi]}_{fi:05d}_*.tif'))[0]
            movie = ZZ_read_reshape_tiff(fn, movie_info)
            movie = np.squeeze(movie)
            # for a given movie, time traces is a 2-d array, rows are glomeruli_id-1, columns are time-points
            traces_file = np.zeros([num_glo, movie.shape[3]])
            # get time traces and area of each glomerulus
            for glo in range(1, num_glo+1):
                ROI = seg==glo
                trace = np.mean(movie[ROI], axis=0)
                if fi in file_ranges_odorEvoked_single[bi]:
                    pre_puff = pre_puff_single
                else:
                    pre_puff = pre_puff_markes
                # convert raw traces to df/f
                real_time = np.arange(len(trace))/imaging_freq - pre_puff
                # at what time point valve open
                time_zero = int(imaging_freq * pre_puff)
                # calculate baseline fluorescence, excluding the 1st and last 3 volumes
                baseline = np.mean(trace[3:(time_zero-3)])
                dff = trace/baseline - 1
                traces_file[glo-1,:] = dff
                dff = ndimage.gaussian_filter1d(dff, 3)
                # calculate area, for single odorants use peakDrop method; for Markes integrate over fixed interval
                if fi in file_ranges_odorEvoked_single[bi]:
                    area = ZZ_calculate_area_peakDrop(dff, real_time, interval=search_interval_single, markes=False, plot_out=False)
                else:
                    left = np.where(real_time>=interval_peak_markes[0])
                    right = np.where(real_time<=interval_peak_markes[1])
                    area_range = np.intersect1d(left, right)
                    # integrate area, unit is df/f * sec
                    area = np.trapz(dff[area_range], real_time[area_range])
                area_brain[count_fi, glo-1] = area
            traces_brain.append(traces_file)
        area_matrix.append(area_brain)
        time_traces.append(traces_brain)
        
    # average the response for each odorant
    area_averaged = []
    channels = list(channel_odorant_single.keys()) + list(channel_odorant_markes.keys())
    odorants = []
    for bi in range(len(brains)):
        area_raw = area_matrix[bi]
        area_this = OrderedDict()
        for ch in channels:
            file_idx =  sorted([v for v in file_ranges_odorEvoked_single[bi]] + [v for v in file_ranges_odorEvoked_markes[bi]])
            # get the file index for the channel
            if ch in channel_odorant_single:
                fidx = ZZ_parse_odor_stringNew(odor_pattern_strings_single[bi], file_ranges_odorEvoked_single[bi], [ch])[ch]
                if bi==0:
                    odorants.append(channel_odorant_single[ch])
            else:
                fidx = ZZ_parse_odor_stringNew(odor_pattern_strings_markes[bi], file_ranges_odorEvoked_markes[bi], [ch])[ch]
                if bi==0:
                    odorants.append(channel_odorant_markes[ch])
            fidx_idx = [file_idx.index(a) for a in fidx]
            ave = np.mean(area_raw[fidx_idx, ], axis=0)
            area_this[ch] = ave
        area_averaged.append(area_this)    
    # stack into one 3d array [brain, channel, glomeruli]
    area_ave_list = [np.stack([area_averaged[bi][ch] for ch in channels]) for bi in range(len(brains))]   
    # convert to pandas dataframes
    ave_area = []
    for bi in range(len(brains)):
        area_pd = pd.DataFrame(area_ave_list[bi])
        # use odorant name as index
        area_pd.index = odorants
        # use glomeruli id as column names
        area_pd.columns = [f'glo_{a+1}' for a in range(area_pd.shape[1])]
        ave_area.append(area_pd)
    
    return(area_matrix, ave_area, time_traces)


def ZZ_subtract_traces_shrink(trace1, time1, trace2, time2, search_interval_single):
    area1, points1 = ZZ_calculate_area_peakDrop(trace1, time1, interval=search_interval_single, markes=False, plot_out=False, return_points=True)
    area2, points2 = ZZ_calculate_area_peakDrop(trace2, time2, interval=search_interval_single, markes=False, plot_out=False, return_points=True)
    # how much to shrink
    if area1 == 0:
        shrink = 1
    else:
        shrink = (area1 - area2) / area1
    shrinked = trace1.copy()
    shrinked[points1[1]:points1[2]] = shrinked[points1[1]:points1[2]] * shrink
    # gaussian smooth
    smoothed = ndimage.gaussian_filter1d(shrinked, 2)
    return(smoothed)