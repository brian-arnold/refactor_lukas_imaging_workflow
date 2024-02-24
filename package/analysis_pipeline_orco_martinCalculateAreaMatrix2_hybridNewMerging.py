## Zhilei Zhao
## run Martin's convex cone algorithm to segment AL based on correlated activity

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
import scipy.io as sio

# import the custom module
# functions in this module are shared between different pipelines
# use the reload function since may change those functions and want to see effects immediately
import ZZ_tiff_file
from ZZ_tiff_file import *
import ZZ_CMTK
from ZZ_CMTK import *
import ZZ_calcium_responseNew
from ZZ_calcium_responseNew import *


def analysis_pipeline_orco_martinCalculateAreaMatrix2_hybridNewMerging(option_ds, refbrain=[], norm_range=[], at_least=3, pre_puff_single=7, pre_puff_markes=30, imaging_freq=3.76, search_interval_single=[0,25], interval_peak_markes=[0,140], solvent_single='paraffin', solvent_markes='conditioned', exception_single=[]):
    """calculate df/f area for segmented glomeruli
    option_ds: Dict, dataset-specific parameters
    refbrain: String, what brain is the target brain (could be brain in other datasets)
    norm_range: List, normalize response to the highest value among a list of odorants
    at_least: when calculate response for matched glomeruli, ignore glomeruli that are matched in fewer brains
    pre_puff_single: Int, how many seconds recorded before the puff for single odorants
    pre_puff_markes: Int, how many seconds recorded before trap fire for Markes puffs
    imaging_freq: Float, volumetric imaging rate
    search_interval_single: List, for single-odorant puffs, where to expect the highest/lowest point, unit is seconds
    solvent_single: String, what to subtract for single odorants, could be empty
    solvent_markes: String, what to subtract for Markes, could be empty
    exception_single: List, what single odorant to use fixed interval method to calculate area
    interval_peak_markes: List, what time interval to integrate df/f for Markes puffs
    return: (area_matrix, time_traces, raw_area_matrix, raw_time_traces, area_matched_ave, time_traces_matched) if refbrain is not empty
    return: (area_matrix, time_traces, raw_area_matrix, raw_time_traces) if refbrain is empty
    area_matrix: List (brains) of pandas dataframe, response matrix, where rows are odorant, columns are glomeruli id; averaged over multiple puffs of the same odorant, null and solvent tube subtracted
    time_traces: List (brains) of Dicts (odorants) of 2-d array (glo*timepoints), time traces of each glomerulus of each brain for given odor; averaged over multiple puffs of the same odorant, null and solvent tube subtracted
    raw_area_matrix: Lists of Dicts, raw area for each puff before subtracting null and solvent response
    raw_time_traces: List of Dicts of 2-d array, time traces of each glomerulus of each brain for given puff
    area_matched_ave: pandas dataframe, averaged response for matched glomeruli with at_least brains
    time_traces_matched: List of Dicts of 2-d array, time traces of matched glomeruli, if not matched in one brain fill in with nan
    """

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
    folder_home = '/jukebox/mcbride/bjarnold/refactor_lukas_imaging_workflow/data'
    movie_info = {'xRes': 128, 'yRes': 128, 'nChannel': 1, 'zStacks': 24, 'dtype': np.dtype('float32'), 'zoomin': 10}

    # read the segmentation results from Convex cone algorithm
    folder_seg = os.path.join(folder_home, 'MartinStrauch_hybridNewMerging', 'SegmentationResults', experiment_name)
    # put all segmentation into one list
    seg_all = []
    for brain in brains:
        # for refbrain need to read the new segmentation where some glomeruli are merged
    #         if brain == refbrain:
    #    fn = os.path.join(folder_seg, f'{brain}_segmentRemoved.tif')
    #         else:
        fn = os.path.join(folder_seg, f'{brain}_segment.tif')
        vol = ZZ_read_reshape_tiff(fn, movie_info)
        vol = np.squeeze(vol)
        seg_all.append(vol)# save results in list

    area_matrix = []
    time_traces = []
    raw_area_matrix = []
    raw_time_traces = []
    # loop through each brain    
    for bi in range(len(brains)):
    #     for bi in range(1):
        seg = seg_all[bi]
        num_glo = int(np.amax(seg))
        area_brain = {}
        traces_brain = {}
        all_file_range = sorted([v for v in file_ranges_odorEvoked_single[bi]] + [v for v in file_ranges_odorEvoked_markes[bi]])

        # define a local function to calculate area
        def ZZ_glo_area_local(fi):
            """return the area of all glomeruli for this file
            input: fi, Int, file index
            return: area for this file, should be 1-d array"""
            # this folder originally grabbed from registered movies?
            # fn = glob.glob(os.path.join(folder_home, 'MartinStrauch_hybridNewMerging', 'MotionCorrectedData', experiment_name, brains[bi], 'odorEvoked', f'{brains[bi]}_{fi:05d}_*.tif'))[0]            
            fn = glob.glob(os.path.join(folder_home, 'Registration', 'RegisteredMovies', experiment_name, brains[bi], f'{brains[bi]}_{fi:05d}_*.tif'))[0]
            movie = ZZ_read_reshape_tiff(fn, movie_info)


            #movie = np.reshape(movie, (movie_info['xRes'], movie_info['yRes'], movie_info['nChannel'], 113, -1), order='F')
            
            movie = np.squeeze(movie)
            if fi in file_ranges_odorEvoked_single[bi]:
                pre_puff = pre_puff_single
                ch = ZZ_reverse_parse_odorString(odor_pattern_strings_single[bi], file_ranges_odorEvoked_single[bi], fi)
                odorant = channel_odorant_single[ch]
            else:
                pre_puff = pre_puff_markes
                ch = ZZ_reverse_parse_odorString(odor_pattern_strings_markes[bi], file_ranges_odorEvoked_markes[bi], fi)
                odorant = channel_odorant_markes[ch]
            # for a given file, time traces is a 2-d array, rows are glomeruli_id-1, columns are time-points
            traces_file = np.zeros([num_glo, movie.shape[3]])
            # get area of each glomerulus
            area_this = []
            for glo in range(1, num_glo+1):
                ROI = seg==glo
                trace = np.mean(movie[ROI], axis=0)
                # convert raw traces to df/f
                real_time = np.arange(len(trace))/imaging_freq - pre_puff
                # at what time point valve open
                time_zero = int(imaging_freq * pre_puff)
                # calculate baseline fluorescence, excluding the 1st and last 3 volumes
                baseline = np.mean(trace[3:(time_zero-3)])
                dff = trace/baseline - 1
                traces_file[glo-1,:] = ndimage.gaussian_filter1d(dff,1)
                dff = ndimage.gaussian_filter1d(dff, 3)
                # calculate area, for single odorants and fakeMarkes use peakDrop method; for Markes integrate over fixed interval
                if (fi in file_ranges_odorEvoked_single[bi] or ch in channel_fakeMarkes) and odorant not in exception_single:
                    #print('peakDrop for ', ch, fi) 
                    area = ZZ_calculate_area_peakDrop(dff, real_time, interval=search_interval_single, markes=False, plot_out=False)
                    #area = ZZ_calculate_area_volume_fixedInterval(fn, movie_info, channel=0, imaging_freq=3.76, pre_puff=7, interval_baseline_markes=[-6,-1], interval_peak_markes=[0,35], threshold=10)
                else:
                    left = np.where(real_time>=interval_peak_markes[0])
                    right = np.where(real_time<=interval_peak_markes[1])
                    area_range = np.intersect1d(left, right)
                    # integrate area, unit is df/f * sec
                    area = np.trapz(dff[area_range], real_time[area_range])
                area_this.append(area)
            return(np.array(area_this), traces_file, odorant)

        # run the function in pararell
        with Parallel(n_jobs=22, verbose=5) as parallel:
            res = parallel(delayed(ZZ_glo_area_local)(fi) for fi in all_file_range)
        for a in res:
            # for some odorant, there are many puffs
            if a[2] in area_brain:
                area_brain[a[2]].append(a[0])
                traces_brain[a[2]].append(a[1])
            else:
                area_brain[a[2]] = [a[0]]
                traces_brain[a[2]] = [a[1]]

        # loop through each odorant, average the multiple puffs
        all_odorant = area_brain.keys()
        # what odorants are fake Markes
        odorant_fakeMarkes = [channel_odorant_markes[a] for a in channel_fakeMarkes]
        ave_trace = {}  
        ave_area = {}
        for odorant in all_odorant:
            if odorant in channel_odorant_single.values():
                pre_puff = pre_puff_single 
            else:
                pre_puff = pre_puff_markes

            tt = traces_brain[odorant]
            # the total length of recording may differ by 1, choose the shortest
            record_len = min([a.shape[1] for a in tt])
            # average the two puffs
            tt_temp = np.stack([a[:,0:record_len] for a in tt])
            tt_ave_odorant = np.mean(tt_temp, axis=0)
            # calculate a vector of real time
            real_time = np.arange(tt_ave_odorant.shape[1])/imaging_freq - pre_puff
            # append the real time to the averaged trace
            tt_ave_odorant = np.vstack([tt_ave_odorant, real_time])
            ave_trace[odorant] = tt_ave_odorant
            # also average the area
            area_temp = np.stack(area_brain[odorant])
            ave_area[odorant] = np.mean(area_temp, axis=0)

        # subtract the null response
        null_subtracted_area = {}
        null_subtracted = {}
        for odorant in all_odorant:
            tt_ave_odorant = ave_trace[odorant]
            if odorant+'Null' in all_odorant:
                null = odorant + 'Null'
                tt_ave_null = ave_trace[null]
                # align the average traces by time 0
                # get where time zero is
                zero_odorant = np.argmin(np.abs(tt_ave_odorant[-1, :]))
                zero_null = np.argmin(np.abs(tt_ave_null[-1, :]))
                # calculate which is longer on each side
                left_long = np.max([zero_odorant, zero_null])
                right_long = np.max([tt_ave_odorant.shape[1]-zero_odorant, tt_ave_null.shape[1]-zero_null])
                # put two response into new numpy arrays
                sub_odorant = np.zeros((tt_ave_odorant.shape[0]-1, left_long+right_long))
                sub_null = np.zeros((tt_ave_null.shape[0]-1, left_long+right_long))
                offset_odorant = left_long - zero_odorant
                offset_null = left_long - zero_null
                sub_odorant[:, offset_odorant:(offset_odorant+tt_ave_odorant.shape[1])] = tt_ave_odorant[0:-1,:]
                sub_null[:, offset_null:(offset_null+tt_ave_null.shape[1])] = tt_ave_null[0:-1,:]
                sub = sub_odorant - sub_null
                real_time = (np.arange(sub.shape[1]) - left_long) / imaging_freq
                sub = np.vstack([sub, real_time])
                null_subtracted[odorant] = sub
                # also subtract the area
                null_subtracted_area[odorant] = ave_area[odorant] - ave_area[null]
            else:
                null_subtracted[odorant] = tt_ave_odorant
                null_subtracted_area[odorant] = ave_area[odorant]

        # subtract the solvent response
        solvent_subtracted = {}
        solvent_subtracted_area = {}
        for odorant in all_odorant:
            tt_ave_odorant = null_subtracted[odorant]
            if 'Null' not in odorant:
                if odorant in channel_odorant_single.values() or odorant in odorant_fakeMarkes:
                    solvent = solvent_single
                else:
                    solvent = solvent_markes
                if solvent:
                    tt_ave_solvent = null_subtracted[solvent]
                    # align the average traces by time 0
                    # get where time zero is
                    zero_odorant = np.argmin(np.abs(tt_ave_odorant[-1, :]))
                    zero_solvent = np.argmin(np.abs(tt_ave_solvent[-1, :]))
                    # calculate which is longer on each side
                    left_long = np.max([zero_odorant, zero_solvent])
                    right_long = np.max([tt_ave_odorant.shape[1]-zero_odorant, tt_ave_solvent.shape[1]-zero_solvent])
                    # put two response into new numpy arrays
                    sub_odorant = np.zeros((tt_ave_odorant.shape[0]-1, left_long+right_long))
                    sub_solvent = np.zeros((tt_ave_solvent.shape[0]-1, left_long+right_long))
                    offset_odorant = left_long - zero_odorant
                    offset_solvent = left_long - zero_solvent
                    sub_odorant[:, offset_odorant:(offset_odorant+tt_ave_odorant.shape[1])] = tt_ave_odorant[0:-1,:]
                    sub_solvent[:, offset_solvent:(offset_solvent+tt_ave_solvent.shape[1])] = tt_ave_solvent[0:-1,:]
                    sub = sub_odorant - sub_solvent
                    real_time = (np.arange(sub.shape[1]) - left_long) / imaging_freq
                    sub = np.vstack([sub, real_time])
                    solvent_subtracted[odorant] = sub
                    # also subtract the area
                    solvent_subtracted_area[odorant] = null_subtracted_area[odorant] - null_subtracted_area[solvent]
                else:
                    solvent_subtracted[odorant] = tt_ave_odorant
                    solvent_subtracted_area[odorant] = null_subtracted_area[odorant]
            else:
                solvent_subtracted[odorant] = tt_ave_odorant
                solvent_subtracted_area[odorant] = null_subtracted_area[odorant]

        # convert to pandas dataframe
        solvent_subtracted_area = pd.DataFrame(solvent_subtracted_area)
        solvent_subtracted_area = solvent_subtracted_area.T
        # append to the mast list
        area_matrix.append(solvent_subtracted_area)
        raw_area_matrix.append(area_brain)
        time_traces.append(solvent_subtracted)
        raw_time_traces.append(traces_brain)

    # extract response for matched glomeruli if specified refbrain
    if refbrain:
        if refbrain in brains:
            fn_matching = os.path.join(folder_home, 'MartinStrauch_hybridNewMerging', 'MatchingResults', experiment_name, f'{experiment_name}_{refbrain}_asRefMatchingIndex.mat')
        else:
            fn_matching = os.path.join(folder_home, 'MartinStrauch_hybridNewMerging', 'MatchingResults', experiment_name, f'{experiment_name}_MatchingIndex.mat')
        match_idx = sio.loadmat(fn_matching)['matching_index']
        # use what reference brain? If matched to brain in a different dataset, set ref_brain_idx=-1
        if refbrain in brains:
            ref_brain_idx = brains.index(refbrain)
            ref_idx = match_idx[ref_brain_idx, :]
        else:
            ref_idx = match_idx[0, :]
            match_idx = match_idx[1:match_idx.shape[0],:]

        # extract the response and time traces
        area_matched = []
        time_traces_matched = []
        for bi in range(len(brains)): 
            area_data = area_matrix[bi]
            trace_data = time_traces[bi]
            # normalize if specified
            norm_factor = 1
            norm_factor2 = 1
            if norm_range:
                norm_factor = np.amax(np.array(area_data.loc[norm_range,:]))
                norm_factor2 = np.amax([np.amax(trace_data[odorant][0:-1,:]) for odorant in norm_range])
            area_data = area_data / norm_factor
            # all brains use the same index order
            area_data = area_data.loc[area_matrix[0].index,:]
            idx_this = match_idx[bi,:]
            area_this = []
            for ii in idx_this:
                if ii > 0:
                    area_this.append(area_data.loc[:, ii-1])
                else:
                    area_this.append(np.array([np.nan]*area_data.shape[0]))
            trace_this = {}
            for odorant in trace_data.keys():
                trace_this_odorant = trace_data[odorant] / norm_factor2
                # add nan if not matched
                trace_temp = np.zeros([len(idx_this), trace_this_odorant.shape[1]])
                for idx, ii in enumerate(idx_this):
                    if ii > 0:
                        trace_temp[idx,:] = trace_this_odorant[ii-1,:]
                    else:
                        trace_temp[idx,:] = np.array([np.nan]*trace_this_odorant.shape[1])
                # append the real time
                if odorant in channel_odorant_single.values():
                    pre_puff = pre_puff_single
                else:
                    pre_puff = pre_puff_markes
                real_time = np.arange(trace_temp.shape[1])/imaging_freq - pre_puff
                trace_temp = np.vstack([trace_temp, real_time])
                trace_this[odorant] = trace_temp
            area_this = np.stack(area_this, axis=1)
            area_matched.append(area_this)
            time_traces_matched.append(trace_this)
        area_matched = np.stack(area_matched)
        # calculate averaged response for matched glomeruli
        # ignore glomeruli that don't have many matched glomeruli
        area_matched_ave = []
        passed_idx = []
        for glo in range(area_matched.shape[2]):
            temp = area_matched[:,:,glo]
            nan_count = np.isnan(temp).astype('int')
            nan_count = np.sum(nan_count, axis=0)
            if nan_count[0] <= len(brains) - at_least:
                area_matched_ave.append(np.nanmean(temp, axis=0))
                passed_idx.append(ref_idx[glo])
        area_matched_ave = np.stack(area_matched_ave, axis=1)
        # convert to pandas dataframe
        area_matched_ave = pd.DataFrame(area_matched_ave)
        area_matched_ave.index = area_matrix[0].index
        area_matched_ave.columns = passed_idx
        #return(area_matrix, time_traces, raw_area_matrix, raw_time_traces, area_matched_ave, time_traces_matched)

        folder_save = os.path.join(folder_home, 'MartinStrauch_hybridNewMerging', 'ResponseResults', experiment_name)
        if not os.path.exists(folder_save):
            os.makedirs(folder_save)
        for bi in range(len(brains)):
            ave_temp = area_matrix[bi]
            fn_csv = os.path.join(folder_save, f'{brains[bi]}_aveArea.csv')
            ave_temp.to_csv(fn_csv)
    
        # save matched response to csv files
        fn_csv = os.path.join(folder_save, f'{experiment_name}_matchedArea.csv')
        area_matched_ave.to_csv(fn_csv)
    else:
        #return(area_matrix, time_traces, raw_area_matrix, raw_time_traces)

        folder_save = os.path.join(folder_home, 'MartinStrauch_hybridNewMerging', 'ResponseResults', experiment_name)
        if not os.path.exists(folder_save):
            os.makedirs(folder_save)
        for bi in range(len(brains)):
            ave_temp = area_matrix[bi]
            fn_csv = os.path.join(folder_save, f'{brains[bi]}_aveArea.csv')
            ave_temp.to_csv(fn_csv)
            
            ave_ttemp = time_traces[bi]
            fn_t_csv = os.path.join(folder_save, f'{brains[bi]}_aveTraces.pkl')
            

            import pickle
            
            with open(fn_t_csv, "wb") as tf:
            	pickle.dump(ave_ttemp,tf)
            	
            ave_ttemp_raw = raw_time_traces[bi]
            fn_t_raw_csv = os.path.join(folder_save, f'{brains[bi]}_aveTracesRaw.pkl')
            

            import pickle
            
            with open(fn_t_raw_csv, "wb") as tf:
            	pickle.dump(ave_ttemp_raw,tf)


            
            



def main():
    
    option_ds = {'experiment_name': '240111', 
                'brains': ['240111_Camphor_1U_F1', '240111_Camphor_1U_F2', '240111_Camphor_1U_F3', '240111_Camphor_1U_F4'],
                'datasets': ['240111', '240111', '240111', '240111'], 
                'meta_file': '240111_meta_data.mat', 
                'AL_loc': ['lower_AL', 'lower_AL', 'upper_AL', 'lower_AL'], 
                'odor_pattern_strings_single': ['100A_50A_10A_100B_10B_100C_10C_100D_100F_100G_100H_100J_100L_100M_100N_10N_100O_100P_100Q_100R_100S_100T_100U_100A_50A_10A_100B_10B_100C_10C_100D_100F_100G_100H_100J_100L_100M_100N_10N_100O_100P_100Q_100R_100S_100T_100U_', '100A_50A_10A_100B_10B_100C_10C_100D_100F_100G_100H_100J_100L_100M_100N_10N_100O_100P_100Q_100R_100S_100T_100U_100A_50A_10A_100B_10B_100C_10C_100D_100F_100G_100H_100J_100L_100M_100N_10N_100O_100P_100Q_100R_100S_100T_100U_', '100A_50A_10A_100B_10B_100C_10C_100D_100F_100G_100H_100J_100L_100M_100N_10N_100O_100P_100Q_100R_100S_100T_100U_100A_50A_10A_100B_10B_100C_10C_100D_100F_100G_100H_100J_100L_100M_100N_10N_100O_100P_100Q_100R_100S_100T_100U_', '100A_50A_10A_100B_10B_100C_10C_100D_100F_100G_100H_100J_100L_100M_100N_10N_100O_100P_100Q_100R_100S_100T_100U_100A_50A_10A_100B_10B_100C_10C_100D_100F_100G_100H_100J_100L_100M_100N_10N_100O_100P_100Q_100R_100S_100T_100U_'], 
                'file_ranges_odorEvoked_single': [[2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47], [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47], [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47], [2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46, 47]], 
                'odor_pattern_strings_markes': ['', '', '', ''], 
                'file_ranges_odorEvoked_markes': [[], [], [], []], 
                'channel_odorant_single': {'100A': 'camphor-1', '100B': 'camphor-3', '100C': '1-8-cineole-2', '100D': 'sulcatone', '100F': 'alpha-terpinene', '100G': 'beta-ocimene', '100H': 'p-cymene', '100J': 'methy-butyl-acetate', '100L': 'paraffin', '100M': 'beta-myrcene', '100N': '4-terpineol-2', '100O': '3-carene', '100P': 'benzaldehyde-methyloctane', '100Q': 'phenol-cresol', '100R': 'nonanal-furfural', '100S': '1-octen-3-ol', '100T': 'heptyl-acetate-ethylhexanol', '100U': 'oxoisophorone-acetophenone', '10A': 'camphor-2', '10B': 'camphor-4', '10C': '1-8-cineole-3', '10N': '4-terpineol-3', '50A': 'camphor-1.5'}, 
                'channel_solvent_single': {'100A': 'paraffin', '100B': 'paraffin', '100C': 'paraffin', '100D': 'paraffin', '100F': 'paraffin', '100G': 'paraffin', '100H': 'paraffin', '100J': 'paraffin', '100L': 'paraffin', '100M': 'paraffin', '100N': 'paraffin', '100O': 'paraffin', '100P': 'paraffin', '100Q': 'paraffin', '100R': 'paraffin', '100S': 'paraffin', '100T': 'paraffin', '100U': 'paraffin', '10A': 'paraffin', '10B': 'paraffin', '10C': 'paraffin', '10N': 'paraffin', '50A': 'paraffin'}, 
                'channel_odorant_markes': {}, 
                'channel_solvent_markes': {}, 
                'remove_odor_brain': [], 
                'remove_odor_odor': [], 
                'single_baseline': '', 
                'single_peak': '', 
                'markes_baseline': '', 
                'markes_peak': '', 
                'channel_fakeMarkes': {}}

    analysis_pipeline_orco_martinCalculateAreaMatrix2_hybridNewMerging(option_ds)


if __name__ == '__main__':
    main()
            

            
            
            
     