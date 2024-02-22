# Zhilei Zhao, McBride Lab, Princeton University
# ## 0. Goal
# this is part 1 of the analysis pipeline for the Orco strain, which organize raw data into subfolders, do motion correction for time-lapse and high-resolution fastZ.

# general modules
import os
import shutil
import glob
import numpy as np
import subprocess
from collections import OrderedDict
import nrrd
import importlib
import datetime

# import the custom module
# functions in this module are shared between different pipelines
# use the reload function since may change those functions and want to see effects immediately
import ZZ_tiff_file
from ZZ_tiff_file import *
importlib.reload(ZZ_tiff_file)
import ZZ_CMTK
from ZZ_CMTK import *
importlib.reload(ZZ_CMTK)
import ZZ_calcium_response
from ZZ_calcium_response import *
importlib.reload(ZZ_calcium_response)


def analysis_pipeline_orco_part1(option_ds, have_markes=False):
    ## inputs:
    # option_ds: dataset-specific parameters
    # have_markes: whether this brain has Markes puffs

    ## 1. Parsing inputs
    # unpack the dataset-specific parameters
    name_dataset = option_ds['name_dataset']
    brain_names = option_ds['brain_names']
    file_ranges_spontaneous = option_ds['file_ranges_spontaneous']
    file_ranges_fastZ = option_ds['file_ranges_fastZ']
    num_stacks = option_ds['num_stacks']
    odor_pattern_strings_single = option_ds['odor_pattern_strings_single']
    channel_odorant_single = option_ds['channel_odorant_single']
    channel_solvent_single = option_ds['channel_solvent_single']
    file_ranges_odorEvoked_single = option_ds['file_ranges_odorEvoked_single']
    odor_pattern_strings_markes = option_ds['odor_pattern_strings_markes']
    channel_odorant_markes = option_ds['channel_odorant_markes']
    channel_solvent_markes = option_ds['channel_solvent_markes']
    file_ranges_odorEvoked_markes = option_ds['file_ranges_odorEvoked_markes']
    fname_markes_timestamp = option_ds['fname_markes_timestamp']
    datetime_format = option_ds['datetime_format']
    date_experiment = option_ds['date_experiment']
    AL_loc = option_ds['AL_loc']

    # folder setting
    # home folder that serves as the base folder for all analysis
    folder_home = '/jukebox/mcbride/bjarnold/refactor_lukas_imaging_workflow/data'
    # folder with additonal scripts/code
    folder_code = '/jukebox/mcbride/bjarnold/refactor_lukas_imaging_workflow/package/matlab'
    # sub-foler in MotionCorrection to store the data before motion correction
    folder_mc_raw = "BeforeCorrection"
    # sub-foler in MotionCorrection to store the with-in movie motion correction results
    folder_mc_within = 'CorrectedDataWithin'
    # sub-foler in MotionCorrection to store the between movie motion correction results
    folder_mc_between = 'CorrectedDataBetween'
    # sub-foler in MotionCorrection to store the temporary results for motion correction
    folder_mc_temp = 'Temp'

    folder_matlab_templates = '/jukebox/mcbride/bjarnold/refactor_lukas_imaging_workflow/package'

    # movie information
    # information for time-lapse and fastZ movies
    # don't specify the time points, since it may be different for different movies. It will be inferred from the file dimensions. 
    # the 'zoomin' value may not be exactly 10 for every brain, but most brains are around 10
    # registration algorithm includes scaling, so this shouldn't be a problem
    movie_info_timelapse = {'xRes': 128, 'yRes': 128, 'nChannel': 1,
                            'zStacks': num_stacks[0], 'dtype': np.dtype('int16'), 'zoomin': 10}
    # pad one zero-stack before and after the actual volume to allow movement around boundaries, so add 2 to num_stacks
    movie_info_timelapse_mc = {'xRes': 128, 'yRes': 128, 'nChannel': 1,
                               'zStacks': num_stacks[0] + 2, 'dtype': np.dtype('int16'), 'zoomin': 10}
    movie_info_timelapse_mc_template = {'xRes': 128, 'yRes': 128, 'nChannel': 1,
                                        'zStacks': num_stacks[0] + 2, 'dtype': np.dtype('float32'), 'zoomin': 10}
    movie_info_fastZ = {'xRes': 256, 'yRes': 256, 'nChannel': 1,
                        'zStacks': num_stacks[1], 'dtype': np.dtype('int16'), 'zoomin': 10}

    # header file for nrrd files, this is needed to convert Tiff files to Nrrd files that CMTK uses
    # the resolution is measured with a grid at zoom-in 10X
    nrrd_header_timelapse = OrderedDict([('type', 'float'), ('encoding', 'raw'),
                                         ('endian', 'big'), ('dimension', 3),
                                         ('sizes', np.array(
                                             [128, 128, num_stacks[0] + 2])), ('space dimension', 3),
                                         ('space directions', np.array(
                                             [[0.908, 0, 0], [0, 0.807, 0], [0, 0, 4]])),
                                         ('space units', ['microns', 'microns', 'microns'])])
    nrrd_header_fastz = OrderedDict([('type', 'float'), ('encoding', 'raw'),
                                     ('endian', 'big'), ('dimension', 3),
                                     ('sizes', np.array(
                                         [256, 256, num_stacks[1] + 2])), ('space dimension', 3),
                                     ('space directions', np.array(
                                         [[0.454, 0, 0], [0, 0.4035, 0], [0, 0, 1]])),
                                     ('space units', ['microns', 'microns', 'microns'])])
    
    # algorithm parameters
    # options for CMTK, see munger --help for explanation
    cmtk_option_warp_mc = OrderedDict([('overwrite', True), ('energy_weight', 5e-1),
                                       ('rigidity_weight', 0),
                                       ('exploration', 16), ('accuracy', 0.4),
                                       ('threads', 22), ('reformat', '01'),
                                       ('coarsest', 4), ('grid', 32),
                                       ('refine', 3)])

    # remove intermediate files or not (CMTK regitration files)
    # note this automatically removing intermediate files function is not implementated yet
    # need to remove manually 
    rm_interm = False

    
    ## 2. Motion Correction

    # ## 2.1 Organize raw files into subfolders under MotionCorrection
    # Raw data should be in the folder_raw defined above. May have several different brains in one folder. <br>
    # Separate the files into different brains. For each brain, separate into spontaneous, odorEvoked, Markes and fastZ subfolders. <br>
    # Since we don't touch the files in RawData folder, we can re-run the function over and over again.

    # create the log file
    # note the each dataset may have several log files if the script is run multiple times at different days
    datetime_now = datetime.datetime.now().strftime("%Y%m%d")
    fname_log = os.path.join(folder_home, 'Logs', f'log_{name_dataset}_{datetime_now}.txt')
    os.makedirs(os.path.join(folder_home, 'Logs'), exist_ok=True)
    log = open(fname_log, 'a')

    # make timelapse and fastZ subfolder for each brain
    log.write(f'Dataset name: {name_dataset}\n')
    log.write(f'Run pipeline at: {datetime_now}\n')
    log.write('1. Organize raw files into subfolders under MotionCorrection...\n')
    for idx, brain in enumerate(brain_names):
        folder_brain_old = os.path.join(folder_home, 'RawData', name_dataset)
        folder_brain_new = os.path.join(folder_home, 'MotionCorrection',
                                        folder_mc_raw, name_dataset, brain)
        folder_raw_odorEvoked = os.path.join(folder_brain_new, 'odorEvoked')
        folder_raw_spontaneous = os.path.join(folder_brain_new, 'spontaneous')
        folder_raw_fastZ = os.path.join(folder_brain_new, 'fastZ')
        # check if folder already exists, if so remove completely and make new
        # don't remove the entire folder, since Markes files already there after running part0 of the pipeline
        if os.path.exists(folder_raw_odorEvoked):
            shutil.rmtree(folder_raw_odorEvoked)
        if os.path.exists(folder_raw_spontaneous):
            shutil.rmtree(folder_raw_spontaneous)
        if os.path.exists(folder_raw_fastZ):
            shutil.rmtree(folder_raw_fastZ)
        os.makedirs(folder_raw_spontaneous)
        os.makedirs(folder_raw_odorEvoked)
        os.makedirs(folder_raw_fastZ)

        # copy the spontaneous data
        for sn in file_ranges_spontaneous[idx]:
            fn = brain + '_' + f'{sn:05d}' + '_00001.tif'
            fn_full = os.path.join(folder_brain_old, fn)
            fn_new = os.path.join(folder_raw_spontaneous, fn)
            shutil.copyfile(fn_full, fn_new)
        # copy the time-lapse data
        for sn in file_ranges_odorEvoked_single[idx]:
            fn = brain + '_' + f'{sn:05d}' + '_00001.tif'
            fn_full = os.path.join(folder_brain_old, fn)
            fn_new = os.path.join(folder_raw_odorEvoked, fn)
            shutil.copyfile(fn_full, fn_new)
        # copy the fastZ data
        for sn in file_ranges_fastZ[idx]:
            fn = brain + '_' + f'{sn:05d}' + '_00001.tif'
            fn_full = os.path.join(folder_brain_old, fn)
            fn_new = os.path.join(folder_raw_fastZ, fn)
            shutil.copyfile(fn_full, fn_new)

    ## 2.2. Motion correction separately within each movie
    # Process each movie separately. Use NoRMCorre package in Matlab to do 3D piece-wise rigid motion correction.
    # Seems that iteration number matters, before I only use 3, but now increase to 5

    # process one brain at a time, loop through all brains
    for brain in brain_names:
        # read-in the template script for motion correction of timelapse movies
        normcorre_timelapse = os.path.join(
            folder_matlab_templates, 'ZZ_normcorre_pwrigid_python_timelapse.m')
        script_lines = open(normcorre_timelapse, 'r').readlines()

        # first do motion correction for the spontaneous activity movie without template input
        # create a new script file for this brain
        fname_mc_run_spon = f'ZZ_normcorre_pwrigid_python_spon_{name_dataset}_{brain}.m'
        # replace the input parameters in the template script
        os.makedirs(os.path.join(folder_home, 'MotionCorrection/Scripts'), exist_ok=True)
        normcorre_timelapse_run = open(os.path.join(
            folder_home, 'MotionCorrection/Scripts', fname_mc_run_spon), 'w')
        for line in script_lines:
            line = line.replace('FOLDER_HOME', folder_home)
            line = line.replace('FOLDER_CODE', folder_code)
            line = line.replace('MOVIETYPE', 'spontaneous')
            line = line.replace('RAWDATAFOLDER', folder_mc_raw)
            line = line.replace('DATASETNAME', name_dataset)
            line = line.replace('BRAINNAME', brain)
            line = line.replace('SAVEFOLDER', folder_mc_within)
            line = line.replace('NUMSTACKS', str(num_stacks[0]))
            normcorre_timelapse_run.write(line)
        normcorre_timelapse_run.close()
        # run the normcorre algorithm for spontaneous activity, save the logs
        datetime_now = datetime.datetime.now()
        log.write(f'{datetime_now}: running normcorre for spontaneous movie...\n')
        cmd_normcorre_timelapse = f'cd {folder_home}/MotionCorrection && module load matlab/R2016b && matlab -nodisplay -nosplash < {folder_home}/MotionCorrection/Scripts/{fname_mc_run_spon} &> {folder_home}/Logs/NoRMCorre_{fname_mc_run_spon}.log'
        cmd_run_timelapse = subprocess.run(cmd_normcorre_timelapse, shell=True, check=True)

        # if have movies for Markes, run NoRMCorre for Markes movies first
        if have_markes:
            # create a new script file for this dataset
            fname_mc_run_markes = f'ZZ_normcorre_pwrigid_python_markes_{name_dataset}_{brain}.m'
            normcorre_timelapse_run = open(os.path.join(
                folder_home, 'MotionCorrection/Scripts', fname_mc_run_markes), 'w')
            for line in script_lines:
                line = line.replace('FOLDER_HOME', folder_home)
                line = line.replace('FOLDER_CODE', folder_code)
                line = line.replace('MOVIETYPE', 'Markes')
                line = line.replace('RAWDATAFOLDER', folder_mc_raw)
                line = line.replace('DATASETNAME', name_dataset)
                line = line.replace('BRAINNAME', brain)
                line = line.replace('SAVEFOLDER', folder_mc_within)
                line = line.replace('NUMSTACKS', str(num_stacks[0]))
                normcorre_timelapse_run.write(line)
            normcorre_timelapse_run.close()
            # run the normcorre algorithm for Markes movie, save the logs
            datetime_now = datetime.datetime.now()
            log.write(f'{datetime_now}: running normcorre for Markes movie...\n')
            cmd_normcorre_timelapse = f'cd {folder_home}/MotionCorrection && module load matlab/R2016b && matlab -nodisplay -nosplash < {folder_home}/MotionCorrection/Scripts/{fname_mc_run_markes} &> {folder_home}/Logs/NoRMCorre_{fname_mc_run_markes}.log'
            cmd_run_timelapse = subprocess.run(cmd_normcorre_timelapse, shell=True, check=True)

        # then do motion correction for odorEvoked movies
        # create a new script file for this dataset
        fname_mc_run_odor = f'ZZ_normcorre_pwrigid_python_odor_{name_dataset}_{brain}.m'
        normcorre_timelapse_run = open(os.path.join(
            folder_home, 'MotionCorrection/Scripts', fname_mc_run_odor), 'w')
        for line in script_lines:
            line = line.replace('FOLDER_HOME', folder_home)
            line = line.replace('FOLDER_CODE', folder_code)            
            line = line.replace('MOVIETYPE', 'odorEvoked')
            line = line.replace('RAWDATAFOLDER', folder_mc_raw)
            line = line.replace('DATASETNAME', name_dataset)
            line = line.replace('BRAINNAME', brain)
            line = line.replace('SAVEFOLDER', folder_mc_within)
            line = line.replace('NUMSTACKS', str(num_stacks[0]))
            normcorre_timelapse_run.write(line)
        normcorre_timelapse_run.close()
        # run the normcorre algorithm for odor movie, save the logs
        datetime_now = datetime.datetime.now()
        log.write(f'{datetime_now}: running normcorre for odorEvoked movie...\n')
        cmd_normcorre_timelapse = f'cd {folder_home}/MotionCorrection && module load matlab/R2016b && matlab -nodisplay -nosplash < {folder_home}/MotionCorrection/Scripts/{fname_mc_run_odor} &> {folder_home}/Logs/NoRMCorre_{fname_mc_run_odor}.log'
        cmd_run_timelapse = subprocess.run(cmd_normcorre_timelapse, shell=True, check=True)

    ## 2.3. Motion correction between movies with CMTK
    # Between movies, there may be slow movement or slight deformation. Use CMTK to register the templates of each odor movie from NoRMCorre to the template of spontaneous movie, then apply the transformation to all timepoints in the movie.
    # run the normcorre algorithm for odor movie, save the logs
    datetime_now = datetime.datetime.now()
    log.write(f'{datetime_now}: running between-movies correction with CMTK...\n')
    for brain in brain_names:
        # subfolder to storage temporary files for this brain
        folder_CMTK_mc_brain = os.path.join(
            folder_home, 'MotionCorrection', folder_mc_temp, name_dataset, brain)
        # remove if already exists
        if os.path.exists(folder_CMTK_mc_brain):
            shutil.rmtree(folder_CMTK_mc_brain)
        os.makedirs(folder_CMTK_mc_brain)
        # folder to save the corrected movie
        folder_mc_between_save = os.path.join(
            folder_home, 'MotionCorrection', folder_mc_between, name_dataset, brain)
        if os.path.exists(folder_mc_between_save):
            shutil.rmtree(folder_mc_between_save)
        os.makedirs(folder_mc_between_save)

        # use the template generated from the spontaneous movie motion correction as the master template
        # use the template2 file, since the template file generated by NoRMCorre itself has some boundary artifacts
        fname_template = sorted(glob.glob(os.path.join(folder_home, 'MotionCorrection',
                                                       folder_mc_within, name_dataset, brain, 'spontaneous', '*template2.tif')))[0]
        # read the template
        template = ZZ_read_reshape_tiff(fname_template, movie_info_timelapse_mc_template)
        # squeeze the other dimensions to a 3D array
        template = np.squeeze(template)
        # what movie to register, glob all movies for this brain
        fnames_movie_template = sorted(glob.glob(os.path.join(
            folder_home, 'MotionCorrection', folder_mc_within, name_dataset, brain, 'odorEvoked', '*template2.tif'))) + sorted(glob.glob(os.path.join(folder_home, 'MotionCorrection', folder_mc_within, name_dataset, brain, 'Markes', '*template2.tif')))
        fnames_movie = sorted(glob.glob(os.path.join(
            folder_home, 'MotionCorrection', folder_mc_within, name_dataset, brain, 'odorEvoked', '*-mc.tif'))) + sorted(glob.glob(os.path.join(folder_home, 'MotionCorrection', folder_mc_within, name_dataset, brain, 'Markes', '*-mc.tif')))

        # template_header and movie_header are the same here
        template_header = nrrd_header_timelapse
        movie_header = nrrd_header_timelapse
        # cmtk option
        cmtk_option = cmtk_option_warp_mc

        # loop through each movie
        for movie_idx in range(len(fnames_movie)):
            # first read the template for each movie
            fname_movie_template = fnames_movie_template[movie_idx]
            movie_template = ZZ_read_reshape_tiff(
                fname_movie_template, movie_info_timelapse_mc_template)
            movie_template = np.squeeze(movie_template)
            # then read the mc movie, convert into float32 type, since int16 may be problematic
            fname_movie = fnames_movie[movie_idx]
            movie = ZZ_read_reshape_tiff(fname_movie, movie_info_timelapse_mc_template)
            movie = np.squeeze(movie)
            # save folder for intermediate files
            save_folder = os.path.join(folder_CMTK_mc_brain, f'odorEvoked_{movie_idx:05}')
            # run the registration
            movie_registered = ZZ_cmtk_register_movie(
                template, movie_template, movie, template_header, movie_header, cmtk_option, save_folder, rm_interm=rm_interm)
            # save corrected movie
            folder_mc_between_this = os.path.join(
                folder_mc_between_save, os.path.basename(fname_movie) + '.between.tif')
            ZZ_save_as_tiff(folder_mc_between_this, movie_registered)

    # 2.4. Motion correction for fastZ high-res movies
    for brain in brain_names:
                
        normcorre_fastZ = os.path.join(
            folder_matlab_templates, 'ZZ_normcorre_pwrigid_python_fastZ.m')
        # create a new script file for this brain
        fname_mc_run_fastZ = f'ZZ_normcorre_pwrigid_python_fastZ_{name_dataset}_{brain}.m'
        normcorre_fastZ_run = open(os.path.join(
            folder_home, 'MotionCorrection/Scripts', fname_mc_run_fastZ), 'w')
        for line in open(normcorre_fastZ, 'r').readlines():
            line = line.replace('FOLDER_HOME', folder_home)
            line = line.replace('FOLDER_CODE', folder_code)            
            line = line.replace('RAWDATAFOLDER', folder_mc_raw)
            line = line.replace('DATASETNAME', name_dataset)
            line = line.replace('BRAINNAME', brain)
            line = line.replace('SAVEFOLDER', folder_mc_within)
            line = line.replace('NUMSTACKS', str(num_stacks[1]))
            normcorre_fastZ_run.write(line)
        normcorre_fastZ_run.close()

        # run the normcorre algorithm for spontaneous activity, save the logs
        datetime_now = datetime.datetime.now()
        log.write(f'{datetime_now}: running NoRMCorre for fastZ high-resolution...\n')
        cmd_normcorre_fastZ = f'cd {folder_home}/MotionCorrection && module load matlab/R2016b && matlab -nodisplay -nosplash < {folder_home}/MotionCorrection/Scripts/{fname_mc_run_fastZ} &> {folder_home}/Logs/{fname_mc_run_fastZ}.log'
        cmd_run_fastZ = subprocess.run(cmd_normcorre_fastZ, shell=True, check=True)

    # close the log files

    log.close()


def main():
    option_ds = {'name_dataset': '240111', 
                 'brain_names': ['240111_Camphor_1U_F1'], 
                 'file_ranges_spontaneous': [[1]], 
                 'file_ranges_fastZ': [range(48, 49)], 
                 'num_stacks': [22, 161], 
                 'odor_pattern_strings_single': '100A_50A_10A_100B_10B_100C_10C_100D_100F_100G_100H_100J_100L_100M_100N_10N_100O_100P_100Q_100R_100S_100T_100U_100A_50A_10A_100B_10B_100C_10C_100D_100F_100G_100H_100J_100L_100M_100N_10N_100O_100P_100Q_100R_100S_100T_100U_', 
                 'channel_odorant_single': {'100A': 'camphor-1', '100B': 'camphor-3', '100C': '1-8-cineole-2', '100D': 'sulcatone', '100F': 'alpha-terpinene', '100G': 'beta-ocimene', '100H': 'p-cymene', '100J': 'methy-butyl-acetate', '100L': 'paraffin_oil', '100M': 'beta-myrcene', '100N': '4-terpineol-2', '100O': '3-carene', '100P': 'benzaldehyde-methyloctane', '100Q': 'phenol-cresol', '100R': 'nonanal-furfural', '100S': '1-octen-3-ol', '100T': 'heptyl-acetate-ethylhexanol', '100U': 'oxoisophorone-acetophenone', '10A': 'camphor-2', '10B': 'camphor-4', '10C': '1-8-cineole-3', '10N': '4-terpineol-3', '50A': 'camphor-1.5'}, 
                 'channel_solvent_single': {'100A': 'paraffin', '100B': 'paraffin', '100C': 'paraffin', '100D': 'paraffin', '100F': 'paraffin', '100G': 'paraffin', '100H': 'paraffin', '100J': 'paraffin', '100L': 'paraffin', '100M': 'paraffin', '100N': 'paraffin', '100O': 'paraffin', '100P': 'paraffin', '100Q': 'paraffin', '100R': 'paraffin', '100S': 'paraffin', '100T': 'paraffin', '100U': 'paraffin', '10A': 'paraffin', '10B': 'paraffin', '10C': 'paraffin', '10N': 'paraffin', '50A': 'paraffin'}, 
                 'file_ranges_odorEvoked_single': [range(2, 48)], 
                 'odor_pattern_strings_markes': [''], 
                 'channel_odorant_markes': {}, 
                 'channel_solvent_markes': {}, 
                 'file_ranges_odorEvoked_markes': [[]], 
                 'fname_markes_timestamp': '20190827_HumanDosageBrp_MarkesTimestamp.txt', 
                 'datetime_format': '%Y-%m-%d %H:%M:%S.%f', 
                 'date_experiment': '2022-07-13', 
                 'AL_loc': ['lower_AL']}
    
    analysis_pipeline_orco_part1(option_ds, have_markes=False)

if __name__ == '__main__':
    main()