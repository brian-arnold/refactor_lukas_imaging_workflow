# Zhilei Zhao, McBride Lab, Princeton University
# ## 0. Goal
# this is part 3 of the pipeline: registered motion corrected movies to a template
# note that the 2P template is of a 'upper_AL'

# import general modules
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
import ZZ_calcium_responseNew
from ZZ_calcium_responseNew import *
importlib.reload(ZZ_calcium_responseNew)


def analysis_pipeline_orco_partRegisterMovie(option_ds, have_markes=False):
    ## inputs
    # option_ds: dataset-specific parameters
    # have_markes: if we have Markes puffs for this brain 
    
    # unpack the dataset-specific parameters
    name_dataset = option_ds['name_dataset']
    brain_names = option_ds['brain_names']
    AL_loc = option_ds['AL_loc']

   # folder setting
    # home folder that serves as the base folder for all analysis
    folder_home = '/jukebox/mcbride/bjarnold/refactor_lukas_imaging_workflow/data'
    # sub-foler in MotionCorrection to store the data before motion correction
    folder_mc_raw = "BeforeCorrection"
    # sub-foler in MotionCorrection to store the with-in movie motion correction results
    folder_mc_within = 'CorrectedDataWithin'
    # sub-foler in MotionCorrection to store the between movie motion correction results
    folder_mc_between = 'CorrectedDataBetween'
    # sub-foler in MotionCorrection to store the temporary results for motion correction
    folder_mc_temp = 'Temp'
    # folder to store the response
    folder_response = 'Response'
    # sub-folder to store the raw peak df/f
    folder_response_rawpeak = 'RawPeakDff'
    # sub-folder in 'Response' to store the registered dff
    folder_response_registered = 'RegisteredDff'
    # what sub-folders in 'RawPeakDff' to apply CMTK transformation
    folders_to_apply_CMTK = ['MovieDff', 'AveragedDff', 'SubstractedDff']
    # sub-folder to store the Registration results
    folder_registration = 'Registration'
    # sub-folder in Registration to store the templates
    folder_registration_templates = 'Templates'
    # sub-folder in Registration to store registered movies
    folder_registration_movies = 'RegisteredMovies'
    # sub-folder in Template to store the 2P template
    folder_registration_templates_2p = 'Make2PtemplateNew'
    # what specific 'make' of 2P template to use (since 2P template keeps updated with more data come in)
    template2P_make = '20200123_humanRepresentation'
    # sub-folder in Registration to store CMTK results
    folder_registration_cmtk = 'CMTKfiles'

    # movie information
    # information for time-lapse and fastZ movies
    # don't specify the time points, since it may be slightly different for different movies
    movie_info_timelapse_mc_template = {'xRes': 128, 'yRes': 128, 'nChannel': 1,
                                        'zStacks': 22 + 2, 'dtype': np.dtype('float32'), 'zoomin': 10}
    # header file for nrrd files
    nrrd_header_timelapse = OrderedDict([('type', 'float'), ('encoding', 'raw'),
                                         ('endian', 'big'), ('dimension', 3),
                                         ('sizes', np.array(
                                             [128, 128, 22 + 2])), ('space dimension', 3),
                                         ('space directions', np.array(
                                             [[0.908, 0, 0], [0, 0.807, 0], [0, 0, 4]])),
                                         ('space units', ["microns", "microns", "microns"])])

    # algorithm parameters
    # options for CMTK, see munger --help for explanation
    cmtk_option_warp_mc = OrderedDict([('overwrite', True), ('energy_weight', 5e-1),
                                       ('rigidity_weight', 0),
                                       ('exploration', 16), ('accuracy', 0.4),
                                       ('threads', 22), ('reformat', '01'),
                                       ('coarsest', 4), ('grid', 32),
                                       ('refine', 3)])


    cmtk_option_warp = OrderedDict([('overwrite', True), ('energy_weight', 3e-1),
                                                      ('rigidity_weight', 0),
                                                      ('exploration', 16), ('accuracy', 0.4),
                                                      ('threads', 22), ('reformat', '01'),
                                                      ('coarsest', 4), ('grid', 32),
                                                      ('refine', 3)])

    # remove intermediate files or not (CMTK regitration files)
    # this automatic removing function is not implementated yet, need to remove manually
    rm_interm = False

    ## read in the 2P template
    # templates formed by registrating spontaneous activity templates
    # default side is upper AL, need to flip if imaged lower AL
    fname_2ptemplate = sorted(glob.glob(os.path.join(folder_home, 'Registration', folder_registration_templates, folder_registration_templates_2p, template2P_make, 'refbrain', '*.nrrd')))[-1]
    twop_template, twop_header = nrrd.read(fname_2ptemplate)
    
    # loop through brains
    for bi, brain in enumerate(brain_names):
        # folder to store the CMTK files for this brain
        folder_CMTK_brain = os.path.join(folder_home, 'Registration', folder_registration_cmtk, name_dataset, brain)
        # remove if already exists
        if os.path.exists(folder_CMTK_brain):
            shutil.rmtree(folder_CMTK_brain)
        os.makedirs(folder_CMTK_brain)
        
        # read the spontaneous template of the brain
        fname_spon = sorted(glob.glob(os.path.join(folder_home, 'MotionCorrection',
                                                   folder_mc_within, name_dataset, brain, 'spontaneous', '*template2.tif')))[0]
        spon = ZZ_read_reshape_tiff(fname_spon, movie_info_timelapse_mc_template)
        # squeeze the other dimensions to a 3D array
        spon = np.squeeze(spon)
        # if not default side (upper AL), needs to flip vertically
        if AL_loc[bi]!='upper_AL':
            spon = np.flip(spon, 1)
        
        # save folder for intermediate files
        save_folder = os.path.join(folder_CMTK_brain, 'spontaneousToTemplate')
        if os.path.exists(save_folder):
            shutil.rmtree(save_folder)
        os.makedirs(save_folder)

        # run the registration
        # header is the same 
        spon_registered_volume, spon_registered_folder = ZZ_cmtk_register_volumes(twop_template, spon, nrrd_header_timelapse, nrrd_header_timelapse, cmtk_option_warp, save_folder)
        
        ## apply to motion corrected movies
        # folder to save registered movies 
        folder_final = os.path.join(folder_home, 'Registration', folder_registration_movies, name_dataset, brain)
        if os.path.exists(folder_final):
            shutil.rmtree(folder_final)
        os.makedirs(folder_final)

        # loop through the motion corrected movies
        folder_mc = os.path.join(folder_home, 'MotionCorrection', folder_mc_between, name_dataset, brain)
        fns_mc = sorted(glob.glob(os.path.join(folder_mc, f'{brain}*between.tif')))
        for fidx in range(len(fns_mc)):
            fn = fns_mc[fidx]
            movie = ZZ_read_reshape_tiff(fn, movie_info_timelapse_mc_template)
            movie = np.squeeze(movie)
            # folder to save intermediate files
            folder_interm = os.path.join(save_folder, 'odorEvoked', f'odorEvoked_{fidx:04d}')
            if os.path.exists(folder_interm):
                shutil.rmtree(folder_interm)
            os.makedirs(folder_interm)
            folder_interm_rf = os.path.join(save_folder, 'odorEvoked', f'odorEvoked_{fidx:04d}_rf')
            if os.path.exists(folder_interm_rf):
                shutil.rmtree(folder_interm_rf)
            os.makedirs(folder_interm_rf)
            # loop through each time point of movie, save individual volume as nrrd in the images2 folder
            for tp in range(movie.shape[3]):
                volume_now = movie[:,:,:,tp]
                # need to flip if not default side
                if AL_loc[bi]!='upper_AL':
                    volume_now = np.flip(volume_now, 1)
                nrrd.write(os.path.join(folder_interm, f'timePoint{tp:05}_01.nrrd'), volume_now, nrrd_header_timelapse)
            # apply the transformation to all time points
            fname_ref = glob.glob(os.path.join(save_folder, 'refbrain', '*.nrrd'))[0]
            ZZ_CMTK_reformat(fname_ref, spon_registered_folder, folder_interm, folder_interm_rf, cmtk_option_warp) 
            # stitch reformated volumes back to a movie
            output_movie = ZZ_stitch_reformatted(folder_interm_rf, movie.shape[0:3])
            # save registered movie
            fn_final = os.path.basename(fn).replace('between', 'registered')
            ZZ_save_as_tiff(os.path.join(folder_final, fn_final), output_movie)

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

    analysis_pipeline_orco_partRegisterMovie(option_ds, have_markes=False)

if __name__ == '__main__':
    main()