# a set of functions that warp CMTK for registration
import numpy as np
import tifffile
import os
import shutil
import subprocess
import nrrd
import datetime
from ZZ_tiff_file import *
from collections import OrderedDict
from joblib import Parallel, delayed
import glob

def ZZ_CMTK_rigid(refbrain, float_images, CMTK_folder, options):
    """use the CTMK registration tool
    refbrain is full path to the ref brain
    float_images is a list containing full paths to the floating images
    CMTK_folder is where to perform the registration
    options include parameters for registration
    return a list of registration output folders"""
    # folder for cmtk bin
    cmtk_bin_folder = '/mnt/cup/people/lw1397/Fiji.app/bin/cmtk'
    # make sub-folders for output if not exists
    if not os.path.exists(CMTK_folder + '/' + 'Registration'):
        os.mkdir(CMTK_folder + '/' + 'Registration')
    output_folders = []
    # loop through each image
    refbrain_fname = refbrain.split('/')[-1]
    for fn in float_images:
        fn_fname = fn.split('/')[-1]
        # name for output registration
        out_name = CMTK_folder + '/' + 'Registration/' + refbrain_fname.replace('.nrrd','') + '_' + fn_fname.replace('.nrrd','') 
    #     cmtk_command= f"cd {CMTK_folder} && {cmtk_bin_folder}/munger -b {cmtk_bin_folder} -i -a -r {options['reformat']} -T {options['threads']} -l af  -C 4 -G 32 -R 3 -A '--exploration {options['exploration']} --accuracy {options['accuracy']}' -s refbrain/{refbrain_fname} images"
        # use registration and reformatx directly
        cmtk_command= f"cd {CMTK_folder} && export CMTK_NUM_THREADS={options['threads']} && {cmtk_bin_folder}/registration --initxlate --dofs {options['dofs']} --coarsest 4 --exploration {options['exploration']} --accuracy {options['accuracy']} --outlist {out_name} {refbrain} {fn}"
        # run the command, report error is exist status is non-zero
        print(f'Doing registration for {fn_fname}...\ncommands: {cmtk_command}')
        cmd = subprocess.run(cmtk_command, shell=True, check=True)
        output_folders.append(out_name)
    return output_folders


def ZZ_CMTK_warp(refbrain, CMTK_folder, options):
    """use the CTMK registration tool
    refbrain is full path to the ref brain
    CMTK_folder is where to perform the registration, should have a subfolder named images
    that contains all the floating images need to be registered
    options include parameters for registration
    return a list of registration output folders"""
    # folder for cmtk bin
    cmtk_bin_folder = '/mnt/cup/people/lw1397/Fiji.app/bin/cmtk'
    # make sub-folders for output if not exists
    folder_register = os.path.join(CMTK_folder, 'Registration')
    if not os.path.exists(folder_register):
        os.mkdir(folder_register)
    # loop through each image
    # munger to do registration and reformat
    cmtk_command= f"cd {CMTK_folder} && {cmtk_bin_folder}/munger -b {cmtk_bin_folder} -i -a -w -r {options['reformat']} -T {options['threads']} -l af -X {options['exploration']} -C {options['coarsest']} -G {options['grid']} -R {options['refine']} -A '--accuracy {options['accuracy']}' -W '--accuracy {options['accuracy']} --energy-weight {options['energy_weight']} --rigidity-weight {options['rigidity_weight']}' -s {refbrain} images"
    # run the command, report error is exist status is non-zero
    print(f'Doing registration...\ncommands: {cmtk_command}')
    cmd = subprocess.run(cmtk_command, shell=True, check=True)
    return folder_register

    
def ZZ_movie_to_nrrd(movie, save_prefix, save_folder, header):
    """split 4D (xyzt) movie into separate nrrd volume files"""
    if os.path.exists(save_folder):
        shutil.rmtree(save_folder)
    os.makedirs(save_folder)
    # loop-through each time-point
    for t in range(movie.shape[3]):
        v = movie[:,:,:,t].astype('single')
        file_name_temp = save_folder + '/' + save_prefix + 't'+'{:04}'.format(t+1)+'_01.nrrd'
#         print(file_name_temp)
        nrrd.write(file_name_temp, v, header)
    return save_folder


def ZZ_CMTK_reformat(refbrain_name, trans_folder, float_folder, output_folder, options, interpolation='cubic'):
    """use the reformatx tool to apply transformation to floating images
    trans_folder is a list, support concatenation of multiple transformation"""
    # folder for cmtk bin
    cmtk_bin_folder = '/mnt/cup/people/lw1397/Fiji.app/bin/cmtk'
    # make the output_folder
    if os.path.exists(output_folder):
        shutil.rmtree(output_folder)
    os.makedirs(output_folder)
    # make a string to represent transformation cascade
    trans_cascade = ' '.join(trans_folder)
    # loop-through each float image
    float_images_all = sorted(os.listdir(float_folder))
    # exclude the rgb volume
    float_images = [fn for fn in float_images_all if 'rgb' not in fn]
#     float_images = sorted(glob.glob(os.path.join(float_folder, '*.nrrd')
    counter = 0
    def ZZ_CMTK_reformat_single(fi):
        float_image = float_folder + '/' + fi
        out_image_name = output_folder + '/' + fi.replace('.nrrd','rf.nrrd')
        # use cubic interpolation 
        cmtk_command= f"export CMTK_NUM_THREADS={options['threads']} && {cmtk_bin_folder}/reformatx --pad-out 0 --{interpolation} -o {out_image_name} --floating {float_image} {refbrain_name} {trans_cascade}"
    #         print(cmtk_command)
        print(f'Doing reformatx for {fi}...\n')
        print(cmtk_command)
        cmd = subprocess.run(cmtk_command, shell=True, check=True)
#     for fi in float_images:
#         ZZ_CMTK_reformat_single(fi)
    # run reformat parallelly
    with Parallel(n_jobs=22) as parallel:
        parallel(delayed(ZZ_CMTK_reformat_single)(fi) for fi in float_images)
        
        
def ZZ_stitch_reformatted(source_folder, volume_size):
    """stitch reformatted time-lapse back to a numpy 4D xyzt array"""
    registered_images = sorted(os.listdir(source_folder))
    registered_timelapse = np.zeros([*volume_size, len(registered_images)])
    # loop through each file
    for idx, ri in enumerate(registered_images):
        file_name_temp = source_folder + '/' + ri
        v_temp, header = nrrd.read(file_name_temp)
        registered_timelapse[:,:,:,idx] = np.squeeze(v_temp)
    return registered_timelapse


def ZZ_movie_register_highres_template(mc_folder, file_name_base, mc_fastZ_folder, file_name_fastZ, trans_folder_template, refbrain_name, CMTK_folder_base):
    """register moview to high-res template"""
    brain_name = file_name_base.split('_')[0]
    session_name = file_name_base.split('_')[0] + 'f' + file_name_base.split('_')[1]
    file_name_timelapse_template = mc_folder + '/' +  file_name_base + '-template.tif'
    file_name_timelapse = mc_folder + '/' +  file_name_base + '-mc.tif'
    brain_folder = CMTK_folder_base + str(datetime.date.today()).replace('-','') + '_' + brain_name
    CMTK_folder = brain_folder + '/' + session_name
    # make folder for registration
    if not os.path.exists(brain_folder):
        os.mkdir(brain_folder)
    if os.path.exists(CMTK_folder):
        shutil.rmtree(CMTK_folder)
    os.mkdir(CMTK_folder)
    # read tiff files, then save as nrrd in subfolders in CMTK_folder
    # save fastZ in nrrd format in fastZ subfolder if it doesn't exist
    fastZ_nrrd = brain_folder + '/fastZ/' + file_name_fastZ.split('_')[0] + 'fastZ_01.nrrd'
    nrrd_header_highres = OrderedDict([('type', 'uint8'), ('encoding', 'raw'), 
                                   ('endian', 'big'), ('dimension', 3), 
                                   ('sizes', np.array([256, 256,  83])), ('space dimension', 3), 
                                   ('space directions', np.array([[0.4880, 0, 0],[0, 0.4880, 0],[0, 0, 1.25]])), 
                                   ('space units', ['"microns"', '"microns"', '"microns"'])])
    if not os.path.isfile(fastZ_nrrd):
        volume_info_fastZ = {'xRes':256, 'yRes':256, 'nChannel':1, 'zStacks': 83, 'dtype':np.dtype('uint8')}
        fastZ = ZZ_read_reshape_tiff(mc_fastZ_folder + '/' + file_name_fastZ, volume_info_fastZ)
        fastZ = np.squeeze(fastZ)
        os.mkdir(brain_folder + '/fastZ')
        nrrd.write(fastZ_nrrd, fastZ, nrrd_header_highres)
    # save the normcorre template in nrrd format in template subfolder
    os.mkdir(CMTK_folder + '/template')
    volume_info = {'xRes':128, 'yRes':128, 'nChannel':1, 'zStacks': 22, 'dtype':np.dtype('float32')}
    template = ZZ_read_reshape_tiff(file_name_timelapse_template, volume_info)
    template = np.squeeze(template)
    template_nrrd = CMTK_folder + '/template/' + file_name_base.split('_')[0] + 'f' + file_name_base.split('_')[1] +'template'+ '_01.nrrd'
    nrrd_header_timelapse_template = OrderedDict([('type', 'double'), ('encoding', 'raw'), 
                                       ('endian', 'big'), ('dimension', 3), 
                                       ('sizes', np.array([128, 128,  22])), ('space dimension', 3), 
                                       ('space directions', np.array([[0.9760, 0, 0],[0, 0.9760, 0],[0, 0, 5]])), 
                                       ('space units', ['"microns"', '"microns"', '"microns"'])])
    nrrd.write(template_nrrd, template, nrrd_header_timelapse_template)
    # read the motion-correcte movie then split the into separate nrrd files
    volume_info_timelapse = {'xRes':128, 'yRes':128, 'nChannel':1, 'zStacks': 22, 'dtype':np.dtype('int16')}
    v2 = ZZ_read_reshape_tiff(file_name_timelapse, volume_info_timelapse)
    movie = np.squeeze(v2)
    nrrd_header_timelapse = OrderedDict([('type', 'float'), ('encoding', 'raw'), 
                                   ('endian', 'big'), ('dimension', 3), 
                                   ('sizes', np.array([128, 128,  22])), ('space dimension', 3), 
                                   ('space directions', np.array([[0.9760, 0, 0],[0, 0.9760, 0],[0, 0, 5]])), 
                                   ('space units', ['"microns"', '"microns"', '"microns"'])])
    folder_name_timelapse = ZZ_movie_to_nrrd(movie, file_name_timelapse, CMTK_folder, nrrd_header_timelapse)
    # register time-lapse template to the fastZ high-res 
    # options for registration
    options = OrderedDict([('overwrite', True), ('dofs', '6,9'),
                          ('exploration', 16), ('accuracy', 0.4),
                          ('threads', 22), ('reformat','01')])
    trans_folders = ZZ_CMTK_rigid(fastZ_nrrd, [template_nrrd], CMTK_folder, options)
    # apply transformation to each time-point including the registration from high-res to high-res template
    trans_folder2 = [trans_folders[0], trans_folder_template]
    output_folder2 = folder_name_timelapse + '_reformatted2'
    ZZ_CMTK_reformat(refbrain_name, trans_folder2, folder_name_timelapse, output_folder2, options)
    # stitch the registered time-lapse together
    registered_timelapse2 = ZZ_stitch_reformatted(output_folder2, nrrd_header_highres['sizes'])
    registered_timelapse2 = registered_timelapse2.astype('int16')
    final_file = CMTK_folder + '/' + file_name_base + '-reg.tif'
    xRes = registered_timelapse2.shape[0]
    data = np.reshape(registered_timelapse2, (xRes, xRes, -1), order='F')
    data = np.swapaxes(data, 0, -1)
    tif.imsave(final_file, data, compress=1)
    # return the filename of output movie
    return final_file


def ZZ_cmtk_movie_highres_template(brain_name, file_timelapse_mc, movie_info, nrrd_header_timelapse, file_fastZ, nrrd_header_fastZ, folder_cmtk_fastZ, refbrain, folder_save, options_cmtk):
    """register moview to high-res template"""
    # get the name for the odor session
    session_name =  file_timelapse_mc.split('/')[-1] 
    session_name = session_name.split('_')[0] + 'f' + session_name.split('_')[1]
    # make folder for saving the registration
    if os.path.exists(folder_save):
        shutil.rmtree(folder_save)
    os.makedirs(folder_save)
    # read the motion-correcte movie then split the into separate nrrd files
    m = ZZ_read_reshape_tiff(file_timelapse_mc, movie_info)
    movie = np.squeeze(m)
    folder_name_timelapse = ZZ_movie_to_nrrd(movie, session_name, folder_save + '/' + session_name, nrrd_header_timelapse)
    # save timelapse template as nrrd
    movie_info_timelapse_template = movie_info.copy()
    movie_info_timelapse_template['dtype'] = 'float32'
    volume_timelapse_template = ZZ_read_reshape_tiff(file_timelapse_mc.replace('-mc','-template'), movie_info_timelapse_template)
    volume_timelapse_template = ZZ_scanimage_to_uint8(volume_timelapse_template, 0.001, 0.05)
    if os.path.exists(folder_save + '/template/'):
        shutil.rmtree(folder_save + '/template/')
    os.makedirs(folder_save + '/template/')
    file_timelapse_template = folder_save + '/template/' + session_name + 'template_01.nrrd'
    nrrd.write(file_timelapse_template, np.squeeze(volume_timelapse_template), nrrd_header_timelapse)
    # register time-lapse template to the fastZ 
    trans_folders = ZZ_CMTK_rigid(file_fastZ, [file_timelapse_template], folder_save, options_cmtk)
    # apply transformation to each time-point including the registration from high-res to high-res template
    trans_folder2 = [trans_folders[0], folder_cmtk_fastZ]
    output_folder2 = folder_save + '/' + session_name + '_reformatted2'
    options_cmtk_reformat = options_cmtk
    ZZ_CMTK_reformat(refbrain, trans_folder2, folder_name_timelapse, output_folder2, options_cmtk_reformat)
    # stitch the registered time-lapse together
    registered_timelapse2 = ZZ_stitch_reformatted(output_folder2, nrrd_header_fastZ['sizes'])
    registered_timelapse2 = registered_timelapse2.astype('int16')
    final_file = folder_save + '/' + session_name + '-reg.tif'
    xRes = registered_timelapse2.shape[0]
    data = np.reshape(registered_timelapse2, (xRes, xRes, -1), order='F')
    data = np.swapaxes(data, 0, -1)
    tif.imsave(final_file, data, compress=1)
    return final_file


def ZZ_cmtk_register_movie(template, movie_template, movie, template_header, movie_header, cmtk_option, save_folder, rm_interm = True):
# a general function to register movie_template to a template, then apply the transformation to movie
# template, movie_template should be 3D array, movie should be 4D array, whose xyz dimension should be same as movie_template, i.e. they have the same header info
# two header contains the volume and movie information
# cmtk_option is parameters passed to cmtk algorithm
# save_folder is where to store the intermediate files
# rm_interm is to remove the intermediate files or not
# return the registered movie

    # check the save folder, delete if already exists
    if os.path.exists(save_folder):
        shutil.rmtree(save_folder)
    os.makedirs(save_folder)
    os.mkdir(os.path.join(save_folder, 'refbrain'))
    os.mkdir(os.path.join(save_folder, 'images'))
    os.mkdir(os.path.join(save_folder, 'images2'))
    os.mkdir(os.path.join(save_folder, 'images2_rf'))

    ## convert file formats and put into corresponding subfolders in save_folder
    # convert template into nrrd, put it into the refbrain folder
    nrrd.write(os.path.join(save_folder,'refbrain', 'template_01.nrrd'), template, template_header)
    # convert movie_template into nrrd, put it into the images folder
    nrrd.write(os.path.join(save_folder,'images', 'movieTemplate_01.nrrd'), movie_template, movie_header)
    # loop through each time point of movie, save individual volume as nrrd in the images2 folder
    for tp in range(movie.shape[3]):
        volume_now = movie[:,:,:,tp]
        nrrd.write(os.path.join(save_folder,'images2', f'timePoint{tp:05}_01.nrrd'), volume_now, movie_header)

    # register movie_template to template
    fname_fbran = os.path.join(save_folder,'refbrain', 'template_01.nrrd')
    folder_register = ZZ_CMTK_warp(fname_fbran, save_folder, cmtk_option)

    # apply the transformation to all time points
    folder_transform = glob.glob(os.path.join(folder_register, 'warp','*list'))
    ZZ_CMTK_reformat(fname_fbran, folder_transform, os.path.join(save_folder,'images2'), os.path.join(save_folder, 'images2_rf'), cmtk_option) 

    # stitch reformated volumes back to a movie
    output_movie = ZZ_stitch_reformatted(os.path.join(save_folder, 'images2_rf'), movie.shape[0:3])

    # clean up if needed 
    if rm_interm:
        shutil.rmtree(save_folder)
        
    # return the corrected movie
    return output_movie


def ZZ_cmtk_mkfolder(folder_home):
    """make the typical folder structure for CMTK munger
    folder_home is where to perform the registration"""
    # check existence
    if os.path.exists(folder_home):
        shutil.rmtree(folder_home)
    os.makedirs(folder_home)
    # make sub-folders for CMTK registration
    folder_images = os.path.join(folder_home, 'images')
    folder_refbrain = os.path.join(folder_home, 'refbrain')
    os.makedirs(folder_images)
    os.makedirs(folder_refbrain)
    return [folder_images, folder_refbrain]


def ZZ_cmtk_reformat_movie(movie, movie_header, fname_refbrain, trans_folders, save_folder, cmtk_options, rm_interm = True):
    """apply the transformation in trans_folder to a movie
   movie is a 4D array, movie_header is the nrrd header info for movie
   fname_refbrain is the name for the final refbrain
   trans_folders is a list containing the CMTK transformation"""
    # check the save folder, delete if already exists
    if os.path.exists(save_folder):
        shutil.rmtree(save_folder)
    # make folder structure for cmtk
    [folder_images_working, folder_refbrain_working] = ZZ_cmtk_mkfolder(save_folder)
    folder_images_rf_working = os.path.join(save_folder, 'images_rf')
    os.makedirs(folder_images_rf_working)
    fname_refbrain_working = os.path.join(folder_refbrain_working, os.path.basename(fname_refbrain))
    shutil.copy(fname_refbrain, fname_refbrain_working)
    # loop through each time point of movie, save individual volume as nrrd in the images2 folder
    for tp in range(movie.shape[3]):
        volume_now = movie[:,:,:,tp]
        nrrd.write(os.path.join(save_folder,'images', f'timePoint{tp:05}_01.nrrd'), volume_now, movie_header)
    # apply the transformation
    ZZ_CMTK_reformat(fname_refbrain_working, trans_folders, folder_images_working, folder_images_rf_working, cmtk_options)
    refbrain = nrrd.read(fname_refbrain)
    # stitch reformated volumes back to a movie
    output_movie = ZZ_stitch_reformatted(folder_images_rf_working, refbrain[0].shape)
    # clean up if needed 
    if rm_interm:
        shutil.rmtree(save_folder)
    # return the corrected movie
    return output_movie


def ZZ_cmtk_register_volumes(template_volume, float_volume, template_header, float_header, cmtk_option, save_folder):
# a general function to float volume to the template volume
# where both template and float are numpy arrays

    # check the save folder, delete if already exists
    if os.path.exists(save_folder):
        shutil.rmtree(save_folder)
    os.makedirs(save_folder)
    os.mkdir(os.path.join(save_folder, 'refbrain'))
    os.mkdir(os.path.join(save_folder, 'images'))

    ## convert file formats and put into corresponding subfolders in save_folder
    # convert template into nrrd, put it into the refbrain folder
    nrrd.write(os.path.join(save_folder,'refbrain', 'template_01.nrrd'), template_volume, template_header)
    # convert movie_template into nrrd, put it into the images folder
    nrrd.write(os.path.join(save_folder,'images', 'float_01.nrrd'), float_volume, float_header)

    # register float_volume to template_volume
    fname_fbran = os.path.join(save_folder,'refbrain', 'template_01.nrrd')
    folder_register = ZZ_CMTK_warp(fname_fbran, save_folder, cmtk_option)
    folder_transform = glob.glob(os.path.join(folder_register, 'warp','*list'))

    # read the registered volume
    fname_registered = glob.glob(os.path.join(save_folder, 'reformatted', '*.nrrd'))[0]
    registered_volume = nrrd.read(fname_registered)
        
    # return the registered movie and folder of registration
    return [registered_volume, folder_transform]