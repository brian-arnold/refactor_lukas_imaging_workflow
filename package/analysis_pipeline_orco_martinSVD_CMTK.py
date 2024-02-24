## Zhilei Zhao
## run SVD as a preprocessing step 

import os
import shutil
import glob
import subprocess

def analysis_pipeline_orco_martinSVD_CMTK(option_ds):
    # unpack the dataset-specific parameters
    experiment_name = option_ds['experiment_name']
    brains = option_ds['brains']
    datasets = option_ds['datasets']
    file_ranges_odorEvoked_single = option_ds['file_ranges_odorEvoked_single']
    file_ranges_odorEvoked_markes = option_ds['file_ranges_odorEvoked_markes']
    
    # shared parameters
    folder_home = '/jukebox/mcbride/bjarnold/refactor_lukas_imaging_workflow/data'
    folder_code = '/jukebox/mcbride/bjarnold/refactor_lukas_imaging_workflow/package'
    folder_code_matlab = '/jukebox/mcbride/bjarnold/refactor_lukas_imaging_workflow/package/matlab/part5'
    folder_martin = f'{folder_home}/MartinStrauch_hybridNewMerging'
    
    # copy motion corrected movies to the working folder
    # for bi in range(len(brains)):
    #     folder_mc = os.path.join(folder_home, 'Registration', 'RegisteredMovies', datasets[bi], brains[bi])
    #     folder_save = os.path.join(folder_home, 'MartinStrauch', 'MotionCorrectedData', experiment_name, brains[bi], 'odorEvoked')
    #     if os.path.exists(folder_save):
    #         shutil.rmtree(folder_save)
    #     os.makedirs(folder_save)
    #     # loop through each file index
    #     for fi in list(file_ranges_odorEvoked_single[bi]) + list(file_ranges_odorEvoked_markes[bi]):
    #         fn = glob.glob(os.path.join(folder_mc, f'{brains[bi]}_{fi:05d}_*.registered.tif'))[0]
    #         fn_save = os.path.join(folder_save, os.path.basename(fn))
    #         shutil.copy(fn, fn_save)

    script_template = os.path.join(folder_code, 'ZZ_runSVD_scotty_template_CMTK.m')
    # save the script of this experiment to 
    script_save = os.path.join(folder_martin, 'ZZ_Scripts', f'{experiment_name}_runSVD.m')
    # remove the script if it exists
    if os.path.exists(script_save):
        os.remove(script_save)
    f_out = open(script_save, 'w')
    f_in = open(script_template, 'r')
    lines = f_in.readlines()
    # replace the dataset_specific parameter
    line_brains = ' ,'.join([f'\'{a}\'' for a in brains])
    line_datasets = ' ,'.join([f'\'{a}\'' for a in datasets])
    for line in lines:
        line = line.replace('FOLDER_HOME', folder_home)
        line = line.replace('FOLDER_CODE_MATLAB', folder_code_matlab)
        line = line.replace('EXPERIMENTNAME', experiment_name)
        line = line.replace('BRAINS', line_brains)
        line = line.replace('DATASETS', line_datasets)
        line = line.replace('FOLDERWORKING', folder_martin)
        f_out.write(line)
    f_in.close()
    f_out.close()
    # make the command
    cmd = f'cd {folder_martin} && module load matlab/R2018a && matlab -nodisplay -nosplash < {script_save}'
    cmd_run = subprocess.run(cmd, shell=True, check=True)





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
                 'channel_odorant_single': {'100A': 'camphor-1', '100B': 'camphor-3', '100C': '1-8-cineole-2', '100D': 'sulcatone', '100F': 'alpha-terpinene', '100G': 'beta-ocimene', '100H': 'p-cymene', '100J': 'methy-butyl-acetate', '100L': 'paraffin_oil', '100M': 'beta-myrcene', '100N': '4-terpineol-2', '100O': '3-carene', '100P': 'benzaldehyde-methyloctane', '100Q': 'phenol-cresol', '100R': 'nonanal-furfural', '100S': '1-octen-3-ol', '100T': 'heptyl-acetate-ethylhexanol', '100U': 'oxoisophorone-acetophenone', '10A': 'camphor-2', '10B': 'camphor-4', '10C': '1-8-cineole-3', '10N': '4-terpineol-3', '50A': 'camphor-1.5'}, 
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


    analysis_pipeline_orco_martinSVD_CMTK(option_ds)


if __name__ == '__main__':
    main()