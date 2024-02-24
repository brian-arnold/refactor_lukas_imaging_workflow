clear meta_data;
% default side for AL
default_AL = 'upper_AL';

experiment_name = 'EXPERIMENTNAME';
brains = {BRAINS};
odor_names = ODOR_NAMES;
file_names = FILE_NAMES;
flip_AL = {FLIP_AL};

meta_data.animal_names = brains';
meta_data.odor_names = odor_names;
meta_data.file_names = file_names;
meta_data.flip_AL = flip_AL';

% save meta_data file
folder_home = 'FOLDER_HOME';
meta_file = strcat(experiment_name, '_meta_data.mat');
fn_meta = fullfile(folder_home, 'MartinStrauch_hybridNewMerging', 'MetaData', meta_file);
save(fn_meta, 'meta_data');