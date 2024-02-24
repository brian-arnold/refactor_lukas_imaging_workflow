%code by: Martin.Strauch@lfb.rwth-aachen.de
%
%input path structure: <path>\animal_name\<subfolder>\
%
%<subfolder>    : selects one from several folders (experiments) for the current animal, e.g. 'odorEvoked' 
%<files_to_read>: the file names of the image stacks we want to use end like this: '*-mc.tif.between.tif' (there can be variants with/without motion correction etc.)
%<output_folder>: save the output .mat file to this folder
%<width>, <height>, <z>    : dimensions of the image stack
%<reference_interval>      : an interval of n1,...,n2 time points before odor stimulation (used for normalization)
%<gauss_kernel>, <gauss_sd>: parameters for Gaussian smoothing
%
%example:
%process_folder_orco('D:\ShareWithMartin20190906\humanNonhumanOrco\', 'odorEvoked',  '*-mc.tif.between.tif', 'D:\test\', 128, 128, 24, 1:20, 13, 5); 
%
function process_folder_orco(path, subfolder, files_to_read, output_folder, width, height, z,reference_interval, gauss_kernel, gauss_sd)

d        = dir([path, '*.*']);
d_subset = [d(:).isdir];
folders  = {d(d_subset).name}';
folders(ismember(folders,{'.','..'})) = [];

for(i=1:length(folders))
    
    f = strcat(path, folders{i}, '\', subfolder, '\');
    
    [U, S, V, baselines, names, parameters] = read_and_preprocess_orco(f, files_to_read, width, height, z, reference_interval, gauss_kernel, gauss_sd);
    parameters.measurement_name  = folders(i);
    
    max_number_of_singular_vectors = 500;
    %currently, processing of a second channel is deactivated
    %for(j=1:2)
    for(j=1:1)
        if(size(U{j},2)>max_number_of_singular_vectors)
            U{j} = U{j}(:, 1:max_number_of_singular_vectors);
            V{j} = V{j}(:, 1:max_number_of_singular_vectors);
            S{j} = S{j}(1:max_number_of_singular_vectors, 1:max_number_of_singular_vectors);
        end
    end
    
    output_file = strcat(output_folder, '\', folders{i}, '_', subfolder, '.mat');
    save(output_file, 'U', 'S', 'V', 'baselines', 'names', 'parameters', '-v7.3');
end


