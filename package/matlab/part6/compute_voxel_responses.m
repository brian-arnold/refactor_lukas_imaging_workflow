function [F_time_series, F_mean_responses, masks] = compute_voxel_responses(data_path, meta_data_file, intervals)

mat_files   = dir(strcat(data_path,'*.mat'));
meta_data   = load_meta_data(meta_data_file);
num_animals = length(mat_files);

F_mean_responses = cell(1,num_animals);
F_time_series    = cell(1,num_animals);
masks            = cell(1,num_animals);

for(i=1:num_animals)
    disp(strcat('processing data set', {' '}, num2str(i), {' '}, 'of', {' '}, num2str(num_animals)));
    F = load_animal(strcat(data_path,mat_files(i).name));
    
    A.C          = reconstruct_movie2(F,1);
    A.parameters = F.parameters;
    clear F;
    T = reorder_time_series(A, meta_data);
   
    regime1 = grep_subset(T.odor_names, intervals.regime1_odors);
    F_mean_responses{i} = to_mean_response2(T.time_series, T.measurement_lengths, regime1, intervals.regime1_reference, intervals.regime1_signal, intervals.regime2_reference, intervals.regime2_signal);
    F_time_series{i}.measurement_lengths = T.measurement_lengths;
    F_time_series{i}.odor_names = T.odor_names;
    clear T; clear A;
end

for(i=1:num_animals)
     F   = load_animal(strcat(data_path, mat_files(i).name));
     M = F.baselines{1}(:,:,1:F.parameters.number_of_z_slices);
     for(j=1:size(M,3))
        M(:,:,j) = imbinarize(M(:,:,j) - mean(mean(M(:,:,j))));
        M(:,:,j) = imfilter(single(M(:,:,j)), fspecial('Gaussian',13,5),'replicate');
     end
     masks{i} = M;
end
 