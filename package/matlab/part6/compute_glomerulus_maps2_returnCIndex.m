% Code by: Martin.Strauch@lfb.rwth-aachen.de
% Computes 3D glomerulus/AL maps for all .mat files in <mat_path>
%
% <c>: Number of clusters/glomeruli
% <meta_data>: meta data file (.mat); needed if ALs should be flipped
%
% if meta_data=[], ALs will not be flipped 
% 
% regex: regular expression specifying the files that should be used for clustering, e.g.: '[A-Z][0-9]_m1' (first measurement "m1" for all reference odorants)
%
function [F_cc, C_index] = compute_glomerulus_maps2_returnCIndex(mat_path, c, meta_data, regex)

mat_files   = dir(strcat(mat_path,'*.mat'));
num_animals = length(mat_files); 
F           = cell(num_animals,1);
F_cc        = cell(num_animals,1);
C_index = cell(num_animals,1);
for(i=1:num_animals)
    F{i}    = load_animal(strcat(mat_path,mat_files(i).name));
    
    %construct an AL mask:
    M = F{i}.baselines{1}(:,:,1:24);
    for(j=1:size(M,3))
        M(:,:,j) = imbinarize(M(:,:,j) - mean(mean(M(:,:,j))));
        M(:,:,j) = imfilter(single(M(:,:,j)), fspecial('Gaussian',13,5),'replicate');
    end
    mask = M;
    %apply mask to cut out the areas outside of the AL:
    U_size=size(F{i}.U{1});
    U4 = reshape(F{i}.U{1}, F{i}.parameters.width, F{i}.parameters.height, F{i}.parameters.number_of_z_slices, size(F{i}.U{1},2));
    for(j=1:size(U4,4))
        U4(:,:,:,j) = U4(:,:,:,j).*mask; 
    end
    F{i}.U{1} = reshape(U4,U_size);
    
    [F_cc{i}, C_index{i}] = run_convex_cone2_returnCIndex(F{i}, 1, 50, c, 'NA', 0.5, meta_data.odor_names{i}, regex); 
end

%if meta_data.flip_AL exists, use the information to map between left and right ALs:
F_cc = flip_ALs(F_cc, meta_data);

