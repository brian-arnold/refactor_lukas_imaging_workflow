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
function F_cc = ZZ_compute_glomerulus_maps2(mat_path, c, meta_data, regex, masks)

mat_files   = dir(strcat(mat_path,'*.mat'));
num_animals = length(mat_files); 
F           = cell(num_animals,1);
F_cc        = cell(num_animals,1);
for(i=1:num_animals)
    F{i}    = load_animal(strcat(mat_path,mat_files(i).name));
    F_cc{i} = ZZ_run_convex_cone2(F{i}, 1, 50, c, 'NA', 0.5, meta_data.odor_names{i}, regex, masks{i}); 
end

%if meta_data.flip_AL exists, use the information to map between left and right ALs:
F_cc = flip_ALs(F_cc, meta_data);

%try: adaptive functional/spatial weight: the lower the maximum response to any odor, the lower should the functional weight (and the higher the spatial weight) be
%then: no response --> only spatial matching, strong response: mostly functional matching