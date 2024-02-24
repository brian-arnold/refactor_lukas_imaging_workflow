% Code by: Martin.Strauch@lfb.rwth-aachen.de
% Computes 3D glomerulus/AL maps for all .mat files in <mat_path>
%
% <c>: Number of clusters/glomeruli
% <meta_data>: meta data file (.mat); needed if ALs should be flipped
% <border_width>: cut out <border_width.x/.y> voxels at the borders of the volume (to remove motion artifacts)
%
% if meta_data=[], ALs will not be flipped 
% if border_with=[], borders will not be cut
%
function F_cc = compute_glomerulus_maps(mat_path, c, meta_data, border_width)

mat_files   = dir(strcat(mat_path,'*.mat'));
num_animals = length(mat_files); 
F           = cell(num_animals,1);
F_cc        = cell(num_animals,1);
for(i=1:num_animals)
    F{i}      = load_animal(strcat(mat_path,mat_files(i).name));
    
    if(isfield(border_width,'x') && isfield(border_width,'y'))
        F{i}.U{1} = mask_borders(F{i}.U{1}, border_width.x, border_width.y, F{i}.parameters.width, F{i}.parameters.height, F{i}.parameters.number_of_z_slices); 
    end
    
    F_cc{i} = run_convex_cone(F{i}, 1, 50, c, 'NA', 0.5);
end

%if meta_data.flip_AL exists, use the information to map between left and right ALs:
F_cc = flip_ALs(F_cc, meta_data);


