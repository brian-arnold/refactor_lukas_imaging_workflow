% <U>: animal.U from animal = load_animal(<mat_file>)
% <x,y,z>: stack dimensions
% <border_width_x> and <border_width_y> pixels are cut out (to remove motion artifacts at the borders of the volume)
%
function U = mask_borders(U, border_width_x, border_width_y, x, y, z)

num = size(U,2);
U   = reshape(U, x, y, z, num);
U(1:border_width_x,:,:,:)       = 0;
U((x-border_width_x):x,:,:,:)   = 0;
U(:, 1:border_width_y,:,:)      = 0;
U(:, (y-border_width_y):y,:,:)  = 0;
U = reshape(U,x*y*z, num);
