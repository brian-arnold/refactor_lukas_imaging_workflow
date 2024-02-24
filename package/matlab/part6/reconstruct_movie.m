%code by: Martin.Strauch@lfb.rwth-aachen.de
%
%F = load_animal(<mat_file>)
%
function A = reconstruct_movie(F, channel)

U = F.U{channel};
V = F.V{channel};
S = F.S{channel};
num_odors = size(F.parameters.odor_names,1);
baselines = reshape(F.baselines{channel}, F.parameters.width, F.parameters.height, F.parameters.number_of_z_slices, num_odors);

%reconstruct the data from SVD:
A = U*S*V';

%add the baselines that have been subtracted before SVD:
num_time_points = size(A,2);
A = reshape(A, F.parameters.width, F.parameters.height, F.parameters.number_of_z_slices, num_time_points);

tp      = 1;
counter = 1; 
for(i=1:num_odors)
    B = baselines(:,:,:,counter);
      
    for(j=1:F.parameters.measurement_lengths(i))
        
        A(:,:,:,tp) = A(:,:,:,tp) + B;
        tp = tp + 1;
    end
    
    counter = counter+1; 
end

A = reshape(A, F.parameters.width * F.parameters.height * F.parameters.number_of_z_slices, num_time_points); 


