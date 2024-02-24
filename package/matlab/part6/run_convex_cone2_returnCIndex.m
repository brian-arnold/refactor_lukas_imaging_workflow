%code by: Martin.Strauch@lfb.rwth-aachen.de
%
% <F>           : input data from load_animal(<mat_file>)
% <channel>     : channel number (always 1 if only one dye has been used) 
% <k>           : number of PCs/singular vectors
% <c>           : number of columns to choose with the ConvexCone algorithm = number of glomeruli/clusters
% <cc_threshold>: parameter for smoothness-constrained variant of ConvexCone (currently we do not use this variant, i.e. the paramter is irrelevant)
% <x_threshold> : threshold to obtain binary glomerulus clusters from the continuous-valued clusters computed by ConvexCone
% <regex>       : regular expression specifying the files that should be used for clustering, e.g.: '[A-Z][0-9]_m1' (first measurement "m1" for all reference odorants)
%
function [glomeruli, C_indices] = run_convex_cone2_returnCIndex(F, channel, k, c, cc_threshold, x_threshold, odor_names, regex)

U = F.U{channel};
V = F.V{channel};
S = F.S{channel};

A = reconstruct_movie2(F,channel);

%work only on the measurements selected through <regex>:
%m1 = grep_subset(odor_names, regex);
%[m1_time_points, measurements] = extract_time_points(F.parameters.measurement_lengths, m1);
%R = A(:,m1_time_points);
%SVD of the selected submatrix:
%[U, S, V] = svd(R, 'econ');

%run ConvexCone:
[C, X, C_indices, reconstruction_accuracy, norm_error] = scca3d((U(:,1:k)*S(1:k,1:k))', F.parameters.width, F.parameters.height, F.parameters.number_of_z_slices, c, cc_threshold);

for(i=1:size(X,1))
     X(i,:) = local_smoothness_filter_3d(X(i,:), F.parameters.width, F.parameters.height, F.parameters.number_of_z_slices, C_indices(i), x_threshold);
     X(i,:) = X(i,:)./sum(X(i,:));
end
C     = (A(C_indices,:));
X_vis = binarize_X(X,0);

glomeruli.C = C;
glomeruli.X = X;
glomeruli.X_vis = X_vis;
glomeruli.U = U;
glomeruli.V = V;
glomeruli.S = S;
glomeruli.parameters = F.parameters;
glomeruli.names = F.names;
glomeruli.baselines = F.baselines;


clear A;
clear R;

%output:
%SVD results: left/right singular vectors (U,V) and singular values (S)
%ConvexCone results (NNCX factorization) results: matrices C (time series) and X (spatial components: glomerulus clusters)
%X_vis is a binarized/thresholded version of X for visualization purposes