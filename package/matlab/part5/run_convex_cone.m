%code by: Martin.Strauch@lfb.rwth-aachen.de
%
% <F>           : input data from load_animal(<mat_file>)
% <channel>     : channel number (always 1 if only one dye has been used) 
% <k>           : number of PCs/singular vectors
% <c>           : number of columns to choose with the ConvexCone algorithm = number of glomeruli/clusters
% <cc_threshold>: parameter for smoothness-constrained variant of ConvexCone (currently we do not use this variant, i.e. the paramter is irrelevant)
% <x_threshold> : threshold to obtain binary glomerulus clusters from the continuous-valued clusters computed by ConvexCone
%
function glomeruli = run_convex_cone(F, channel, k, c, cc_threshold, x_threshold)

U = F.U{channel};
V = F.V{channel};
S = F.S{channel};

[C, X, C_indices, reconstruction_accuracy, norm_error] = scca3d(U(:,1:k)', F.parameters.width, F.parameters.height, F.parameters.number_of_z_slices, c, cc_threshold);

for(i=1:size(X,1))
    X(i,:) = local_smoothness_filter_3d(X(i,:), F.parameters.width, F.parameters.height, F.parameters.number_of_z_slices, C_indices(i), x_threshold);
end
X_vis = binarize_X(X, 0);

A = reconstruct_movie(F,channel);
C = pinv(X_vis.*X') * A;

glomeruli.C = C;
glomeruli.X = X;
glomeruli.X_vis = X_vis;
glomeruli.U = U;
glomeruli.V = V;
glomeruli.S = S;
glomeruli.parameters = F.parameters;
glomeruli.names = F.names;
glomeruli.baselines = F.baselines;

%output:
%SVD results: left/right singular vectors (U,V) and singular values (S)
%ConvexCone results (NNCX factorization) results: matrices C (time series) and X (spatial components: glomerulus clusters)
%X_vis is a binarized/thresholded version of X for visualization purposes