% F_cc = compute_glomerulus_maps(out_path, num_clusters, meta_data, border_width);
% F_mean_responses = compute_mean_responses(F_time_series, regime1, regime1_intervals, regime2_intervals);
%
% 3D heatmap for the first odor from the first animal:
% show_3D_heatmap(F_cc{1}, F_mean_responses{1}, 1)
%
% example with graphics options:
% show_3D_heatmap(F_cc{1}, F_mean_responses{1}, 1, {'transparent', 60, 30})
% 'transparent' --> glomeruli appear the more transparent the weaker their responses are
% view options: azimuth=60, elevation = 30
%
function show_3D_heatmap(convex_cone_result, mean_responses, odor_nr, graphics_options)

num_colors = 4096;
w = convex_cone_result.parameters.width;
h = convex_cone_result.parameters.height;
z = convex_cone_result.parameters.number_of_z_slices;

if exist('graphics_options','var')
    transparency = graphics_options{1};
    azimuth = graphics_options{2};
    elevation = graphics_options{3};
else
    transparency = [];
    azimuth =  -37.5;
    elevation = 30;
end

%global min-max color scale for all mean odor responses from the current animal:
scaled = global_scaling(mean_responses, num_colors);
figure();
response_pattern_3d(scaled(:,odor_nr), convex_cone_result.X_vis, w, h, z, transparency, azimuth, elevation); 
