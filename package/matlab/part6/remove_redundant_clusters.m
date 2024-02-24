function [F_time_series, F_mean_responses, F_cc] = remove_redundant_clusters(F_time_series, F_mean_responses, F_cc, threshold)

[unique_clusters, redundant, non_redundant, representative, to_merge] = find_redundant_clusters(F_mean_responses, threshold);

merged = merge_clusters(F_cc.X_vis, F_time_series, F_mean_responses, to_merge);
F_cc.X_vis       = cat(2, merged.X_vis, F_cc.X_vis(:, non_redundant));
F_mean_responses = cat(1, merged.F_mean_responses,  F_mean_responses(non_redundant,:));
F_time_series.time_series = cat(1, merged.time_series, F_time_series.time_series(non_redundant,:));

