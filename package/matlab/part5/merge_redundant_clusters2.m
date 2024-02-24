function [F_time_series, F_mean_responses, F_cc] = merge_redundant_clusters2(F_time_series, F_mean_responses, F_cc, correlation_threshold, dice_threshold)

[redundant, non_redundant, to_merge] = merge_X2(F_cc.X_vis, F_mean_responses, correlation_threshold, dice_threshold);

if(length(to_merge)>0)
    merged = merge_clusters(F_cc.X_vis, F_time_series, F_mean_responses, to_merge);
    F_cc.X_vis = merged.X_vis;
    F_mean_responses = merged.F_mean_responses;
    F_time_series.time_series = merged.time_series;
end


