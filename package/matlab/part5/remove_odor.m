% [F_time_series{1}, F_mean_responses{1}] = remove_odor(F_time_series{1}, F_mean_responses{1}, 37);
%
function [F_time_series, F_mean_responses] = remove_odor(F_time_series, F_mean_responses, odor_nr)

if(isfield(F_time_series, 'time_series'))
    cut_start = sum(F_time_series.measurement_lengths(1:odor_nr-1)) +1;
    cut_end   = cut_start + F_time_series.measurement_lengths(odor_nr);
    l = size(F_time_series.time_series,2);

    F_time_series.time_series = F_time_series.time_series(:, [1:cut_start, cut_end:l]);
end

F_mean_responses(:,odor_nr)  = [];

F_time_series.odor_names(odor_nr)          = [];
F_time_series.measurement_lengths(odor_nr) = [];
