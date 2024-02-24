% [F_time_series{1}, F_mean_responses{1}] = remove_odor(F_time_series{1}, F_mean_responses{1}, 37);
%Zhilei Zhao
% modified from the remove_odor function
function [F_time_series] = ZZ_remove_odor(F_time_series, odor_nr)


cut_start = sum(F_time_series.measurement_lengths(1:odor_nr-1)) +1;
cut_end   = cut_start + F_time_series.measurement_lengths(odor_nr);
l = size(F_time_series.time_series,2);

F_time_series.time_series = F_time_series.time_series(:, [1:cut_start, cut_end:l]);


%F_mean_responses(:,odor_nr)  = [];

F_time_series.odor_names(odor_nr)          = [];
F_time_series.measurement_lengths(odor_nr) = [];
