function F_mean_responses_log = logscale_mean_responses(F_mean_responses)

F_mean_responses_log = cell(length(F_mean_responses),1);

for(i=1:length(F_mean_responses_log))
    F_mean_responses_log{i} = logscale(F_mean_responses{i});
end