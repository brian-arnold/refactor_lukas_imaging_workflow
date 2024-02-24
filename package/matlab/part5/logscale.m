function F_mean_responses_log = logscale(F_mean_responses)

F_mean_responses_log = log(min_max_normalize(F_mean_responses)+1);
