function indices = sort_by_max_response(F_mean_responses)

max_responses = max(F_mean_responses,[],2);
[sorted, indices] = sort(max_responses,'descend');