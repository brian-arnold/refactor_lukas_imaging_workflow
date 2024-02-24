function [mapping, mapped_indices, cost] = register_all_to_target(F_mean_responses, F_cc, constraints, target_index)

num_animals    = length(F_cc);
mapping        = cell(num_animals, 1);
mapped_indices = cell(num_animals,1);
cost = zeros(num_animals,1);

for(i=1:num_animals)
   
    [subject, target, this_mapping, this_mapped_indices, this_cost] = register_animals(F_mean_responses, F_cc, constraints, target_index, i);

    mapping{i}        = this_mapping;
    mapped_indices{i} = this_mapped_indices;
    cost(i) = this_cost;
end