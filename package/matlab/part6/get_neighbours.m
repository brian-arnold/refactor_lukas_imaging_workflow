function [neighbours, already_found] = get_neighbours(C, nr, already_found)

n   = find(C(nr,:)==1);
idx = find(ismember(n, already_found)==0);
if(length(idx)>0)
    neighbours = n(idx);
    already_found = [already_found, neighbours];
else 
    neighbours=[];
    return;
end

for(i = 1:length(neighbours))
    this_n = find(C(neighbours(i),:)==1);
    idx    = find(ismember(this_n, already_found)==0);
    %recursively visit also neighbours of neighbours:
    for(j = 1:length(idx))
        [current_neighbours, already_found] = get_neighbours(C, this_n(idx(j)), already_found);
        if(length(current_neighbours)>0)
            neighbours = [neighbours, current_neighbours];
        end
    end
    %if(length(idx)>0)
    %    neighbours = [neighbours, this_n(idx)];
    %    already_found = [already_found, this_n(idx)];
    %end
end

