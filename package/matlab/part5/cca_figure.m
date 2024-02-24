function [C, names] = cca_figure(F_time_series, F_mean_responses, masks, odor_subset, out_file)

[xdim, ydim, zdim] = size(masks{1});
num_animals = length(F_mean_responses);

X = []; 
for(i=1:num_animals)
    s2 = size(F_mean_responses{i}(:,odor_subset),2);
    R4 = reshape(F_mean_responses{i}(:,odor_subset), xdim, ydim, zdim, s2);
    for(j=1:s2)
        R4(:,:,:,j) = medfilt3(R4(:,:,:,j),[5 5 1]);
        R4(:,:,:,j) = R4(:,:,:,j) .* masks{i};
    end
    R = reshape(R4, xdim*ydim*zdim, s2);
    X = cat(2, X, whiten_matrix(R, 100)');
end
X = X';
[U, S, V] = svd(X, 'econ');

C = cell(num_animals,1);
s = xdim*ydim*zdim; 
k = 3;
for(i=1:num_animals)  
    start = (i-1)*s+1;
    stop  = i*s;
    S = U(start:stop,1:k);
    C{i} = pinv(S)*X(start:stop,:);
    C{i} = C{i}';
end

names = F_time_series{i}.odor_names(odor_subset);
le = length(odor_subset);
figure();
for(i=1:length(C))
    scatter3(C{i}(1:le,1), C{i}(1:le,2), C{i}(1:le,3));
    text(C{i}(1:le,1), C{i}(1:le,2), C{i}(1:le,3), names, 'Interpreter', 'none', 'FontSize',8);
    hold on;
end

write_to_csv(C, names, out_file);



