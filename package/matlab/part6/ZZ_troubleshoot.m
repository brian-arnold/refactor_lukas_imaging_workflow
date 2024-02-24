folder_home = '/Volumes/mcbride/zhileizhao/Scotty/Analysis/CalciumImaging/AnalysisPipelineOrco/MartinStrauch';
out_path  = fullfile(folder_home, 'SVDResults', 'humanNonhumanOrco/');
meta_data = load_meta_data('humanNonhumanOrco_meta_data.mat');
num_clusters = 40; 
selected_odors = '[A-Z][0-9]_m1'; %

mat_path=out_path;
c=40;
regex = selected_odors;

mat_files   = dir(strcat(mat_path,'*.mat'));
num_animals = length(mat_files); 
F           = cell(num_animals,1);
F_cc        = cell(num_animals,1);

i = 3;
F{i}    = load_animal(strcat(mat_path,mat_files(i).name));

F=F{i};channel=1; k=50;, cc_threshold='NA'; x_threshold=0.5; odor_names=meta_data.odor_names{i}

A = reconstruct_movie(F,channel);

m1 = grep_subset(odor_names, regex);
[m1_time_points, measurements] = extract_time_points(F.parameters.measurement_lengths, m1);
R = A(:,m1_time_points);

R = zscore(R,[],2);

%construct an AL mask:
M = max(F.baselines{1},[],3);
M = imbinarize(M-mean(mean(M)));
M = imfilter(single(M), fspecial('Gaussian',29,11),'replicate');

%figure, imagesc(M);
mask = zeros(F.parameters.width, F.parameters.height, F.parameters.number_of_z_slices);
for(i=1:size(mask,3))
    mask(:,:,i) = M;
end

% save the mask
options.overwrite = true;
fname = fullfile(folder_home, 'brainF3_mask.tif');
saveastiff(single(mask), fname, options);