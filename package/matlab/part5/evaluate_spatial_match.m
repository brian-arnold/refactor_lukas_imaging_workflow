% visualizes distances between matchted glomeruli in animal F1 (this should be the target; the one with the larger number of glomeruli) and animal F2
%
function distances = evaluate_spatial_match(F1, F2)

X1 = F1.X_vis;
X2 = F2.X_vis;
w1 = F1.parameters.width;
h1 = F1.parameters.height;
num_slices1 = F1.parameters.number_of_z_slices;
w2 = F2.parameters.width;
h2 = F2.parameters.height;
num_slices2 = F2.parameters.number_of_z_slices;

num_clusters = size(X1,2);

distances = zeros(num_clusters,1);
means_X1  = zeros(num_clusters,3);
means_X2  = zeros(num_clusters,3);

for(i=1:num_clusters)

    X1_component = reshape(X1(:,i),w1,h1,num_slices1);
    [x1,y1,z1] = ind2sub(size(X1_component),find(X1_component==1));
    means_X1(i,:) = [mean(x1), mean(y1), mean(z1)];
 
    X2_component = reshape(X2(:,i),w2,h2,num_slices2);
    [x2,y2,z2] = ind2sub(size(X2_component),find(X2_component==1));
    means_X2(i,:) = [mean(x2), mean(y2), mean(z2)];
end

figure();
scatter3(means_X1(:,1), means_X1(:,2), means_X1(:,3),'filled'); hold on;
scatter3(means_X2(:,1), means_X2(:,2), means_X2(:,3), 'filled', 'MarkerFaceColor','red'); hold on;
for(i=1:num_clusters)
  line([means_X1(i,1), means_X2(i,1)], [means_X1(i,2), means_X2(i,2)], [means_X1(i,3), means_X2(i,3)]); hold on;  
  distances(i) = sqrt(sum((means_X1(i,:) - means_X2(i,:)).^2));
end

axis vis3d tight, box on, grid on    
camproj perspective    
daspect([1 1 1])           
xlim([1,w1]);
ylim([1,h1]);
zlim([1,num_slices1]);

