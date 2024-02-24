%code by: Martin.Strauch@lfb.rwth-aachen.de
%
function response_pattern_3d(pattern, X, w, h, num_slices, transparency, azimuth, elevation)

num_colors = 4096;
map        = CubeHelix(num_colors,0.5,-1.5,1.2,1.0);

for i=1:size(X,2)
    temp = reshape(X(:,i),w,h,num_slices);
    [x,y,z] = ind2sub(size(temp),find(temp==1));
   
    p=plot(alphaShape(x,y,z)); set(gca,'Color',[0.75 0.75 0.75]);
    
    alpha = 0.75;
    if(strcmp(transparency, 'transparent'))
        alpha = pattern(i)/num_colors;
    end
      
    set(p, 'FaceColor',map(floor(pattern(i)),:), 'EdgeColor','none', 'FaceAlpha', alpha)  
    hold on;
end

daspect([1 1 1])                             
view(azimuth,elevation);
axis vis3d tight, box on, grid on    
camproj perspective                           
camlight 
lighting phong  
 
xlim([1,w]);
ylim([1,h]);
zlim([1,num_slices]);
colormap(map);
colorbar;
set(gcf,'color','w');



