function [correlation, cripness] = ZZ_plotpw_normcore(Y, M, bnd, fig_name)
%ZZ_PLOT_NORMCORE plot shifts and image similarity metrics
% Y is raw data
% M is motion corrected data
% since both Y and M are 3D volume, need to loop through all slices when
% calculate cripness
% bnd is how many pixels to remove at the boundaries
% return the correlation and cripness before and after motion correction
% purpose, save the figure in .fig format

zStacks = size(Y,3);
for z=1:zStacks
    [cYz,mYz,vYz] = motion_metrics(squeeze(Y(:,:,z,:)),bnd);
    [cMz,mMz,vMz] = motion_metrics(squeeze(M(:,:,z,:)),bnd);
    cripnessY(z) = vYz;
    cripnessM(z) = vMz;
end

% plot the cripness
f = figure('visible','off');
subplot(2,1,1);
cripness(:,1) = cripnessY;
cripness(:,2) = cripnessM;
bar(cripness); 
xlabel('slice #','fontsize',14,'fontweight','bold');
ylabel('cripness','fontsize',14,'fontweight','bold');
legend('raw data', 'corrected');


% plot correlation
[cY, mY, vY] = motion_metrics(Y, bnd);
[cM1, mM1, vM1] = motion_metrics(M, bnd);

T = length(cY);
subplot(2,1,2);
plot(1:T,cY,1:T,cM1); 
xlabel('time point #','fontsize',14,'fontweight','bold');
ylabel('correlation coefficient','fontsize',14,'fontweight','bold');
legend('raw data','corrected'); 

savefig(f, fig_name);
correlation = [(1:T)',cY,cM1];
