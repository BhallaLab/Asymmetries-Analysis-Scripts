function frameStretcher(cameraImage,peakMap)

%% Rotation correction: This corrects the rotation of the grid peaks w.r.t the camera frame
% When overlaying two images: camera and heatmap, the heatmap pixel#1
% should correspond to left bottom and pixel#29 should correspond to left
% top. This rotation takes care of that.
% The transformation is essentially flipping in vertical and horizontal both and then taking a transpose.
% This approach is better than changing the original peakMap because the
% original peakMap preserves the pixel identity also as 'ij' address.

peakMapFT = (rot90(peakMap,2))';

%%
offset = 50;
scale = 0.5;

figA = figure;
imA = imagesc(cameraImage);
colormap('gray')
axis off
figA.Position = [0 0 822.5 612];
axA = gca;
axA.Units = 'pixels';
axA.Position = [132.5 50 640 512];
% figA.PaperPositionMode = 'auto';
% figA.PaperUnits = 'centimeters';
% figA.PaperPosition = [0 0 8.225 6.12];
saveas(gcf, 'cameraImage', 'bmp')


figB = figure;
imB = imagesc(peakMapFT);
colormap gray
axis off
figB.Position = [0 0 822.5 612];
axB = gca;
axB.Units = 'pixels';
axB.Position = [50 100 712 400];
% figB.PaperPositionMode = 'auto';
% figB.PaperUnits = 'centimeters';
% figB.PaperPosition = [0 0 8.225 6.12];
saveas(gcf, 'heatMap', 'bmp')


close all
