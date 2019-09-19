function makeHeatPlots(peakMap, AUCMap, timeofpeakMap,cellIR,exptSpecs,scalebarSize)
objMag=exptSpecs.objMag;
gridSize=exptSpecs.gridSize;
peakThres=exptSpecs.peakThres;
cellID=exptSpecs.cellID;

[a,b,c,d,e]=gridScaleBar(scalebarSize,objMag,gridSize);
cmp = 'jet';
%% Rotation correction: This corrects the rotation of the grid peaks w.r.t the camera frame
% When overlaying two images: camera and heatmap, the heatmap pixel#1
% should correspond to left bottom and pixel#29 should correspond to left
% top. This rotation takes care of that.
% The transformation is essentially flipping in vertical and horizontal both and then taking a transpose.
% This approach is better than changing the original peakMap because the
% original peakMap preserves the pixel identity also as 'ij' address.

peakMapFT = (rot90(peakMap,2))';
AUCMapFT = (rot90(AUCMap,2))';
timeofpeakMapFT = (rot90(timeofpeakMap,2))';

%% Peak heatmap
figure
gridPeakMap = imagesc(peakMapFT);
daspect([0.4,0.7,1])
axis on
colormap(cmp)
h = colorbar();
h.Label.String = exptSpecs.unit;
%titleString = sprintf('Peak Response from baseline(Spikes clipped at %s %s)', int2str(peakThres), exptSpecs.unit);
titleString = sprintf('Peak Response from baseline (%s)',exptSpecs.unit);
title(titleString)

%scale bar
hold on;
line(a,b,'Color','w','LineWidth',5);
text(c,d,e,'Color','w','FontWeight','bold','FontSize',14)

plotFile = strcat(cellID,'_gridPeakMap_',num2str(gridSize),'x');
print(plotFile,'-dpng')

%% AuC heatmap
figure
gridAUCMap = imagesc(AUCMapFT);
daspect([0.4,0.7,1])
axis off
colormap(cmp)
h = colorbar();
h.Label.String = 'a.u.';
title('AuC of Response')

%scale bar
hold on;
line(a,b,'Color','w','LineWidth',5);
text(c,d,e,'Color','w','FontWeight','bold','FontSize',14)

plotFile = strcat(cellID,'_gridAUCMap_',num2str(gridSize),'x');
print(plotFile,'-dpng')

%% Time to Peak heatmap
figure
timeofpeakMapFT = timeofpeakMapFT./20;
gridtimetopeakMap = imagesc(timeofpeakMapFT);
daspect([0.4,0.7,1])
axis off
colormap(cmp)
h = colorbar();
h.Label.String = 'ms';
title('Time to Peak of Responses')

% scale bar
hold on;
line(a,b,'Color','w','LineWidth',5);
text(c,d,e,'Color','w','FontWeight','bold','FontSize',14)

plotFile = strcat(cellID,'_gridtimetoPeakMap_',num2str(gridSize),'x');
print(plotFile,'-dpng')

%% IR Trend
if isnan(cellIR)
    disp('No IR Trend Found')
else
    figure
    plot(1:length(cellIR),cellIR,'b','LineWidth',1);hold on;
    plot(1:length(cellIR),movmean(cellIR,[10 0]),'r','LineWidth',3)
    xlabel('Trials')
    ylabel('IR (M\Omega)')
    title('Input Resistance during the grid experiment')
    plotFile = strcat(cellID,'_IRtrend');
    print(plotFile,'-dpng')
end

%close all
end

function [lineposX, lineposY, textposX, textposY, scaleSize] = gridScaleBar(scalebarSize,objMag,gridSize)
%SCALEBAR gives out specifications to create a scale bar for the heatmaps
%of the responses. The scale factor depends on the model of the Polygon and
%objective magnification used in the experiments.
%
%Resolution data from Mightex website for Olympus system:
%Full frame = 13900 x 7800 µm (Diagonal = 16000 µm)
%Pixel size = 16.2 µm

%---Variables--------------------------------------------------------------
frameSize = 16000; %diagonal in µm, fixed for the Polygon 400E model
res = frameSize/(objMag*gridSize); %pixel size in µm/square for grid29
pxinBar = scalebarSize/res; %number of pixels in the scale bar

%----Scale Bar-------------------------------------------------------------
%ceil function with added 0.5 is to ensure that scale bar start point aligns
%with a square edge.
scalebarPos = ceil(0.8*gridSize)+0.5; 
lineposX = [scalebarPos scalebarPos+pxinBar]; %Line start and end points in X-axis
lineposY = [scalebarPos+0.5 scalebarPos+0.5]; %Line start and end points in Y-axis
textposX = scalebarPos; %text should begin where scale bar begins on X-axis
textposY = scalebarPos+1.5; %text should be 1.5 square below the scale bar in Y-axis
scaleSize = strcat(num2str(scalebarSize),' µm'); %text content

end



