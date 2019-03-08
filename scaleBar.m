function [lineposX, lineposY, textposX, textposY, scaleSize] = scaleBar(scalebarSize,objMag,gridSize)
% Resolution data from Mightex website for Olympus system
% Full frame = 13900 x 7800 µm (Diagonal = 16000 µm)
% Pixel size = 16.2 µm

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
