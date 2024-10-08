function heatmapAligner(heatMap,fluoImage)

dimFI = size(fluoImage); %1280x1024
dimHM = [1400,800]; % Aspect Ratio = 7/4 = 1.75

%the coordinates where the centre of the polygon frame lies in the
%fluoImage
offset = [562,490];



figure
% imshow(fluoImage,'XData',[0,dimFI(2)],'YData',[0,dimFI(1)]);
hold on
x = [offset(1)-(dimHM(1)/2),offset(1)+(dimHM(1)/2)];
y = [offset(2)-(dimHM(2)/2),offset(2)+(dimHM(2)/2)];
imagesc(heatMap,'XData',[x(1),x(2)],'YData',[y(1),y(2)]);
axis([x(1),max([x(2),dimFI(2)]),0,dimFI(1)])

