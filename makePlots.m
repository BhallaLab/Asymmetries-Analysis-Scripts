function makePlots(peakMap, AUCMap, timeofpeakMap,scalebarSize,objMag,gridSize,ExptID)
[a,b,c,d,e]=scaleBar(scalebarSize,objMag,gridSize);
cmp = 'jet';

% Peak heatmap
figure
gridPeakMap = imagesc(peakMap);

colormap(cmp)
h = colorbar();
h.Label.String = 'mV';
title('Peak Response from baseline(Spikes clipped at 30)')

%scale bar
hold on;
line(a,b,'Color','w','LineWidth',5);
text(c,d,e,'Color','w','FontWeight','bold','FontSize',14)

plotFile = strcat(ExptID,'_gridPeakMap_',num2str(gridSize),'x');
print(plotFile,'-dpng')

% AuC heatmap
gridAUCMap = imagesc(AUCMap);
colormap(cmp)
h = colorbar();
h.Label.String = 'mV';
title('AuC of Response')

%scale bar
hold on;
line(a,b,'Color','w','LineWidth',5);
text(c,d,e,'Color','w','FontWeight','bold','FontSize',14)

plotFile = strcat(ExptID,'_gridAUCMap_',num2str(gridSize),'x');
print(plotFile,'-dpng')

% Time to Peak heatmap
gridtimetopeakMap = imagesc(timeofpeakMap);
colormap(cmp)
h = colorbar();
h.Label.String = 'mV';
title('Time to Peak of Responses')

% scale bar
hold on;
line(a,b,'Color','w','LineWidth',5);
text(c,d,e,'Color','w','FontWeight','bold','FontSize',14)

plotFile = strcat(ExptID,'_gridtimetoPeakMap_',num2str(gridSize),'x');
print(plotFile,'-dpng')
 
close all
end
