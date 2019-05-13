function [LM_thresholded, LM_weighted, CoM,asym,aD,LRAsym] = AsymmetryCalc(cell,cellID)
%% This function gives out structural/functional long axis of a neuron from the heatmap

%% Close figures
close all
clc

%% Get Matrices

% weighted matrix
rescaleGridWts = rescale(cell);
wts = [];
coordsWts = [];
for i=1:29
    for j = 1:29 
        coordsWts = [coordsWts; i,j];
        wts = [wts;rescaleGridWts(i,j)];
    end
end
xWts = coordsWts(:,2);
yWts = coordsWts(:,1);


% Thresholded matrix without weights (threshold = median heat value)

%Rescale the grid
binThres = median(wts);
rescaleGrid = rescale(cell);
rescaleGrid(rescaleGrid<binThres)=0;
rescaleGrid(rescaleGrid>binThres)=1;
coords = [];
for i=1:29
    for j = 1:29 
        if rescaleGrid(i,j)>0
            coords = [coords; i,j];        
        end
    end
end

% coordinates of points above threshold
x = coords(:,2);
y = coords(:,1);


%% getting regression
% first draw the figure
figure
plotregression(x,y)
print(strcat(cellID,'_reg'),'-dpng')

figure
gridAucMap = imagesc(rescaleGridWts);
colormap('jet')
h = colorbar();

% weighted regression
LM_weighted = fitlm(xWts,yWts,'Weights',wts);
c = LM_weighted.Coefficients.Estimate(1);
m = LM_weighted.Coefficients.Estimate(2);


%non-weighted (thresholded) regression
LM_thresholded = fitlm(x,y);
c2 = LM_thresholded.Coefficients.Estimate(1);
m2 = LM_thresholded.Coefficients.Estimate(2);

% Get Left/Right asymmetry above and below the line
LRAsym = getLRAsym(cell,c2,m2);
ann = strcat('LRAsym=',num2str(LRAsym));
text(1.5,24,strcat('\bf ',ann),'Color','w')

%plot the weighted regression line
hold on
lineX = linspace(0.5,29.5,100);
lineY = m.*lineX+c;

% Plot the non weighted (thresholded) regression line
hold on
lineY2 = m2.*lineX+c2;
plot(lineX,lineY,'k',lineX,lineY2,'m','LineWidth',1.0)
legend('Weighted','Thresholded','Location','SouthOutside')

x1 = lineX(1); x2=lineX(end);
y1 = lineY2(1); y2=lineY2(end);

% calculate centre of mass
tot_mass = sum(rescaleGridWts(:));
[ii,jj] = ndgrid(1:size(rescaleGridWts,1),1:size(rescaleGridWts,2));
R = sum(ii(:).*rescaleGridWts(:))/tot_mass;
C = sum(jj(:).*rescaleGridWts(:))/tot_mass;
CoM = [tot_mass,R,C];

% plot centre of mass on the same plot
hold on
plot(R,C,'.w','MarkerSize',20)
hold off

% distance between CoM and the line
asym = point_to_line([R,C,0],[x1,y1,0],[x2,y2,0]);
ann = num2str(asym);
text(1.5,28,strcat('\bf ','dist=',ann),'Color','w')

% Angle difference between the lines
X = [x1 x2];Y=[y1 y2];
a = [lineX(1) lineX(100)];b=[lineY(1) lineY(100)];
% angleDiff = (atan((Y(2)-X(2))/(Y(1)-X(1))) - atan((b(2)-a(2))/(b(1)-a(1)))) * 180/pi;
aD = abs((180/pi)*(atan(m2)-atan(m)));
ann = strcat(num2str(aD),'°');
text(1.5,26,strcat('\bf ','Angle Diff=',ann),'Color','w')



% save figure
print(cellID,'-dpng')


end

function d = point_to_line(pt, v1, v2)
      a = v1 - v2;
      b = pt - v2;
      d = norm(cross(a,b)) / norm(a);
end

function LRAsym = getLRAsym(currentCell,c,m)
heatAbove = 0;
heatBelow = 0; 

for x=1:29
    for y=1:29
        lineY = m*x+c;
        if y<lineY
            heatAbove = heatAbove+currentCell(x,y);
        else
            heatBelow = heatBelow+currentCell(x,y);
        end
    end
end

LRAsym = heatBelow/heatAbove;
end

