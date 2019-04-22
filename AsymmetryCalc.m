%% example grid matrix


%% This is an analysis script for measuring asymmetry in the grid recording data
figure
gridAucMap = imagesc(exGrid);
colormap('jet')
h = colorbar();
h.Label.String = 'a.u.';
title('Area Under the Curve for Responses after the stimulation')

% %scale bar
% hold on;
% line(lineX,lineY,'Color','w','LineWidth',5);
% text(textX,textY,scaleSize,'Color','w','FontWeight','bold','FontSize',14)

%% calculate centre of mass

tot_mass = sum(exGrid(:));
[ii,jj] = ndgrid(1:size(exGrid,1),1:size(exGrid,2));
R = sum(ii(:).*exGrid(:))/tot_mass;
C = sum(jj(:).*exGrid(:))/tot_mass;
out = [tot_mass,R,C];

%% plot centre of mass on the same plot
hold on
plot(R,C,'.w','MarkerSize',20)


