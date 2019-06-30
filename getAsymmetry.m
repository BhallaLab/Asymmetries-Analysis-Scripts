clear all
clc

%% Load cells

load('Grid_response_all_cells.mat')
allCells = who('cell*');

%% Get Asymmetry regression plots
tic
cellData = struct([]);
for i=1:numel(allCells)
    currentCell = eval(allCells{i});
    cellData(i).cellID = allCells{i};
    cellData(i).heatMap=currentCell;
    [cellData(i).regFit,~,cellData(i).CoM,cellData(i).CoMDistance,cellData(i).angleDiff,cellData(i).heatDist] = AsymmetryCalc(currentCell,allCells{i});
    
end
toc

%%
clear currentCell i LM*

save('Asymmetry_all_cells.mat','cellData')
