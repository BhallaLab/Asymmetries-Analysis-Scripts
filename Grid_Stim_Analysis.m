clear variables
clc

%% Traceparsing

[recFileName,recPathName]=uigetfile('*.atf','Pick the Data File'); %Opens a file selection box
traceFile = strcat(recPathName,recFileName); %Make it the current working directory

ExptID = strsplit(recFileName,'.'); %Separate file name from extension
ExptFile = ExptID{1};
mkdir(recPathName,ExptFile)
exportFileDir = strcat(recPathName,ExptFile);
cd(exportFileDir)

[coordFileName, coordFilePathName]=uigetfile('*.txt','Select the file with Grid Order');
coordFile = strcat(coordFilePathName,coordFileName);
tic;
[patchTrace, polygonTrace, timeTrace, metadata] = traceParse(traceFile,coordFile);


%% Input resistance trend during grid stimulation
pulseCurrent = -50; %pA, hyperpolarizing current pulse given to measure IR
[avgIR, cellIR] =  IRTrend(patchTrace,pulseCurrent);
plot(1:881,cellIR,'b','LineWidth',1);hold on;
plot(1:881,movmean(cellIR,[10 0]),'r','LineWidth',3)
title('Input Resistance during the experiment')
plotFile = strcat(ExptFile,'_IRtrend_');
savefig(plotFile)

%% Grid Analysis
peakThres = 30 ; %above the baseline
traceDur = 100; %time in millisecond for area under the curve calculation
[peakMap, AUCMap, timeofpeakMap] = heatmap(patchTrace, peakThres,traceDur);

%% Generate traces figure
tracePlot(patchTrace,polygonTrace,timeTrace,ExptFile)


%% plot heatmaps
scalebarSize = 50; 
gridSize = 29;
objMag = 40;
gridSiz = 29;
makePlots(peakMap, AUCMap, timeofpeakMap,scalebarSize,objMag,gridSize,ExptFile)


%% Saving File

clear exportFileDir
clear ExptID
clear PathName
clear TraceFile

save(ExptFile)

toc