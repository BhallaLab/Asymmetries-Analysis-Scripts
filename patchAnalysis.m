% This file allows the analysis of the grid stimulation experiment in
% current clamp. Different parts do different functions. A directory with a
% single cell's data should be loaded. Rest all the things are taken care
% of. The file names in the directory should be correctly named and must
% contain keywords 'IR', 'IF', 'Grid' and 'Coordinates' for the program to
% indentify them and process. This scrip will analyse and produce results
% for passive properties, IF relationship and grid stimulation heat maps.
% Check the parameters before running.
% check the function files also: GETPASSIVEPROP, IFPLOTTER, GRIDANALYSIS,
% MAKEHEATPLOTS, ASYMMETRYCALC.
%
% GNU GPL: Aditya Asopa, Bhalla Lab, NCBS, 2019

%% Clear the workspace
clear all
clc

%% Get the directory path containing the recordings
recDir = uigetdir('','Select the directory with recording.');
[~,cellID] = fileparts(recDir);  % get the directory from the end of the path
cd(recDir)
fileList = dir('*.*t*'); %list all the files containing data *.atf, *.txt

for i=1:size(fileList,1)
    fil = fileList(i).name;
    if contains(fil,'grid')
        gridRecFile = fileList(i).name;
    elseif contains(fil,'IR')
        IRFile = fileList(i).name;
    elseif contains(fil,'IF')
        IFFile = fileList(i).name;
    elseif contains(fil,'coordinates')
        coordFile = fileList(i).name;
    else
        error('No ePhys files or coordinates file found in the folder')
    end
    i=i+1;
end


%% Passive Properties
pulseStart = 300;
pulseEnd = 600;
currentPulse = -100; %current pulse in pA

[IR,Cm,tau]=getPassiveProp(IRFile,pulseStart,pulseEnd,currentPulse);

%% IF Relationship
currentSteps = -50:10:140;
stepStart = 100;
stepEnd = 600;

numSpikes = IFPlotter(IFFile,stepStart,stepEnd,currentSteps);

%% Grid Response Analysis
blankFrames = 0; % number of blank frames on either side of grid stimulation 
lightPulseDur = 20; %Light pulse duration, ms
lightIntensity = 100; %LED brightness
IRpulse = -50; %hyperpolarizing pulse amplitude for IR measurement
gridSize = 29;

[IRtrend, peakMap, AUCMap, timetopeakMap, gridMetadata] = gridAnalysis(gridRecFile,coordFile,blankFrames,lightPulseDur,lightIntensity,gridSize,cellID);

%% Make Heat Map Plots
objMag = 40; % objective used during grid stim
scaleBar = 50; % 50um scale bar
makeHeatPlots(peakMap, AUCMap, timetopeakMap,IRtrend,scaleBar,objMag,gridSize,cellID)

%% Save workspace
save(cellID)
close all
