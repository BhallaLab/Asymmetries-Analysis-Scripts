% This file allows the analysis of the grid stimulation experiments.
% Different parts do different functions. A directory with a
% single cell's data should be loaded. Rest all the things are taken care
% of. The file names in the directory should be correctly named and must
% contain keywords 'IR', 'IF', 'Grid' and 'Coordinates' for the program to
% indentify them and process. This script will analyse and produce results
% for passive properties, IF relationship and grid stimulation heat maps.
% Check the parameters before running.
% check the function files also: GETPASSIVEPROP, IFPLOTTER, GRIDANALYSIS,
% MAKEHEATPLOTS, ASYMMETRYCALC.
%
% GNU GPL: Aditya Asopa, Bhalla Lab, NCBS, 2019

%% Clear the workspace
clear variables
clc

%% Experiment Specifications

% A metadata container to hold all the experimental parameters
exptSpecs = struct();

% User inputs
prompt = {'Voltage or Current Clamp',...                   % User specified parameters
    'Grid Size',...
    'Number of Flanking Frames:'...
    'Light Intensity in %',...
    'Pulse Duration in ms',...
    'Grid Inter Stim Interval in Seconds'};
dlgtitle = 'Enter Experiment Parameters';
dims = 1;
definput = {'Current','29','20','100','10','1'};           %default values
exptParam = inputdlg(prompt,dlgtitle,dims,definput);

exptSpecs(:).clamp = lower(string(exptParam{1}));          % voltage/current clamp
exptSpecs(:).gridSize = str2num(exptParam{2});             % size of the grid
exptSpecs(:).blankFrames = str2num(exptParam{3});          % flanking blank frames
exptSpecs(:).lightPulseDur = str2double(exptParam{4});     % Light pulse duration
exptSpecs(:).lightIntensity = str2double(exptParam{5});    % Light intensity
exptSpecs(:).gridISI = 1000*str2double(exptParam{6});      % Inter frame interval

clear definput dims dlgtitle exptParam prompt              % remove useless variables

%% Get the directory path containing the recordings

recDir = uigetdir('','Select the directory with recording.');
[~,cellID] = fileparts(recDir);                            % get the directory from the end of the path
[exptSpecs(:).cellID] = cellID;                            % Store cell ID

cd(recDir)
fileList = dir('*.*t*'); %list all the files containing data *.atf, *.txt

for i=1:size(fileList,1)
    fil = lower(fileList(i).name);                         % make all filenames lower case
    if contains(fil,'grid')
        gridRecFile = fileList(i).name; 
    elseif contains(fil,'ir')
        IRFile = fileList(i).name;
    elseif contains(fil,'if')
        IFFile = fileList(i).name;
    elseif contains(fil,'coordinates')
        coordFile = fileList(i).name;
    else
        error('No ePhys files or coordinates file found in the folder')
    end
end


%% Passive Properties

if exist('IRFile','var')
    %Change following values if needed
    [exptSpecs(:).IRpulseStart] = 300;                     % Refer to axon protocol
    [exptSpecs(:).IRpulseEND] = 600;                       % Refer to axon protocol
    [exptSpecs(:).IRpulse] = -100;                         % Refer to axon protocol

    [IRPre,CmPre,tauPre]=getPassiveProp(IRFile,exptSpecs);
end


%% IF Relationship

if exist('IFFile','var')
    %Change following values if needed
    [exptSpecs(:).IFStepStart] = 100;                      % Refer to axon protocol
    [exptSpecs(:).IFStepEnd] = 600;                        % Refer to axon protocol
    [exptSpecs(:).IFcurrentSteps] = [-50:10:140];          % Refer to axon protocol
 
    numSpikesPre = IFPlotter(IFFile,exptSpecs);
end

%% Grid Response Analysis

if exist('gridRecFile','var') && exist('coordFile','var')
    exptSpecs(:).IRPulse = -20;                            % Refer to axon protocol
    exptSpecs(:).IFStepStart = 100;                        % Refer to axon protocol
    exptSpecs(:).gridBaseline = 1:4000;                    % Refer to axon protocol
    exptSpecs(:).gridPre = 200;                            % Refer to axon protocol
    
    recFileInfo = dir(gridRecFile);
    exptSpecs(:).dateofExpt = recFileInfo.date;            % Save experiment data from file
    clear recFileInfo;
    
    if exptSpecs.clamp =='voltage'                         % Create unit based on clamp type
        exptSpecs(:).unit = 'pA';
    elseif exptSpecs.clamp == 'current'
        exptSpecs(:).unit = 'mV';
    else
        exptSpecs(:).unit = '';
    end
            
    % Main Analysis to generate heatmaps, tracelets,and update metadata
    [IRtrend, peakMap, AUCMap, timetopeakMap,peakThres,polygonTracelet,patchTracelets,exptSpecs] = gridAnalysis(gridRecFile,coordFile,exptSpecs);
    exptSpecs(:).peakThres = peakThres;clear peakThres     % Refer to gridAnalysis function 
                                                           
end

%% Plots
% Make Heat Map Plots

exptSpecs(:).objMag = 40;                                  % Objective used
scaleBar = 50;                                             % 50um scale bar
makeHeatPlots(peakMap, AUCMap, timetopeakMap,IRtrend,exptSpecs,scaleBar)

% Make distribution histogram of peaks of responses
responsePeaks = reshape(peakMap,[1,exptSpecs.gridSize^2]); % get a row vector
figure;
hist(responsePeaks,100);                                   % 100 bins
title('Distribution of Response Amplitudes')
plotFile = strcat(cellID,'_responseDist_',num2str(exptSpecs.gridSize),'x');
print(plotFile,'-dpng')

%% Heatmap images for alignment
[cameraImageFile,~] = uigetfile('*.bmp','Pick the Camera Image File');
cameraImage = imread(cameraImageFile);clear cameraImageFile;
frameStretcher(cameraImage,peakMap);


%% All response grid (Optional and Heavy)
% figure;
% figureResponses=gcf;
% figureResponses.Units='normalized';
% figureResponses.OuterPosition=[0 0 1 1];
% gridSize = exptSpecs.gridSize;
% 
% plotLim = [min(min(patchTracelets)) max(max(patchTracelets))];
% for i=1:gridSize^2
%     subplot(gridSize,gridSize,i)
%     plot(patchTracelets(i,:))
%     ylim(plotLim)
%     axis off
% end
% 
% title('All Responses to Grid Stimulation')
% allResponses = strcat(cellID,'allResponses_',num2str(exptSpecs.gridSize),'x');
% print(allResponses,'-dpng')
% clear fig*



%% Save workspace

save(cellID)

% close all
