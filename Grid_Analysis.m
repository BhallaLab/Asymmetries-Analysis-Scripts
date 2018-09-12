clearvars
clc

%% Load .mat file with parsed tracelets

% These tracelets are going to be ordered whether external order or
% internal. This script will analyse and plot the data. Run this script
% after Patch_Tracelet_Parser.mat
% PatchTracelets is not shuffeled, orderedTracelets is always shuffled
% according to the stimluation square order

[FileName,PathName]=uigetfile('*.mat','Pick the Parsed Patch tracelet File'); %Opens a file selection box
TraceletFile = strcat(PathName,FileName);
cd(PathName) %Change the working directory to the path
load(FileName)


%% calculate peak of responses
gridPeak=zeros(gridSize);

for i=1:length(locs)
    gridPeak(i)=max(orderedPatchTracelets(i,:));
end
gridPeak = gridPeak';

% Clip any response that grows above 30 mV so that
% the resolution of heat map can be improved.
gridPeak(gridPeak>30)=30; % this is better way of capping values than doing it inside the for loop

% Generate and save map of peak responses
figure
gridPeakMap = imagesc(gridPeak);
colormap('jet')
h = colorbar();
title('Peak Response from baseline(Spikes clipped at 30)')
PeakImageFile = strcat(ExptID,'_gridPeakMap_',num2str(gridSize),'x');
print(PeakImageFile,'-dpng')


%% calculate AuC of responses

gridAuc=zeros(gridSize);
% aucDuration = time period for determination of AuC
aucDuration = pre*acqRate:(pre+post)*acqRate;

for i=1:length(locs)
    gridAuc(i)=trapz(orderedPatchTracelets(i,aucDuration));
end
gridAuc = gridAuc';

% Generate and save map of AUC values
figure
gridAucMap = imagesc(gridAuc);
colormap('jet')
h = colorbar();
title('Area Under the Curve for Responses after the stimulation')
AucImageFile = strcat(ExptID,'_gridAucMap_',num2str(gridSize),'x');
print(AucImageFile,'-dpng')

%% Generate traces figure
figure;
%Y-axis values are set to be 10% more than the min and max of all the
%tracelets
axis([-1*pre post 1.1*min(min(orderedPatchTracelets)) 1.1*max(max(orderedPatchTracelets))])

%drawing a figure and making it fullscreen size
figurePSTH=gcf;
figurePSTH.Units='normalized';
figurePSTH.OuterPosition=[0 0 1 1];

%plotting only one polygon TTL trace for reference
for row=1:length(locs)
    hold on
    plot(timeTracelet,orderedPatchTracelets(row,:),'k')        
end

hold on;
orderedPolTracelets(orderedPolTracelets>10)=10; %cap the TTL signal at 10
TTLGraph = area(timeTracelet,orderedPolTracelets(1,:));
TTLGraph.FaceColor = 'red'; %setting colour to red
TTLGraph.FaceAlpha = 0.5; 

hold on;
stimGraph = area(timeTracelet,exStimTrace); %plotting the stimulus trace for reference
stimGraph.FaceColor = 'blue' ; %setting colour to be blue
stimGraph.FaceAlpha = 0.5 ; %setting transparency

title('Response traces from baseline')
response_traces = strcat(ExptID,'_response_traces_',num2str(gridSize),'x');
print(response_traces,'-dpng')

%% calculate time of peak for responses

% for PSTH we need to get the timings of the peaks in all the traces. As
% there are going to be 841 traces for each trial, we can generate a
% distribution of peaks.

% 1. Find out the time points of peaks
% the function max finds maxima along columns, therefore the matrix of
% responses needs to be transposed.

[TraceletMax, TraceletPeakTime]= max(orderedPatchTracelets,[],2);

% 2. Get a histogram of timings

% Create figure

figure;
figurePSTH=gcf;
figurePSTH.Units='normalized'; %making the figure size in normalized units
figurePSTH.OuterPosition=[0 0 1 1]; %size of the figure -> full screen

% Create axes
axesPSTH = axes('Parent',figurePSTH);
hold(axesPSTH,'on');
box(axesPSTH,'on');
set(axesPSTH, 'fontsize', 10)

% Set the remaining axes properties
set(axesPSTH,'XTick',...
    [0 100 200 300 400 500 600 700 800 900 1000 1100 1200 1300 1400 1500 1600 1700 1800 1900 2000 2100 2200 2300 2400 2500 2600 2700 2800 2900 3000 3100 3200 3300 3400 3500 3600 3700 3800 3900 4000],...
    'XTickLabel',...
    {'-50','-45','-40','-35','-30','-25','-20','-15','-10','-5','0','5','10','15','20','25','30','35','40','45','50','55','60','65','70','75', '80','85','90','95','100','105','110','115','120','125','130','135','140','145','150'});

% Create histogram
histogram(TraceletPeakTime,'Parent',axesPSTH,'BinWidth',100);
title('PSTH')
xlim(axesPSTH,[0 4000])

% 3. Save histogram

PSTH = strcat(ExptID,'_PSTH_',num2str(gridSize),'x');
print(PSTH,'-dpng')

%% Time of Peak heatmap

% Plotting Time of Peak on the grid

%  For that create a matrix denoting the timing of the peak of the cell in each row
%  corresponding to the EPSP peak or spike

gridTimeofPeak=zeros(gridSize);

for i=1:length(locs)
    % The remapping of squares is done here using a variable 'j' which maps
    % the index i onto the correct coordinate of the square from the coord
    % data
    gridTimeofPeak(i)= TraceletPeakTime(i);
end
gridTimeofPeak = gridTimeofPeak';

% Generate and save map of Time of Peak values

figure
gridTimeofPeakMap = imagesc(gridTimeofPeak);
colormap('jet')
h = colorbar();
title('Time of Peak of Responses for the stimulation')
TimeofPeakImageFile = strcat(ExptID,'_gridTimeofPeakMap_',num2str(gridSize),'x');
print(TimeofPeakImageFile,'-dpng')

