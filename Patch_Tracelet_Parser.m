clearvars
clc

%% This script will load a file with two channel gap free data Ch#1 Patch trace
%  and Ch#2 Polygon trace, and use the grid coordinates file (sequentially or
%  randomly ordered) and generate a matrix of tracelets for all the polygon
%  squares in order.

acqRate = 20;
pre = 50; 
post = 150;
points = acqRate*(pre+post)+1; %Number of datapoints

%getting a pulse duration value for the experiment
pulseDur = acqRate*input('how long is the pulse in ms?');

%% Load files

[FileName,PathName]=uigetfile('*.mat;*.txt','Pick the Data File'); %Opens a file selection box
TraceFile = strcat(PathName,FileName);
cd(PathName) %Change the working directory to the path
load(FileName) %Load the file
ExptID = strsplit(FileName,'.'); %Extract the filename ignoring the extension
ExptID = ExptID(1); ExptID = ExptID{1};

% To import the external grid coordinates

%Get the Coordinates file using GUI window
[coordFileName, coordFilePathName]=uigetfile('*.txt','Select the file with Grid Order');

%Open and scan it for the fourth column as the coordinates data is always
%in the fourth column. In textscan, %u = unsigned integer format, an
%asterisk in %*u means that that column is to be skipped

fid = fopen(coordFileName);
coord = textscan(fid,'%*u%*u%*u%u');
fclose(fid);
coord = coord{1};

%Estimating grid size from the coordinates file
gridSize = sqrt(length(coord));

%% Create Traces and Separate out Triggered Responses

% All data files are parsed and therefore have two variables
% 1)PatchTrace
% 2)PolygonTrace

%First get a thesholded Polygon trace and binarize it.
% The expression "PolygonTrace>50" would generate a logical array which
% won't be useful for calculation or findpeaks function. To convert it into
% a double array, we can multiply it with 1.00. Therefore the logical
% expression binarizes the Polygon trace and the multiplication generates
% an array of numbers.

PolygonTraceThres = 1.00*(PolygonTrace>50);

% Locations of the TTLs in channel 2 i.e. Polygon
% In the line below the original output of the function was to [peaks, locs]
% MATLAB suggested I use '~' instead of peaks because I was not using the
% variable peaks anywhere. Using a tilde instead of a variable ignores that
% particular function output and saves computing
%
%The thresholded trace has sharpe rise and fall and therefore is better
%in locating the TTL onset. Here a threshold of 50 is enough to
%generate precise TTL location.In TTL standard, the TTL is counted with the
%threshold of 0.8 V (whether a TTL is of 3.3 V or 5.0 V).

% MeanPeakDistance depends upon the inter stimulus interval and is going to
% be 1 sec (20000 points) for 29x29 and 3 sec (60000 points) for 10x10
% grid.

%Number of peaks is just an added check to get only as many TTLs deteced as
%there are stimuli.
[~, locs] = findpeaks(PolygonTraceThres,'MinPeakDistance',18000,'Npeaks',gridSize^2);

% Create a matrix in which each row corresponds to a section of 
% patch trace around the stimulus
PatchTracelets=zeros(length(locs),points);

%Fill in the matrix using locations of TTL peaks
for i=1:length(locs)
    PatchTracelets(i,:)=PatchTrace((locs(i)-pre*acqRate):(locs(i)+post*acqRate));
    % Baseline subtraction
    % A mean value is calculated between datapoints 100 and 400
    % this mean is then subtracted from the entire traceline
    % thus shifting the trace to zero.
    baseline = mean(PatchTracelets(i,100:400));
    PatchTracelets(i,:) = PatchTracelets(i,:)-baseline;
end

PolygonTracelets=zeros(length(locs),points);

%Fill in the matrix using locations of TTL peaks
for i=1:length(locs)
    PolygonTracelets(i,:)=PolygonTrace((locs(i)-pre*acqRate):(locs(i)+post*acqRate));
    % Baseline subtraction
    % A mean value is calculated between datapoints 100 and 400
    % this mean is then subtracted from the entire traceline
    % thus shifting the trace to zero.
    baseline = mean(PolygonTracelets(i,100:400));
    PolygonTracelets(i,:) = PolygonTracelets(i,:)-baseline;
end


%% Reshuffeling the tracelets according to the external order

orderedPatchTracelets = zeros(size(PatchTracelets));
orderedPolTracelets = zeros(size(PolygonTracelets));

for i=1:length(locs)
    % The remapping of squares is done here using a variable 'j' which maps
    % the index i onto the correct coordinate of the square from the coord
    % data
    j = coord(i);
    orderedPatchTracelets(j,:)= PatchTracelets(i,:);
    orderedPolTracelets(j,:)= PolygonTracelets(i,:);
end

%% Creating an example trace for the stimulus
timeTracelet = linspace(pre*(-1),post,points); % time tracelet

%creating an example stimulus trace
exStimTrace = zeros(1,points);
exStimTrace(pre*acqRate:pre*acqRate+pulseDur) = 5; % 5 is a safe value for the stim



%% Save the data

mkdir(ExptID)
ParsedFilePath = strcat(PathName,ExptID,'\');
cd(ParsedFilePath) %Change the working directory to the path
ParsedFile = strcat(ParsedFilePath,ExptID,'_Parsed_Tracelets_',num2str(gridSize),'.mat');

clear ans baseline coord* fid i j max* Patch* Path* Poly* Trace*

save(ParsedFile)

%% Generate images
figure;
subplot(1,2,1)
imagesc(orderedPatchTracelets)
subplot(1,2,2)
imagesc(orderedPolTracelets)

print('Tracelets.png','-dpng')

%% Run the grid analysis file
% Running the next script in the pipeline

run('Grid_Analysis.m')

