function [patchTrace, polygonTrace, timeTrace, metadata]= traceParse(traceFile,coordFile)
%TRACEPARSE Parses the data in Axon Text Format to matrices for each
%channel and also creates a time vector denoting the time period of each
%sweep.
%Rows    -> Number of sweeps i.e. row vectors = sweeps
%Columns -> Sample points per sweep i.e. sweep duration*sampling frequency
%
%   
%   See also IRTREND, HEATMAP, SCALEBAR, MAKEPLOTS, TRACEPLOT.

%   Aditya Asopa, GNU GPL 2019

%#ok<*EMCLS>
%#ok<*EMCA>
%#codegen

delimiter = '\t';
formatSpec = '%q%[^\n\r]';
metaDataRows = 11;
acqRate = 20;
pre = 200; 
post = 300;
baselineWindow = 1:4000;
% points = acqRate*(pre+post)+1; %Number of datapoints
% pulseDur = acqRate*10;

%The data in Tracefiles is tab deliminted and starts at row#12 and column#1
data = dlmread(traceFile,delimiter,metaDataRows,0);
data = data';
recSize = size(data);

%Meta Data
fileID = fopen(traceFile);
metadata = textscan(fileID, formatSpec, metaDataRows, 'Delimiter', delimiter, 'ReturnOnError', false, 'EndOfLine', '\r\n');
acqMode = strsplit(metadata{1}{3},'=');acqMode = string(acqMode{2});
% comments = strsplit(metadata{1}{4}, '=');%any comments recorded
% comments = string(comments{2});
fclose(fileID);

%coordiantes
fid = fopen(coordFile);
coord = textscan(fid,'%*u%*u%*u%u');
fclose(fid);
coord = coord{1};
gridSize = sqrt(length(coord)-40);


%PatchTrace is every even column starting from 2
%PolygonTrace is every odd column starting from 2

if acqMode == 'Episodic Stimulation'
    timeTrace = data(1,:); %first column of the data
    patchTrace = data((2:2:size(data,1)),:); %channel 1 (IN0), all even rows
    polygonTrace = data((3:2:size(data,1)),:); %channel 2 (IN1), all odd rows
    polygonTrace(polygonTrace<75)=0;
    polygonTrace(polygonTrace>75)=5;
    
    patchTrace(21:861,:) = matrixReorder(patchTrace(21:861,:),coord(21:861));
    polygonTrace(21:861,:) = matrixReorder(polygonTrace(21:861,:),coord(21:861));
    
    for i=1:size(patchTrace,1)
        patchTrace(i,:) = baselineCorrect(patchTrace(i,:),baselineWindow);
        polygonTrace(i,:) = baselineCorrect(polygonTrace(i,:),baselineWindow);
    end
    
elseif acqMode == 'Gap Free'
    timeTrace = linspace(pre*(-1),pre*(-1)+20000,20000);
    PatchTrace = data(2,:);    
    PolygonTrace = data(3,:);
    
    PolygonTraceThres = 1.00*(PolygonTrace>50);
    [~, locs] = findpeaks(PolygonTraceThres,'MinPeakDistance',18000,'Npeaks',gridSize^2);
    
    patchTracelets=zeros(length(locs),20000);
    polygonTracelets=zeros(length(locs),20000);
    
    for i=1:length(locs)
        point = locs(i)-pre*acqRate;
        j = coord(i);
        patchTracelets(j,:)=PatchTrace(point:point+19999);
        polygonTracelets(j,:)=PolygonTrace(point:point+19999);
        % Baseline subtraction
        baseline = mean(patchTracelets(j,1:4000));
        patchTracelets(j,:) = patchTracelets(j,:)-baseline;
        baseline = mean(polygonTracelets(j,1:4000));
        polygonTracelets(j,:) = polygonTracelets(j,:)-baseline;        
    end    
end

function trace = baselineCorrect(row,baselineWindow)
baseline = mean(row(baselineWindow));
trace = row-baseline;

function orderedMatrix = matrixReorder(matrix,coord)
orderedMatrix = zeros(size(matrix));
for i=1:size(matrix,1)
    j = coord(i);
    orderedMatrix(j,:)=matrix(i,:);
end


    