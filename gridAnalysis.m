function [IRtrend, peakMap, AUCMap, timetopeakMap, metadata] = gridAnalysis(recFile,coordFile,flank,pulseDur,intensity,gridSize,responseThres,cellID)
% GRIDANALYSIS takes:
% recFile   = File where data is saved (*.atf)
% coordFile = File with grid square coordinates for polygon frames (*.txt)
% flank     = Number of flanking frames
% pulseDur  = Duration of light pulse, stored in metadata
% intensity = Intensity of light pulse, stored in metadata
% gridSize  = Size of grid (29x29, 10x10...)
%
%  see also GETPASSIVEPROP, IFPLOTTER, MAKEHEATPLOTS, ASYMMETRYCALC.

%% Stage 1: Parsing recording file into data
delimiter = '\t';
formatSpec = '%q%[^\n\r]';
metaDataRows = 11;
isi = 1000; % inter-stim interval in ms, constant
baselineWindow = 1:4000;
pre = 200; %ms before the stim

disp('Parsing Data...')
disp('Reading Files...')

% Data
data = dlmread(recFile,delimiter,11,0);
data = data';
sampFreq = 0.001/(data(1,2)-data(1,1));

% Metadata
metadata = metadataParser(recFile);
acqMode = metadata{1,2};

% Coordinates
fid = fopen(coordFile);
coord = textscan(fid,'%*u%*u%*u%u');
fclose(fid);
coord = coord{1};
gridSizeCalc = sqrt(length(coord)-2*flank);
if ~gridSizeCalc == gridSize
    error('Grid size given does not look right. Aborting!')
end

disp('File reading complete.')

%% Stage 2: Traceparse
disp('Extracting channel data')

[patchTrace, polygonTrace, timeTrace]= traceParse(data,coord,flank);
disp('Parsing Done.')

%% Stage 3: Heatmaps
disp('Getting response maps...')
%responseThres = 30; not hardcoding the response threshold
[peakMap, AUCMap, timetopeakMap] = heatmap(patchTrace,responseThres);

disp('Response maps made.')

%% Stage 4: IR Trend
if acqMode == "Episodic Stimulation"
    pulseCurrent = -50;
    IRtrend = IRchange(patchTrace,pulseCurrent);
else
    IRtrend = NaN;
end

%% Trace plots
post = 100;
tracePlot(patchTrace, polygonTrace, timeTrace)

%% Local Functions

    function [patchTrace, polygonTrace, timeTrace]= traceParse(data,coord,flank)
    %PatchTrace is every even column starting from 2
    %PolygonTrace is every odd column starting from 2

    if acqMode == "Episodic Stimulation"
        timeTrace = data(1,:); %first column of the data
        patchTrace = data((2:2:size(data,1)),:); %channel 1 (IN0), all even rows
        polygonTrace = data((3:2:size(data,1)),:); %channel 2 (IN1), all odd rows
        
        polygonTrace(polygonTrace<75)=0; polygonTrace(polygonTrace>75)=5;

        stimFrames = 1+flank:size(patchTrace,1)-flank;
        patchTrace(stimFrames,:) = matrixReorder(patchTrace(stimFrames,:),coord(stimFrames));
        polygonTrace(stimFrames,:) = matrixReorder(polygonTrace(stimFrames,:),coord(stimFrames));

        for i=1:size(patchTrace,1)
            patchTrace(i,:) = baselineCorrect(patchTrace(i,:),baselineWindow);
            polygonTrace(i,:) = baselineCorrect(polygonTrace(i,:),baselineWindow);
        end

    elseif acqMode == "Gap Free"
        isiDatapoints = isi*sampFreq; %sample points per stim, usually 20000
        
        timeTrace = linspace(pre*(-1),pre*(-1)+isi,isiDatapoints);
        PatchTrace = data(2,:);    
        PolygonTrace = data(3,:);

        PolygonTraceThres = 1.00*(PolygonTrace>50);
        [~, locs] = findpeaks(PolygonTraceThres,'MinPeakDistance',18000,'Npeaks',gridSize^2);

        patchTracelets=zeros(length(locs),isiDatapoints);
        polygonTracelets=zeros(length(locs),isiDatapoints);

        for i=1:length(locs)
            point = locs(i)-pre*sampFreq;
            j = coord(i);
            
            patchTracelets(j,:)=PatchTrace(point:point+isiDatapoints-1);            
            patchTracelets(j,:) = baselineCorrect(patchTracelets(j,:),baselineWindow);
            
            polygonTracelets(j,:)=PolygonTrace(point:point+isiDatapoints-1);
            polygonTracelets(j,:) = baselineCorrect(polygonTracelets(j,:),baselineWindow);        
        end
        patchTrace = patchTracelets;
        polygonTrace = polygonTracelets;
    end
end

    function trace = baselineCorrect(row,baselineWindow)
        baseline = mean(row(baselineWindow));
        trace = row-baseline;
    end

    function orderedMatrix = matrixReorder(matrix,coord)
        orderedMatrix = zeros(size(matrix));
        for i=1:size(matrix,1)
            j = coord(i);
            orderedMatrix(j,:)=matrix(i,:);
        end
    end

    function [peakMap, AUCMap, timeofpeakMap] = heatmap(patchTrace,peakThres)
    % Update 28Jun19: I found out that the polygon frame is rotated
    % therefore there is no need to take transpose of the heatmaps. In
    % future updates, all the camera rotations will be addressed.
        %---variables--------------------------------------------------------------
        TTLstart = pre*sampFreq; % number of points before TTL comes
        Window = 100*sampFreq; % 100 ms window
        stimFrames = 1+flank:size(coord)-flank;
        gridTraces = patchTrace(stimFrames,:);


        %---Grid Peak HeatMap------------------------------------------------------
        peakMap=zeros(gridSize);
        for i=1:size(gridTraces,1)
            peakMap(i)=max(gridTraces(i,TTLstart:TTLstart+Window));
        end
        % Update 28Jun19: I found out that the polygon frame is rotated
        % therefore there is no need to take transpose.
%         peakMap = peakMap';
        peakMap(peakMap>peakThres)=peakThres;

        %----Grid AUC Heatmap------------------------------------------------------
        AUCMap=zeros(gridSize);
        for i=1:size(gridTraces,1)
            AUCMap(i)=trapz(gridTraces(i,TTLstart:TTLstart+Window));
        end
        % Update 28Jun19: I found out that the polygon frame is rotated
        % therefore there is no need to take transpose.
%         AUCMap = AUCMap';

        %----Slope Heatmap-(Future update)-----------------------------------------

        %----Time of Peak Heatmap--------------------------------------------------
        timeofpeakMap=zeros(gridSize);
        [~, TraceletPeakTime]= max(gridTraces(:,TTLstart:TTLstart+Window),[],2);
        for i=1:size(gridTraces,1)
            timeofpeakMap(i)= TraceletPeakTime(i);
        end
        % Update 28Jun19: I found out that the polygon frame is rotated
        % therefore there is no need to take transpose.
%         timeofpeakMap = timeofpeakMap';
    end

    function trend = IRchange(patchTrace,pulseCurrent)
    %IRTREND takes the patch clamp traces and gets the trend in input resistance
    %of the cell during the time period of the grid stimulation experiment.
    %patchTrace file is a matrix with row vectors representing sweeps.
    %
    %For measurement of input resistance, a current pulse is given to
    %hyperpolarize the neuron. The voltage change (in mV) is divided by the
    %pulse current (in pA) to obtain IR in GOhm. The function generates output in
    %MOhm by multiplying the values with 1000.
    %Input arguments -> patchTrace, pulse current
    %Output values   -> Average input resistance during the experiment and
    %                   A vector with input resistance values for every sweep

        sampFreq = 20;
        IRWindow = sampFreq*500:sampFreq*800;
        IRTracelets = patchTrace(:,IRWindow);
        trend = zeros(size(patchTrace,1),1);
        for i=1:size(patchTrace,1)
            level2 = mean(IRTracelets(i,sampFreq*250:sampFreq*300));
            level1 = mean(IRTracelets(i,sampFreq*50:sampFreq*100));
            trend(i) = 1000*(level2-level1)/pulseCurrent; % input resistance in MOhm (mV/pA = GOhm)
        end
    end

    function metadata = metadataParser(recFile)
        fileID = fopen(recFile);
        meta = textscan(fileID, formatSpec, metaDataRows, 'Delimiter', delimiter, 'ReturnOnError', false, 'EndOfLine', '\r\n');
        acqMode = strsplit(meta{1}{3},'=');
        acqMode = string(acqMode{2});
        comments = strsplit(meta{1}{4}, '=');%any comments recorded
        comments = string(comments{2});
        comments = strsplit(comments,','); %comments are comma separated
        
        metadata = {'Acquisition Mode' acqMode;'Comments' comments;'Intensity of Optical Stimulation' intensity;'Pulse Duration of Optical Stimulation' pulseDur};
        
        fclose(fileID);
    end

    function tracePlot(patchTrace,polygonTrace,timeTrace)
    traceStart = 50;
    traceEnd = 150;
    tracePoints = traceStart+traceEnd;
    windowStart = ((pre-traceStart)*sampFreq)+1; % from 3001st point
    windowEnd = (pre+traceEnd)*sampFreq; % to 7000th point
    
    traceWindow = windowStart:windowEnd;
    
    figure;
    axis([-1*traceStart traceEnd -5 1.1*max(max(patchTrace))])
    figurePSTH=gcf;
    figurePSTH.Units='normalized';
    figurePSTH.OuterPosition=[0 0 1 1];

    for row=1:size(patchTrace,1)
        hold on
        plot(linspace(-1*traceStart,traceEnd,sampFreq*tracePoints),patchTrace(row,traceWindow),'k','LineWidth',1)       
    end
    hold on;
    plot(linspace(-1*traceStart,traceEnd,sampFreq*tracePoints),polygonTrace(1,traceWindow),'r','LineWidth',1)

    title('Response traces from baseline')
    xlabel('Time (ms)');
    ylabel('mV');
    response_traces = strcat(cellID,'_response_traces_',num2str(gridSize),'x');
    print(response_traces,'-dpng')

    % %% PSTH
    % % for PSTH we need to get the timings of the peaks in all the traces. As
    % % there are going to be 841 traces for each trial, we can generate a
    % % distribution of peaks.
    % 
    % % 1. Find out the time points of peaks
    % % the function max finds maxima along columns, therefore the matrix of
    % % responses needs to be transposed.
    % 
    % [TraceletMax, TraceletPeakTime]= max(orderedPatchTracelets,[],2);
    % 
    % % 2. Get a histogram of timings
    % 
    % % Create figure
    % 
    % figure;
    % figurePSTH=gcf;
    % figurePSTH.Units='normalized'; %making the figure size in normalized units
    % figurePSTH.OuterPosition=[0 0 1 1]; %size of the figure -> full screen
    % 
    % % Create axes
    % axesPSTH = axes('Parent',figurePSTH);
    % hold(axesPSTH,'on');
    % box(axesPSTH,'on');
    % set(axesPSTH, 'fontsize', 10)
    % 
    % % Set the remaining axes properties
    % set(axesPSTH,'XTick',...
    %     [0 100 200 300 400 500 600 700 800 900 1000 1100 1200 1300 1400 1500 1600 1700 1800 1900 2000 2100 2200 2300 2400 2500 2600 2700 2800 2900 3000 3100 3200 3300 3400 3500 3600 3700 3800 3900 4000],...
    %     'XTickLabel',...
    %     {'-50','-45','-40','-35','-30','-25','-20','-15','-10','-5','0','5','10','15','20','25','30','35','40','45','50','55','60','65','70','75', '80','85','90','95','100','105','110','115','120','125','130','135','140','145','150'});
    % 
    % % Create histogram
    % histogram(TraceletPeakTime,'Parent',axesPSTH,'BinWidth',100);
    % title('PSTH')
    % xlabel('Time after Stimulus (ms)');
    % ylabel('Frequency');
    % xlim(axesPSTH,[0 4000])
    % 
    % % 3. Save histogram
    % 
    % PSTH = strcat(ExptID,'_PSTH_',num2str(gridSize),'x');
    % savefig(PSTH,'-dpng')

    end

end