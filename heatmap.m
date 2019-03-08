function [peakMap, AUCMap, timeofpeakMap] = heatmap(patchTrace,peakThres,traceDur)

%---variables--------------------------------------------------------------
gridSize = 29;
acqRate = 20; %timepoints per millisecond
pre = 200; % 200 comes from the value of pre used in traceParse function
TTLstart = pre*acqRate; % number of points before TTL comes
Window = traceDur*acqRate;
gridTraces = patchTrace(21:861,:);



%---Grid Peak HeatMap------------------------------------------------------
peakMap=zeros(gridSize);
for i=1:size(gridTraces,1)
    peakMap(i)=max(gridTraces(i,TTLstart:TTLstart+Window));
end
peakMap = peakMap';
peakMap(peakMap>peakThres)=peakThres;

%----Grid AUC Heatmap------------------------------------------------------
AUCMap=zeros(gridSize);
for i=1:size(gridTraces,1)
    AUCMap(i)=trapz(gridTraces(i,TTLstart:TTLstart+Window));
end
AUCMap = AUCMap';

%----Slope Heatmap-(Future update)-----------------------------------------

%----Time of Peak Heatmap--------------------------------------------------
timeofpeakMap=zeros(gridSize);
[~, TraceletPeakTime]= max(gridTraces(:,TTLstart:TTLstart+Window),[],2);
for i=1:size(gridTraces,1)
    timeofpeakMap(i)= TraceletPeakTime(i);
end
timeofpeakMap = timeofpeakMap';
end
%--------------------------------------------------------------------------