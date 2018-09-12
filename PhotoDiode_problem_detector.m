
acqRate = 20;
pre = 25; 
post = 15;
points = acqRate*(pre+post)+1; %Number of datapoints


%% Normal Way

TimeTrace = linspace(0,length(PatchTrace)/(1000*acqRate),length(PatchTrace));

maxPolygon = max(PolygonTrace);

[~, locs] = findpeaks(PolygonTrace,'MinPeakHeight',0.95*maxPolygon,'MinPeakDistance',800);

PatchTracelets=zeros(length(locs),points);

for i=1:length(locs)
    PatchTracelets(i,:)=PatchTrace((locs(i)-pre*acqRate):(locs(i)+post*acqRate));
    baseline = mean(PatchTracelets(i,100:400));
    PatchTracelets(i,:) = PatchTracelets(i,:)-baseline;
end

PolygonTracelets=zeros(length(locs),points);
for i=1:length(locs)
    PolygonTracelets(i,:)=PolygonTrace((locs(i)-pre*acqRate):(locs(i)+post*acqRate));
    baseline = mean(PolygonTracelets(i,100:400));
    PolygonTracelets(i,:) = PolygonTracelets(i,:)-baseline;
end

figure;subplot(1,2,1);imagesc(PatchTracelets);subplot(1,2,2);imagesc(PolygonTracelets)
print('Comparison_of_traces.png','-dpng')

% New way for TTL detection

PolygonTraceDiff = diff(PolygonTrace);

maxPolygonDiff = max(PolygonTraceDiff);

[~, locsDiff] = findpeaks(PolygonTraceDiff,'MinPeakHeight',0.2*maxPolygonDiff,'MinPeakDistance',800);

PatchTraceletsDiff=zeros(length(locsDiff),points);

for i=1:length(locsDiff)
    PatchTraceletsDiff(i,:)=PatchTrace((locsDiff(i)-pre*acqRate):(locsDiff(i)+post*acqRate));
    baseline = mean(PatchTraceletsDiff(i,100:400));
    PatchTraceletsDiff(i,:) = PatchTraceletsDiff(i,:)-baseline;
end

PolygonTraceletsDiff =zeros(length(locsDiff),points);

for i=1:length(locsDiff)
    PolygonTraceletsDiff(i,:)=PolygonTrace((locsDiff(i)-pre*acqRate):(locsDiff(i)+post*acqRate));
    baseline = mean(PolygonTraceletsDiff(i,100:400));
    PolygonTraceletsDiff(i,:) = PolygonTraceletsDiff(i,:)-baseline;
end

figure;subplot(1,2,1);imagesc(PatchTraceletsDiff);subplot(1,2,2);imagesc(PolygonTraceletsDiff)
print('Comparison_of_traces_Diff.png','-dpng')

%% The new method of using the differentiation of polygon trace to detect the TTL, works.
%  I want to see how different are the locations of TTL in both the cases.

gap = locs - locsDiff;

% combined matrix of polygon and Photodiode Data
combNew = [10*PatchTraceletsDiff ;PolygonTraceletsDiff]; %multiplying photodiode matrix to show it better on the heatmap
imagesc(combNew)
print('Time_of_TTL_and_PhotodiodeTrace.png','-dpng')

combOld = [10*PatchTracelets(:,400:800) ;PolygonTracelets(:,400:800)];
imagesc(combOld)
print('Time_of_TTL_and_PhotodiodeTrace_Old.png','-dpng')

%% Trace plots for the two methods
figure;

for i = 1:500
plot(20*PatchTracelets(i,400:800),'k')
hold on
plot(PolygonTracelets(i,400:800),'b')
end

figure;
for i = 1:500
plot(20*PatchTraceletsDiff(i,400:800),'k')
hold on
plot(PolygonTraceletsDiff(i,400:800),'b')
end



