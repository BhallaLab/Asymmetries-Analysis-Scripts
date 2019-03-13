function tracePlot(patchTrace,polygonTrace,timeTrace,ExptID)
pre = 200;
post = 300;
pulseDur = 20;
acqRate = 20;
gridSize = 29;

figure;
axis([-1*pre post 1.1*min(min(patchTrace)) 1.1*max(max(patchTrace))])
figurePSTH=gcf;
figurePSTH.Units='normalized';
figurePSTH.OuterPosition=[0 0 1 1];

for row=1:size(patchTrace,1)
    hold on
    plot(linspace(-pre,post,10000),patchTrace(row,1:acqRate*(pre+post)),'k','LineWidth',1)        
end
hold on;
TTLGraph = area(polygonTrace(1,:));
TTLGraph.FaceColor = 'red'; %setting colour to red
TTLGraph.FaceAlpha = 0.5;
TTLGraph.LineWidth = 1;
TTLGraph.EdgeColor = 'r';

exStimTrace = zeros(size(timeTrace));
exStimTrace(pre*acqRate:pre*acqRate+pulseDur*acqRate)=2;

hold on;
stimGraph = area(exStimTrace); %plotting the stimulus trace for reference
stimGraph.FaceColor = 'blue' ; %setting colour to be blue
stimGraph.FaceAlpha = 0.5 ; %setting transparency
stimGraph.LineWidth = 1;
stimGraph.EdgeColor = 'b';

title('Response traces from baseline')
xlabel('Time (ms)');
ylabel('mV');
response_traces = strcat(ExptID,'_response_traces_',num2str(gridSize),'x');
savefig(response_traces)

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