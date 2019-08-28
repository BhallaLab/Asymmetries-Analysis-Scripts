function numSpikes = IFPlotter(recFile,exptSpecs)

tStart = exptSpecs(:).IFStepStart;
tEnd = exptSpecs(:).IFStepEnd;
iSteps = exptSpecs(:).IFcurrentSteps;

data = dlmread(recFile,'\t',11,0);
traceTime = data(:,1);
trials = data(:,2:end);
sampFreq = 0.001/(traceTime(2)-traceTime(1));
tStart = tStart*sampFreq;
tEnd = tEnd*sampFreq;
pulseDur = tEnd-tStart;

numSpikes = zeros(size(trials,2),1);
for i=1:size(trials,2)
    x = trials(tStart:tEnd,i);
    if max(x)>0
        [~,spikeLocs] = findpeaks(x,'MinPeakHeight',0.8*max(x),'MinPeakDistance',20);
        numSpikes(i) = length(spikeLocs);
    end
end

figure
plot(iSteps,numSpikes)
axis([min(iSteps)-20 max(iSteps)+20 0 max(numSpikes)+2])
xlabel('Time (ms)')
ylabel('Number of spikes')
print('IF_plot.png','-dpng')

% inR = zeros(size(trials,2),1);
% for i=1:size(trials,2)
%     x =trials(:,i);
%     if max(x)<0
%         inR(i) = getPassiveProp(x,tStart,tEnd,iSteps(i),sampFreq);        
%     end
% end

