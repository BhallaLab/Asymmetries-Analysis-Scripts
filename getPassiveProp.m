function [IR, Cm, Tau] = getPassiveProp(recFile,exptSpecs)
% GETPASSIVEPROP takes data file of a current clamp recording and the pulse
% parameters (start time, end time and pulse amplitude) to compute the
% input resistance (IR), capacitance (Cm) and time constant (Tau) of the
% cell. IR is calculated from the change in membrane potential due to the
% current pulse given. Cm and Tau are calcualted from exponential charging
% of the the membrane.
%   
%   See also TRACEPARSE, HEATMAP, SCALEBAR, MAKEPLOTS, TRACEPLOT, IRTREND, IFPLOTTER.

tStart = exptSpecs.IRpulseStart;
tEnd = exptSpecs.IRpulseEND;
delI = exptSpecs.IRpulse;

data = dlmread(recFile,'\t',11,0);

rec = data(:,2:end);
timeTrace = data(:,1);
sampFreq = 0.001/(timeTrace(2)-timeTrace(1));
numPoints = sampFreq*100; %100 ms average

window1 = ((tStart*sampFreq)-numPoints):tStart*sampFreq; %how big the window of mean is? (datapoints)
window2 = ((tEnd*sampFreq)-numPoints):tEnd*sampFreq; %how big the window of mean is? (datapoints)
prePulseAvg = mean(rec(window1,:)); % Average membrane potential before the pulse starts
ssAvg = mean(rec(window2,:)); % Average membrane potential at stead state near the end of the pulse

delV = ssAvg-prePulseAvg; % change in membrane potential as a result of current pulse
tauV = (1-(1/exp(1)))*delV; % 63% of delV

tauThres = prePulseAvg+tauV;
tauTime = zeros(size(tauThres));

if numel(delI)==1
    iSteps = repmat(delI,1,size(rec,2));
else
    iSteps = delI;
end


for i = 1:size(tauThres,2)
    if iSteps(i)>0
        tauTime(i) = find(rec(:,i)>tauThres(i),1)/sampFreq;
    else
        tauTime(i) = find(rec(:,i)<tauThres(i),1)/sampFreq;
    end
end

IR = 1000*delV./delI; %Input Resistance in MOhm
Tau = tauTime-tStart; %Time Constant in ms
Cm = 1000*Tau./IR; %Cm in picoFarads

x = linspace(0,size(rec,1)/sampFreq, size(rec,1));

figure('units','normalized','outerposition',[0 0 1 1])
subplot(2,2,1)
plot(IR,'b')
title('Cell Input Resistance')
xlabel('trials')
ylabel('IR (M\Omega)')
axis([1 size(rec,2) 0 1.2*max(IR)])
str = strcat('Avg. IR = ', num2str(mean(IR)),' M\Omega');
text(2,20,str)


subplot(2,2,2)
plot(Tau,'r')
title('Cell Time Constant')
xlabel('Trials')
ylabel('Tau (ms)')
axis([1 size(rec,2) 0 1.2*max(Tau)])
str = strcat('Avg. \tau = ', num2str(mean(Tau)),' ms');
text(2,10,str)

subplot(2,2,3)
plot(Cm,'g')
title('Cell Capacitance')
xlabel('Trials')
ylabel('Cm (pF)')
axis([1 size(rec,2) 0.5*min(Cm) 1.5*max(Cm)])
str = strcat('Avg. Cm = ', num2str(mean(Cm)),' pF');
text(2,-20+min(Cm),str)

subplot(2,2,4)
hold on
plot(x,rec(:,1),'b')
plot(tauThres(1)*ones(size(rec,1),1),'-k')
plot(x(window1),rec(window1,1),'r','lineWidth',2)
plot(x(window2),rec(window2,1),'r','lineWidth',2)
title('Example Trace with Sampling Time Marked in Red')
xlabel('Time (ms)')
ylabel('Membrane Voltage')
axis([1 x(end) -2+min(rec(:,1)) 2+max(rec(:,1))])

print('passProp.png','-dpng')

end