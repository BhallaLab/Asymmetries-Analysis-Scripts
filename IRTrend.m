function [meanIR, IR] = IRTrend(patchTrace,pulseCurrent)
acqRate = 20;
IRWindow = acqRate*500:acqRate*800;
IRTracelets = patchTrace(:,IRWindow);
IR = zeros(size(patchTrace,1),1);
for i=1:size(patchTrace,1)
    level2 = mean(IRTracelets(i,acqRate*250:acqRate*300));
    level1 = mean(IRTracelets(i,acqRate*50:acqRate*100));
    IR(i) = 1000*(level2-level1)/pulseCurrent; % input resistance in M? (mV/pA = G?)
end
meanIR = mean(IR);

