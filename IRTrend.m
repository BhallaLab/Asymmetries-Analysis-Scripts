function [meanIR, IR] = IRTrend(patchTrace,pulseCurrent)
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

