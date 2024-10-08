import numpy as np
from scipy.optimize  import curve_fit
from eidynamics.utils import PSP_start_time
from eidynamics import ephys_classes
from eidynamics import pattern_index
import pickle

def get_expected_response(cellFile):
    with open(cellFile,'rb') as fin:
        x = pickle.load(fin)

def alphafunc(t,A,tau,delta):
    T   = delta+t
    a   = 1/tau
    y   = A*(a*T)*np.exp(1-a*T)
    return y

def delayed_alpha(t,A,tau,delta):
    tdel = np.zeros(int(2e4*delta))
    T   = np.append(tdel,t)
    T = T[:len(t)]
    a   = 1/tau
    y   = A*(a*T)*np.exp(1-a*T)
    return y

def make_spot_profile(exptObj1sq):
    pd                      = exptObj1sq.extract_trial_averaged_data(channels=[2])[2] # 45 x 40000
    cell                    = exptObj1sq.extract_trial_averaged_data(channels=[0])[0] # 45 x 40000
    Fs                      = exptObj1sq.Fs
    numSweeps               = exptObj1sq.numSweeps
    spotCoords              = dict([(k, exptObj1sq.stimCoords[k]) for k in range(1,1+int(exptObj1sq.numSweeps/exptObj1sq.numRepeats))])
    IPI                     = exptObj1sq.IPI # 0.05 seconds
    firstPulseTime          = int(Fs*(exptObj1sq.stimStart)) # 4628 sample points
    secondPulseTime         = int(Fs*(exptObj1sq.stimStart + IPI)) # 5628 sample points

    avgResponseStartTime,_  = PSP_start_time(cell)   # 0.2365 seconds
    avgSecondResponseStartTime = avgResponseStartTime + IPI # 0.2865 seconds
    avgSynapticDelay        = avgResponseStartTime[0]-exptObj1sq.stimStart # ~0.0055 seconds

    spotExpectedDict        = {}

    for i in range(len(spotCoords)):
        spotPD_trialAvg     = pd[i,int(Fs*avgResponseStartTime):int(Fs*avgSecondResponseStartTime)] # 1 x 1000
        # spotCell_trialAvg   = cell[0][i,int(Fs*avgResponseStartTime):int(Fs*avgSecondResponseStartTime)] # 1 x 1000
        spotCell_trialAVG_pulse2pulse = cell[i,firstPulseTime:secondPulseTime+200]

        t                   = np.linspace(0,IPI+0.01,len(spotCell_trialAVG_pulse2pulse)) # seconds at Fs sampling
        T                   = np.linspace(0,0.4,int(0.4*Fs)) # seconds at Fs sampling
        popt,_              = curve_fit(delayed_alpha,t,spotCell_trialAVG_pulse2pulse,p0=([0.5,0.05,0.005])) #p0 are initial guesses A=0.5 mV, tau=50ms,delta=5ms
        A,tau,delta         = popt
        fittedSpotRes       = delayed_alpha(T,*popt)
        
        spotExpectedDict[spotCoords[i+1][0]] = [avgSynapticDelay,A,tau,delta,spotCell_trialAVG_pulse2pulse,fittedSpotRes]
    
    all1sqAvg           = np.mean(cell[:,firstPulseTime:secondPulseTime+200],axis=0)
    popt,_              = curve_fit(delayed_alpha,t,all1sqAvg,p0=([0.5,0.05,0.005])) #p0 are initial guesses A=0.5 mV, tau=50ms,delta=5ms
    A,tau,delta         = popt
    all1sqAvg_fittedSpotRes       = delayed_alpha(T,*popt)
    spotExpectedDict['1sqAvg']           = [avgSynapticDelay,A,tau,delta,all1sqAvg,all1sqAvg_fittedSpotRes]

    return spotExpectedDict

def find_frame_expected(exptObj,spotExpectedDict_1sq):
    stimFreq        = exptObj.stimFreq
    IPI             = exptObj.IPI # IPI of current freq sweep experiment
    numPulses       = exptObj.numPulses
    numSweeps       = exptObj.numSweeps
    numRepeats      = exptObj.numRepeats
    cell            = exptObj.extract_trial_averaged_data(channels=[0])[0][:,:20000] # 8 x 40000
    stimCoords      = dict([(k, exptObj.stimCoords[k]) for k in range(1,1+int(numSweeps/numRepeats))]) # {8 key dict}
    stimStart       = 0.2314
    Fs              = 2e4
    
    frameExpected   = {}
    for i in range(len(stimCoords)):
        coordsTemp = stimCoords[i+1] # nd array of spot coords
        frameID    = pattern_index.get_patternID(coordsTemp)
        numSq      = len(coordsTemp)
        firstPulseExpected = np.zeros((int(Fs*(IPI+0.01))))
        print(len(firstPulseExpected))
        firstPulseFitted   = np.zeros((int(Fs*(0.4))))# added 0.01 second = 10 ms to IPI, check line 41,43,68,69
        for spot in coordsTemp:            
            spotExpected        = spotExpectedDict_1sq[spot][4]
            spotFitted          = spotExpectedDict_1sq['1sqAvg'][5]
            firstPulseExpected  += spotExpected[:len(firstPulseExpected)]
            firstPulseFitted    += spotFitted[:len(firstPulseFitted)]
        avgSynapticDelay    = spotExpectedDict_1sq[spot][0]
        expectedResToPulses = np.zeros(len(cell[0,:]))
        fittedResToPulses = np.zeros(len(cell[0,:]))
        t1 = int(Fs*stimStart)
        for k in range(numPulses):
            
            t2 = t1+int(Fs*IPI+avgSynapticDelay)
            T2 = t1+int(0.4*Fs)
            window1 = range(t1,t2)
            window2 = range(t1,T2)
            expectedResToPulses[window1] += firstPulseExpected[:len(window1)]
            fittedResToPulses[window2]   += firstPulseFitted[:len(window2)]
            t1 = t1+int(Fs*IPI)

        frameExpected[frameID] = [numSq,stimFreq,exptObj.stimIntensity,exptObj.pulseWidth,expectedResToPulses,fittedResToPulses,firstPulseFitted,firstPulseExpected]

    return frameExpected