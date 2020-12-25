import iminuit
from icecube import icetray, dataio, dataclasses, simclasses, clsim
from icecube.icetray import I3Units, OMKey, I3Frame
from icecube.dataclasses import ModuleKey
import numpy as np
import matplotlib.pylab as plt
from likelihoodHelpers import log_likelihood_biGauss, log_likelihood_doublePeak, likelihood_ratio_doublePeak, likelihood_ratio_biGauss, biGauss, double_peak
import corner
import sys

def numHits(frame):
    recoPulseMap = frame['I3RecoPulses']
    numHits = 0
    for omkey in recoPulseMap.keys():
        recoPulseList = recoPulseMap[omkey]
        recoPulse_chargeList = np.array([recoPulse.charge for recoPulse in recoPulseList])
        numHits = numHits+sum(recoPulse_chargeList)
    return numHits
