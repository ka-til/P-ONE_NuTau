from icecube import dataclasses, dataio, simclasses
from icecube.icetray import I3Units, OMKey, I3Frame
from icecube.dataclasses import ModuleKey
import matplotlib.pyplot as plt
import numpy as np
import random

##########################################
#
#   Producing the list of OMs in Geometry
#
#########################################

gcdPath = '/home/users/akatil/P-ONE/GCD_files/PONE_Phase1.i3.gz'
geofile = dataio.I3File(gcdPath)
cframe = geofile.pop_frame(I3Frame.Calibration)
geometry = cframe["I3Geometry"]
geoMap = geometry.omgeo

#################################################
#
# Adding Noise data to simualtion
#
#################################################

file = dataio.I3File('/data/p-one/akatil/test/genHitsNoiseTrial.i3.gz')
#timestamps = loadFile('/data/p-one/gaertner/akanksha/20200424_222144_UTC_SDOM1_MUON_EMBLA_DARK_RUN1_60s_2020-04-24_2222_20115222154.corrected.txt')

frame=file.pop_frame()

for f in range(0, 1):
    mcpeMap = frame['MCPESeriesMap']
    noiseMap = frame['NoiseSeriesMap']
    mcpeOMKeys = mcpeMap.keys()
    print(mcpeOMKeys)
    mcpeTime_ph = ([])
    for omkey in mcpeMap.keys():
        mcpeList = mcpeMap[omkey]
        noiseList = noiseMap[omkey]
        timeList_ph = [mcpe.time for mcpe in mcpeList]
        mcpeTime_ph = np.append(mcpeTime_ph, timeList_ph)
        medVal = np.mean(mcpeTime_ph)
        timeList_nh = [mcpe.time for mcpe in noiseList]

        print(timeList_nh, 'noiseHits')
        print(timeList_ph, 'physicsHits')

        if omkey == OMKey(8,8,0):
            bins = np.linspace(medVal-(72e2/2), medVal+(72e2/2), 11)
            plt.hist(timeList_nh, bins = bins, density = False, log=False, histtype='step', label = 'Noise Hits')
            plt.hist(timeList_ph, bins = bins, density = False, log=False, histtype='step', label = 'Physics Hits')
            plt.legend()
            plt.xlabel("Time(ns)")
            plt.title('Hits in a single DOM')
            plt.savefig('/home/users/akatil/P-ONE/git/PONE_NuTau/Noise/Noise_Physics_Hits_correct.pdf', dpi=200)
            plt.clf()
    print ('median calculated')
