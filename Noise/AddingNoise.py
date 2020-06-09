from icecube import dataclasses, dataio, simclasses
from icecube.icetray import I3Units, OMKey, I3Frame
from icecube.dataclasses import ModuleKey
import matplotlib.pyplot as plt
import numpy as np
import random
import decimal

file = dataio.I3File('/data/p-one/akatil/step_3_medium_water/Custom/step_3_1000_medium_water_custom.i3.gz')

frames = []
while(file.more()):
    frames.append(file.pop_daq())

count = 0

for frame in frames:
    mcpeMap = frame['MCPESeriesMap']
    mcpeTime = ([])
    for mcpeList in mcpeMap.values():
        timeList = [mcpe.time for mcpe in mcpeList]
        mcpeTime = np.append(mcpeTime, timeList)

    if count == 14:
        mcpeList1=mcpeMap[OMKey(6,5,0)]
        mcpeList2=mcpeMap[OMKey(0,5,0)]
        mcpeHitTime1=[mcpe.time for mcpe in mcpeList1]
        mcpeHitTime2=[mcpe.time for mcpe in mcpeList2]

    count+=1

data = np.loadtxt('/data/p-one/gaertner/akanksha/20200424_222144_UTC_SDOM1_MUON_EMBLA_DARK_RUN1_60s_2020-04-24_2222_20115222154.corrected.txt')

channel = data[:, 0]
upperPMT = data[(channel==1)]
lowerPMT = data[(channel==5)]
risingEdge = upperPMT[upperPMT[:, 1] == 0]
risingEdgeL = lowerPMT[lowerPMT[:, 1] == 0]
timestamps = risingEdge[:, 2]

ten_us_window = np.arange(0, 60*1e9, 72e2)
seperateWindows = ([])
sWNoCorrection = ([])

for i in range(len(ten_us_window)):
    if i == 48:
        break

    vals = timestamps[(timestamps>ten_us_window[i]) & (timestamps<ten_us_window[i+1])]
    sWNoCorrection = np.append(sWNoCorrection, vals)
    if len(vals) == 0:
        continue
    else:
        vals = vals - ten_us_window[i]
    seperateWindows = np.append(seperateWindows,vals)

NoiseTime1 = min(mcpeHitTime1) + seperateWindows[0:23] - random.uniform(0, (72e2 - (max(mcpeHitTime1)-min(mcpeHitTime1))))
NoiseTime2 = min(mcpeHitTime2) + seperateWindows[23:48] - random.uniform(0, (72e2 - (max(mcpeHitTime2)-min(mcpeHitTime2))))

bins = np.linspace(min(NoiseTime1), max(NoiseTime2), 11)
plt.hist(NoiseTime1, bins = bins, density = False, log=False, histtype='step', label = 'Noise Hits')
plt.hist(mcpeHitTime1, bins = bins, density = False, log=False, histtype='step', label = 'Physics Hits')
plt.legend()
plt.xlabel("Time(ns)")
plt.title('Hits in a single DOM')
plt.savefig('/home/users/akatil/P-ONE/git/PONE_NuTau/Noise/Noise_Physics_Hits1.pdf', dpi=200)
plt.clf()

plt.hist(NoiseTime2, bins = bins, density = False, log=True, histtype='step', label = 'Noise Hits')
plt.hist(mcpeHitTime2, bins = bins, density = False, log=True, histtype='step', label = 'Physics Hits')
plt.legend()
plt.xlabel("Time(ns)")
plt.title('Hits in a single DOM')
plt.savefig('/home/users/akatil/P-ONE/git/PONE_NuTau/Noise/Noise_Physics_Hits2.pdf', dpi=200)
plt.clf()
