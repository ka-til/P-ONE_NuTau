#Testing mcpes in I3Files

from icecube import dataclasses, dataio, simclasses
from icecube.icetray import I3Units, OMKey, I3Frame
from icecube.dataclasses import ModuleKey
import numpy as np
import random, re, os
from optparse import OptionParser
from os.path import expandvars
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm

#######
'''Geometry'''
######

#gcdPath = '/home/users/akatil/P-ONE/GCD_files/PONE_Phase1.i3.gz'
geofile = dataio.I3File('/home/users/akatil/P-ONE/GCD_files/PONE_Phase1.i3.gz')
cframe = geofile.pop_frame(I3Frame.Calibration)
geometry = cframe["I3Geometry"]
geoMap = geometry.omgeo
#######

#######
'''I3 File'''
#######

directory='/data/p-one/akatil/step_4_medium_water/'

def choose_file():
    filename = random.choice(os.listdir("/data/p-one/akatil/step_4_medium_water/"))
    if os.path.getsize(directory+str(filename)) < 100000:
        return choose_file()
    else:
        return filename

filename=choose_file()#'step_4_58_medium_water_custom_mDOM_noise.i3.gz'
file = dataio.I3File(directory+filename)

frames = []
while(file.more()):
    frames.append(file.pop_daq())

noiseFrame = np.empty(len(frames), dtype=object)
physicsFrame = np.empty(len(frames), dtype=object)

randFrame = random.choice(frames)
print(randFrame)

j = 0
noisePE = ([])
physPE = ([])
totPE = ([])
bNoise = ([])
bPhys = ([])
bTot  = ([])

physHits = ([])
noiseHits = ([])
physDOMs = ([])
noiseDOMs = ([])

grad = np.linspace(1, 0.1, 200)
for frame in range(0, 1):
    print(j)
    primary = randFrame["NuGPrimary"]
    nuEnergy = primary.energy
    mcpeMap = randFrame['MCPESeriesMap']
    noiseMap = randFrame['NoiseSeriesMap']
    mcpeOMKeys = mcpeMap.keys()

    i = 1
    mcpeTime = ([])

    for omkey in mcpeMap.keys():
        mcpeList = mcpeMap[omkey]
        timeList_ph = [mcpe.time for mcpe in mcpeList]
        mcpeTime = np.append(mcpeTime, timeList_ph)
    medVal = np.mean(mcpeTime)

    s, sn = 0, 0
    for modkey in geoMap.keys():
        if modkey in mcpeOMKeys:
            mcpeList = mcpeMap[modkey]
            noiseList = noiseMap[modkey]
            timeList_ph = [mcpe.time for mcpe in mcpeList]
            timeList_nh = [mcpe.time for mcpe in noiseList]
            physHits = np.append(physHits, timeList_ph)
            noiseHits = np.append(noiseHits, timeList_nh)
            physDOMs = np.append(physDOMs, np.ones(len(timeList_ph))*i)
            noiseDOMs = np.append(noiseDOMs, np.ones(len(timeList_nh))*i)
            #if len(noiseDOMs[i]) == 0:
                #print('true')
        else:
            noiseList = noiseMap[modkey]
            timeList_nh = [mcpe.time for mcpe in noiseList]
            noiseHits = np.append(noiseHits, timeList_nh)
            noiseDOMs = np.append(noiseDOMs, np.ones(len(timeList_nh))*i)
            #if len(physicsDOMs[i]) == 0:
                #print('true in else')
                #print(modkey)

        #print(totalDOMs[i].shape)

        i=i+1
    j=j+1

#######
'''Plot'''
#######

fig, axs = plt.subplots(3, 1, sharex=True, gridspec_kw={'height_ratios': [1, 1, 1]})
fig.set_figheight(10)
fig.set_figwidth(15)
#fig.subplots_adjust(hspace=0)

a_string = filename.replace(".i3.gz", "")
a_string = a_string.replace("step_4_", "")
a_string = a_string.replace("_medium_water_custom_mDOM_noise", "")

binsDOMs = np.arange(1, 201, 1)
bins = np.linspace(medVal-(92e2/2), medVal+(92e2/2), 11)
a0 = axs[0].hist2d(noiseDOMs, noiseHits, bins=[binsDOMs, bins], cmap=plt.cm.Reds)
axs[0].set_ylabel('Noise Hits Time(us)', fontsize=15)
axs[0].tick_params(axis="x", labelsize=15)
axs[0].tick_params(axis="y", labelsize=15)
fig.colorbar(a0[3], ax=axs[0])

a1 = axs[1].hist2d(physDOMs, physHits, bins=[binsDOMs, bins], cmap=plt.cm.Reds, norm = LogNorm())
axs[1].set_ylabel('Physics Hits Time(us)', fontsize=15)
axs[1].tick_params(axis="x", labelsize=15)
axs[1].tick_params(axis="y", labelsize=15)
fig.colorbar(a1[3], ax=axs[1])

totHits = np.append(physHits, noiseHits)
totDOMs = np.append(physDOMs, noiseDOMs)

a2 = axs[2].hist2d(totDOMs, totHits, bins=[binsDOMs, bins], cmap=plt.cm.Reds, norm = LogNorm())
axs[2].set_ylabel('Physics+Noise Time(us)', fontsize=15)
axs[2].tick_params(axis="x", labelsize=15)
axs[2].tick_params(axis="y", labelsize=15)
fig.colorbar(a2[3], ax=axs[2])

axs[2].set_xlabel('DOMS - File No: '+a_string, fontsize=25)
plt.tight_layout()
plt.savefig("NoiseDistribution_"+a_string+"_2D.pdf", dpi = 200)
plt.clf()
