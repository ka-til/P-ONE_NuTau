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
'''Plot'''
#######

fig, axs = plt.subplots(3, 1, sharex=True, gridspec_kw={'height_ratios': [1, 1, 1]})
#fig.set_figheight(10)
#fig.set_figwidth(30)
fig.subplots_adjust(hspace=0)

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

grad = np.linspace(1, 0.1, 200)
for frame in range(0, 1):
    print(j)
    primary = randFrame["NuGPrimary"]
    nuEnergy = primary.energy
    noiseDOMs = np.empty(200, dtype=object)
    physicsDOMs = np.empty(200, dtype=object)
    totalDOMs = np.empty(200, dtype=object)
    mcpeMap = randFrame['MCPESeriesMap']
    noiseMap = randFrame['NoiseSeriesMap']
    mcpeOMKeys = mcpeMap.keys()

    i = 0
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
            physicsDOMs[i] = np.array(timeList_ph)
            noiseDOMs[i] = np.array(timeList_nh)
            mcpeTime = np.append(mcpeTime, timeList_ph)
            #if len(noiseDOMs[i]) == 0:
                #print('true')
        else:
            noiseList = noiseMap[modkey]
            timeList_nh = [mcpe.time for mcpe in noiseList]
            noiseDOMs[i] = np.array(timeList_nh)
            physicsDOMs[i] = np.array([])
            #if len(physicsDOMs[i]) == 0:
                #print('true in else')
                #print(modkey)

        totalDOMs[i] = physicsDOMs[i]
        #rint(totalDOMs[i].shape)
        totalDOMs[i] = np.append(totalDOMs[i], noiseDOMs[i])
        #print(totalDOMs[i].shape)

        bins = np.linspace(medVal-(72e2/2), medVal+(72e2/2), 11)
        numHitsPhys, binsPhys = np.histogram(physicsDOMs[i], bins=bins)
        numHitsNoise, binsNoise = np.histogram(noiseDOMs[i], bins=bins)
        totHits, binsTot = np.histogram(totalDOMs[i], bins=bins)

        noisePE = np.append(noisePE, numHitsNoise)
        physPE = np.append(physPE, numHitsPhys)
        totPE = np.append(totPE, totHits)
        bNoise = np.append(bNoise, (binsNoise[1:]-binsNoise[:-1])/2)
        bPhys = np.append(bPhys, (binsPhys[1:]-binsNoise[:-1])/2)
        bTot = np.append(bTot, (binsTot[1:]-binsTot[:-1])/2)
        print('x-', len(bNoise), 'y-', len(noisePE))


        s = s+len(physicsDOMs[i])
        sn = sn+len(noiseDOMs[i])
        #print(s, sn)
        axs[0].hist(physicsDOMs[i], bins=bins, facecolor='deepskyblue', alpha = 0.2, log=True, label = 'Physics Hits')
        axs[1].hist(noiseDOMs[i], bins=bins, facecolor='deepskyblue', alpha=0.2, label='Noise Hits')
        axs[2].hist(totalDOMs[i], bins=bins, facecolor='deepskyblue', alpha = 0.2,log=True, label='Phys+Noise')
        #axs[0].plot(binsPhys[:-1], np.log10(numHitsPhys), '-')
        #axs[1].plot(binsNoise[:-1], numHitsNoise, '-')
        #axs[2].plot(binsTot[:-1], np.log10(totHits), '-')

        i=i+1

    #medVal = np.mean(mcpeTime)
    #print(medVal)

    noiseFrame[j] = noiseDOMs
    physicsFrame[j] = physicsDOMs

    j=j+1


a_string = filename.replace(".i3.gz", "")
a_string = a_string.replace("step_4_", "")
a_string = a_string.replace("_medium_water_custom_mDOM_noise", "")

axs[2].set_xlabel('TimeWindow (7.2 us)', fontsize=14)
axs[0].text(medVal+1000, max(physPE)/2, 'Physics', fontsize=10)
axs[1].text(medVal+1000, max(noisePE)/2, 'Noise', fontsize=10)
axs[2].text(medVal+1000, max(totPE)/2, 'Physics+Noise', fontsize=10)
axs[0].text(medVal+700, 1, 'Neutrino Energy='+str(int(nuEnergy))+'GeV', style='italic',
bbox={'facecolor': 'red', 'alpha': 0.5, 'pad': 1})

fig.suptitle('PE Hits '+a_string, fontsize=20)
#plt.colorbar()
plt.savefig("NoiseDistribution_"+a_string+".pdf", dpi = 200)
plt.clf()
