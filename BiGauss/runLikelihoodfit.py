from icecube import icetray, dataio, dataclasses, simclasses, clsim
from icecube.icetray import I3Units, OMKey, I3Frame
from icecube.dataclasses import ModuleKey
from os.path import expandvars
import numpy as np
#import matplotlib.pylab as plt
import argparse
from scipy import stats
from scipy.optimize import minimize
from scipy.stats.distributions import chi2
import scipy
from likelihoodRatioAnalysis import likelihoodfit

#file = dataio.I3File(str(infile))
parser = argparse.ArgumentParser(description = "Likelihood analysis")
parser.add_argument('-cmin', '--min_charge', dest = 'cmin', help= 'Minimum charge of the events')
parser.add_argument('-cmax', '--max_charge', dest = 'cmax', help= 'Maximum charge of the events')
args = parser.parse_args()

gcd_file = dataio.I3File('/home/users/akatil/P-ONE/GCD_files/PONE_Phase1.i3.gz')
cframe = gcd_file.pop_frame()
geometry = cframe["I3Geometry"]
omgeo = geometry.omgeo
print('loaded geometry')

tauLLR = ([])
eLLR = ([])

tauTimeDiff = ([])
eTimeDiff = ([])

rejected = ([])
for i in range(0, 2000):
    print(i)
    file = dataio.I3File('/data/p-one/akatil/step_4_medium_water/NuTau_NuE_20Events/step_4_'+str(i)+'_medium_water_custom_mDOM_noise.i3.gz')

    for frame in file:
        print('Starting')
        numHitsinDOM = ([])
        numHits = ([])
        mctree = frame["I3MCTree"]
        primary = mctree.primaries
        lepton = dataclasses.I3MCTree.first_child(mctree, primary[0].id)

        mcpeMap = frame['MCPESeriesMap']
        noiseMap = frame['NoiseSeriesMap']

        #looping through doms that have physics hits
        for omkey in mcpeMap.keys():
            oKey = omgeo.get(omkey)

            '''
            Obtaining the timeList
            '''
            noise_mcpeList = noiseMap[omkey]
            noise_timeList = np.array([mcpe.time for mcpe in noise_mcpeList])
            mcpeList = mcpeMap[omkey]
            timeList = np.array([mcpe.time for mcpe in mcpeList])
            tot_timeList = np.append(timeList, noise_timeList)
            numHitsinDOM = np.append(tot_timeList, noise_timeList)

            #print(numHitsinDOM)

        numHits = np.append(numHits, sum(numHitsinDOM))
        log_numHits = np.log10(numHits[numHits > 0])
        print(log_numHits)
        if log_numHits >= float(args.cmin) and log_numHits < float(args.cmax):
            val, timeDiff, string = likelihoodfit(frame, omgeo)
            if string == 'tau':
                tauLLR = np.append(tauLLR, val)
                tauTimeDiff = np.append(tauTimeDiff, timeDiff)
            if string == 'e':
                eLLR = np.append(eLLR, val)
                eTimeDiff = np.append(eTimeDiff, timeDiff)
            if string == 'rejected events':
                rejected = np.append(rejected, val)

np.savetxt('/data/p-one/akatil/analysis/NuTau_NuE_20Events_allCharge_2DelLLR_cmin'+str(args.cmin)+'_cmax_'+str(args.cmax)+ '_tau.csv', tauLLR, delimiter=',')
np.savetxt('/data/p-one/akatil/analysis/NuTau_NuE_20Events_allCharge_2DelLLR_cmin'+str(args.cmin)+'_cmax_'+str(args.cmax)+ '_e.csv', eLLR, delimiter=',')

np.savetxt('/data/p-one/akatil/analysis/NuTau_NuE_20Events_allCharge_timeDiff_cmin'+str(args.cmin)+'_cmax_'+str(args.cmax)+ '_tau.csv', tauTimeDiff, delimiter=',')
np.savetxt('/data/p-one/akatil/analysis/NuTau_NuE_20Events_allCharge_timeDiff_cmin'+str(args.cmin)+'_cmax_'+str(args.cmax)+ '_e.csv', eTimeDiff, delimiter=',')

np.savetxt('/data/p-one/akatil/analysis/NuTau_NuE_20Events_allCharge_cmin'+str(args.cmin)+'_cmax_'+str(args.cmax)+ '_rejected_.csv', rejected, delimiter=',')
