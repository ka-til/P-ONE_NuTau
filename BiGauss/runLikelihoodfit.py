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

tauLLR2 = ([])
eLLR2 = ([])

tauTimeDiff = ([])
eTimeDiff = ([])

tauTimeDiff2 = ([])
eTimeDiff2 = ([])

rejected = ([])

file_num = ([])
frame_num = ([])
select_dom = ([])
e_in_amp = ([])
min_fail = ([])
weird_timeDiff = ([])

for i in range(0, 2000):
    print(i)
    file = dataio.I3File('/data/p-one/akatil/step_4_medium_water/NuTau_NuE_20Events/step_4_'+str(i)+'_medium_water_custom_mDOM_noise.i3.gz')

    f = 0
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
            numHitsinDOM = np.append(numHitsinDOM, tot_timeList) #this should be len(tot_timeList) ----> this could have affected the analysis!!!!

            #print(numHitsinDOM)

        numHits = np.append(numHits, sum(numHitsinDOM))
        log_numHits = np.log10(numHits[numHits > 0])
        print(log_numHits)
        if log_numHits >= float(args.cmin) and log_numHits < float(args.cmax):
            num_doms_selected, error_in_amp, minimizer_fail, weirdTimeDifferences, val, timeDiff, val2, timeDiff2, string = likelihoodfit(frame, omgeo)

            file_num = np.append(file_num, i)
            frame_num = np.append(frame_num, f)
            select_dom = np.append(select_dom, len(num_doms_selected))
            e_in_amp = np.append(e_in_amp, len(error_in_amp))
            min_fail = np.append(min_fail, len(minimizer_fail))
            weird_timeDiff = np.append(weird_timeDiff, weirdTimeDifferences)

            if string == 'tau':
                tauLLR = np.append(tauLLR, val)
                tauLLR2 = np.append(tauLLR2, val2)
                tauTimeDiff = np.append(tauTimeDiff, timeDiff)
                tauTimeDiff2 = np.append(tauTimeDiff2, timeDiff2)
            if string == 'e':
                eLLR = np.append(eLLR, val)
                eLLR2 = np.append(eLLR2, val2)
                eTimeDiff = np.append(eTimeDiff, timeDiff)
                eTimeDiff2 = np.append(eTimeDiff2, timeDiff2)
            if string == 'rejected events':
                rejected = np.append(rejected, val)
        f += 1

'''
np.savetxt('/data/p-one/akatil/analysis/NuTau_NuE_20Events_allCharge_2DelLLR_100_cmin'+str(args.cmin)+'_cmax_'+str(args.cmax)+ '_tau.csv', tauLLR, delimiter=',')
np.savetxt('/data/p-one/akatil/analysis/NuTau_NuE_20Events_allCharge_2DelLLR_100_cmin'+str(args.cmin)+'_cmax_'+str(args.cmax)+ '_e.csv', eLLR, delimiter=',')

np.savetxt('/data/p-one/akatil/analysis/NuTau_NuE_20Events_allCharge_2DelLLRall_100_cmin'+str(args.cmin)+'_cmax_'+str(args.cmax)+ '_tau.csv', tauLLR2, delimiter=',')
np.savetxt('/data/p-one/akatil/analysis/NuTau_NuE_20Events_allCharge_2DelLLRall_100_cmin'+str(args.cmin)+'_cmax_'+str(args.cmax)+ '_e.csv', eLLR2, delimiter=',')

np.savetxt('/data/p-one/akatil/analysis/NuTau_NuE_20Events_allCharge_timeDiff_100_cmin'+str(args.cmin)+'_cmax_'+str(args.cmax)+ '_tau.csv', tauTimeDiff, delimiter=',')
np.savetxt('/data/p-one/akatil/analysis/NuTau_NuE_20Events_allCharge_timeDiff_100_cmin'+str(args.cmin)+'_cmax_'+str(args.cmax)+ '_e.csv', eTimeDiff, delimiter=',')

np.savetxt('/data/p-one/akatil/analysis/NuTau_NuE_20Events_allCharge_timeDiff_100_cmin'+str(args.cmin)+'_cmax_'+str(args.cmax)+ '_tau.csv', tauTimeDiff2, delimiter=',')
np.savetxt('/data/p-one/akatil/analysis/NuTau_NuE_20Events_allCharge_timeDiff_100_cmin'+str(args.cmin)+'_cmax_'+str(args.cmax)+ '_e.csv', eTimeDiff2, delimiter=',')

np.savetxt('/data/p-one/akatil/analysis/NuTau_NuE_20Events_allCharge_100_cmin'+str(args.cmin)+'_cmax_'+str(args.cmax)+ '_rejected_.csv', rejected, delimiter=',')
'''

np.savetxt('/data/p-one/akatil/analysis/NuTau_NuE_20Events_allCharge_'+str(args.cmin)+'_cmax_'+str(args.cmax)+ '_fileNum_.csv', file_num, delimiter=',')
np.savetxt('/data/p-one/akatil/analysis/NuTau_NuE_20Events_allCharge_'+str(args.cmin)+'_cmax_'+str(args.cmax)+ '_frameNum_.csv', frame_num, delimiter=',')
np.savetxt('/data/p-one/akatil/analysis/NuTau_NuE_20Events_allCharge_'+str(args.cmin)+'_cmax_'+str(args.cmax)+ '_numDOMSselected_.csv', select_dom, delimiter=',')
np.savetxt('/data/p-one/akatil/analysis/NuTau_NuE_20Events_allCharge_'+str(args.cmin)+'_cmax_'+str(args.cmax)+ '_errorInAmp_.csv', e_in_amp, delimiter=',')
np.savetxt('/data/p-one/akatil/analysis/NuTau_NuE_20Events_allCharge_'+str(args.cmin)+'_cmax_'+str(args.cmax)+ '_minimizerFail_.csv', min_fail, delimiter=',')
np.savetxt('/data/p-one/akatil/analysis/NuTau_NuE_20Events_allCharge_'+str(args.cmin)+'_cmax_'+str(args.cmax)+ '_weirdTimeDiff_.csv', weird_timeDiff, delimiter=',')
