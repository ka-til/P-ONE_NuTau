'''
Using Liklihood Ratio to seperate between the Taus and the electrons
'''

from icecube import icetray, dataio, dataclasses, simclasses, clsim
from icecube.icetray import I3Units, OMKey, I3Frame
from icecube.dataclasses import ModuleKey
from os.path import expandvars
import numpy as np
#import matplotlib.pylab as plt
from scipy import stats
from scipy.optimize import minimize
from scipy.stats.distributions import chi2
from likelihoodHelpers import log_likelihood_biGauss, log_likelihood_doublePeak, likelihood_ratio_doublePeak, likelihood_ratio_biGauss
import scipy

def likelihoodfit(frame, omgeo):
    print('Likelihood fit function called')

    tau_timeDiff = ([])
    #tau_pVal = ([])
    tau_LRR = ([])

    e_timeDiff = ([])
    #e_pVal = ([])
    e_LRR = ([])

    minimizer_fail = ([])
    error_in_amp = ([])
    weirdTimeDifferences = ([])
    num_doms_selected = ([])
    #file = dataio.I3File(str(infile))

    mctree = frame["I3MCTree"]
    primary = mctree.primaries
    lepton = dataclasses.I3MCTree.first_child(mctree, primary[0].id)

    mcpeMap = frame['MCPESeriesMap']
    noiseMap = frame['NoiseSeriesMap']
    recoPulseMap = frame['I3RecoPulses']

    #for omkey in mcpeMap.keys():

    for omkey in recoPulseMap.keys():
        oKey = omgeo.get(omkey)

        '''
        Obtaining the timeList

        noise_mcpeList = noiseMap[omkey]
        noise_timeList = np.array([mcpe.time for mcpe in noise_mcpeList])
        mcpeList = mcpeMap[omkey]
        timeList = np.array([mcpe.time for mcpe in mcpeList])
        tot_timeList = np.append(timeList, noise_timeList)

        '''

        recoPulseList = recoPulseMap[omkey]
        recoPulse_timeList = np.array([recoPulse.time for recoPulse in recoPulseList])
        recoPulse_chargeList = np.array([recoPulse.charge for recoPulse in recoPulseList])


        '''
        Removing DOMs with hits less than 250 Hits
        '''
        if len(recoPulse_timeList) < 100:
            #print('less than 250 hits - 1')
            continue


        '''
        Calculating the mean and removing the tails
        '''

        mean = recoPulse_timeList.mean()

        select_time = recoPulse_timeList[(recoPulse_timeList > mean-50) & (recoPulse_timeList < mean+50)]
        new_mean = select_time.mean()

        #if len(select_time) < 50:
            #print('less than 250 hits - 1')
            #continue

        bins = np.arange(mean - 100, mean + 100, 1)
        #bins = np.arange(min(select_time), max(select_time), 1)
        max_hitTimes = select_time[(select_time > (new_mean-40))&(select_time < (new_mean+40))]


        #[using zscore to remove the effect of outliers from the analysis]

        z = stats.zscore(max_hitTimes)

        max_hitTimes = max_hitTimes[(z>-1.6)&(z < 1.2)]
        new_mean = max_hitTimes.mean()

        #Shifting mean to zero

        timestamps = max_hitTimes - new_mean
        final_mean = timestamps.mean()

        if len(max_hitTimes) < 50:
            #print('less than 250 hits - 2')
            continue

        num_doms_selected = np.append(num_doms_selected, omkey)

        '''
        Histogramming the data from simulation
        '''

        bins = np.arange(min(timestamps), max(timestamps), 1)
        num, bin_edges = np.histogram(timestamps, bins=bins)
        bin_centers = (bin_edges[:-1]+bin_edges[1:])/2

        '''
        Fitting bifurcated Gaussian and double bifurcated gaussian to
        the mcpe hit time distributions for both tau and electron.
        '''

        #Single Peak

        nll = lambda *args: log_likelihood_biGauss(*args)
        initial_biGauss = np.array([final_mean, 50, 5, max(num)])
        bnds_biGauss = ((min(bin_centers), max(bin_centers)), (0, 100), (0, 10), (0, 1e6))
        soln_biGauss = minimize(log_likelihood_biGauss, initial_biGauss,
                                args=(num, bin_centers),
                                method='Powell',
                                bounds = bnds_biGauss)

        #Double Peak

        nll = lambda *args: log_likelihood_doublePeak(*args)
        initial_doublePeak = np.array([min(bin_centers)+10, 20, 1, max(num), final_mean, 20, 1, max(num)])
        bnds_doublePeak = ((min(bin_centers), final_mean), (0, 100), (0, 10), (0, 1e6),
                            (final_mean, max(bin_centers)), (0, 100), (0, 10), (0,1e6))
        soln_doublePeak = minimize(log_likelihood_doublePeak, initial_doublePeak,
                                    args=(num, bin_centers),
                                    method='Powell',
                                    bounds=bnds_doublePeak)


        '''
        Checking if there are any terrible fits
        '''
        amp1 = soln_doublePeak.x[3]
        amp2 = soln_doublePeak.x[7]

        #if amp1/amp2 < 1/4 and amp1/amp2 > 4:
            #print('Removing terrible fits')

        if amp1 < 0 or amp2 < 0:
            error_in_amp = np.append(error_in_amp, omkey)
            print('Error in amp')
            continue

        '''
        Removing DOMs whose minimization is not successful
        '''
        if soln_biGauss.success == False or soln_doublePeak.success == False:
            minimizer_fail = np.append(minimizer_fail, omkey)
            print('Removing DOMs whose minimization is not successful')
            continue

        '''
        Calculating the Likelihood ratio for bifurcated gaussian
        and double double bifurcated gaussian
        '''
        LR_biGauss = likelihood_ratio_biGauss(bin_centers[num>0], num[num>0], soln_biGauss.x[0],
                                              soln_biGauss.x[1], soln_biGauss.x[2], soln_biGauss.x[3])
        LR_doublePeak = likelihood_ratio_doublePeak(bin_centers[num>0], num[num>0], soln_doublePeak.x[0],
                                                    soln_doublePeak.x[1],soln_doublePeak.x[2],
                                                    soln_doublePeak.x[3], soln_doublePeak.x[4],
                                                    soln_doublePeak.x[5], soln_doublePeak.x[6],
                                                    soln_doublePeak.x[7])

        '''
        Calculating the p-value using the likelihood ratio
        '''
        #pVal_biGauss = chi2.sf(LR_biGauss, len(num) - 4)
        #pVal_doublePeak = chi2.sf(LR_doublePeak, len(num) - 8)

        '''
        Are there any messed up p-values?

        if pVal_biGauss != pVal_biGauss:
            print('BiGauss is not well defined - ', str(lepton.type))
            print('Minimisation - ', soln_biGauss.success)
            print('Degrees of Freedom - ', len(num) - 4)
            print('Log Likelihood - ', LR_biGauss)
        if pVal_doublePeak != pVal_doublePeak:
            print('double peak is not well defined - ', str(lepton.type))
            print('Minimisation - ', soln_doublePeak.success)
            print('Degrees of Freedom - ', len(num) - 8)
            print('Log Likelihood - ', LR_doublePeak)

        '''

        '''
        (x, y) values for the fit

        x = np.linspace(min(bin_centers), max(bin_centers), 1000)
        #x = np.linspace(0, max(bin_centers)+1e5, 1000)
        y_biGauss = biGauss(x, soln_biGauss.x[0],
                                soln_biGauss.x[1], soln_biGauss.x[2], soln_biGauss.x[3])
        y_doublePeak = double_peak(x, soln_doublePeak.x[0], soln_doublePeak.x[1],
                                            soln_doublePeak.x[2], soln_doublePeak.x[3], soln_doublePeak.x[4],
                                            soln_doublePeak.x[5], soln_doublePeak.x[6], soln_doublePeak.x[7])

        '''

        '''
        Calculating the time difference and p-value ratio of bigauss and double peak
        '''
        timeDifference_doublePeak = soln_doublePeak.x[4] - soln_doublePeak.x[0]
        #pVal_ratio = pVal_doublePeak/pVal_biGauss
        LRR = LR_doublePeak/LR_biGauss
        twoDelLRR = 2*abs(LR_doublePeak - LR_biGauss)
        #print(twoDelLRR)

        '''
        Note - Minimize for double peak can sometimes return a large
        time difference despite the fact the data actually contains only one peak
        Accounting for that below using Likelihood ratios. If the Likehihood ratios of
        the single peak and double peak are similar the timedifference will be forced to
        be zero
        '''

        if LRR >= 0.85 and LRR <=1.3:
            print('Time Difference before -', timeDifference_doublePeak)
            weirdTimeDifferences = np.append(weirdTimeDifferences, timeDifference_doublePeak)
            timeDifference_doublePeak = 0.
            #if abs(timeDifference_doublePeak) < 100:
                #plt.figure(figsize=(10,9))
                #_ = plt.hist(timestamps, bins=bins, histtype='step')
                #plt.title(str(lepton.type)+' timeDifference < 100 ' + str(abs(LR_doublePeak/LR_biGauss)))
                #plt.plot(x, y_biGauss, '--', c = 'r')
                #plt.plot(x, y_doublePeak, '--', c = 'k')
                #plt.axvline(final_mean, c = 'b')



        '''
        Checking if there are any terrible fits

        amp1 = soln_doublePeak.x[3]
        amp2 = soln_doublePeak.x[7]

        #if amp1/amp2 < 1/4 and amp1/amp2 > 4:
            #print('Removing terrible fits')

        if amp1 < 0 or amp2 < 0:
            print('Error in amp')

        '''

        '''
        Separating the time difference calculated above and appending the values
        '''

        '''
        Tau
        '''
        if lepton.type == 15 or lepton.type == -15:
            tau_timeDiff = np.append(tau_timeDiff, timeDifference_doublePeak)
            #tau_pVal = np.append(tau_pVal, pVal_ratio)
            tau_LRR = np.append(tau_LRR, twoDelLRR)
            #plt.title('E')

        '''
        Electron and Neutral Current
        '''

        if lepton.type == 11 or lepton.type == -11 or lepton.type == 12 or lepton.type == -12 or lepton.type == 16 or lepton.type == -16:

            e_timeDiff = np.append(e_timeDiff, timeDifference_doublePeak)
            #e_pVal = np.append(e_pVal, pVal_ratio)
            e_LRR = np.append(e_LRR, twoDelLRR)


    if len(tau_LRR) == 0 and len(e_LRR) == 0:
        print('event rejected')
        return num_doms_selected, error_in_amp, minimizer_fail, weirdTimeDifferences, 0, 0, 0, 0, 'rejected events'

    if len(tau_LRR) > 0:
        tau_max_val = max(tau_LRR)
    else:
        tau_max_val = -1e-9

    if len(e_LRR) > 0:
        e_max_val = max(e_LRR)
    else:
        e_max_val = -1e-9

    if tau_max_val > e_max_val:
        tau_LRR_2 = tau_LRR[tau_LRR != tau_max_val]
        #tau_max_val2 = max(tau_LRR[tau_LRR != tau_max_val])
        #tau_LRR_3 = tau_LRR_2[tau_LRR_2 != tau_max_val2]
        #tau_max_val3 = max(tau_LRR_2[tau_LRR_2 != tau_max_val2])

        tau_numMaxLLR = len(tau_LRR[tau_LRR == tau_max_val])
        tau_timeDiff_LR = tau_timeDiff[tau_LRR == tau_max_val]
        tau_max2DelLLR = ([])
        for j in range(0, tau_numMaxLLR):
            tau_max2DelLLR = np.append(tau_max2DelLLR, tau_max_val)
            return num_doms_selected, error_in_amp, minimizer_fail, weirdTimeDifferences, tau_max2DelLLR, tau_timeDiff_LR, tau_LRR, tau_timeDiff, 'tau'

    if e_max_val > tau_max_val:
        e_LRR_2 = e_LRR[e_LRR != e_max_val]
        #e_max_val2 = max(e_LRR[e_LRR != e_max_val])
        #e_LRR_3 = e_LRR_2[e_LRR_2 != e_max_val2]
        #e_max_val3 = max(e_LRR_2[e_LRR_2 != e_max_val2])

        e_numMaxLLR = len(e_LRR[e_LRR == e_max_val])
        e_timeDiff_LR = e_timeDiff[e_LRR == e_max_val]
        e_max2DelLLR = ([])
        for j in range(0, e_numMaxLLR):
            e_max2DelLLR = np.append(e_max2DelLLR, e_max_val)
            return num_doms_selected, error_in_amp, minimizer_fail, weirdTimeDifferences, e_max2DelLLR, e_timeDiff_LR, e_LRR, e_timeDiff, 'e'
    else:
        print('they are equal?!')
