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
from likelihoodHelpers import log_likelihood_biGauss, log_likelihood_doublePeak, likelihood_ratio_doublePeak, likelihood_ratio_biGauss, biGauss, double_peak
import scipy, csv

def likelihoodfit(frame, omgeo, file_num, frame_num, csv_writer):
    print('Likelihood fit function called')

    mctree = frame["I3MCTree"]
    primary = mctree.primaries
    lepton = dataclasses.I3MCTree.first_child(mctree, primary[0].id)

    recoPulseMap = frame['I3RecoPulses']

    tot_infoList = []

    for omkey in recoPulseMap.keys():
        recoPulseList = recoPulseMap[omkey]
        recoPulse_timeList = np.array([recoPulse.time for recoPulse in recoPulseList])
        recoPulse_chargeList = np.array([recoPulse.charge for recoPulse in recoPulseList])

        '''
        Removing DOMs with hits less than 150 Hits
        '''
        if sum(recoPulse_chargeList) < 150:
            '''
            info_list = [file_num, frame_num, lepton.type, omkey[1], omkey[0], 0, 0, 0,
                        0, 0, 0, 0,
                        0, 0,0, 0,
                        0, 0,0, 0,
                        0, 0, 0, 0, 0, 0, 0]

            csv_writer.writerow({'file': info_list[0], 'frame': info_list[1], 'lepton_type': info_list[2], 'DOM': info_list[3], 'string': info_list[4],
                                'binEntries_mean': info_list[5], 'success_biGauss': info_list[6], 'success_doublePeak': info_list[7],
                                'biGauss_pos': info_list[8], 'biGauss_wid': info_list[9], 'biGauss_rat': info_list[10], 'biGauss_amp': info_list[11],
                                'doublePeak_pos1': info_list[12], 'doublePeak_wid1': info_list[13], 'doublePeak_rat1': info_list[14], 'doublePeak_amp1': info_list[15],
                                'doublePeak_pos2': info_list[16], 'doublePeak_wid2': info_list[17], 'doublePeak_rat2': info_list[18], 'doublePeak_amp2': info_list[19],
                                'area_data': info_list[20], 'area_biGauss_fit': info_list[21], 'area_doublePeak_fit': info_list[22],
                                'gof_biGauss': info_list[23], 'gof_doublePeak': info_list[24], 'dof_biGauss' : info_list[25], 'dof_doublePeak' : info_list[26]})
            '''

            continue

        '''
        Calculating the mean and removing the tails
        '''

        #mean = recoPulse_timeList.mean()
        mean = sum(recoPulse_timeList*recoPulse_chargeList)/sum(recoPulse_chargeList) #mean is weighted
        select_time = recoPulse_timeList[(recoPulse_timeList > mean-50) & (recoPulse_timeList < mean+50)]
        select_charge = recoPulse_chargeList[(recoPulse_timeList > mean-50) & (recoPulse_timeList < mean+50)]
        #print('SELECT CHARGE', select_charge, select_time, mean, recoPulse_timeList, recoPulse_chargeList)

        if len(select_time) < 10:
            continue

        mean_select_time = sum(select_time*select_charge)/sum(select_charge)
        max_hitTimes = recoPulse_timeList[(recoPulse_timeList > (mean_select_time-100))&(recoPulse_timeList < (mean_select_time+100))]
        max_charge = recoPulse_chargeList[(recoPulse_timeList > (mean_select_time-100))&(recoPulse_timeList < (mean_select_time+100))]

        #[using zscore to remove the effect of outliers from the analysis]
        z = stats.zscore(max_hitTimes)
        max_hitTimes = max_hitTimes[(z>-1.2)&(z < 1.2)]
        max_charge = max_charge[(z>-1.2)&(z < 1.2)]

        if len(max_hitTimes) < 10:
            continue

        #Shifting mean to zero
        max_hitTimes_mean = sum(max_hitTimes*max_charge)/sum(max_charge)
        timestamps = max_hitTimes - max_hitTimes_mean
        final_mean = timestamps.mean()

        '''
        Histogramming the data from simulation
        '''
        print('Now Histogramming')
        bins = np.arange(min(timestamps), max(timestamps), 1)
        num, bin_edges = np.histogram(timestamps, bins=bins, weights=max_charge)
        bin_centers = (bin_edges[:-1]+bin_edges[1:])/2

        #removing bins with < 0 entries ---> Recommended for binned likelihood
        entries_in_bins = num[num > 0]
        bin_centers = bin_centers[num > 0]

        #Getting data for the chi2 fit
        chi2_entries_in_bins = entries_in_bins[entries_in_bins > 10]
        chi2_bin_centers = bin_centers[entries_in_bins > 10]

        #Degrees of freedom should be greater than zero!
        if len(chi2_entries_in_bins) < 8:
            '''
            info_list = [file_num, frame_num, lepton.type, omkey[1], omkey[0], -1, -1, -1,
                        -1, -1, -1, -1,
                        -1, -1,-1, -1,
                        -1, -1,-1, -1,
                        -1, -1, -1, -1, -1, -1, -1]

            csv_writer.writerow({'file': info_list[0], 'frame': info_list[1], 'lepton_type': info_list[2], 'DOM': info_list[3], 'string': info_list[4],
                                'binEntries_mean': info_list[5], 'success_biGauss': info_list[6], 'success_doublePeak': info_list[7],
                                'biGauss_pos': info_list[8], 'biGauss_wid': info_list[9], 'biGauss_rat': info_list[10], 'biGauss_amp': info_list[11],
                                'doublePeak_pos1': info_list[12], 'doublePeak_wid1': info_list[13], 'doublePeak_rat1': info_list[14], 'doublePeak_amp1': info_list[15],
                                'doublePeak_pos2': info_list[16], 'doublePeak_wid2': info_list[17], 'doublePeak_rat2': info_list[18], 'doublePeak_amp2': info_list[19],
                                'area_data': info_list[20], 'area_biGauss_fit': info_list[21], 'area_doublePeak_fit': info_list[22],
                                'gof_biGauss': info_list[23], 'gof_doublePeak': info_list[24], 'dof_biGauss' : info_list[25], 'dof_doublePeak' : info_list[26]})

            '''
            continue

        num_dataPoints = len(chi2_entries_in_bins)
        area_data = sum(entries_in_bins)
        mean_entries = entries_in_bins.mean()

        '''
        Fitting bifurcated Gaussian and double bifurcated gaussian to
        the mcpe hit time distributions for both tau and electron.
        '''

        #Single Peak

        nll = lambda *args: log_likelihood_biGauss(*args)
        initial_biGauss = np.array([final_mean, 50, 5, max(entries_in_bins)])
        bnds_biGauss = ((min(bin_centers), max(bin_centers)), (0, 500), (0, 10), (0, 1e6))
        soln_biGauss = minimize(log_likelihood_biGauss, initial_biGauss,
                                args=(entries_in_bins, bin_centers),
                                method='Powell',
                                bounds = bnds_biGauss)

        #Double Peak

        nll = lambda *args: log_likelihood_doublePeak(*args)
        initial_doublePeak = np.array([min(bin_centers)+10, 20, 1, max(entries_in_bins), final_mean, 20, 1, max(entries_in_bins)])
        bnds_doublePeak = ((min(bin_centers), final_mean), (0, 500), (0, 10), (0, 1e6),
                            (final_mean, max(bin_centers)), (0, 500), (0, 10), (0,1e6))
        soln_doublePeak = minimize(log_likelihood_doublePeak, initial_doublePeak,
                                    args=(entries_in_bins, bin_centers),
                                    method='Powell',
                                    bounds=bnds_doublePeak)

        dof_biGauss = num_dataPoints - 4
        dof_doublePeak = num_dataPoints - 8

        success_biGauss = soln_biGauss.success
        success_doublePeak = soln_doublePeak.success
        '''
        Calculating the Likelihood ratio for bifurcated gaussian
        and double double bifurcated gaussian
        '''
        LR_biGauss = likelihood_ratio_biGauss(bin_centers, entries_in_bins, soln_biGauss.x[0],
                                              soln_biGauss.x[1], soln_biGauss.x[2], soln_biGauss.x[3])
        LR_doublePeak = likelihood_ratio_doublePeak(bin_centers, entries_in_bins, soln_doublePeak.x[0],
                                                    soln_doublePeak.x[1],soln_doublePeak.x[2],
                                                    soln_doublePeak.x[3], soln_doublePeak.x[4],
                                                    soln_doublePeak.x[5], soln_doublePeak.x[6],
                                                    soln_doublePeak.x[7])


        '''
        (x, y) values for the fit
        '''
        x = bin_centers
        #x = np.linspace(0, max(bin_centers)+1e5, 1000)
        y_biGauss = biGauss(x, soln_biGauss.x[0],
                                soln_biGauss.x[1], soln_biGauss.x[2], soln_biGauss.x[3])
        y_doublePeak = double_peak(x, soln_doublePeak.x[0], soln_doublePeak.x[1],
                                            soln_doublePeak.x[2], soln_doublePeak.x[3], soln_doublePeak.x[4],
                                            soln_doublePeak.x[5], soln_doublePeak.x[6], soln_doublePeak.x[7])

        area_biGauss_fit = sum(y_biGauss)
        area_doublePeak_fit = sum(y_doublePeak)


        '''
        Goodness of fit
        '''

        #ChiSquare

        #chi2_biGauss = 2*LR_biGauss
        #chi2_doublePeak = 2*LR_doublePeak

        x = chi2_bin_centers
        #x = np.linspace(0, max(bin_centers)+1e5, 1000)
        y_biGauss = biGauss(x, soln_biGauss.x[0],
                                soln_biGauss.x[1], soln_biGauss.x[2], soln_biGauss.x[3])
        y_doublePeak = double_peak(x, soln_doublePeak.x[0], soln_doublePeak.x[1],
                                            soln_doublePeak.x[2], soln_doublePeak.x[3], soln_doublePeak.x[4],
                                            soln_doublePeak.x[5], soln_doublePeak.x[6], soln_doublePeak.x[7])

        chi2_biGauss = np.sum((chi2_bin_centers - y_biGauss)**2/y_biGauss)
        chi2_doublePeak = np.sum((chi2_bin_centers - y_doublePeak)**2/y_doublePeak)

        gof_biGauss = chi2_biGauss/dof_biGauss
        gof_doublePeak = chi2_doublePeak/dof_doublePeak

        info_list = [file_num, frame_num, lepton.type, omkey[1], omkey[0], mean_entries, success_biGauss, success_doublePeak,
                    soln_biGauss.x[0], soln_biGauss.x[1], soln_biGauss.x[2], soln_biGauss.x[3],
                    soln_doublePeak.x[0], soln_doublePeak.x[1],soln_doublePeak.x[2], soln_doublePeak.x[3],
                    soln_doublePeak.x[4], soln_doublePeak.x[5],soln_doublePeak.x[6], soln_doublePeak.x[7],
                    area_data, area_biGauss_fit, area_doublePeak_fit, gof_biGauss, gof_doublePeak, dof_biGauss, dof_doublePeak]

        #print('Declared InfoList - ', info_list)
        csv_writer.writerow({'file': info_list[0], 'frame': info_list[1], 'lepton_type': info_list[2], 'DOM': info_list[3], 'string': info_list[4],
                            'binEntries_mean': info_list[5], 'success_biGauss': info_list[6], 'success_doublePeak': info_list[7],
                            'biGauss_pos': info_list[8], 'biGauss_wid': info_list[9], 'biGauss_rat': info_list[10], 'biGauss_amp': info_list[11],
                            'doublePeak_pos1': info_list[12], 'doublePeak_wid1': info_list[13], 'doublePeak_rat1': info_list[14], 'doublePeak_amp1': info_list[15],
                            'doublePeak_pos2': info_list[16], 'doublePeak_wid2': info_list[17], 'doublePeak_rat2': info_list[18], 'doublePeak_amp2': info_list[19],
                            'area_data': info_list[20], 'area_biGauss_fit': info_list[21], 'area_doublePeak_fit': info_list[22],
                            'gof_biGauss': info_list[23], 'gof_doublePeak': info_list[24], 'dof_biGauss' : info_list[25], 'dof_doublePeak' : info_list[26]})

    return
