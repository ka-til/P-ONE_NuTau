from icecube import icetray, dataio, dataclasses, simclasses, clsim
from icecube.icetray import I3Units, OMKey, I3Frame
from icecube.dataclasses import ModuleKey
from os.path import expandvars
import numpy as np
from scipy import stats
from iminuit import minimize
from scipy.stats.distributions import chi2
from likelihoodHelpers import log_likelihood_biGauss, log_likelihood_doublePeak
from likelihoodHelpers import likelihood_ratio_doublePeak, likelihood_ratio_biGauss, biGauss, double_peak
from likelihoodHelpers import log_likelihood_expGauss, log_likelihood_expDoublePeak, expGauss, expDoublePeak
import scipy, csv
from tabulate import tabulate

class curveFit(icetray.I3ConditionalModule):
    """
    Fitting single and double bifurcated gaussian
    """

    def __init__(self, context):
        icetray.I3ConditionalModule.__init__(self, context)
        self.AddParameter("omgeo",
                        "geometry map given for the analysis",
                        "omgeo")

        self.AddParameter("InputMCPETree",
                         "Input MCPETree name for analysis",
                         "MCPESeriesMap")

        self.AddParameter("OutputMCPETree",
                         "Output MCPETree name",
                         "Curvefit Parameters")

        self.AddParameter("FrameList",
                          "List of frame numbers to debug",
                          [])

        self.AddParameter("StringList",
                          "List of string numbers to debug",
                          [])

        self.AddParameter("DOMList",
                          "List of DOM numbers to debug",
                          [])

        self.AddOutBox("OutBox")

    def Configure(self):

        self.omgeo = self.GetParameter("omgeo")
        self.input = self.GetParameter("InputMCPETree")
        self.output = self.GetParameter("OutputMCPETree")
        self.frames = self.GetParameter("FrameList")
        self.strings = self.GetParameter("StringList")
        self.doms = self.GetParameter("DOMList")

        self.frame_counter = 0

    def DAQ(self, frame):

        debug_mode = False

        recoPulseMap = frame[self.input]

        biGauss_valuesMap = dataclasses.I3MapKeyVectorDouble()
        doublePeak_valuesMap = dataclasses.I3MapKeyVectorDouble()

        for omkey in recoPulseMap.keys():

            # Check if I want to debug this frame
            if self.frame_counter in self.frames and omkey[0] in self.strings and omkey[1] in self.doms:
                print('Frame number - '+ str(self.frame_counter), 'String number - ' + str(omkey[0]), 'DOM number - '+ str(omkey[1]))
                debug_mode = True
            else:
                debug_mode = False

            recoPulseList = recoPulseMap[omkey]
            recoPulse_timeList = np.array([recoPulse.time for recoPulse in recoPulseList])
            recoPulse_chargeList = np.array([recoPulse.charge for recoPulse in recoPulseList])

            '''
            Removing DOMs with hits less than 200 Hits
            '''
            if sum(recoPulse_chargeList) < 200:
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

            if len(max_hitTimes) < 10:
                continue

            #Shifting mean to zero
            max_hitTimes_mean = sum(max_hitTimes*max_charge)/sum(max_charge)
            timestamps = max_hitTimes - max_hitTimes_mean
            final_mean = sum(timestamps*max_charge)/sum(max_charge)

            '''
            Histogramming the data from simulation
            '''
            #print('Now Histogramming')
            bins = np.arange(min(timestamps), max(timestamps), 3)
            num, bin_edges = np.histogram(timestamps, bins=bins, weights=max_charge)
            bin_centers = (bin_edges[:-1]+bin_edges[1:])/2

            num_ampRatio = num/max(num)

            #removing bins which are <1/5 the max(num), removing the tails this way.
            num_select = num[num_ampRatio > 0.2]
            bin_centers_select = bin_centers[num_ampRatio > 0.2]

            '''
            Including continuity in the bins
            '''

            #considering two extra bins on both sides
            bin_center_bool = (bin_centers >= min(bin_centers_select) - 6)&(bin_centers <= max(bin_centers_select) + 6)
            entries_in_bins = num[bin_center_bool]
            bin_centers = bin_centers[bin_center_bool]

            '''
            Removing DOMs which don't have enough hits
            '''

            if len(entries_in_bins) <= 9:
                continue


            if max(bin_centers) <= 0:
                maxBinCenter = max(bin_centers) + abs(max(bin_centers)) + 3
            else:
                maxBinCenter = max(bin_centers)


            time_window = max(bin_centers) - min(bin_centers)


            '''
            Fitting bifurcated Gaussian and double bifurcated gaussian to
            the mcpe hit time distributions for both tau and electron.
            '''

            #Single Peak

            nll = lambda *args: log_likelihood_biGauss(*args)
            initial_biGauss = np.array([final_mean, time_window/2, 2, max(entries_in_bins)])
            bnds_biGauss = [[min(bin_centers), maxBinCenter], 
                            [-time_window, time_window], # Let the width be negative
                            [0.1, 20], # Restrict k to be positive, but only up to 20 
                            [0.1, 1E10]] # Don't restrict the amplitude, it will vary greatly with K
            if debug_mode == True:
                print('Bounds on single peak')
                # Don't repeat code
                print(tabulate(bnds_biGauss,
                               tablefmt=u'fancy_grid'))
                headers = [["llh","pos1", "wid1", "k1", "amp1"]]
                print(tabulate(headers))

            soln_biGauss = minimize(log_likelihood_expGauss, initial_biGauss,
                                        args=(entries_in_bins, bin_centers, debug_mode),
                                        #method='TNC',
                                        bounds = bnds_biGauss)

            #Double Peak

            nll = lambda *args: log_likelihood_doublePeak(*args)
            peak_time_boundary = final_mean-6.
            initial_doublePeak = np.array([peak_time_boundary-1, 
                                           -time_window/2, 
                                           2, 
                                           max(entries_in_bins), 
                                           peak_time_boundary+1,
                                           time_window/2, 
                                           2, 
                                           max(entries_in_bins)])

            # Re-use the single peak bounds, no need to reinvent stuff
            #peak_time_boundary = final_mean-6.
            bnds_doublePeak = [[min(bin_centers), peak_time_boundary],
                               bnds_biGauss[1],
                               bnds_biGauss[2],
                               bnds_biGauss[3],
                               [peak_time_boundary, maxBinCenter], 
                               bnds_biGauss[1],
                               bnds_biGauss[2],
                               bnds_biGauss[3]]

            # JP: start the loop here
            best_fcn = 1E9
            # JP: Define initial values down here
            initial_doublePeak = np.zeros(8)
            # Get random values within the bounds for each parameter

            if debug_mode == True:
                print('Bounds on double peak')
                # Don't repeat existing stuff
                print(tabulate(bnds_doublePeak,
                               tablefmt=u'fancy_grid'))

            soln_doublePeak = minimize(log_likelihood_expDoublePeak, initial_doublePeak,
                                        args=(entries_in_bins, bin_centers, debug_mode),
                                        #method='TNC',
                                        bounds=bnds_doublePeak)

            # Compare to best_fcn, save results if better
            # JP: End loop. Do something similar for the single fit

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
            Update values
            '''

            biGauss_values = np.array([soln_biGauss.fun, soln_biGauss.x[0],
                                            soln_biGauss.x[1], soln_biGauss.x[2], soln_biGauss.x[3],
                                            LR_biGauss])

            doublePeak_values = np.array([soln_doublePeak.fun, soln_doublePeak.x[0],
                                        soln_doublePeak.x[1],soln_doublePeak.x[2],
                                        soln_doublePeak.x[3], soln_doublePeak.x[4],
                                        soln_doublePeak.x[5], soln_doublePeak.x[6],
                                        soln_doublePeak.x[7], LR_doublePeak])

            # Define omkey:vector dictionary
            biGauss_valuesMap.update({omkey: dataclasses.I3VectorDouble(biGauss_values)})
            doublePeak_valuesMap.update({omkey: dataclasses.I3VectorDouble(doublePeak_values)})

            if debug_mode==True:
                print('Print message single peak -', soln_biGauss.message)
                print('Print message double peak -', soln_doublePeak.message)
                print('Log Likelihood Value single peak -', soln_biGauss.fun)
                print('Log Likelihood Value double peak -', soln_doublePeak.fun)
                import matplotlib.pyplot as plt
                '''
                (x, y) values for the fit
                '''
                #x = bin_centers
                x = np.linspace(min(bin_centers)-30, max(bin_centers)+30, 1000)
                #y_biGauss = biGauss(x, vals_single[1], vals_single[2], vals_single[3], vals_single[4])
                #y_doublePeak = double_peak(x, vals[1], vals[2], vals[3], vals[4],
                #                           vals[5], vals[6], vals[7], vals[8])

                y_biGauss = expGauss(x, soln_biGauss.x[0],
                                                soln_biGauss.x[1], soln_biGauss.x[2], soln_biGauss.x[3])
                y_doublePeak = expDoublePeak(x, soln_doublePeak.x[0],
                                            soln_doublePeak.x[1],soln_doublePeak.x[2],
                                            soln_doublePeak.x[3], soln_doublePeak.x[4],
                                            soln_doublePeak.x[5], soln_doublePeak.x[6],
                                            soln_doublePeak.x[7])

                plt.figure(figsize=(10,9))
                _ = plt.hist(timestamps, bins=bins, weights=max_charge, histtype='step', linewidth = 5)
                plt.plot(bin_centers, entries_in_bins, '*', c='k', label = 'Bins for fit', markersize=12, linewidth=6)
                plt.plot(x, y_biGauss, '--', c = 'k', label = 'expGauss', linewidth=3)
                plt.plot(x, y_doublePeak, '--', c = 'r', label = 'double expGauss', linewidth=3)
                plt.axvline(x=soln_doublePeak.x[0], c='orange', label='Postion of First Gaussian')
                plt.axvline(x=soln_doublePeak.x[4], c='g', label='Postion of Second Gaussian')
                plt.axvline(x=peak_time_boundary,c='k')
                #plt.axhline(y=soln_doublePeak.x[3], c='orange', label='Amplitude of First Gaussian')
                #plt.axhline(y=soln_doublePeak.x[7], c='g', label='Amplitude of Second gaussian')
                plt.legend()
                plt.xlabel('Time(ns)', fontsize = 16)
                plt.title(str(omkey), fontsize=14)

        frame[self.output+'_biGauss'] = biGauss_valuesMap
        frame[self.output+ '_doublePeak'] = doublePeak_valuesMap
        # Increase the frame counter
        self.frame_counter += 1

        self.PushFrame(frame)
