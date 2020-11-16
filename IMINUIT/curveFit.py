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
            print('Now Histogramming')
            bins = np.arange(min(timestamps), max(timestamps), 3)
            num, bin_edges = np.histogram(timestamps, bins=bins, weights=max_charge)
            bin_centers = (bin_edges[:-1]+bin_edges[1:])/2

            num_ampRatio = num/max(num)

            #removing bins which are <1/5 the max(num), removing the tails this way.
            num_select = num[num_ampRatio > 0.1]
            bin_centers_select = bin_centers[num_ampRatio > 0.1]

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
            initial_biGauss = np.array([final_mean, 50, 5, max(entries_in_bins)])
            bnds_biGauss = ((min(bin_centers), maxBinCenter), (0, time_window), (1, 20), (1, 2*max(entries_in_bins)))
            if debug_mode == True:
                print('bounds on single peak')
                print(tabulate([(min(bin_centers), maxBinCenter), (0, time_window), (1, 20), (1, 2*max(entries_in_bins))], tablefmt=u'fancy_grid'))
            soln_biGauss = minimize(log_likelihood_expGauss, initial_biGauss,
                                        args=(entries_in_bins, bin_centers, debug_mode),
                                        #method='TNC',
                                        bounds = bnds_biGauss)

            #Double Peak

            nll = lambda *args: log_likelihood_doublePeak(*args)
            initial_doublePeak = np.array([min(bin_centers)+10, 20, 1, max(entries_in_bins), final_mean, 20, 1, max(entries_in_bins)])
            bnds_doublePeak = ((min(bin_centers), final_mean-6), (0, time_window), (1, 20), (1, 2*max(entries_in_bins)),
                                    (final_mean-6, maxBinCenter), (0, time_window), (1, 20), (1, 2*max(entries_in_bins)))
            if debug_mode == True:
                print('bounds on double peak')
                print(tabulate([(min(bin_centers), final_mean-6), (0, time_window), (1, 20), (1, 2*max(entries_in_bins)),
                                        (final_mean-6, maxBinCenter), (0, time_window), (1, 20), (1, 2*max(entries_in_bins))], tablefmt=u'fancy_grid'))

            soln_doublePeak = minimize(log_likelihood_expDoublePeak, initial_doublePeak,
                                        args=(entries_in_bins, bin_centers, debug_mode),
                                        #method='TNC',
                                        bounds=bnds_doublePeak)

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

        frame[self.output+'_biGauss'] = biGauss_valuesMap
        frame[self.output+ '_doublePeak'] = doublePeak_valuesMap
        # Increase the frame counter
        self.frame_counter += 1

        self.PushFrame(frame)
