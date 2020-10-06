from icecube import icetray, dataio, dataclasses, simclasses, clsim
from icecube.icetray import I3Units, OMKey, I3Frame
from icecube.dataclasses import ModuleKey
from os.path import expandvars
import numpy as np
import scipy.constants as spc
import matplotlib.pylab as plt
import sys
sys.path.insert(0, '/home/users/akatil/P-ONE/git/PONE_NuTau/BiGauss')
from likelihoodHelpers import log_likelihood_biGauss, log_likelihood_doublePeak, likelihood_ratio_doublePeak, likelihood_ratio_biGauss
from scipy.optimize import minimize

class definingDOMs(icetray.I3ConditionalModule):
    """
    Finding DOMs that show clear double peak structures.
    """

    def __init__(self, context):
        icetray.I3ConditionalModule.__init__(self, context)
        self.AddParameter("omgeo",
                        "geometry map given for the analysis",
                        "omgeo")

        self.AddParameter("InputMCPETree",
                         "Input MCPETree name for analysis",
                         "MCPESeriesMap")

        self.AddParameter("NoiseMCPETree",
                         "Noise MCPETree name for analysis",
                         "NoiseSeriesMap")

        self.AddParameter("OutputMCPETree",
                         "Output MCPETree name",
                         "DoublePeakSeriesMap")
        self.AddOutBox("OutBox")

    def Configure(self):

        self.omgeo = self.GetParameter("omgeo")
        self.mcpeSeries = self.GetParameter("InputMCPETree")
        self.noiseSeries = self.GetParameter("NoiseMCPETree")
        self.output = self.GetParameter("OutputMCPETree")

    def DAQ(self, frame):

        mctree = frame["I3MCTree"]
        doublePeakMap = simclasses.I3MCPESeriesMap()
        primary = mctree.primaries
        lepton = dataclasses.I3MCTree.first_child(mctree, primary[0].id)

        if lepton.type == 15 or lepton.type == -15:

            tau_daughters = dataclasses.I3MCTree.get_daughters(mctree, lepton.id)
            tau_pos = lepton.pos
            x_tau_pos = tau_pos.x
            y_tau_pos = tau_pos.y
            z_tau_pos = tau_pos.z

            for td in range(0, len(tau_daughters)):
                if tau_daughters[td].type == 16 or tau_daughters[td].type == -16:
                    tau_daughters_pos = tau_daughters[td].pos
                    x_td_pos = tau_daughters_pos.x
                    y_td_pos = tau_daughters_pos.y
                    z_td_pos = tau_daughters_pos.z

            mcpeMap = frame[self.mcpeSeries]
            noiseMap = frame[self.noiseSeries]

            x_center, y_center, z_center = ([]), ([]), ([])
            x_double_peak, y_double_peak, z_double_peak = ([]), ([]), ([])
            charge = ([])
            double_peak_charge = ([])

            for omkey in noiseMap.keys():
                oKey = self.omgeo.get(omkey)

                domPos = oKey.position
                x_dom = domPos.x
                y_dom = domPos.y
                z_dom = domPos.z

                noise_mcpeList = noiseMap[omkey]
                noise_timeList = np.array([mcpe.time for mcpe in noise_mcpeList])

                if omkey in mcpeMap.keys():
                    mcpeList = mcpeMap[omkey]
                    timeList = np.array([mcpe.time for mcpe in mcpeList])
                    tot_timeList = np.append(timeList, noise_timeList)
                else:
                    tot_timeList = noise_timeList

                charge = np.append(charge, len(tot_timeList))

                x_center = np.append(x_center, x_dom)
                y_center = np.append(y_center, y_dom)
                z_center = np.append(z_center, z_dom)

                firstVertex = np.sqrt((x_dom - x_tau_pos)**2 + (y_dom - y_tau_pos)**2 + (z_dom - z_tau_pos)**2)
                secondVertex = np.sqrt((x_dom - x_td_pos)**2 + (y_dom - y_td_pos)**2 + (z_dom - z_td_pos)**2)
                refractiveIndex = 1.333
                speed_of_light_water = (spc.c)/refractiveIndex #[Units: m/seconds]
                speed_of_light_ns = speed_of_light_water

                tDiff_ns = ((firstVertex - secondVertex)/speed_of_light_water) * 1e9 #[Units: nanoseconds]

                doublePeakList = simclasses.I3MCPESeries()

                if abs(tDiff_ns) > 20 and abs(tDiff_ns) < 400 and firstVertex < 250 and secondVertex < 250 and len(tot_timeList) > 100:

                    double_peak_charge = np.append(double_peak_charge, len(tot_timeList))

                    for time in tot_timeList:
                        mcpe_new = simclasses.I3MCPE()
                        mcpe_new.npe = 1
                        mcpe_new.time = time
                        doublePeakList.append(mcpe_new)

                    doublePeakMap[omkey] = doublePeakList

                    plt.figure(figsize=(10,9))
                    mean_timestamps = tot_timeList.mean()
                    #bins = np.arange(mean_timestamps - 100, mean_timestamps + 100, 1)
                    num, bins, _ = plt.hist(timestamps, bins = bins, histtype='step')
                    plt.title('OMKey:' + str(omkey) + ' TimeDiff:' + str(timeDifference_doublePeak))

                    x_double_peak = np.append(x_double_peak, x_dom)
                    y_double_peak = np.append(y_double_peak, y_dom)
                    z_double_peak = np.append(z_double_peak, z_dom)

            if len(doublePeakMap) > 0:
                xCoG = sum(x_center*charge)/sum(charge)
                yCoG = sum(y_center*charge)/sum(charge)
                zCoG = sum(z_center*charge)/sum(charge)

                xDoublePeak_CoG = sum(x_double_peak*double_peak_charge)/sum(double_peak_charge)
                yDoublePeak_CoG = sum(y_double_peak*double_peak_charge)/sum(double_peak_charge)
                zDoublePeak_CoG = sum(z_double_peak*double_peak_charge)/sum(double_peak_charge)

                CoG = dataclasses.I3Position(xCoG, yCoG, zCoG)
                CoG_doublePeak = dataclasses.I3Position(xDoublePeak_CoG, yDoublePeak_CoG, zDoublePeak_CoG)

                frame['CoG'] = CoG
                frmae['CoG_doublePeak'] = CoG_doublePeak

                #distCoG = np.sqrt((x_double_peak - xCoG)**2 + (y_double_peak - yCoG)**2 + (z_double_peak - zCoG)**2)
                #print('Center of Mass - ', xCoG, yCoG, zCoG)
                #print('Distance to CoG -', distCoG)
                frame[self.output] = doublePeakMap
                self.PushFrame(frame)
