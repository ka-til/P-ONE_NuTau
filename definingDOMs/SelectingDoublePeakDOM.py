from icecube import icetray, dataio, dataclasses, simclasses, clsim
from icecube.phys_services import I3Calculator
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

        doms = ([])

        mctree = frame["I3MCTree"]
        doublePeakMap = simclasses.I3MCPESeriesMap()
        primary = mctree.primaries
        lepton = dataclasses.I3MCTree.first_child(mctree, primary[0].id)

        if lepton.type == 15 or lepton.type == -15:
            #print(lepton.type)

            tau_daughters = dataclasses.I3MCTree.get_daughters(mctree, lepton.id)
            tau_pos = lepton.pos
            x_tau_pos = tau_pos.x
            y_tau_pos = tau_pos.y
            z_tau_pos = tau_pos.z

            lepton.shape = dataclasses.I3Particle.InfiniteTrack

            #print(lepton.shape)

            if abs(x_tau_pos) < 150 and abs(y_tau_pos) < 150 and lepton.length > 20:
                print('condition met - 1')

                for td in range(0, len(tau_daughters)):
                    #print('Daughters', tau_daughters[td].type)
                    if tau_daughters[td].type == 16 or tau_daughters[td].type == -16:
                        tau_daughters_pos = tau_daughters[td].pos
                        x_td_pos = tau_daughters_pos.x
                        y_td_pos = tau_daughters_pos.y
                        z_td_pos = tau_daughters_pos.z


                #rho =
                physics_vector = [x_td_pos - x_tau_pos, y_td_pos - y_tau_pos, z_td_pos - z_tau_pos]

                phys_magnitude = np.sqrt(physics_vector[0]**2 + physics_vector[1]**2 + physics_vector[2]**2)

                if abs(x_td_pos) < 150 and abs(y_td_pos) < 150 and abs(z_td_pos) < 500:
                    print('condition met - 2')

                    mcpeMap = frame[self.mcpeSeries]
                    noiseMap = frame[self.noiseSeries]

                    x_center, y_center, z_center = ([]), ([]), ([])
                    x_physics, y_physics, z_physics = ([]), ([]), ([])
                    x_double_peak, y_double_peak, z_double_peak = ([]), ([]), ([])
                    charge = ([])
                    charge_physics = ([])
                    double_peak_charge = ([])

                    dom_number = 0

                    for omkey in noiseMap.keys():
                        oKey = self.omgeo.get(omkey)

                        domPos = oKey.position
                        x_dom = domPos.x
                        y_dom = domPos.y
                        z_dom = domPos.z

                        tau_daughter_vector = [x_td_pos - x_dom, y_td_pos - y_dom, z_td_pos - z_dom]
                        tau_vector = [x_tau_pos - x_dom, y_tau_pos - y_dom, z_tau_pos - z_dom]

                        td_magnitude = np.sqrt(tau_daughter_vector[0]**2 + tau_daughter_vector[1]**2 + tau_daughter_vector[2]**2)
                        tau_magnitude = np.sqrt(tau_vector[0]**2 + tau_vector[1]**2 + tau_vector[2]**2)


                        dotProduct = tau_daughter_vector[0]*tau_vector[0] + tau_daughter_vector[1]*tau_vector[1] + tau_daughter_vector[2]*tau_vector[2]

                        theta = np.arccos(dotProduct/(td_magnitude*tau_magnitude))
                        theta_deg = np.degrees(theta)

                        dist1 = np.sqrt((x_tau_pos - x_dom)**2 + (y_tau_pos - y_dom)**2 + (z_tau_pos - z_dom)**2)
                        dist2 = np.sqrt((x_td_pos - x_dom)**2 + (y_td_pos - y_dom)**2 + (z_td_pos - z_dom)**2)


                        perp_dist = I3Calculator.distance_along_track(lepton, domPos)
                        two_vertex_distances = np.sqrt((x_td_pos - x_tau_pos)**2 + (y_td_pos - y_tau_pos)**2 + (z_td_pos - z_tau_pos)**2)



                        #dotProduct = x_dom*physics_vector[0] + y_dom*physics_vector[1] + z_dom*physics_vector[2]
                        #dom_pos_magnitude = np.sqrt(x_dom**2+y_dom**2+z_dom**2)

                        #theta = np.arccos(dotProduct/(phys_magnitude*dom_pos_magnitude))
                        #theta_deg = np.degrees(theta)

                        noise_mcpeList = noiseMap[omkey]
                        noise_timeList = np.array([mcpe.time for mcpe in noise_mcpeList])

                        if omkey in mcpeMap.keys():
                            x_physics = np.append(x_physics, x_dom)
                            y_physics = np.append(y_physics, y_dom)
                            z_physics = np.append(z_physics, z_dom)

                            mcpeList = mcpeMap[omkey]
                            timeList = np.array([mcpe.time for mcpe in mcpeList])
                            tot_timeList = np.append(timeList, noise_timeList)
                            charge_physics = np.append(charge_physics, len(timeList))
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

                        if len(tot_timeList) > 100: # and firstVertex < 250 and secondVertex < 250:  abs(tDiff_ns) > 0 and abs(tDiff_ns) < 200 and

                            double_peak_charge = np.append(double_peak_charge, len(tot_timeList))

                            doms = np.append(doms, omkey)


                            for time in tot_timeList:
                                mcpe_new = simclasses.I3MCPE()
                                mcpe_new.npe = 1
                                mcpe_new.time = time
                                doublePeakList.append(mcpe_new)

                            doublePeakMap[omkey] = doublePeakList

                            plt.figure(figsize=(10,9))
                            mean_timestamps = tot_timeList.mean()
                            bins = np.arange(mean_timestamps - 100, mean_timestamps + 100, 1)
                            num, bins, _ = plt.hist(tot_timeList, bins = bins, histtype='step')
                            #plt.title('OMKey:' + str(omkey) + ' DOM_NUM: ' + str(dom_number) + ' Degrees: ' + str(theta_deg), fontsize = 18)
                            plt.title('OMKey:' + str(omkey) + ' DOM_NUM: ' + str(dom_number) + ' dist1: ' + str(round(perp_dist)) + ' dist2: ' + str(round(two_vertex_distances)), fontsize = 18)

                            x_double_peak = np.append(x_double_peak, x_dom)
                            y_double_peak = np.append(y_double_peak, y_dom)
                            z_double_peak = np.append(z_double_peak, z_dom)

                            dom_number += 1

                    print('------------------------NUMBER OF DOMS--------------------------', len(doms))
                    print('OMKEYS - ', doms)
                    if len(doublePeakMap) > 0:
                        xCoG = sum(x_center*charge)/sum(charge)
                        yCoG = sum(y_center*charge)/sum(charge)
                        zCoG = sum(z_center*charge)/sum(charge)

                        xCoG_physics = sum(x_physics*charge_physics)/sum(charge_physics)
                        yCoG_physics = sum(y_physics*charge_physics)/sum(charge_physics)
                        zCoG_physics = sum(z_physics*charge_physics)/sum(charge_physics)

                        xDoublePeak_CoG = sum(x_double_peak*double_peak_charge)/sum(double_peak_charge)
                        yDoublePeak_CoG = sum(y_double_peak*double_peak_charge)/sum(double_peak_charge)
                        zDoublePeak_CoG = sum(z_double_peak*double_peak_charge)/sum(double_peak_charge)

                        distCoG = np.sqrt((x_double_peak - xCoG)**2 + (y_double_peak - yCoG)**2 + (z_double_peak - zCoG)**2)
                        min_dist = min(distCoG)

                        CoG = dataclasses.I3Position(xCoG, yCoG, zCoG)
                        CoG_physics = dataclasses.I3Position(xCoG_physics, yCoG_physics, zCoG_physics)
                        CoG_doublePeak = dataclasses.I3Position(xDoublePeak_CoG, yDoublePeak_CoG, zDoublePeak_CoG)

                        frame['CoG'] = CoG
                        frame['CoG_physics'] = CoG_physics
                        frame['CoG_doublePeak'] = CoG_doublePeak

                        distCoG = np.sqrt((x_double_peak - xCoG)**2 + (y_double_peak - yCoG)**2 + (z_double_peak - zCoG)**2)
                        #print('Center of Mass - ', xCoG, yCoG, zCoG)
                        #print('Distance to CoG -', distCoG)
                        frame[self.output] = doublePeakMap
                        self.PushFrame(frame)
