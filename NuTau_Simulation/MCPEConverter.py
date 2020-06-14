from icecube import icetray, dataio, dataclasses, simclasses, clsim
from icecube.icetray import I3Units, OMKey, I3Frame
from icecube.dataclasses import ModuleKey
from os.path import expandvars
import numpy as np
import argparse

class makeHits(icetray.I3ConditionalModule):
    """
    Generate hits from I3Photons for P-ONE
    """

    def __init__(self, context):
        icetray.I3ConditionalModule.__init__(self, context)
        self.AddParameter("WavelengthAcceptance",
                          "Path to the wavelength acceptance of the DOM",
                          "/home/users/akatil/P-ONE/git/PONE_NuTau/DOM/MDOM/DOMEfficiency.dat")
        self.AddParameter("AngularAcceptance",
                          "Input the flat angular acceptance for the mDOM",
                          0.811397114255)
        self.AddParameter("GCDFile",
                          "GCD to be simulated",
                          dataio.I3File('/home/users/akatil/P-ONE/GCD_files/PONE_Phase1.i3.gz'))
        self.AddParameter("PropagatedPhotons",
                          "I3MCTree containing the list of propagated I3Photons",
                          "I3Photons")
        self.AddParameter("OutputMCTreeName",
                          "Name of the out put I3MCTree series name",
                          "MCPESeriesMap")
        self.AddParameter("DOMOversizeFactor",
                          "Specifiy the \"oversize factor\" (i.e. DOM radius scaling factor) you used during the CLSim run.\n"
	                      "The photon arrival times will be corrected. In practice this means your large spherical DOMs will\n"
	                      "become ellipsoids.",
                          1.)
        self.AddParameter("DefaultRelativeDOMEfficiency",
	                      "Default relative efficiency.",
	                      1.0);
        self.AddOutBox("OutBox")

    def Configure(self):

        self.wavelengthAcceptance = self.GetParameter("WavelengthAcceptance")
        self.angularAcceptance = self.GetParameter("AngularAcceptance")
        self.gcdFile = self.GetParameter("GCDFile")
        self.propagatedPhotons = self.GetParameter("PropagatedPhotons")
        self.outputMCTreeName = self.GetParameter("OutputMCTreeName")
        self.domOversizeFactor = self.GetParameter("DOMOversizeFactor")
        self.relativeEfficieny = self.GetParameter("DefaultRelativeDOMEfficiency")


    def DAQ(self, frame):

        gcd_file = dataio.I3File(self.gcdFile)
        cframe = gcd_file.pop_frame()
        geometry = cframe["I3Geometry"]
        geoMap = geometry.omgeo
        #calibration = cframe["I3Calibration"]
        #calMap = calibration.dom_cal

        DOMRadius = 0.16510*I3Units.m
        domAcceptance = clsim.GetIceCubeDOMAcceptance(
            	                                domRadius = DOMRadius*self.domOversizeFactor,
            	                                efficiency=self.relativeEfficieny)

        def getSurvivalProbability(photon, omkey):

            domGeo = geoMap[omkey]
            domDirection = domGeo.direction
            photonDirection = photon.dir
            dotProduct = photonDirection.x*domDirection.x + photonDirection.y*domDirection.y + photonDirection.z*domDirection.z
            # photon coming in, direction coming out. Sign on dot product should be flipped
            # directions are unit vectors already so cos_theta = -dotProduct (due to sign flip)

            probDOMAcc = domAcceptance.GetValue(photon.wavelength)

            probAngAcc = self.angularAcceptance

            return probAngAcc*probDOMAcc*photon.weight*self.relativeEfficieny


        def survived(photon,omkey):
            probability = getSurvivalProbability(photon, omkey)
            randomNumber = np.random.uniform()

            if(probability > 1):

                raise ValueError("Probability = " + str(probability) + " > 1. Most likely photon weights are too high")

            if(probability > randomNumber):
                return True

            return False

        def generateMCPEList(photons, modkey):
            omkey = OMKey(modkey.string, modkey.om, 0)
            mcpeList = simclasses.I3MCPESeries()
            photonList = []
            for photon in photons:
                if survived(photon,omkey):
                    mcpe = simclasses.I3MCPE()
                    #mcpe.id = dataclasses.I3ParticleID(photon.particleMajorID, photon.particleMinorID)
                    mcpe.npe = 1
                    mcpe.time = photon.time #TODO: change to corrected time
                    mcpeList.append(mcpe)
                    photonList.append(photon)

            return mcpeList, photonList

        photonDOMMap = frame[self.propagatedPhotons]
        mcpeMap = simclasses.I3MCPESeriesMap()
        succPhotonMap = simclasses.I3CompressedPhotonSeriesMap()
        for modkey in photonDOMMap.keys():
            mcpeList, photonList = generateMCPEList(photonDOMMap[modkey], modkey)
            omkey = OMKey(modkey.string, modkey.om, 0)
            if len(mcpeList) > 0:
                mcpeMap[omkey] = mcpeList
                succPhotonMap[modkey] = photonList

        if len(mcpeMap) > 0:
            frame[self.outputMCTreeName] = mcpeMap

            frame["SuccPhotonMap"] = succPhotonMap
            #frame.Delete("I3Photons")

            self.PushFrame(frame)
