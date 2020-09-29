from icecube import icetray, dataclasses, dataio, simclasses
from icecube.icetray import I3Units, OMKey, I3Frame
from icecube.dataclasses import ModuleKey
import numpy as np
import random, re, os
from optparse import OptionParser
from os.path import expandvars

class combiningHits(icetray.I3ConditionalModule):
    """
    Merging physics and noise hits
    """

    def __init__(self, context):
        icetray.I3ConditionalModule.__init__(self, context)

        self.AddParameter("PhysicsMCPETreeName",
                          "Name of the Physics I3MCTree name",
                          "MCPESeriesMap")
        self.AddParameter("NoiseMCPEMCTreeName",
                          "Name of the noise I3MCTree name",
                          "NoiseSeriesMap")
        self.AddParameter("MergedMCPETreeName",
                          "Name of the Merged MCPE tree name",
                          "MergedSeriesMap")
        self.AddOutBox("OutBox")

    def Configure(self):

        self.mcpeSeriesName = self.GetParameter("PhysicsMCPETreeName")
        self.noiseSeriesName = self.GetParameter("NoiseMCPEMCTreeName")
        self.mergedSeriesName = self.GetParameter("MergedMCPETreeName")

    def DAQ(self, frame):

        mcpeMap = frame[self.mcpeSeriesMap]
        mcpeOMKeys = mcpeMap.keys()

        noiseMap = frame[self.noiseSeriesName]
        noiseOMKeys = noiseMap.keys()

        mergedHitsMap = simclasses.I3MCPESeriesMap()

        for omkey in noiseOMKeys:

            newMCPEList = simclasses.I3MCPESeries()
            noise_mcpeList = noiseMap[omkey]

            for mcpe in noise_mcpeList:
                newMCPEList.append(mcpe)

            if omkey in mcpeOMKeys:
                mcpeList = mcpeMap[omkey]
                for mcpe in mcpeList:
                    newMCPEList.append(mcpe)

            mergedHitsMap[omkey] = newMCPEList

        frame[self.mergedSeriesName] = mergedHitsMap
        self.PushFrame(frame)

        #for omkey in noiseOMKeys:
            #noise_mcpeList = noiseMap[omkey]
            #noiseList = mcpe for mcpe in noise_mcpeList

            #if omkey in mcpeOMKeys:
            #    mcpeList = mcpeMap[omkey]
            #    phsicsList = mcpe for mcpe in mcpeList
            #    totList = noiseList.append(PhysicsList)
            #else:
            #    totList = noiseList

            #mergedHitsMap[omkey] = totList

        #frame[self.mergedSeriesName] = mergedHitsMap

        #self.PushFrame(frame)
