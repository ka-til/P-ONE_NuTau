from icecube import icetray, dataclasses, dataio, simclasses
from icecube.icetray import I3Units, OMKey, I3Frame
from icecube.dataclasses import ModuleKey
import numpy as np
import random, re, os
from optparse import OptionParser
from os.path import expandvars

class addNoise(icetray.I3ConditionalModule):
    """
    Generate hits from I3Photons for P-ONE
    """

    def __init__(self, context):
        icetray.I3ConditionalModule.__init__(self, context)
        self.AddParameter("GCDFile",
                          "GCD to be simulated",
                          '/home/users/akatil/P-ONE/GCD_files/PONE_Phase1.i3.gz')
        self.AddParameter("PhysicsMCPETreeName",
                          "Name of the Physics I3MCTree name",
                          "MCPESeriesMap")
        self.AddParameter("NoiseMCPEMCTreeName",
                          "Name of the noise I3MCTree name",
                          "NoiseSeriesMap")
        self.AddParameter("STRAWTimestamps",
                          "Noise Hits for mDOMs to be injected")
        self.AddOutBox("OutBox")

    def Configure(self):

        self.gcdFile = self.GetParameter("GCDFile")
        self.mcpeSeriesName = self.GetParameter("PhysicsMCPETreeName")
        self.noiseSeriesName = self.GetParameter("NoiseMCPEMCTreeName")
        #self.mDOM_Noise = self.GetParameter("mDOM_noise")
        self.timestamps = self.GetParameter("STRAWTimestamps")


    def DAQ(self, frame):

        def selectStartTime(timestamps):
            randStartTime = random.uniform(0, max(timestamps))
            noiseWindow = timestamps[(timestamps >= randStartTime) & (timestamps <= randStartTime+72e2)]
            return randStartTime, noiseWindow

        def addHits(timestamps):
            print "adding Hits"
            mDOM_noise = ([])
            startTimeVals = ([])
            for segment in range(0, 24):
                startTime, noiseWindow = selectStartTime(timestamps)
                movingVals = (noiseWindow - startTime)-(72e2/2)
                startTimeVals = np.append(startTimeVals, startTime)
                mDOM_noise = np.append(mDOM_noise, movingVals)
            #print "mDOM_noise", mDOM_noise
            return mDOM_noise

        gcd_file = dataio.I3File(self.gcdFile)
        cframe = gcd_file.pop_frame()
        geometry = cframe["I3Geometry"]
        geoMap = geometry.omgeo

        def generateNoiseMCPEList(modkey, medVal):
            omkey = OMKey(modkey.string, modkey.om, 0)
            noiseHits = addHits(self.timestamps)
            noiseList = simclasses.I3MCPESeries()

            for noiseHit in noiseHits:
                mcpeNoise = simclasses.I3MCPE()
                #mcpe.id = dataclasses.I3ParticleID(photon.particleMajorID, photon.particleMinorID)
                mcpeNoise.npe = 1
                mcpeNoise.time = noiseHit + medVal
                #print "mcpeNoise.time", mcpeNoise.time
                noiseList.append(mcpeNoise)

            #print 'NoiseList', noiseList
            return noiseList

        mcpeMap = frame[self.mcpeSeriesName]
        mcpeOMKeys = mcpeMap.keys()
        mcpeTime = ([])
        for mcpeList in mcpeMap.values():
            timeList = [mcpe.time for mcpe in mcpeList]
            mcpeTime = np.append(mcpeTime, timeList)
        medVal = np.mean(mcpeTime)
        #print 'mcpeTime', mcpeTime
        #print 'mean', medVal
        print 'median calculated'

        mcpeNoiseMap = simclasses.I3MCPESeriesMap()
        for modkey in geoMap.keys():
            #print 'generating noise in DOM'
            noiseList = generateNoiseMCPEList(modkey, medVal)
            omkey = OMKey(modkey.string, modkey.om, 0)
            print omkey
            mcpeNoiseMap[omkey] = noiseList

        frame[self.noiseSeriesName] = mcpeNoiseMap

        self.PushFrame(frame)
