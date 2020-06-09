from icecube import dataclasses, dataio, simclasses
from icecube.icetray import I3Units, OMKey, I3Frame
from icecube.dataclasses import ModuleKey
import numpy as np
import random, re, os
from optparse import OptionParser
from os.path import expandvars

usage = "usage: %prog [options]"
parser = OptionParser(usage)
parser.add_option("-o", "--outfile",default="./test_output.i3",
                  dest="OUTFILE", help="Write output to OUTFILE (.i3{.gz} format)")
parser.add_option("-i", "--infile",default="./test_input.i3",
                  dest="INFILE", help="Read input from INFILE (.i3{.gz} format)")
parser.add_option("-g", "--gcdfile", default=os.getenv('GCDfile'),
		          dest="GCDFILE", help="Read in GCD file")

(options,args) = parser.parse_args()


######
"""Noise"""
######

def selectStartTime(timestamps):
    randStartTime = random.uniform(0, max(timestamps))
    noiseWindow = timestamps[(timestamps >= randStartTime) & (timestamps <= randStartTime+72e2)]
    return randStartTime, noiseWindow

def loadFile(File):
    print 'loading file'
    data = np.loadtxt(File)

    channel = data[:, 0]
    upperPMT = data[(channel==1)]
    lowerPMT = data[(channel==5)]
    risingEdge = upperPMT[upperPMT[:, 1] == 0]
    #risingEdgeL = lowerPMT[lowerPMT[:, 1] == 0]
    timestamps = risingEdge[:, 2]
    print 'timestamps generated'
    return timestamps

def addHits(timestamps):
    mDOM_noise = ([])
    startTimeVals = ([])
    for segment in range(0, 24):
        startTime, noiseWindow = selectStartTime(timestamps)
        movingVals = (noiseWindow - startTime)-(72e2/2)
        startTimeVals = np.append(startTimeVals, startTime)
        mDOM_noise = np.append(mDOM_noise, movingVals)
    #print "mDOM_noise", mDOM_noise
    return mDOM_noise


#######
'''Geometry'''
######

#gcdPath = '/home/users/akatil/P-ONE/GCD_files/PONE_Phase1.i3.gz'
geofile = dataio.I3File(options.GCDFILE)
cframe = geofile.pop_frame(I3Frame.Calibration)
geometry = cframe["I3Geometry"]
geoMap = geometry.omgeo
#######

#######
'''I3 File'''
#######

outPath = '/data/p-one/akatil/test/' + 'genHitsNoiseTrial.i3.gz'
outfile = dataio.I3File(options.OUTFILE, 'w')

def choose_file():
    filename = random.choice(os.listdir("/data/p-one/gaertner/2004_embla/unpacked/"))
    FLASH_re = re.compile('FLASH')
    if len(FLASH_re.findall(filename)) > 0:
        choose_file()
        #print 'true', choose_file()
        return choose_file()
    else:
        return filename

#file = dataio.I3File('/data/p-one/akatil/step_3_medium_water/Custom/step_3_1000_medium_water_custom.i3.gz')
file = dataio.I3File(options.INFILE)

filename = choose_file()

timestamps = loadFile('/data/p-one/gaertner/2004_embla/unpacked/' + filename)

frames = []
while(file.more()):
    frames.append(file.pop_daq())

def generateNoiseMCPEList(modkey, medVal):
    omkey = OMKey(modkey.string, modkey.om, 0)
    noiseHits = addHits(timestamps)
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

for frame in frames:
    #frame = frames[f]
    #print frame
    mcpeMap = frame['MCPESeriesMap']
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
        print 'generating noise in DOM'
        noiseList = generateNoiseMCPEList(modkey, medVal)
        omkey = OMKey(modkey.string, modkey.om, 0)
        print omkey
        mcpeNoiseMap[omkey] = noiseList

    frame["NoiseSeriesMap"] = mcpeNoiseMap

    # only add frame to file if a hit was generated
    #if passFrame(frame, mcpeMap.keys(), int(args.hitThresh), int(args.domThresh), ):
        #frame.Delete("I3Photons")
    outfile.push(frame)

outfile.close()
