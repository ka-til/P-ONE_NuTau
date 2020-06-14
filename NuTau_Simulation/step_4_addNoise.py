from I3Tray import *
from icecube import icetray, dataclasses, dataio, simclasses
from icecube.icetray import I3Units, OMKey, I3Frame
from icecube.dataclasses import ModuleKey
import numpy as np
import random, re, os
from optparse import OptionParser
from os.path import expandvars
import NoiseGenerator
import time

start_time = time.time()
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

#######
'''I3Module'''
#######

def choose_file():
    filename = random.choice(os.listdir("/data/p-one/gaertner/2004_embla/unpacked/"))
    FLASH_re = re.compile('FLASH')
    if len(FLASH_re.findall(filename)) > 0:
        return choose_file()
    else:
        return filename

#file = dataio.I3File('/data/p-one/akatil/step_3_medium_water/Custom/step_3_1000_medium_water_custom.i3.gz')
filename = choose_file()
timestamps = loadFile('/data/p-one/gaertner/2004_embla/unpacked/' + filename)

tray = I3Tray()

tray.AddModule('I3Reader', 'reader',
            FilenameList = [options.GCDFILE, options.INFILE]
            )

tray.AddModule(NoiseGenerator.addNoise, "InjectSTRAWNoise",
               GCDFile = options.GCDFILE,
               PhysicsMCPETreeName = "MCPESeriesMap",
               NoiseMCPEMCTreeName = "NoiseSeriesMap",
               STRAWTimestamps = timestamps)

tray.AddModule("I3Writer","writer",
               #SkipKeys=SkipKeys,
               Filename = options.OUTFILE,
               Streams = [icetray.I3Frame.DAQ, icetray.I3Frame.Physics],
              )

tray.AddModule("TrashCan","adios")
tray.Execute()
tray.Finish()

print "--- %s seconds ---" % (time.time() - start_time)
