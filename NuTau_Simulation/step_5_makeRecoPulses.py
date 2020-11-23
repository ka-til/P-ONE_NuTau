from I3Tray import *
from icecube import icetray, dataclasses, dataio, simclasses
from icecube.icetray import I3Units, OMKey, I3Frame
from icecube.dataclasses import ModuleKey
import numpy as np
import random, re, os
from optparse import OptionParser
from os.path import expandvars
import RecoPulseGenerator, TimeShifted, NoisePhysicsHitsMerger
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

#######
'''I3Modules'''
#######


tray = I3Tray()

tray.AddModule('I3Reader', 'reader',
            FilenameList = [options.GCDFILE, options.INFILE]
            )

tray.AddModule(NoisePhysicsHitsMerger.combiningHits, "Merge Physics and Noise hits",
               PhysicsMCPETreeName = "MCPESeriesMap",
               NoiseMCPEMCTreeName = "NoiseSeriesMap",
               MergedMCPETreeName = "MergedSeriesMap")

tray.AddModule(TimeShifted.timeShift, "TimeShifter",
               MergedMCPETreeName = "MergedSeriesMap",
               TimeShiftedMCPE = "TimeShiftedMCPEMap")

tray.AddModule(RecoPulseGenerator.merge_recoPulse, "RecoPulseGenerator",
               TimeShiftedMCPE = "TimeShiftedMCPEMap",
               RecoPulseTreeName = "I3RecoPulses",
               MergeTime = 3)

SkipKeys = ["MergedSeriesMap"]

tray.AddModule("I3Writer","writer",
               SkipKeys=SkipKeys,
               Filename = options.OUTFILE,
               Streams = [icetray.I3Frame.DAQ, icetray.I3Frame.Physics],
              )

tray.AddModule("TrashCan","adios")
tray.Execute()
tray.Finish()

print "--- %s seconds ---" % (time.time() - start_time)
