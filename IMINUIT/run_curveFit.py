from icecube import icetray, dataio, dataclasses, simclasses, clsim
from icecube.icetray import I3Units, OMKey, I3Frame
from I3Tray import *
from icecube.dataclasses import ModuleKey
from optparse import OptionParser
from os.path import expandvars
import numpy as np
from scipy import stats
from iminuit import minimize
from likelihoodHelpers import log_likelihood_biGauss, log_likelihood_doublePeak
from likelihoodHelpers import likelihood_ratio_doublePeak, likelihood_ratio_biGauss, biGauss, double_peak
import scipy, csv
import curveFit

usage = "usage: %prog [options]"
parser = OptionParser(usage)
parser.add_option("-o", "--outfile",default="./test_output.i3",
                  dest="OUTFILE", help="Write output to OUTFILE (.i3{.gz} format)")
parser.add_option("-i", "--infile",default="./test_input.i3",
                  dest="INFILE", help="Read input from INFILE (.i3{.gz} format)")

(options,args) = parser.parse_args()


'''
Loading Geometry
'''

gcd_file = '/home/users/akatil/P-ONE/GCD_files/PONE_Phase1.i3.gz'
gcd = dataio.I3File(gcd_file)
cframe = gcd.pop_frame()
geometry = cframe["I3Geometry"]
omgeo = geometry.omgeo
print('loaded geometry')

infile = '/data/p-one/akatil/step_5_medium_water/NuTau_NuE_20Events/step_5_713_medium_water_custom_mDOM_recoPulse.i3.gz'
outfile = '/data/p-one/akatil/test/step_6_713_parameters.i3.gz'

tray = I3Tray()

tray.AddModule('I3Reader', 'reader',
            FilenameList = [gcd_file, options.INFILE]
            )

tray.AddModule(curveFit.curveFit, "Double Peak Selector",
               omgeo = omgeo,
               InputMCPETree = "I3RecoPulses",
               OutputMCPETree = "Parameters")

tray.AddModule("I3Writer","writer",
               #SkipKeys=SkipKeys,
               Filename = options.OUTFILE,
               Streams = [icetray.I3Frame.DAQ, icetray.I3Frame.Physics],
              )

tray.AddModule("TrashCan","adios")
tray.Execute()
tray.Finish()
