from icecube import icetray, dataio, dataclasses, simclasses, clsim
from icecube.icetray import I3Units, OMKey, I3Frame
from I3Tray import *
from icecube.dataclasses import ModuleKey
import numpy as np
from scipy import stats
from scipy.optimize import minimize
from scipy.stats.distributions import chi
import scipy.constants as spc
import SelectingDoublePeakDOM

'''

Loading the geometry

'''

gcd_file = '/home/users/akatil/P-ONE/GCD_files/PONE_Phase1.i3.gz'
gcd = dataio.I3File(gcd_file)
cframe = gcd.pop_frame()
geometry = cframe["I3Geometry"]
omgeo = geometry.omgeo
print('loaded geometry')

infile = '/data/p-one/akatil/step_4_medium_water/step_4_506_medium_water_custom_mDOM_noise.i3.gz'
outfile = '/data/p-one/akatil/test/step_4_506_medium_water_doublePeak.i3.gz'

tray = I3Tray()

tray.AddModule('I3Reader', 'reader',
            FilenameList = [gcd_file, infile]
            )

tray.AddModule(SelectingDoublePeakDOM.definingDOMs, "Double Peak Selector",
               omgeo = omgeo,
               InputMCPETree = "MCPESeriesMap",
               NoiseMCPETree = "NoiseSeriesMap",
               OutputMCPETree = "DoublePeakSeriesMap")

tray.AddModule("I3Writer","writer",
               #SkipKeys=SkipKeys,
               Filename = outfile,
               Streams = [icetray.I3Frame.DAQ, icetray.I3Frame.Physics],
              )

tray.AddModule("TrashCan","adios")
tray.Execute()
tray.Finish()
