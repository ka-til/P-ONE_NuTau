from icecube import icetray, dataio, dataclasses, simclasses, clsim
from icecube.icetray import I3Units, OMKey, I3Frame
from icecube.dataclasses import ModuleKey
from os.path import expandvars
import numpy as np
#import matplotlib.pylab as plt
import argparse
from scipy import stats
from scipy.optimize import minimize
from scipy.stats.distributions import chi2
import scipy
from RecoPulseLikelihoodRatio import likelihoodfit
import csv

gcd_file = dataio.I3File('/home/users/akatil/P-ONE/GCD_files/PONE_Phase1.i3.gz')
cframe = gcd_file.pop_frame()
geometry = cframe["I3Geometry"]
omgeo = geometry.omgeo
print('loaded geometry')

f = open('/data/p-one/akatil/analysis/RecoPulses/RecoPulseFitInfo_correcectedchi2_20201022.csv', 'w')

with f:
    fnames = ['file', 'frame', 'lepton_type', 'DOM', 'string', 'binEntries_mean', 'success_biGauss', 'success_doublePeak',
             'biGauss_pos', 'biGauss_wid', 'biGauss_rat', 'biGauss_amp',
             'doublePeak_pos1', 'doublePeak_wid1', 'doublePeak_rat1', 'doublePeak_amp1',
             'doublePeak_pos2', 'doublePeak_wid2', 'doublePeak_rat2', 'doublePeak_amp2',
             'area_data', 'area_biGauss_fit', 'area_doublePeak_fit', 'gof_biGauss', 'gof_doublePeak', 'dof_biGauss', 'dof_doublePeak']
    writer = csv.DictWriter(f, fieldnames=fnames)
    writer.writeheader()

    for i in range(0, 2000):
        print(i)
        file = dataio.I3File('/data/p-one/akatil/step_5_medium_water/NuTau_NuE_20Events/step_5_'+str(i)+'_medium_water_custom_mDOM_recoPulse.i3.gz')

        frame_num = 0
        for frame in file:
            likelihoodfit(frame, omgeo, i, frame_num, writer)
            #print(info_list)
            frame_num += 1
