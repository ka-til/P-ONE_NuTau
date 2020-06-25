from icecube import icetray, dataio, dataclasses, simclasses, clsim
from icecube.icetray import I3Units, OMKey, I3Frame
from icecube.dataclasses import ModuleKey
from os.path import expandvars
import scipy.constants as spc
import numpy as np
import matplotlib.pylab as plt

file = dataio.I3File('/data/p-one/akatil/step_3_medium_water/Custom/step_3_525_medium_water_custom_mDOM_katil.i3.gz')

gcd_file = dataio.I3File('/home/users/akatil/P-ONE/GCD_files/PONE_Phase1.i3.gz')
cframe = gcd_file.pop_frame()
geometry = cframe["I3Geometry"]
geoMap = geometry.omgeo
print('loaded geometry')

maxMCPE = 0
for frame in file:
    print('new Frame')
    mctree = frame["I3MCTree"]
    primary = mctree.primaries
    tau = dataclasses.I3MCTree.first_child(mctree, primary[0].id)
    tauDirection = tau.dir
    tauDecayProd = dataclasses.I3MCTree.get_daughters(mctree, tau.id)
    daughterTauPos = tau.pos
    daughterTauEnergies = tau.energy

    '''
    Note: The MCTree has a weird book keeping and a tau lepton has daughter tau leptons. There is no violation of
    physics here but rather the MCTree keeps track of all the stochastic losses of the tau. The final "daughter" tau will
    be responsible for production of the signature double pulse signal.
    '''

    print('finding the right tau')
    for i in range(0, len(tauDecayProd)):
        if tauDecayProd[i].type == 15 or tauDecayProd[i].type == -15:
            print('Has daughter Taus')
            daughterTauPos = tauDecayProd[i].pos                           #the last of the tau lepton will be the first interaction vertex in the double bang
            daughterTauEnergies = tauDecayProd[i].energy



    refractiveIndex = 1.333
    speed_of_light_water = (spc.c*1e-9)/refractiveIndex
    print('SPEED OF LIGHT IN WATER - ', speed_of_light_water)

    mcpeMap = frame['MCPESeriesMap']

    print('Finding OM Positions and time residuals')
    for omkey in mcpeMap.keys():
        oKey = geoMap.get(omkey)
        domPos = oKey.position
        dist_to_vertex = np.sqrt((domPos.x - daughterTauPos.x)**2+(domPos.y - daughterTauPos.y)**2+(domPos.z - daughterTauPos.z)**2)
        timeTaken = dist_to_vertex/speed_of_light_water

        mcpeList = mcpeMap[omkey]
        timeList = [mcpe.time for mcpe in mcpeList]
        timeList = np.array(timeList)
        timeRes = timeList - timeTaken
        if len(timeList) > maxMCPE and len(timeList) < 500:
            print('-----TIME LIST ------', len(timeList))
            maxMCPE = len(timeList)
            timeRes_single_DOM = timeList

print('plotting')
bins = np.arange(min(timeRes_single_DOM), min(timeRes_single_DOM)+181, 1)
print(bins[1:]-bins[:-1])
plt.hist(timeRes_single_DOM, bins=bins)
plt.savefig('testTimeRes1.pdf', dpi=200)
