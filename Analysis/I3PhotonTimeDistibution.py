from icecube import dataclasses, dataio, icetray, simclasses
from icecube.icetray import I3Units, I3Frame, OMKey
import matplotlib.pyplot as plt
import numpy as np
import argparse, matplotlib, csv
from statistics import mean
from icecube.phys_services import I3Calculator
import os

directory='/data/p-one/akatil/step_2_medium_water/'
file = ([])
for filename in os.listdir(directory):
    if filename.startswith("step_2"):
        if os.path.getsize(directory+str(filename)) > 0 and os.path.getsize(directory+str(filename)) < 2000000:
            file.append(filename)

time = ([])
timeDiffinFrame = ([])
for i in range(0, 50):
    pfile = dataio.I3File(directory+str(file[i]))
    frame = pfile.pop_frame()
    for pframe in pfile:
        photonMap = pframe['I3Photons']
        photonTime = ([])
        for modKey,photons in photonMap:
            for photon in photons:
                photonTime = np.append(photonTime, photon.time)

        timeDiffinFrame = np.append(timeDiffinFrame, (max(photonTime) - min(photonTime))*1e-3)


print(len(timeDiffinFrame,))
plt.hist(timeDiffinFrame, bins = 100, density = False, log=False, histtype='step')
plt.xlabel("T_last - T_first(us)")
plt.title('I3Photon TimeDifference')
plt.savefig('/home/users/akatil/P-ONE/git/PONE_NuTau/Analysis/I3PhotonsTimeDiff.pdf', dpi=200)
plt.clf()
