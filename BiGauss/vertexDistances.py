from icecube import dataclasses, dataio, icetray
from icecube.icetray import I3Units
from icecube.icetray import OMKey
import matplotlib.pyplot as plt
import scipy.constants as spc
import numpy as np

gcd_file = dataio.I3File('/home/users/akatil/P-ONE/GCD_files/PONE_Phase1.i3.gz')
cframe = gcd_file.pop_frame()
geometry = cframe["I3Geometry"]
omgeo = geometry.omgeo

timeDiff = ([])
fv = ([])
sv = ([])

for file_num in range(0, 2000):
    file = dataio.I3File('/data/p-one/akatil/step_1_medium_water/step_1_'+str(file_num)+'_PONE_Phase1_NuTau_NuE.i3.gz')
    for frame in file:
        mctree = frame["I3MCTree"]
        primary = mctree.primaries
        lepton = dataclasses.I3MCTree.first_child(mctree, primary[0].id)

        if lepton.type == 15 or lepton.type == -15:
            tau_daughters = dataclasses.I3MCTree.get_daughters(mctree, lepton.id)
            tau_pos = lepton.pos
            x_tau_pos = tau_pos.x
            y_tau_pos = tau_pos.y
            z_tau_pos = tau_pos.z

            for td in range(0, len(tau_daughters)):
                if tau_daughters[td].type == 16 or tau_daughters[td].type == -16:
                    #print(tau_daughters[td])
                    tau_daughters_pos = tau_daughters[td].pos
                    x_td_pos = tau_daughters_pos.x
                    y_td_pos = tau_daughters_pos.y
                    z_td_pos = tau_daughters_pos.z

            for i in omgeo.keys():
                oKey = omgeo.get(i)
                domPos = oKey.position
                x_dom = domPos.x
                y_dom = domPos.y
                z_dom = domPos.z

                firstVertex = np.sqrt((x_dom - x_tau_pos)**2 + (y_dom - y_tau_pos)**2 + (z_dom - z_tau_pos)**2)
                secondVertex = np.sqrt((x_dom - x_td_pos)**2 + (y_dom - y_td_pos)**2 + (z_dom - z_td_pos)**2)
                refractiveIndex = 1.333
                speed_of_light_water = (spc.c)/refractiveIndex #[Units: m/seconds]
                speed_of_light_ns = speed_of_light_water
                #print(firstVertex-secondVertex)
                tDiff_ns = ((firstVertex - secondVertex)/speed_of_light_water) * 1e9 #[Units: nanoseconds]
                timeDiff = np.append(timeDiff, tDiff_ns)
                fv = np.append(fv, firstVertex)
                sv = np.append(sv, secondVertex)
                #print(tDiff_ns)

    print('file_number -', file_num)


bins = np.arange(min(timeDiff), max(timeDiff)+1, 1)
plt.hist(timeDiff,  bins = bins, histtype = 'step', log=True)
plt.title('Time Difference between 2 Vertices')
plt.xlabel('Time Difference(ns)')
plt.ylabel('Count')
plt.savefig('/home/users/akatil/P-ONE/git/PONE_NuTau/BiGauss/timeDiff_vertices.pdf', dpi=200)
plt.clf()

plt.hist(fv, bins=500, histtype='step', label = 'First Vertex')
plt.hist(sv, bins=500, histtype='step', label = 'Second Vertex')
plt.title('Distance to DOM')
plt.xlabel('Distance(m)')
plt.ylabel('Count')
plt.legend()
plt.savefig('/home/users/akatil/P-ONE/git/PONE_NuTau/BiGauss/distance_vertices.pdf', dpi=200)
plt.clf()
