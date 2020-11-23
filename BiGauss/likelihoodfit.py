from icecube import icetray, dataio, dataclasses, simclasses, clsim
from icecube.icetray import I3Units, OMKey, I3Frame
from icecube.dataclasses import ModuleKey
from os.path import expandvars
import scipy.constants as spc
import scipy as sc
import numpy as np
import matplotlib.pylab as plt
from scipy import stats
from scipy.optimize import minimize
from scipy.stats.distributions import chi2

import scipy

parser = argparse.ArgumentParser(description = "Takes I3Photons from step2 of the simulations and generates DOM hits")
parser.add_argument('-i', '--infile', dest = 'infile', help= 'input file and path')
parser.add_argument('-o', '--outfile', dest = 'outfile', help= 'output file and path')
args = parser.parse_args()

def gaussian(x, pos, wid, amp):
    y = amp*np.exp(-4*np.log(2)*((x-pos)/(wid))**2)
    return y

def biGauss(x, pos, wid, r, amp):
    mask = x < pos

    y_all = ([])
    for i in range(0, len(mask)):

        if mask[i] == True:
            m = 1
            nm = 0
        else:
            m = 0
            nm = 1
        if r != 0:
            y1 = gaussian(x[i],pos,r*wid/(r+1),amp)*m
            y2 = gaussian(x[i],pos,wid/(r+1),amp)*nm
            y = y1 + y2
        else:
            y = gaussian(x[i],pos,wid, amp)*nm

        y_all = np.append(y_all, y)
    return y_all

def double_peak(x, pos1, wid1, r1, amp1, pos2, wid2, r2, amp2):
    b1 = biGauss(x, pos1, wid1, r1, amp1)
    b2 = biGauss(x, pos2, wid2, r2, amp2)
    b = np.append(b1, b2)
    return b1+b2

def log_likelihood_biGauss(theta, n, x):
    pos, wid, r, amp = theta
    model = biGauss(x, pos, wid, r, amp)
    L = np.log(scipy.special.factorial(n)) + model - (n*np.log(model))
    return np.sum(L)

def log_likelihood_doublePeak(theta, n, x):
    pos1, wid1, r1, amp1, pos2, wid2, r2, amp2 = theta
    model = double_peak(x, pos1, wid1, r1, amp1, pos2, wid2, r2, amp2)
    L = np.log(scipy.special.factorial(n)) + model - (n*np.log(model))
    return np.sum(L)

def likelihood_ratio_doublePeak(x, n, pos1, wid1, r1, amp1, pos2, wid2, r2, amp2):
    model = double_peak(x, pos1, wid1, r1, amp1, pos2, wid2, r2, amp2)
    val = model - n + (n*np.log(n/model))
    #print('log - ', n/model, 'n - ', n)
    return np.sum(val)

def likelihood_ratio_biGauss(x, n, pos, wid, r, amp):
    model = biGauss(x, pos, wid, r, amp)
    val = model - n + (n*np.log(n/model))
    #print('log - ', n/model, 'n - ', n)
    return np.sum(val)

gcd_file = dataio.I3File('/home/users/akatil/P-ONE/GCD_files/PONE_Phase1.i3.gz')
cframe = gcd_file.pop_frame()
geometry = cframe["I3Geometry"]
omgeo = geometry.omgeo
print('loaded geometry')

tau_timeDiff = ([])
tau_pVal = ([])
tau_LRR = ([])
tau_wid1_ratio = ([])
tau_wid2_ratio = ([])
tau_amp1_ratio = ([])
tau_amp2_ratio = ([])

e_timeDiff = ([])
e_pVal = ([])
e_LRR = ([])
e_wid1_ratio = ([])
e_wid2_ratio = ([])
e_amp1_ratio = ([])
e_amp2_ratio = ([])

tau_wid_ratio_dp = ([])
e_wid_ratio_dp = ([])
tau_amp_ratio_dp = ([])
e_amp_ratio_dp = ([])

for i in range(0, 1):
    print('FILE NUMBER - ', i)
    #file = dataio.I3File('/data/p-one/akatil/step_4_medium_water/NuTau_NuE_20Events/step_4_'+str(i)+'_medium_water_custom_mDOM_noise.i3.gz')
    file = dataio.I3File(str(args.infile))

    f = 1
    for frame in file:
        print('frame num - ', f)
        mctree = frame["I3MCTree"]
        primary = mctree.primaries
        lepton = dataclasses.I3MCTree.first_child(mctree, primary[0].id)

        '''
        Removing NCC interations of the neutrino
        '''

        if lepton.type == 12 or lepton.type == -12 or lepton.type == 16 or lepton.type == -16:
            continue

        '''
        Lepton position
        '''
        lepton_pos = lepton.pos
        x_lepton_pos = lepton_pos.x
        y_lepton_pos = lepton_pos.y
        z_lepton_pos = lepton_pos.z

        mcpeMap = frame['MCPESeriesMap']
        noiseMap = frame['NoiseSeriesMap']

        #print('Finding OM Positions and time residuals')


        #looping through doms that have physics hits
        for omkey in mcpeMap.keys():
            oKey = omgeo.get(omkey)

            '''
            Dom Positons
            '''
            domPos = oKey.position
            x_dom = domPos.x
            y_dom = domPos.y
            z_dom = domPos.z

            '''
            Distance between event vertex and DOM
            '''
            distance = np.sqrt((x_dom - x_lepton_pos)**2 + (y_dom - y_lepton_pos)**2 +
                               (z_dom - z_lepton_pos)**2)

            #removing doms with distances > 100m from the event vertex
            if distance > 200:
                continue

            '''
            Obtaining the timeList
            '''
            noise_mcpeList = noiseMap[omkey]
            noise_timeList = np.array([mcpe.time for mcpe in noise_mcpeList])
            mcpeList = mcpeMap[omkey]
            timeList = np.array([mcpe.time for mcpe in mcpeList])
            tot_timeList = np.append(timeList, noise_timeList)


            '''
            Removing DOMs with hits less than 100
            '''
            if len(tot_timeList) < 100:
                continue


            '''
            Calculating the mean and removing the tails
            '''

            timeList = timeList[timeList < min(timeList)+30]

            mean_physicsHits = timeList.mean()
            mean_tot = tot_timeList.mean()

            select_time = tot_timeList[(tot_timeList > mean_physicsHits-50) & (tot_timeList < mean_physicsHits+50)]
            new_mean = select_time.mean()

            bins = np.arange(min(select_time), max(select_time), 1)
            max_hitTimes = select_time[(select_time > (new_mean-40))&(select_time < (new_mean+40))]

            z = stats.zscore(max_hitTimes)
            #using zscore to remove the effect of outliers from the analysis]
            max_hitTimes = max_hitTimes[(z < 1.2)]
            new_mean = max_hitTimes.mean()

            num_photons = len(max_hitTimes[max_hitTimes>0])

            #if len(max_hitTimes) < 100:
                #continue

            if len(max_hitTimes) < 10:
                continue

            if np.log10(num_photons) >= 3.0 or np.log10(num_photons) < 2.5:
                continue

            '''
            Histogramming the data from simulation
            '''

            bins = np.arange(min(max_hitTimes), max(max_hitTimes), 1)
            num, bin_edges = np.histogram(max_hitTimes, bins=bins)
            bin_centers = (bin_edges[:-1]+bin_edges[1:])/2

            '''
            Removing hits for DOMs that have more than 200 hits in 1 second bin.
            '''
            if max(num) > 175 or len(num) == 0:
                continue

            print('LOG LIKELIHOOD')

            '''
            Removing DOMs that don't have less than 8 non zero bins
            '''
            if len(num[num>0]) <= 8:
                continue


            '''
            Fitting bifurcated Gaussian and double bifurcated gaussian to the mcpe hit time distributions
            for both tau and electron.
            '''

            nll = lambda *args: log_likelihood_biGauss(*args)
            initial_biGauss = np.array([new_mean, 20, 1, 10])
            #bnds_biGauss = ((min(bin_centers), mean_timeArrival), (0, 20), (0, 2), (0, max(num)), (mean_timeArrival, max(bin_centers)), (0, 20), (0, 2), (0, max(num)))

            #print(len(num), len(initial_biGauss), initial_biGauss)
            bnds_biGauss = ((min(bin_centers), max(bin_centers)), (0, 1e3), (0, 100), (0, 1e6))
            soln_biGauss = minimize(log_likelihood_biGauss, initial_biGauss, args=(num, bin_centers),
                                    method='Powell', bounds = bnds_biGauss)

            nll = lambda *args: log_likelihood_doublePeak(*args)
            initial_doublePeak = np.array([min(bin_centers), 20, 1, 10, new_mean, 20, 1, 10])
            bnds_doublePeak = ((min(bin_centers), new_mean), (0, 1e3), (0, 100), (0, 1e6),
                               (new_mean, max(bin_centers)), (0, 1e3), (0, 100), (0,1e6))
            soln_doublePeak = minimize(log_likelihood_doublePeak, initial_doublePeak, args=(num, bin_centers),
                                       method='Powell',bounds=bnds_doublePeak)

            '''
            Removing DOMs whose minimization is not successful
            '''
            if soln_biGauss.success == False or soln_doublePeak.success == False:
                continue

            '''
            Calculating the Likelihood ratio for bifurcated gaussian and double double bifurcated gaussian
            '''
            LR_biGauss = likelihood_ratio_biGauss(bin_centers[num>0], num[num>0], soln_biGauss.x[0],
                                             soln_biGauss.x[1], soln_biGauss.x[2], soln_biGauss.x[3])
            LR_doublePeak = likelihood_ratio_doublePeak(bin_centers[num>0], num[num>0], soln_doublePeak.x[0], soln_doublePeak.x[1],
                                             soln_doublePeak.x[2], soln_doublePeak.x[3], soln_doublePeak.x[4],
                                             soln_doublePeak.x[5], soln_doublePeak.x[6], soln_doublePeak.x[7])


            '''
            Calculating the p-value using the likelihood ratio
            '''
            pVal_biGauss = chi2.sf(LR_biGauss, len(num) - 4)
            pVal_doublePeak = chi2.sf(LR_doublePeak, len(num) - 8)

            if pVal_biGauss != pVal_biGauss:
                print('BiGauss gives not well defined - ', str(lepton.type))
                print('Minimisation - ', soln_biGauss.success)
                print('Degrees of Freedom - ', len(num) - 4)
                print('Log Likelihood - ', LR_biGauss)
            if pVal_doublePeak != pVal_doublePeak:
                print('double peak gives not well defined - ', str(lepton.type))
                print('Minimisation - ', soln_doublePeak.success)
                print('Degrees of Freedom - ', len(num) - 8)
                print('Log Likelihood - ', LR_doublePeak)


            '''
            (x, y) values for the fit
            '''
            x = np.linspace(min(bin_centers), max(bin_centers), 1000)
            #x = np.linspace(0, max(bin_centers)+1e5, 1000)
            y_biGauss = biGauss(x, soln_biGauss.x[0],
                                             soln_biGauss.x[1], soln_biGauss.x[2], soln_biGauss.x[3])
            y_doublePeak = double_peak(x, soln_doublePeak.x[0], soln_doublePeak.x[1],
                                             soln_doublePeak.x[2], soln_doublePeak.x[3], soln_doublePeak.x[4],
                                             soln_doublePeak.x[5], soln_doublePeak.x[6], soln_doublePeak.x[7])

            '''
            Calculating the time difference and p-value ratio of bigauss and double peak
            '''
            timeDifference_doublePeak = soln_doublePeak.x[4] - soln_doublePeak.x[0]
            pVal_ratio = pVal_doublePeak/pVal_biGauss
            LRR = LR_doublePeak/LR_biGauss
            wid1_ratio = soln_doublePeak.x[1]/soln_biGauss.x[1]
            wid2_ratio = soln_doublePeak.x[5]/soln_biGauss.x[1]
            amp1_ratio = soln_doublePeak.x[3]/soln_biGauss.x[3]
            amp2_ratio = soln_doublePeak.x[7]/soln_biGauss.x[3]
            wid1_wid2 = soln_doublePeak.x[1]/soln_doublePeak.x[5]
            amp1_amp2 = soln_doublePeak.x[3]/soln_doublePeak.x[7]



            '''
            Removing terrible fits
            '''
            if abs(timeDifference_doublePeak) > 100:
                continue

            amp1 = soln_doublePeak.x[3]
            amp2 = soln_doublePeak.x[7]
            if amp1/amp2 < 1/4 and amp1/amp2 > 4:
                continue

            if amp1 < 0 or amp2 < 0:
                print('Error in amp')
                continue


            '''
            plot mcpe time distributions obtained using simulations and the fits
            '''
            #plt.figure(figsize=(10,9))
            #_ = plt.hist(max_hitTimes, bins=bins, histtype='step')
            #plt.title(str(lepton.type))
            #plt.plot(x, y_biGauss, '--', c = 'r')
            #plt.plot(x, y_doublePeak, '--', c = 'k')
            #plt.axvline(new_mean, c = 'b')



            '''
            Separating the time difference calculated above and appending the values
            '''

            '''
            Tau
            '''
            if lepton.type == 15 or lepton.type == -15:
                tau_timeDiff = np.append(tau_timeDiff, timeDifference_doublePeak)
                tau_pVal = np.append(tau_pVal, pVal_ratio)
                tau_LRR = np.append(tau_LRR, LRR)
                tau_wid1_ratio = np.append(tau_wid1_ratio, wid1_ratio)
                tau_wid2_ratio = np.append(tau_wid2_ratio, wid2_ratio)
                tau_amp1_ratio = np.append(tau_amp1_ratio, amp1_ratio)
                tau_amp2_ratio = np.append(tau_amp2_ratio, amp2_ratio)

                tau_wid_ratio_dp = np.append(tau_wid_ratio_dp, wid1_wid2)
                tau_amp_ratio_dp = np.append(tau_amp_ratio_dp, amp1_amp2)
                #plt.title('E')

            '''
            Electron
            '''

            if lepton.type == 11 or lepton.type == -11:
                e_timeDiff = np.append(e_timeDiff, timeDifference_doublePeak)
                e_pVal = np.append(e_pVal, pVal_ratio)
                e_LRR = np.append(e_LRR, LRR)
                e_wid1_ratio = np.append(e_wid1_ratio, wid1_ratio)
                e_wid2_ratio = np.append(e_wid2_ratio, wid2_ratio)
                e_amp1_ratio = np.append(e_amp1_ratio, amp1_ratio)
                e_amp2_ratio = np.append(e_amp2_ratio, amp2_ratio)

                e_wid_ratio_dp = np.append(e_wid_ratio_dp, wid1_wid2)
                e_amp_ratio_dp = np.append(e_amp_ratio_dp, amp1_amp2)



                #plt.title('Tau')


            print('P-VAL CALCULATED')


        '''
        print(tot_timeList)

        bins = np.arange(min(tot_timeList), min(tot_timeList)+41, 1)
        num, bin_edges, _ = plt.hist(tot_timeList, bins=bins, histtype='step')
        plt.title('Tau')
        '''
        f = f+1
#print(tot_timeList)

np.savetxt(str(args.outfile) + '_timeDifference_tau.csv', tau_timeDiff, delimiter=',')
np.savetxt(str(args.outfile) + '_timeDifference_e.csv', e_timeDiff, delimiter=',')

np.savetxt(str(args.outfile) + '_pval_e.csv', e_pVal, delimiter=',')
np.savetxt(str(args.outfile) + '_pval_tau.csv', tau_pVal, delimiter=',')

np.savetxt(str(args.outfile) + '_LRR_tau.csv', tau_LRR, delimiter=',')
np.savetxt(str(args.outfile) + '_LRR_e.csv', e_LRR, delimiter=',')

np.savetxt(str(args.outfile) + '_wid1_tau.csv', tau_wid1_ratio, delimiter=',')
np.savetxt(str(args.outfile) + '_wid1_e.csv', e_wid1_ratio, delimiter=',')

np.savetxt(str(args.outfile) + '_wid2_tau.csv', tau_wid2_ratio, delimiter=',')
np.savetxt(str(args.outfile) + '_wid2_e.csv', e_wid2_ratio, delimiter=',')

np.savetxt(str(args.outfile) + '_amp1_tau.csv', tau_amp1_ratio, delimiter=',')
np.savetxt(str(args.outfile) + '_amp1_e.csv', e_amp1_ratio, delimiter=',')

np.savetxt(str(args.outfile) + '_amp2_tau.csv', tau_amp2_ratio, delimiter=',')
np.savetxt(str(args.outfile) + '_amp2_e.csv', e_amp2_ratio, delimiter=',')

np.savetxt(str(args.outfile) + '_wid1_wid2_tau.csv', tau_wid_ratio_dp, delimiter=',')
np.savetxt(str(args.outfile) + '_wid1_wid2_e.csv', e_wid_ratio_dp, delimiter=',')

np.savetxt(str(args.outfile) + '_amp1_amp2_tau.csv', tau_amp_ratio_dp, delimiter=',')
np.savetxt(str(args.outfile) + '_amp1_amp2_e.csv', e_amp_ratio_dp, delimiter=',')
