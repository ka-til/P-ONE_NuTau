import corner
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm
import matplotlib.lines as mlines
from likelihoodHelpers import log_likelihood_biGauss, log_likelihood_doublePeak, likelihood_ratio_doublePeak, likelihood_ratio_biGauss, biGauss, double_peak
from likelihoodHelpers import log_likelihood_expGauss, log_likelihood_expDoublePeak, expGauss, expDoublePeak

class plots2(object):

    def __init__(self, params, params_perEvent):
        print("Initializing parameters")

        #Params
        self.LLHDiff_t = params[0]
        self.energy_t = params[1]
        self.numHits_t = params[2]
        self.aAsym_t = params[3]
        self.tLength =params[4]

        self.LLHDiff_e = params[5]
        self.energy_e = params[6]
        self.numHits_e =params[7]
        self.aAsym_e = params[8]

        #Params_perEvent
        self.largest_LLHDiff_t = params_perEvent[0]
        self.energyPerEvent_t = params_perEvent[1]
        self.numHitsPerEvent_t = params_perEvent[2]
        self.aAsymPerEvent_t = params_perEvent[3]
        self.tLenghtPerEvent = params_perEvent[4]

        self.largest_LLHDiff_e = params_perEvent[5]
        self.energyPerEvent_e = params_perEvent[6]
        self.numHitsPerEvent_e = params_perEvent[7]
        self.aAsymPerEvent_e = params_perEvent[8]

    def setlog(self, t, e):
        return np.log10(t), np.log10(e)

    def setRange(self, t, e, r_min, r_max):
        return t[(t>=r_min)&(t<=r_max)], e[(e>=r_min)&(e<=r_max)]

    def logRange(self, range):
        return [np.log10(range[0]), np.log10(range[1])]

    def bin_range(self, log, t, e, numBins):
        if min(t) > min(e):
            min_val = min(e)
        else:
            min_val = min(t)

        if max(t) > max(e):
            max_val = max(t)
        else:
            max_val = max(e)

        '''
        if log == True:
            step = 1
            bins = np.linspace(min_val, max_val, 20)
            print(min_val, max_val, step)
        else:
            diff = max_val - min_val

            if diff > 0 and diff < 50:
                step = 10
                add = 2
            if diff > 50 and diff < 9e3:
                step = 10
                add = 20
            if diff > 9e3 and diff < 9e4:
                step = 1e2
                add = 200
            if diff > 9e4 and diff < 9e5:
                step = 1e3
                add = 2000
            if diff > 9e5 and diff < 9e6:
                step = 1e4
                add = step*2
            if diff > 9e6 and diff < 1e12:
                step = 1e5
                add = step*2
            '''
        bins = np.linspace(min_val, max_val, numBins)

        return bins

    def scatter(self, log_LLH=False, log_energy=False, log_numHits=False, log_tlength=False, log_asym=False, log_all=False,
                range_energy=[1e-12, 1e12], range_numHits=[1e-12, 1e12], range_tauLength=[1e-12, 1e12], range_LLH=[1e-12, 1e12], range_asym=[1e-12, 1e12]):

        '''
        Setting Range
        '''
        l_t, l_e = self.LLHDiff_t, self.LLHDiff_e
        e_t, e_e = self.energy_t, self.energy_e
        n_t, n_e = self.numHits_t, self.numHits_e
        a_t, a_e = self.aAsym_t, self.aAsym_e
        tl = self.tLength

        '''
        Changing to log scale
        '''
        if log_all == True:
            LLHDiff_t, LLHDiff_e = self.setlog(l_t, l_e)
            energy_t, energy_e = self.setlog(e_t, e_e)
            numHits_t, numHits_e = self.setlog(n_t, n_e)
            asym_t, asym_e = self.setlog(a_t, a_e)
            tLength = np.log10(tl)

            range_LLH = self.logRange(range_LLH)
            range_energy = self.logRange(range_energy)
            range_numHits = self.logRange(range_numHits)
            range_asym = self.logRange(range_asym)
            range_tLength = self.logRange(range_tauLength)

            label_LLH ='log 2DelLLH'
            label_energy = 'log energy[GeV]'
            label_numHits = 'log Num Hits'
            label_asym = 'log Amp Asymmetry'
            label_tLength = 'log Tau Length'
        else:
            if log_energy == True:
                energy_t, energy_e = self.setlog(e_t, e_e)
                range_energy = self.logRange(range_energy)
                label_energy ='log energy[GeV]'
            else:
                energy_t, energy_e = e_t, e_e
                label_energy ='energy[GeV]'

            if log_numHits == True:
                numHits_t, numHits_e = self.setlog(n_t, n_e)
                range_numHits = self.logRange(range_numHits)
                label_numHits ='log Num Hits'
            else:
                numHits_t, numHits_e = n_t, n_e
                label_numHits ='Num Hits'

            if log_skew == True:
                tLength = np.log10(tl)
                range_tLength = self.logRange(range_tauLength)
                label_skew ='log Tau Length'
            else:
                tLength = tl
                label_skew ='Tau Length'

            if log_asym == True:
                asym_t, asym_e = self.setlog(a_t, a_e)
                range_asym = self.logRange(range_asym)
                label_asym ='log Amp Asymmetry'
            else:
                asym_t, asym_e = a_t, a_e
                label_asym ='Amp Asymmetry'

            if log_LLH == True:
                LLHDiff_t, LLHDiff_e = self.setlog(l_t, l_e)
                range_LLH = self.logRange(range_LLH)
                label_LLH ='log 2DelLLH'
            else:
                LLHDiff_t, LLHDiff_e = l_t, l_e
                label_LLH ='2DelLLH'

        print(len(numHits_t), len(LLHDiff_t))

        plt.figure(figsize=(10,9))
        plt.scatter(energy_t, LLHDiff_t, label='tau', c ='r', alpha = 0.6)
        plt.legend(fontsize=16)
        plt.xlim(range_energy)
        plt.ylim(range_LLH)
        plt.xlabel(label_energy, fontsize = 18)
        plt.ylabel(label_LLH, fontsize = 18)

        plt.figure(figsize=(10,9))
        plt.scatter(energy_e, LLHDiff_e, label='e&NC', c='skyblue', alpha = 0.6)
        plt.legend(fontsize=16)
        plt.xlim(range_energy)
        plt.ylim(range_LLH)
        plt.xlabel(label_energy, fontsize = 18)
        plt.ylabel(label_LLH, fontsize = 18)

        plt.figure(figsize=(10,9))
        plt.scatter(numHits_t, LLHDiff_t, label='tau', c ='r', alpha = 0.6)
        plt.legend(fontsize=16)
        plt.xlim(range_numHits)
        plt.ylim(range_LLH)
        plt.xlabel(label_numHits, fontsize = 18)
        plt.ylabel(label_LLH, fontsize = 18)

        plt.figure(figsize=(10,9))
        plt.scatter(numHits_e, LLHDiff_e, label='e&NC', c='skyblue', alpha = 0.6)
        plt.legend(fontsize=16)
        plt.xlim(range_numHits)
        plt.ylim(range_LLH)
        plt.xlabel(label_numHits, fontsize = 18)
        plt.ylabel(label_LLH, fontsize = 18)

        plt.figure(figsize=(10,9))
        plt.scatter(tLength, LLHDiff_t, label='tau', c ='r', alpha = 0.6)
        plt.legend(fontsize=16)
        plt.xlim(range_tLength)
        plt.ylim(range_LLH)
        plt.xlabel(label_tLength, fontsize = 18)
        plt.ylabel(label_LLH, fontsize = 18)

        plt.figure(figsize=(10,9))
        plt.scatter(asym_t, LLHDiff_t, label='tau', c ='r', alpha = 0.6)
        plt.legend(fontsize=16)
        plt.xlim(range_asym)
        plt.ylim(range_LLH)
        plt.xlabel(label_asym, fontsize = 18)
        plt.ylabel(label_LLH, fontsize = 18)

        plt.figure(figsize=(10,9))
        plt.scatter(asym_e, LLHDiff_e, label='e&NC', c='skyblue', alpha = 0.6)
        plt.legend(fontsize=16)
        plt.xlim(range_asym)
        plt.ylim(range_LLH)
        plt.xlabel(label_asym, fontsize = 18)
        plt.ylabel(label_LLH, fontsize = 18)


    def scatterPerEvent(self, log_LLH=False, log_energy=False, log_numHits=False, log_tlength=False, log_asym=False, log_all=False,
                range_energy=[1e-12, 1e12], range_numHits=[1e-12, 1e12], range_tauLength=[1e-12, 1e12], range_LLH=[1e-12, 1e12], range_asym=[1e-12, 1e12]):

        '''
        Setting Range
        '''
        l_t, l_e = self.largest_LLHDiff_t, self.largest_LLHDiff_e
        e_t, e_e = self.energyPerEvent_t, self.energyPerEvent_e
        n_t, n_e = self.numHitsPerEvent_t, self.numHitsPerEvent_e
        a_t, a_e = self.aAsymPerEvent_t, self.aAsymPerEvent_e
        tl = self.tLenghtPerEvent

        '''
        Changing to log scale
        '''
        if log_all == True:
            LLHDiff_t, LLHDiff_e = self.setlog(l_t, l_e)
            energy_t, energy_e = self.setlog(e_t, e_e)
            numHits_t, numHits_e = self.setlog(n_t, n_e)
            asym_t, asym_e = self.setlog(a_t, a_e)
            tLength = np.log10(tl)

            range_LLH = self.logRange(range_LLH)
            range_energy = self.logRange(range_energy)
            range_numHits = self.logRange(range_numHits)
            range_asym = self.logRange(range_asym)
            range_tLength = self.logRange(range_tauLength)

            label_LLH ='log 2DelLLH'
            label_energy = 'log energy[GeV]'
            label_numHits = 'log Num Hits'
            label_asym = 'log Amp Asymmetry'
            label_tLength = 'log Tau Length'
        else:
            if log_energy == True:
                energy_t, energy_e = self.setlog(e_t, e_e)
                range_energy = self.logRange(range_energy)
                label_energy ='log energy[GeV]'
            else:
                energy_t, energy_e = e_t, e_e
                label_energy ='energy[GeV]'

            if log_numHits == True:
                numHits_t, numHits_e = self.setlog(n_t, n_e)
                range_numHits = self.logRange(range_numHits)
                label_numHits ='log Num Hits'
            else:
                numHits_t, numHits_e = n_t, n_e
                label_numHits ='Num Hits'

            if log_skew == True:
                tLength = np.log10(tl)
                range_tLength = self.logRange(range_tauLength)
                label_skew ='log Tau Length'
            else:
                tLength = tl
                label_skew ='Tau Length'

            if log_asym == True:
                asym_t, asym_e = self.setlog(a_t, a_e)
                range_asym = self.logRange(range_asym)
                label_asym ='log Amp Asymmetry'
            else:
                asym_t, asym_e = a_t, a_e
                label_asym ='Amp Asymmetry'

            if log_LLH == True:
                LLHDiff_t, LLHDiff_e = self.setlog(l_t, l_e)
                range_LLH = self.logRange(range_LLH)
                label_LLH ='log 2DelLLH'
            else:
                LLHDiff_t, LLHDiff_e = l_t, l_e
                label_LLH ='2DelLLH'

        print(len(numHits_t), len(LLHDiff_t))

        plt.figure(figsize=(10,9))
        plt.scatter(energy_t, LLHDiff_t, label='tau', c ='r', alpha = 0.6)
        plt.legend(fontsize=16)
        plt.xlim(range_energy)
        plt.ylim(range_LLH)
        plt.xlabel(label_energy, fontsize = 18)
        plt.ylabel(label_LLH, fontsize = 18)

        plt.figure(figsize=(10,9))
        plt.scatter(energy_e, LLHDiff_e, label='e&NC', c='skyblue', alpha = 0.6)
        plt.legend(fontsize=16)
        plt.xlim(range_energy)
        plt.ylim(range_LLH)
        plt.xlabel(label_energy, fontsize = 18)
        plt.ylabel(label_LLH, fontsize = 18)


        plt.figure(figsize=(10,9))
        plt.scatter(numHits_t, LLHDiff_t, label='tau', c ='r', alpha = 0.6)
        plt.legend(fontsize=16)
        plt.xlim(range_numHits)
        plt.ylim(range_LLH)
        plt.xlabel(label_numHits, fontsize = 18)
        plt.ylabel(label_LLH, fontsize = 18)

        plt.figure(figsize=(10,9))
        plt.scatter(numHits_e, LLHDiff_e, label='e&NC', c='skyblue', alpha = 0.6)
        plt.legend(fontsize=16)
        plt.xlim(range_numHits)
        plt.ylim(range_LLH)
        plt.xlabel(label_numHits, fontsize = 18)
        plt.ylabel(label_LLH, fontsize = 18)

        plt.figure(figsize=(10,9))
        plt.scatter(tLength, LLHDiff_t, label='tau', c ='r', alpha = 0.6)
        plt.legend(fontsize=16)
        plt.xlim(range_tLength)
        plt.ylim(range_LLH)
        plt.xlabel(label_tLength, fontsize = 18)
        plt.ylabel(label_LLH, fontsize = 18)

        plt.figure(figsize=(10,9))
        plt.scatter(asym_t, LLHDiff_t, label='tau', c ='r', alpha = 0.6)
        plt.legend(fontsize=16)
        plt.xlim(range_asym)
        plt.ylim(range_LLH)
        plt.xlabel(label_asym, fontsize = 18)
        plt.ylabel(label_LLH, fontsize = 18)

        plt.figure(figsize=(10,9))
        plt.scatter(asym_e, LLHDiff_e, label='e&NC', c='skyblue', alpha = 0.6)
        plt.legend(fontsize=16)
        plt.xlim(range_asym)
        plt.ylim(range_LLH)
        plt.xlabel(label_asym, fontsize = 18)
        plt.ylabel(label_LLH, fontsize = 18)
