import corner
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm
import matplotlib.lines as mlines
from likelihoodHelpers import log_likelihood_biGauss, log_likelihood_doublePeak, likelihood_ratio_doublePeak, likelihood_ratio_biGauss, biGauss, double_peak
from likelihoodHelpers import log_likelihood_expGauss, log_likelihood_expDoublePeak, expGauss, expDoublePeak
class plots(object):

    def __init__(self, fitParams):
        print("Will plot some pretty plots using corner and matplotlib.pyplot")

        #Tau Parameters
        self.tDiff_t = fitParams[0]
        self.widDiff_t = fitParams[1]
        self.skewDiff_t = fitParams[2]
        self.ampRatio_t = fitParams[3]
        self.LLHDiff_t = fitParams[4]

        #E&NC parameters
        self.tDiff_e = fitParams[5]
        self.widDiff_e = fitParams[6]
        self.skewDiff_e = fitParams[7]
        self.ampRatio_e = fitParams[8]
        self.LLHDiff_e = fitParams[9]


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

    def scatter(self, log_time=False, log_wid=False, log_skew=False, log_amp=False, log_LLH=False, log_all=False,
                range_time=[1e-12, 200], range_wid=[1e-12, 100], range_skew=[1e-12, 10], range_amp=[1e-12, 1e12], range_LLH=[1e-12, 1e6]):

        '''
        Setting Range
        '''
        t_t, t_e = self.tDiff_t, self.tDiff_e
        w_t, w_e = self.widDiff_t, self.widDiff_e
        s_t, s_e = self.skewDiff_t, self.skewDiff_e
        a_t, a_e = self.ampRatio_t, self.ampRatio_e
        l_t, l_e = self.LLHDiff_t, self.LLHDiff_e

        '''
        Changing to log scale
        '''
        if log_all == True:
            tDiff_t, tDiff_e = self.setlog(t_t, t_e)
            widDiff_t, widDiff_e = self.setlog(w_t, w_e)
            skewDiff_t, skewDiff_e = self.setlog(s_t, s_e)
            ampRatio_t, ampRatio_e = self.setlog(a_t, a_e)
            LLHDiff_t, LLHDiff_e = self.setlog(l_t, l_e)
            range_time = self.logRange(range_time)
            range_wid = self.logRange(range_wid)
            range_skew = self.logRange(range_skew)
            range_amp = self.logRange(range_amp)
            range_LLH = self.logRange(range_LLH)
            label_time ='log time difference'
            label_width ='log width difference'
            label_skew ='log skewness difference'
            label_amp ='log amplitude ratio'
            label_LLH ='log 2DelLLH'
        else:
            if log_time == True:
                tDiff_t, tDiff_e = self.setlog(t_t, t_e)
                range_time = self.logRange(range_time)
                label_time ='log time difference'
            else:
                tDiff_t, tDiff_e = t_t, t_e
                label_time ='time difference[ns]'

            if log_wid == True:
                widDiff_t, widDiff_e = self.setlog(w_t, w_e)
                range_wid = self.logRange(range_wid)
                label_width ='log width difference'
            else:
                widDiff_t, widDiff_e = w_t, w_e
                label_width ='width difference[ns]'

            if log_skew == True:
                skewDiff_t, skewDiff_e = self.setlog(s_t, s_e)
                range_skew = self.logRange(range_skew)
                label_skew ='log skewness difference'
            else:
                skewDiff_t, skewDiff_e = s_t, s_e
                label_skew ='skewness difference'

            if log_amp == True:
                ampRatio_t, ampRatio_e = self.setlog(a_t, a_e)
                range_amp = self.logRange(range_amp)
                label_amp ='log amplitude ratio'
            else:
                ampRatio_t, ampRatio_e = a_t, a_e
                label_amp ='amplitude Ratio'

            if log_LLH == True:
                LLHDiff_t, LLHDiff_e = self.setlog(l_t, l_e)
                range_LLH = self.logRange(range_LLH)
                label_LLH ='log 2DelLLH'
            else:
                LLHDiff_t, LLHDiff_e = l_t, l_e
                label_LLH ='2DelLLH'

        plt.figure(figsize=(10,9))
        plt.scatter(tDiff_t, LLHDiff_t, label='tau', c ='r', alpha = 0.6)
        plt.scatter(tDiff_e, LLHDiff_e, label='e&NC', c='skyblue', alpha = 0.6)
        plt.legend(fontsize=16)
        plt.xlim(range_time)
        plt.ylim(range_LLH)
        plt.xlabel(label_time, fontsize = 18)
        plt.ylabel(label_LLH, fontsize = 18)

        plt.figure(figsize=(10,9))
        plt.scatter(tDiff_t, skewDiff_t, label='tau', c ='r', alpha = 0.6)
        plt.scatter(tDiff_e, skewDiff_e, label='e&NC', c='skyblue', alpha = 0.6)
        plt.legend(fontsize=16)
        plt.xlim(range_time)
        plt.ylim(range_skew)
        plt.xlabel(label_time, fontsize = 18)
        plt.ylabel(label_skew, fontsize = 18)

        plt.figure(figsize=(10,9))
        plt.scatter(tDiff_t, ampRatio_t, label='tau', c ='r', alpha = 0.6)
        plt.scatter(tDiff_e, ampRatio_e, label='e&NC', c='skyblue', alpha = 0.6)
        plt.legend(fontsize=16)
        plt.xlim(range_time)
        plt.ylim(range_amp)
        plt.xlabel(label_time, fontsize = 18)
        plt.ylabel(label_amp, fontsize = 18)

        plt.figure(figsize=(10,9))
        plt.scatter(tDiff_t, widDiff_t, label='tau', c ='r', alpha = 0.6)
        plt.scatter(tDiff_e, widDiff_e, label='e&NC', c='skyblue', alpha = 0.6)
        plt.legend(fontsize=16)
        plt.xlim(range_time)
        plt.ylim(range_wid)
        plt.xlabel(label_time, fontsize = 18)
        plt.ylabel(label_width, fontsize = 18)

        plt.figure(figsize=(10,9))
        plt.scatter(widDiff_t, LLHDiff_t, label='tau', c ='r', alpha = 0.6)
        plt.scatter(widDiff_e, LLHDiff_e, label='e&NC', c='skyblue', alpha = 0.6)
        plt.legend(fontsize=16)
        plt.xlim(range_wid)
        plt.ylim(range_LLH)
        plt.xlabel(label_width, fontsize = 18)
        plt.ylabel(label_LLH, fontsize = 18)

        plt.figure(figsize=(10,9))
        plt.scatter(widDiff_t, skewDiff_t, label='tau', c ='r', alpha = 0.6)
        plt.scatter(widDiff_e, skewDiff_e, label='e&NC', c='skyblue', alpha = 0.6)
        plt.legend(fontsize=16)
        plt.xlim(range_wid)
        plt.ylim(range_skew)
        plt.xlabel(label_width, fontsize = 18)
        plt.ylabel(label_skew, fontsize = 18)

        plt.figure(figsize=(10,9))
        plt.scatter(widDiff_t, ampRatio_t, label='tau', c ='r', alpha = 0.6)
        plt.scatter(widDiff_e, ampRatio_e, label='e&NC', c='skyblue', alpha = 0.6)
        plt.legend(fontsize=16)
        plt.xlim(range_wid)
        plt.ylim(range_amp)
        plt.xlabel(label_width, fontsize = 18)
        plt.ylabel(label_amp, fontsize = 18)

        plt.figure(figsize=(10,9))
        plt.scatter(ampRatio_t, LLHDiff_t, label='tau', c ='r', alpha = 0.6)
        plt.scatter(ampRatio_e, LLHDiff_e, label='e&NC', c='skyblue', alpha = 0.6)
        plt.legend(fontsize=16)
        plt.xlim(range_amp)
        plt.ylim(range_LLH)
        plt.xlabel(label_amp, fontsize = 18)
        plt.ylabel(label_LLH, fontsize = 18)

        plt.figure(figsize=(10,9))
        plt.scatter(ampRatio_t, skewDiff_t, label='tau', c ='r', alpha = 0.6)
        plt.scatter(ampRatio_e, skewDiff_e, label='e&NC', c='skyblue', alpha = 0.6)
        plt.legend(fontsize=16)
        plt.xlim(range_amp)
        plt.ylim(range_skew)
        plt.xlabel(label_amp, fontsize = 18)
        plt.ylabel(label_skew, fontsize = 18)

        plt.figure(figsize=(10,9))
        plt.scatter(skewDiff_t, LLHDiff_t, label='tau', c ='r', alpha = 0.6)
        plt.scatter(skewDiff_e, LLHDiff_e, label='e&NC', c='skyblue', alpha = 0.6)
        plt.legend(fontsize=16)
        plt.xlim(range_skew)
        plt.ylim(range_LLH)
        plt.xlabel(label_skew, fontsize = 18)
        plt.ylabel(label_LLH, fontsize = 18)


    def hist2d(self, log_time=False, log_wid=False, log_skew=False, log_amp=False, log_LLH=False, log_all=False,
               range_time=[1e-12, 200], range_wid=[1e-12, 100], range_skew=[1e-12, 10], range_amp=[1e-12, 1e12], range_LLH=[1e-12, 1e6],
               numBins_time=10, numBins_wid=10, numBins_skew=10, numBins_amp=100, numBins_LLH=100):

        print(range_time, range_wid, range_skew, range_amp, range_LLH)

        '''
        Setting Range
        '''
        t_t, t_e = self.tDiff_t, self.tDiff_e
        w_t, w_e = self.widDiff_t, self.widDiff_e
        s_t, s_e = self.skewDiff_t, self.skewDiff_e
        a_t, a_e = self.ampRatio_t, self.ampRatio_e
        l_t, l_e = self.LLHDiff_t, self.LLHDiff_e

        if log_all == True:
            tDiff_t, tDiff_e = self.setlog(t_t, t_e)
            widDiff_t, widDiff_e = self.setlog(w_t, w_e)
            skewDiff_t, skewDiff_e = self.setlog(s_t, s_e)
            ampRatio_t, ampRatio_e = self.setlog(a_t, a_e)
            LLHDiff_t, LLHDiff_e = self.setlog(l_t, l_e)
            range_time = self.logRange(range_time)
            range_wid = self.logRange(range_wid)
            range_skew = self.logRange(range_skew)
            range_amp = self.logRange(range_amp)
            range_LLH = self.logRange(range_LLH)
            label_time ='log time difference'
            label_width ='log width difference'
            label_skew ='log skewness difference'
            label_amp ='log amplitude ratio'
            label_LLH ='log 2DelLLH'
        else:
            if log_time == True:
                tDiff_t, tDiff_e = self.setlog(t_t, t_e)
                range_time = self.logRange(range_time)
                label_time ='log time difference'
            else:
                tDiff_t, tDiff_e = t_t, t_e
                label_time ='time difference[ns]'

            if log_wid == True:
                widDiff_t, widDiff_e = self.setlog(w_t, w_e)
                range_wid = self.logRange(range_wid)
                label_width ='log width difference'
            else:
                widDiff_t, widDiff_e = w_t, w_e
                label_width ='width difference[ns]'

            if log_skew == True:
                skewDiff_t, skewDiff_e = self.setlog(s_t, s_e)
                range_skew = self.logRange(range_skew)
                label_skew ='log skewness difference'
            else:
                skewDiff_t, skewDiff_e = s_t, s_e
                label_skew ='skewness difference'

            if log_amp == True:
                ampRatio_t, ampRatio_e = self.setlog(a_t, a_e)
                range_amp = self.logRange(range_amp)
                label_amp ='log amplitude ratio'
            else:
                ampRatio_t, ampRatio_e = a_t, a_e
                label_amp ='amplitude ratio'

            if log_LLH == True:
                LLHDiff_t, LLHDiff_e = self.setlog(l_t, l_e)
                range_LLH = self.logRange(range_LLH)
                label_LLH ='log 2DelLLH'
            else:
                LLHDiff_t, LLHDiff_e = l_t, l_e
                label_LLH ='2DelLLH'

        '''
        Setting bin range
        '''

        bins_time  = self.bin_range(log_time, tDiff_t, tDiff_e, numBins_time)
        bins_width  = self.bin_range(log_wid, widDiff_t, widDiff_e, numBins_wid)
        bins_skew  = self.bin_range(log_skew, skewDiff_t, skewDiff_e, numBins_skew)
        bins_amp  = self.bin_range(log_amp, ampRatio_t, ampRatio_e, numBins_amp)
        bins_LLH  = self.bin_range(log_LLH, LLHDiff_t, LLHDiff_e, numBins_LLH)

        print(bins_time.shape, bins_width.shape, bins_skew.shape, bins_amp.shape, bins_LLH.shape)
        print(tDiff_t.shape, LLHDiff_t.shape)

        '''
        plot
        '''

        cmap_tau = plt.cm.OrRd
        cmap_e = plt.cm.BuPu
        #TimeDiff, LLH
        fig, axs = plt.subplots(1, 2)
        fig.set_figheight(8)
        fig.set_figwidth(20)

        a0 = axs[0].hist2d(tDiff_t, LLHDiff_t, cmap=cmap_tau, norm=LogNorm(), normed=True, bins=(bins_time,bins_LLH))
        fig.colorbar(a0[3], ax=axs[0])
        axs[0].set_xlabel(label_time, fontsize = 18)
        axs[0].set_ylabel(label_LLH, fontsize = 18)
        axs[0].set_xlim(range_time)
        axs[0].set_ylim(range_LLH)

        a1 = axs[1].hist2d(tDiff_e, LLHDiff_e, cmap=cmap_e, norm=LogNorm(), normed=True, bins=(bins_time,bins_LLH))
        fig.colorbar(a1[3], ax=axs[1])
        axs[1].set_xlabel(label_time, fontsize = 18)
        axs[1].set_ylabel(label_LLH, fontsize = 18)
        axs[1].set_xlim(range_time)
        axs[1].set_ylim(range_LLH)
        plt.show()

        #TimeDiff, Skew
        fig, axs = plt.subplots(1, 2)
        fig.set_figheight(8)
        fig.set_figwidth(20)

        a0 = axs[0].hist2d(tDiff_t, skewDiff_t, cmap=cmap_tau, norm=LogNorm(), normed=True, bins=(bins_time,bins_skew))
        fig.colorbar(a0[3], ax=axs[0])
        axs[0].set_xlabel(label_time, fontsize = 18)
        axs[0].set_ylabel(label_skew, fontsize = 18)
        axs[0].set_xlim(range_time)
        axs[0].set_ylim(range_skew)

        a1 = axs[1].hist2d(tDiff_e, skewDiff_e, cmap=cmap_e, norm=LogNorm(), normed=True, bins=(bins_time,bins_skew))
        fig.colorbar(a1[3], ax=axs[1])
        axs[1].set_xlabel(label_time, fontsize = 18)
        axs[1].set_ylabel(label_skew, fontsize = 18)
        axs[1].set_xlim(range_time)
        axs[1].set_ylim(range_skew)
        plt.show()

        #time, amp
        fig, axs = plt.subplots(1, 2)
        fig.set_figheight(8)
        fig.set_figwidth(20)

        a0 = axs[0].hist2d(tDiff_t, ampRatio_t, cmap=cmap_tau, norm=LogNorm(), normed=True, bins=(bins_time,bins_amp))
        fig.colorbar(a0[3], ax=axs[0])
        axs[0].set_xlabel(label_time, fontsize = 18)
        axs[0].set_ylabel(label_amp, fontsize = 18)
        axs[0].set_xlim(range_time)
        axs[0].set_ylim(range_amp)

        a1 = axs[1].hist2d(tDiff_e, ampRatio_e, cmap=cmap_e, norm=LogNorm(), normed=True, bins=(bins_time,bins_amp))
        fig.colorbar(a1[3], ax=axs[1])
        axs[1].set_xlabel(label_time, fontsize = 18)
        axs[1].set_ylabel(label_amp, fontsize = 18)
        axs[1].set_xlim(range_time)
        axs[1].set_ylim(range_amp)
        plt.show()

        #time, wid
        fig, axs = plt.subplots(1, 2)
        fig.set_figheight(8)
        fig.set_figwidth(20)

        a0 = axs[0].hist2d(tDiff_t, widDiff_t, cmap=cmap_tau, norm=LogNorm(), normed=True, bins=(bins_time,bins_width))
        fig.colorbar(a0[3], ax=axs[0])
        axs[0].set_xlabel(label_time, fontsize = 18)
        axs[0].set_ylabel(label_width, fontsize = 18)
        axs[0].set_xlim(range_time)
        axs[0].set_ylim(range_wid)

        a1 = axs[1].hist2d(tDiff_e, widDiff_e, cmap=cmap_e, norm=LogNorm(), normed=True, bins=(bins_time,bins_width))
        fig.colorbar(a1[3], ax=axs[1])
        axs[1].set_xlabel(label_time, fontsize = 18)
        axs[1].set_ylabel(label_width, fontsize = 18)
        axs[1].set_xlim(range_time)
        axs[1].set_ylim(range_wid)
        plt.show()

        #Wid, LLH
        fig, axs = plt.subplots(1, 2)
        fig.set_figheight(8)
        fig.set_figwidth(20)

        a0 = axs[0].hist2d(widDiff_t, LLHDiff_t, cmap=cmap_tau, norm=LogNorm(), normed=True, bins=(bins_width,bins_LLH))
        fig.colorbar(a0[3], ax=axs[0])
        axs[0].set_xlabel(label_width, fontsize = 18)
        axs[0].set_ylabel(label_LLH, fontsize = 18)
        axs[0].set_xlim(range_wid)
        axs[0].set_ylim(range_LLH)

        a1 = axs[1].hist2d(widDiff_e, LLHDiff_e, cmap=cmap_e, norm=LogNorm(), normed=True, bins=(bins_width,bins_LLH))
        fig.colorbar(a1[3], ax=axs[1])
        axs[1].set_xlabel(label_width, fontsize = 18)
        axs[1].set_ylabel(label_LLH, fontsize = 18)
        axs[1].set_xlim(range_wid)
        axs[1].set_ylim(range_LLH)
        plt.show()

        #wid, skew
        fig, axs = plt.subplots(1, 2)
        fig.set_figheight(8)
        fig.set_figwidth(20)

        a0 = axs[0].hist2d(widDiff_t, skewDiff_t, cmap=cmap_tau, norm=LogNorm(), normed=True, bins=(bins_width,bins_skew))
        fig.colorbar(a0[3], ax=axs[0])
        axs[0].set_xlabel(label_width, fontsize = 18)
        axs[0].set_ylabel(label_skew, fontsize = 18)
        axs[0].set_xlim(range_wid)
        axs[0].set_ylim(range_skew)

        a1 = axs[1].hist2d(widDiff_e, skewDiff_e, cmap=cmap_e, norm=LogNorm(), normed=True, bins=(bins_width,bins_skew))
        fig.colorbar(a1[3], ax=axs[1])
        axs[1].set_xlabel(label_width, fontsize = 18)
        axs[1].set_ylabel(label_skew, fontsize = 18)
        axs[1].set_xlim(range_wid)
        axs[1].set_ylim(range_skew)
        plt.show()

        #wid,amp
        fig, axs = plt.subplots(1, 2)
        fig.set_figheight(8)
        fig.set_figwidth(20)

        a0 = axs[0].hist2d(widDiff_t, ampRatio_t, cmap=cmap_tau, norm=LogNorm(), normed=True, bins=(bins_width,bins_amp))
        fig.colorbar(a0[3], ax=axs[0])
        axs[0].set_xlabel(label_width, fontsize = 18)
        axs[0].set_ylabel(label_amp, fontsize = 18)
        axs[0].set_xlim(range_wid)
        axs[0].set_ylim(range_amp)

        a1 = axs[1].hist2d(widDiff_e, ampRatio_e, cmap=cmap_e, norm=LogNorm(), normed=True, bins=(bins_width,bins_amp))
        fig.colorbar(a1[3], ax=axs[1])
        axs[1].set_xlabel(label_width, fontsize = 18)
        axs[1].set_ylabel(label_amp, fontsize = 18)
        axs[1].set_xlim(range_wid)
        axs[1].set_ylim(range_amp)
        plt.show()

        #amp, LLH
        fig, axs = plt.subplots(1, 2)
        fig.set_figheight(8)
        fig.set_figwidth(20)

        a0 = axs[0].hist2d(ampRatio_t, LLHDiff_t, cmap=cmap_tau, norm=LogNorm(), normed=True, bins=(bins_amp,bins_LLH))
        fig.colorbar(a0[3], ax=axs[0])
        axs[0].set_xlabel(label_amp, fontsize = 18)
        axs[0].set_ylabel(label_LLH, fontsize = 18)
        axs[0].set_xlim(range_amp)
        axs[0].set_ylim(range_LLH)

        a1 = axs[1].hist2d(ampRatio_e, LLHDiff_e, cmap=cmap_e, norm=LogNorm(), normed=True, bins=(bins_amp,bins_LLH))
        fig.colorbar(a1[3], ax=axs[1])
        axs[1].set_xlabel(label_amp, fontsize = 18)
        axs[1].set_ylabel(label_LLH, fontsize = 18)
        axs[1].set_xlim(range_amp)
        axs[1].set_ylim(range_LLH)
        plt.show()

        #amp, skew
        fig, axs = plt.subplots(1, 2)
        fig.set_figheight(8)
        fig.set_figwidth(20)

        a0 = axs[0].hist2d(ampRatio_t, skewDiff_t, cmap=cmap_tau, norm=LogNorm(), normed=True, bins=(bins_amp,bins_skew))
        fig.colorbar(a0[3], ax=axs[0])
        axs[0].set_xlabel(label_amp, fontsize = 18)
        axs[0].set_ylabel(label_skew, fontsize = 18)
        axs[0].set_xlim(range_amp)
        axs[0].set_ylim(range_skew)

        a1 = axs[1].hist2d(ampRatio_e, skewDiff_e, cmap=cmap_e, norm=LogNorm(), normed=True, bins=(bins_amp,bins_skew))
        fig.colorbar(a1[3], ax=axs[1])
        axs[1].set_xlabel(label_amp, fontsize = 18)
        axs[1].set_ylabel(label_skew, fontsize = 18)
        axs[1].set_xlim(range_amp)
        axs[1].set_ylim(range_skew)
        plt.show()

        #skew, LLH
        fig, axs = plt.subplots(1, 2)
        fig.set_figheight(8)
        fig.set_figwidth(20)

        a0 = axs[0].hist2d(skewDiff_t, LLHDiff_t, cmap=cmap_tau, norm=LogNorm(), normed=True, bins=(bins_skew,bins_LLH))
        fig.colorbar(a0[3], ax=axs[0])
        axs[0].set_xlabel(label_skew, fontsize = 18)
        axs[0].set_ylabel(label_LLH, fontsize = 18)
        axs[0].set_xlim(range_skew)
        axs[0].set_ylim(range_LLH)

        a1 = axs[1].hist2d(skewDiff_e, LLHDiff_e, cmap=cmap_e, norm=LogNorm(), normed=True, bins=(bins_skew,bins_LLH))
        fig.colorbar(a1[3], ax=axs[1])
        axs[1].set_xlabel(label_skew, fontsize = 18)
        axs[1].set_ylabel(label_LLH, fontsize = 18)
        axs[1].set_xlim(range_skew)
        axs[1].set_ylim(range_LLH)
        plt.show()

    def corner(self, log_time=False, log_wid=False, log_skew=False, log_amp=False, log_LLH=False, log_all=False,
                range_time=[1e-12, 200], range_wid=[1e-12, 100], range_skew=[1e-12, 10], range_amp=[1e-12, 1e12], range_LLH=[1e-12, 1e6],
                numBins_time=10, numBins_wid=10, numBins_skew=10, numBins_amp=100, numBins_LLH=100):

        '''
        Setting Range
        '''
        t_t, t_e = self.tDiff_t, self.tDiff_e
        w_t, w_e = self.widDiff_t, self.widDiff_e
        s_t, s_e = self.skewDiff_t, self.skewDiff_e
        a_t, a_e = self.ampRatio_t, self.ampRatio_e
        l_t, l_e = self.LLHDiff_t, self.LLHDiff_e

        bool_time_t = (t_t>=range_time[0])&(t_t<=range_time[1])
        bool_time_e = (t_e>=range_time[0])&(t_e<=range_time[1])
        bool_wid_t = (w_t>=range_wid[0])&(w_t<=range_wid[1])
        bool_wid_e = (w_e>=range_wid[0])&(w_e<=range_wid[1])
        bool_skew_t = (s_t>=range_skew[0])&(s_t<=range_skew[1])
        bool_skew_e = (s_e>=range_skew[0])&(s_e<=range_skew[1])
        bool_amp_t = (a_t>=range_amp[0])&(a_t<=range_amp[1])
        bool_amp_e = (a_e>=range_amp[0])&(a_e<=range_amp[1])
        bool_LLH_t = (l_t>=range_LLH[0])&(l_t<=range_LLH[1])
        bool_LLH_e = (l_e>=range_LLH[0])&(l_e<=range_LLH[1])

        bool_tau = bool_time_t&bool_wid_t&bool_skew_t&bool_amp_t&bool_LLH_t
        bool_e = bool_time_e&bool_wid_e&bool_skew_e&bool_amp_e&bool_LLH_e

        tau = [t_t, w_t, a_t, s_t, l_t]
        e_nc = [t_e, w_e, a_e, s_e, l_e]

        '''
        Changing to log scale
        '''
        if log_time == True:
            tDiff_t, tDiff_e = self.setlog(t_t, t_e)
            range_time = self.logRange(range_time)
            label_time ='log time difference'
        else:
            tDiff_t, tDiff_e = t_t, t_e
            label_time ='time difference'

        if log_wid == True:
            widDiff_t, widDiff_e = self.setlog(w_t, w_e)
            range_wid = self.logRange(range_wid)
            label_width ='log width difference'
        else:
            widDiff_t, widDiff_e = w_t, w_e
            label_width ='width difference'

        if log_skew == True:
            skewDiff_t, skewDiff_e = self.setlog(s_t, s_e)
            range_skew = self.logRange(range_skew)
            label_skew ='log skewness difference'
        else:
            skewDiff_t, skewDiff_e = s_t, s_e
            label_skew ='skewness difference'

        if log_amp == True:
            ampRatio_t, ampRatio_e = self.setlog(a_t, a_e)
            range_amp = self.logRange(range_amp)
            label_amp ='log amplitude ratio'
        else:
            ampRatio_t, ampRatio_e = a_t, a_e
            label_amp ='amplitude ratio'

        if log_LLH == True:
            LLHDiff_t, LLHDiff_e = self.setlog(l_t, l_e)
            range_LLH = self.logRange(range_LLH)
            label_LLH ='log 2DelLLH'
        else:
            LLHDiff_t, LLHDiff_e = l_t, l_e
            label_LLH ='2DelLLH'

        print(len(tDiff_t), len(widDiff_t), len(skewDiff_t), len(ampRatio_t), len(LLHDiff_t))
        print(len(tDiff_e), len(widDiff_e), len(skewDiff_e), len(ampRatio_e), len(LLHDiff_e))


        if log_all == True:
            data_tau = np.log10(np.vstack(tau).T)
            data_e = np.log10(np.vstack(e_nc).T)
            range_time = self.logRange(range_time)
            range_wid = self.logRange(range_wid)
            range_skew = self.logRange(range_skew)
            range_amp = self.logRange(range_amp)
            range_LLH = self.logRange(range_LLH)
            label_time ='log time difference'
            label_width ='log width difference'
            label_skew ='log skewness difference'
            label_amp ='log amplitude ratio'
            label_LLH ='log 2DelLLH'
        else:
            data_tau = np.vstack([tDiff_t, widDiff_t, ampRatio_t, skewDiff_t, LLHDiff_t]).T
            data_e = np.vstack([tDiff_e, widDiff_e, ampRatio_e, skewDiff_e, LLHDiff_e]).T

        #labels = ['Time Difference[ns]', 'Width Difference[ns]', 'Amplitude Ratio', 'Skewness Difference', '2DelLLH']
        labels = [label_time, label_width, label_amp, label_skew, label_LLH]
        range_cust = ((range_time[0], range_time[1]), (range_wid[0], range_wid[1]), (range_amp[0], range_amp[1]), (range_skew[0], range_skew[1]), (range_LLH[0], range_LLH[1]))
        numBins = (numBins_time, numBins_wid, numBins_amp, numBins_skew, numBins_LLH)

        figure = corner.corner(data_tau, plot_contours=False, scale_hist=False, labels=labels,
                        plot_density = False, color='r', hist_kwargs={"normed":True,
                                                                     "log":True}, label='tau', range=range_cust, bins=numBins)

        corner.corner(data_e, fig=figure, plot_contours=False, scale_hist=False, labels=labels,
                        plot_density = False, color='b', hist_kwargs={"normed":True,
                                                                     "log":True}, label='e&NC', range=range_cust, bins=numBins)

        red_line = mlines.Line2D([], [], color='red', label='tau')
        blue_line = mlines.Line2D([], [], color='blue', label='e&NC')

        plt.legend(handles=[red_line,blue_line], bbox_to_anchor=(0., 1.0, 1., .0), loc=4)
        plt.show()

def plot_condition(condition, vals_single, vals, recoPulse_timeList, recoPulse_chargeList, tDiff, file_num, frame_num, omkey, lRatio):
    if condition:
        print(lRatio)

        '''
        Calculating the mean and removing the tails
        '''

        #mean = recoPulse_timeList.mean()
        mean = sum(recoPulse_timeList*recoPulse_chargeList)/sum(recoPulse_chargeList) #mean is weighted
        select_time = recoPulse_timeList[(recoPulse_timeList > mean-50) & (recoPulse_timeList < mean+50)]
        select_charge = recoPulse_chargeList[(recoPulse_timeList > mean-50) & (recoPulse_timeList < mean+50)]
        #print('SELECT CHARGE', select_charge, select_time, mean, recoPulse_timeList, recoPulse_chargeList)

        mean_select_time = sum(select_time*select_charge)/sum(select_charge)
        max_hitTimes = recoPulse_timeList[(recoPulse_timeList > (mean_select_time-100))&(recoPulse_timeList < (mean_select_time+100))]
        max_charge = recoPulse_chargeList[(recoPulse_timeList > (mean_select_time-100))&(recoPulse_timeList < (mean_select_time+100))]

        #Shifting mean to zero
        max_hitTimes_mean = sum(max_hitTimes*max_charge)/sum(max_charge)
        timestamps = max_hitTimes - max_hitTimes_mean
        final_mean = sum(timestamps*max_charge)/sum(max_charge)

        '''
        Histogramming the data from simulation
        '''
        print('Now Histogramming')
        bins = np.arange(min(timestamps), max(timestamps), 3)
        num, bin_edges = np.histogram(timestamps, bins=bins, weights=max_charge)
        bin_centers = (bin_edges[:-1]+bin_edges[1:])/2

        num_ampRatio = num/max(num)

        #removing bins which are <1/5 the max(num), removing the tails this way.
        num_select = num[num_ampRatio > 0.1]
        bin_centers_select = bin_centers[num_ampRatio > 0.1]

        '''
        Including continuity in the bins
        '''

        #considering two extra bins on both sides
        bin_center_bool = (bin_centers >= min(bin_centers_select) - 6)&(bin_centers <= max(bin_centers_select) + 6)
        entries_in_bins = num[bin_center_bool]
        bin_centers = bin_centers[bin_center_bool]

        maxBinCenter = max(bin_centers)

        '''
        (x, y) values for the fit
        '''
        #x = bin_centers
        x = np.linspace(min(bin_centers)-30, max(bin_centers)+30, 1000)
        #y_biGauss = biGauss(x, vals_single[1], vals_single[2], vals_single[3], vals_single[4])
        #y_doublePeak = double_peak(x, vals[1], vals[2], vals[3], vals[4],
        #                           vals[5], vals[6], vals[7], vals[8])

        y_biGauss = expGauss(x, vals_single[1], vals_single[2], vals_single[3], vals_single[4])
        y_doublePeak = expDoublePeak(x, vals[1], vals[2], vals[3], vals[4],
                                   vals[5], vals[6], vals[7], vals[8])

        plt.figure(figsize=(10,9))
        _ = plt.hist(timestamps, bins=bins, weights=max_charge, histtype='step', linewidth = 5)
        plt.plot(bin_centers, entries_in_bins, '*', c='k', label = 'Bins for fit', markersize=12, linewidth=6)
        plt.plot(x, y_biGauss, '--', c = 'k', label = 'expGauss', linewidth=3)
        plt.plot(x, y_doublePeak, '--', c = 'r', label = 'double expGauss', linewidth=3)
        plt.axvline(x=vals[1], c='orange', label='Postion of First Gaussian')
        plt.axvline(x=vals[5], c='g', label='Postion of Second Gaussian')
        plt.axhline(y=vals[4], c='orange', label='Amplitude of First Gaussian')
        plt.axhline(y=vals[8], c='g', label='Amplitude of Second gaussian')
        plt.legend()
        plt.xlabel('Time(ns)', fontsize = 16)
        plt.title('Time Difference - ' +str(round(tDiff)) + ' 2DelLLH - ' +str(lRatio) + ' File num - ' + str(file_num) + ' Frame - ' + str(frame_num) + ' ' + str(omkey), fontsize=14)
