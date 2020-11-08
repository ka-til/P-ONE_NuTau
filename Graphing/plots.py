import corner
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm
import matplotlib.lines as mlines
from likelihoodHelpers import log_likelihood_biGauss, log_likelihood_doublePeak, likelihood_ratio_doublePeak, likelihood_ratio_biGauss, biGauss, double_peak

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

    def bin_range(self, log, t, e):
        if min(t) > min(e):
            min_val = min(e)
        else:
            min_val = min(t)

        if max(t) > max(e):
            max_val = max(t)
        else:
            max_val = max(e)

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

            bins = np.linspace(min_val, max_val+add, step)

        return bins

    def scatter(self, log_time=False, log_wid=False, log_skew=False, log_amp=False, log_LLH=False,
                range_time=[0, 200], range_wid=[0, 100], range_skew=[0, 10], range_amp=[0, 1e12], range_LLH=[0, 1e6]):

        '''
        Setting Range
        '''

        t_t, t_e = self.setRange(self.tDiff_t, self.tDiff_e, range_time[0], range_time[1])
        w_t, w_e = self.setRange(self.widDiff_t, self.widDiff_e, range_wid[0], range_wid[1])
        s_t, s_e = self.setRange(self.skewDiff_t, self.skewDiff_e, range_skew[0], range_skew[1])
        a_t, a_e = self.setRange(self.ampRatio_t, self.ampRatio_e, range_amp[0], range_amp[1])
        l_t, l_e = self.setRange(self.LLHDiff_t, self.LLHDiff_e, range_LLH[0], range_LLH[1])

        '''
        Changing to log scale
        '''
        if log_time == True:
            tDiff_t, tDiff_e = self.setlog(t_t, t_e)
        else:
            tDiff_t, tDiff_e = t_t, t_e

        if log_wid == True:
            widDiff_t, widDiff_e = self.setlog(w_t, w_e)
        else:
            widDiff_t, widDiff_e = w_t, w_e

        if log_skew == True:
            skewDiff_t, skewDiff_e = self.setlog(s_t, s_e)
        else:
            skewDiff_t, skewDiff_e = s_t, s_e

        if log_amp == True:
            ampRatio_t, ampRatio_e = self.setlog(a_t, a_e)
        else:
            ampRatio_t, ampRatio_e = a_t, a_e

        if log_LLH == True:
            LLHDiff_t, LLHDiff_e = self.setlog(l_t, l_e)
        else:
            LLHDiff_t, LLHDiff_e = l_t, l_e

        plt.figure(figsize=(10,9))
        plt.scatter(tDiff_t, LLHDiff_t, label='tau', c ='r', alpha = 0.6)
        plt.scatter(tDiff_e, LLHDiff_e, label='e&NC', c='skyblue', alpha = 0.6)
        plt.legend(fontsize=16)
        plt.xlabel('Time Difference[ns]', fontsize = 18)
        plt.ylabel('2DelLLH', fontsize = 18)

        plt.figure(figsize=(10,9))
        plt.scatter(tDiff_t, skewDiff_t, label='tau', c ='r', alpha = 0.6)
        plt.scatter(tDiff_e, skewDiff_e, label='e&NC', c='skyblue', alpha = 0.6)
        plt.legend(fontsize=16)
        plt.xlabel('Time Difference[ns]', fontsize = 18)
        plt.ylabel('Skewness Difference', fontsize = 18)

        plt.figure(figsize=(10,9))
        plt.scatter(tDiff_t, ampRatio_t, label='tau', c ='r', alpha = 0.6)
        plt.scatter(tDiff_e, ampRatio_e, label='e&NC', c='skyblue', alpha = 0.6)
        plt.legend(fontsize=16)
        plt.xlabel('Time Difference[ns]', fontsize = 18)
        plt.ylabel('Amplitude Ratio', fontsize = 18)

        plt.figure(figsize=(10,9))
        plt.scatter(tDiff_t, widDiff_t, label='tau', c ='r', alpha = 0.6)
        plt.scatter(tDiff_e, widDiff_e, label='e&NC', c='skyblue', alpha = 0.6)
        plt.legend(fontsize=16)
        plt.xlabel('Time Difference[ns]', fontsize = 18)
        plt.ylabel('Width Difference[ns]', fontsize = 18)

        plt.figure(figsize=(10,9))
        plt.scatter(widDiff_t, LLHDiff_t, label='tau', c ='r', alpha = 0.6)
        plt.scatter(widDiff_e, LLHDiff_e, label='e&NC', c='skyblue', alpha = 0.6)
        plt.legend(fontsize=16)
        plt.xlabel('Width Difference[ns]', fontsize = 18)
        plt.ylabel('2DelLLH', fontsize = 18)

        plt.figure(figsize=(10,9))
        plt.scatter(widDiff_t, skewDiff_t, label='tau', c ='r', alpha = 0.6)
        plt.scatter(widDiff_e, skewDiff_e, label='e&NC', c='skyblue', alpha = 0.6)
        plt.legend(fontsize=16)
        plt.xlabel('Width Difference[ns]', fontsize = 18)
        plt.ylabel('Skewness Difference', fontsize = 18)

        plt.figure(figsize=(10,9))
        plt.scatter(widDiff_t, ampRatio_t, label='tau', c ='r', alpha = 0.6)
        plt.scatter(widDiff_e, ampRatio_e, label='e&NC', c='skyblue', alpha = 0.6)
        plt.legend(fontsize=16)
        plt.xlabel('Width Difference[ns]', fontsize = 18)
        plt.ylabel('Amplitude Ratio', fontsize = 18)

        plt.figure(figsize=(10,9))
        plt.scatter(ampRatio_t, LLHDiff_t, label='tau', c ='r', alpha = 0.6)
        plt.scatter(ampRatio_e, LLHDiff_e, label='e&NC', c='skyblue', alpha = 0.6)
        plt.legend(fontsize=16)
        plt.xlabel('Amplitude Ratio', fontsize = 18)
        plt.ylabel('2DelLLH', fontsize = 18)

        plt.figure(figsize=(10,9))
        plt.scatter(ampRatio_t, skewDiff_t, label='tau', c ='r', alpha = 0.6)
        plt.scatter(ampRatio_e, skewDiff_e, label='e&NC', c='skyblue', alpha = 0.6)
        plt.legend(fontsize=16)
        plt.xlabel('Amplitude Ratio', fontsize = 18)
        plt.ylabel('Skewness Difference', fontsize = 18)

        plt.figure(figsize=(10,9))
        plt.scatter(skewDiff_t, LLHDiff_t, label='tau', c ='r', alpha = 0.6)
        plt.scatter(skewDiff_e, LLHDiff_e, label='e&NC', c='skyblue', alpha = 0.6)
        plt.legend(fontsize=16)
        plt.xlabel('Skewness Difference', fontsize = 18)
        plt.ylabel('2DelLLH', fontsize = 18)


    def hist2d(self, log_time=False, log_wid=False, log_skew=False, log_amp=False, log_LLH=False,
                    range_time=[0, 200], range_wid=[0, 100], range_skew=[0, 10], range_amp=[0, 1e12], range_LLH=[0, 1e6],
                    numBins_time=10, numBins_wid=10, numBins_skew=10, numBins_amp=100, numBins_LLH=100):

        print(range_time, range_wid, range_skew, range_amp, range_LLH)

        '''
        Setting Range
        '''
        t_t, t_e = self.setRange(self.tDiff_t, self.tDiff_e, range_time[0], range_time[1])
        w_t, w_e = self.setRange(self.widDiff_t, self.widDiff_e, range_wid[0], range_wid[1])
        s_t, s_e = self.setRange(self.skewDiff_t, self.skewDiff_e, range_skew[0], range_skew[1])
        a_t, a_e = self.setRange(self.ampRatio_t, self.ampRatio_e, range_amp[0], range_amp[1])
        l_t, l_e = self.setRange(self.LLHDiff_t, self.LLHDiff_e, range_LLH[0], range_LLH[1])

        '''
        Changing to log scale
        '''
        if log_time == True:
            tDiff_t, tDiff_e = self.setlog(t_t, t_e)
            range_time = self.logRange(range_time)
        else:
            tDiff_t, tDiff_e = t_t, t_e

        if log_wid == True:
            widDiff_t, widDiff_e = self.setlog(w_t, w_e)
            range_width = self.logRange(range_width)
        else:
            widDiff_t, widDiff_e = w_t, w_e

        if log_skew == True:
            skewDiff_t, skewDiff_e = self.setlog(s_t, s_e)
            range_skew = self.logRange(range_skew)
        else:
            skewDiff_t, skewDiff_e = s_t, s_e

        if log_amp == True:
            ampRatio_t, ampRatio_e = self.setlog(a_t, a_e)
            range_amp = self.logRange(range_amp)
        else:
            ampRatio_t, ampRatio_e = a_t, a_e

        if log_LLH == True:
            LLHDiff_t, LLHDiff_e = self.setlog(l_t, l_e)
            range_LLH = self.logRange(range_LLH)
        else:
            LLHDiff_t, LLHDiff_e = l_t, l_e

        '''
        Setting bin range
        '''

        bins_time  = self.bin_range(log_time, tDiff_t, tDiff_e)
        bins_width  = self.bin_range(log_wid, widDiff_t, widDiff_e)
        bins_skew  = self.bin_range(log_skew, skewDiff_t, skewDiff_e)
        bins_amp  = self.bin_range(log_amp, ampRatio_t, ampRatio_e)
        bins_LLH  = self.bin_range(log_LLH, LLHDiff_t, LLHDiff_e)

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
        axs[0].set_xlabel('Time Difference[ns]', fontsize = 18)
        axs[0].set_ylabel('2DelLLH', fontsize = 18)

        a1 = axs[1].hist2d(tDiff_e, LLHDiff_e, cmap=cmap_e, norm=LogNorm(), normed=True, bins=(bins_time,bins_LLH))
        fig.colorbar(a1[3], ax=axs[1])
        axs[1].set_xlabel('Time Difference[ns]', fontsize = 18)
        axs[1].set_ylabel('2DelLLH', fontsize = 18)
        plt.show()

        #TimeDiff, Skew
        fig, axs = plt.subplots(1, 2)
        fig.set_figheight(8)
        fig.set_figwidth(20)

        a0 = axs[0].hist2d(tDiff_t, skewDiff_t, cmap=cmap_tau, norm=LogNorm(), normed=True, bins=(bins_time,bins_skew))
        fig.colorbar(a0[3], ax=axs[0])
        axs[0].set_xlabel('Time Difference[ns]', fontsize = 18)
        axs[0].set_ylabel('Skewness Difference', fontsize = 18)

        a1 = axs[1].hist2d(tDiff_e, skewDiff_e, cmap=cmap_e, norm=LogNorm(), normed=True, bins=(bins_time,bins_skew))
        fig.colorbar(a1[3], ax=axs[1])
        axs[1].set_xlabel('Time Difference[ns]', fontsize = 18)
        axs[1].set_ylabel('Skewness Difference', fontsize = 18)
        plt.show()

        #time, amp
        fig, axs = plt.subplots(1, 2)
        fig.set_figheight(8)
        fig.set_figwidth(20)

        a0 = axs[0].hist2d(tDiff_t, ampRatio_t, cmap=cmap_tau, norm=LogNorm(), normed=True, bins=(bins_time,bins_amp))
        fig.colorbar(a0[3], ax=axs[0])
        axs[0].set_xlabel('Time Difference[ns]', fontsize = 18)
        axs[0].set_ylabel('Amplitude Ratio', fontsize = 18)

        a1 = axs[1].hist2d(tDiff_e, ampRatio_e, cmap=cmap_e, norm=LogNorm(), normed=True, bins=(bins_time,bins_amp))
        fig.colorbar(a1[3], ax=axs[1])
        axs[1].set_xlabel('Time Difference[ns]', fontsize = 18)
        axs[1].set_ylabel('Amplitude Ratio', fontsize = 18)
        plt.show()

        #time, wid
        fig, axs = plt.subplots(1, 2)
        fig.set_figheight(8)
        fig.set_figwidth(20)

        a0 = axs[0].hist2d(tDiff_t, widDiff_t, cmap=cmap_tau, norm=LogNorm(), normed=True, bins=(bins_time,bins_width))
        fig.colorbar(a0[3], ax=axs[0])
        axs[0].set_xlabel('Time Difference[ns]', fontsize = 18)
        axs[0].set_ylabel('Width Difference[ns]', fontsize = 18)

        a1 = axs[1].hist2d(tDiff_e, widDiff_e, cmap=cmap_e, norm=LogNorm(), normed=True, bins=(bins_time,bins_width))
        fig.colorbar(a1[3], ax=axs[1])
        axs[1].set_xlabel('Time Difference[ns]', fontsize = 18)
        axs[1].set_ylabel('Width Difference[ns]', fontsize = 18)
        plt.show()

        #Wid, LLH
        fig, axs = plt.subplots(1, 2)
        fig.set_figheight(8)
        fig.set_figwidth(20)

        a0 = axs[0].hist2d(widDiff_t, LLHDiff_t, cmap=cmap_tau, norm=LogNorm(), normed=True, bins=(bins_width,bins_LLH))
        fig.colorbar(a0[3], ax=axs[0])
        axs[0].set_xlabel('Width Difference[ns]', fontsize = 18)
        axs[0].set_ylabel('2DelLLH', fontsize = 18)

        a1 = axs[1].hist2d(widDiff_e, LLHDiff_e, cmap=cmap_e, norm=LogNorm(), normed=True, bins=(bins_width,bins_LLH))
        fig.colorbar(a1[3], ax=axs[1])
        axs[1].set_xlabel('Width Difference[ns]', fontsize = 18)
        axs[1].set_ylabel('2DelLLH', fontsize = 18)
        plt.show()

        #wid, skew
        fig, axs = plt.subplots(1, 2)
        fig.set_figheight(8)
        fig.set_figwidth(20)

        a0 = axs[0].hist2d(widDiff_t, skewDiff_t, cmap=cmap_tau, norm=LogNorm(), normed=True, bins=(bins_width,bins_skew))
        fig.colorbar(a0[3], ax=axs[0])
        axs[0].set_xlabel('Width Difference[ns]', fontsize = 18)
        axs[0].set_ylabel('Skewness Difference', fontsize = 18)

        a1 = axs[1].hist2d(widDiff_e, skewDiff_e, cmap=cmap_e, norm=LogNorm(), normed=True, bins=(bins_width,bins_skew))
        fig.colorbar(a1[3], ax=axs[1])
        axs[1].set_xlabel('Width Difference[ns]', fontsize = 18)
        axs[1].set_ylabel('Skewness Difference', fontsize = 18)
        plt.show()

        #wid,amp
        fig, axs = plt.subplots(1, 2)
        fig.set_figheight(8)
        fig.set_figwidth(20)

        a0 = axs[0].hist2d(widDiff_t, ampRatio_t, cmap=cmap_tau, norm=LogNorm(), normed=True, bins=(bins_width,bins_amp))
        fig.colorbar(a0[3], ax=axs[0])
        axs[0].set_xlabel('Width Difference[ns]', fontsize = 18)
        axs[0].set_ylabel('Amplitude Ratio', fontsize = 18)

        a1 = axs[1].hist2d(widDiff_e, ampRatio_e, cmap=cmap_e, norm=LogNorm(), normed=True, bins=(bins_width,bins_amp))
        fig.colorbar(a1[3], ax=axs[1])
        axs[1].set_xlabel('Width Difference[ns]', fontsize = 18)
        axs[1].set_ylabel('Amplitude Ratio', fontsize = 18)
        plt.show()

        #amp, LLH
        fig, axs = plt.subplots(1, 2)
        fig.set_figheight(8)
        fig.set_figwidth(20)

        a0 = axs[0].hist2d(ampRatio_t, LLHDiff_t, cmap=cmap_tau, norm=LogNorm(), normed=True, bins=(bins_amp,bins_LLH))
        fig.colorbar(a0[3], ax=axs[0])
        axs[0].set_xlabel('Amplitude Ratio', fontsize = 18)
        axs[0].set_ylabel('2DelLLH', fontsize = 18)

        a1 = axs[1].hist2d(ampRatio_e, LLHDiff_e, cmap=cmap_e, norm=LogNorm(), normed=True, bins=(bins_amp,bins_LLH))
        fig.colorbar(a1[3], ax=axs[1])
        axs[1].set_xlabel('Amplitude Ratio', fontsize = 18)
        axs[1].set_ylabel('2DelLLH', fontsize = 18)
        plt.show()

        #amp, skew
        fig, axs = plt.subplots(1, 2)
        fig.set_figheight(8)
        fig.set_figwidth(20)

        a0 = axs[0].hist2d(ampRatio_t, skewDiff_t, cmap=cmap_tau, norm=LogNorm(), normed=True, bins=(bins_amp,bins_skew))
        fig.colorbar(a0[3], ax=axs[0])
        axs[0].set_xlabel('Amplitude Ratio', fontsize = 18)
        axs[0].set_ylabel('Skewness Difference', fontsize = 18)

        a1 = axs[1].hist2d(ampRatio_e, skewDiff_e, cmap=cmap_e, norm=LogNorm(), normed=True, bins=(bins_amp,bins_skew))
        fig.colorbar(a1[3], ax=axs[1])
        axs[1].set_xlabel('Amplitude Ratio', fontsize = 18)
        axs[1].set_ylabel('Skewness Difference', fontsize = 18)
        plt.show()

        #skew, LLH
        fig, axs = plt.subplots(1, 2)
        fig.set_figheight(8)
        fig.set_figwidth(20)

        a0 = axs[0].hist2d(skewDiff_t, LLHDiff_t, cmap=cmap_tau, norm=LogNorm(), normed=True, bins=(bins_skew,bins_LLH))
        fig.colorbar(a0[3], ax=axs[0])
        axs[0].set_xlabel('Skewness Difference', fontsize = 18)
        axs[0].set_ylabel('2DelLLH', fontsize = 18)

        a1 = axs[1].hist2d(skewDiff_e, LLHDiff_e, cmap=cmap_e, norm=LogNorm(), normed=True, bins=(bins_skew,bins_LLH))
        fig.colorbar(a1[3], ax=axs[1])
        axs[1].set_xlabel('Skewness Difference', fontsize = 18)
        axs[1].set_ylabel('2DelLLH', fontsize = 18)
        plt.show()

    def corner(self, log_time=False, log_wid=False, log_skew=False, log_amp=False, log_LLH=False, log_all=False,
                range_time=[0, 200], range_wid=[0, 100], range_skew=[0, 10], range_amp=[0, 1e12], range_LLH=[0, 1e6]):

        '''
        Setting Range
        '''

        t_t, t_e = self.setRange(self.tDiff_t, self.tDiff_e, range_time[0], range_time[1])
        w_t, w_e = self.setRange(self.widDiff_t, self.widDiff_e, range_wid[0], range_wid[1])
        s_t, s_e = self.setRange(self.skewDiff_t, self.skewDiff_e, range_skew[0], range_skew[1])
        a_t, a_e = self.setRange(self.ampRatio_t, self.ampRatio_e, range_amp[0], range_amp[1])
        l_t, l_e = self.setRange(self.LLHDiff_t, self.LLHDiff_e, range_LLH[0], range_LLH[1])

        labels = ['Time Difference[ns]', 'Width Difference[ns]', 'Amplitude Ratio', 'Skewness Difference', '2DelLLH']
        tau = [t_t, w_t, a_t, s_t, l_t]
        e_nc = [t_e, w_e, a_e, s_e, l_e]

        '''
        Changing to log scale
        '''
        if log_time == True:
            tDiff_t, tDiff_e = self.setlog(t_t, t_e)
        else:
            tDiff_t, tDiff_e = t_t, t_e

        if log_wid == True:
            widDiff_t, widDiff_e = self.setlog(w_t, w_e)
        else:
            widDiff_t, widDiff_e = w_t, w_e

        if log_skew == True:
            skewDiff_t, skewDiff_e = self.setlog(s_t, s_e)
        else:
            skewDiff_t, skewDiff_e = s_t, s_e

        if log_amp == True:
            ampRatio_t, ampRatio_e = self.setlog(a_t, a_e)
        else:
            ampRatio_t, ampRatio_e = a_t, a_e

        if log_LLH == True:
            LLHDiff_t, LLHDiff_e = self.setlog(l_t, l_e)
        else:
            LLHDiff_t, LLHDiff_e = l_t, l_e

        if log_all == True:
            data_tau = np.log10(np.vstack(tau).T)
            data_e = np.log10(np.vstack(e_nc).T)
        else:
            data_tau = np.vstack([tDiff_t, widDiff_t, ampRatio_t, skewDiff_t, LLHDiff_t]).T
            data_e = np.vstack(tDiff_e, widDiff_e, ampRatio_e, skewDiff_e, LLHDiff_e).T


        figure = corner.corner(data_tau, plot_contours=False, scale_hist=False, labels=labels,
                        plot_density = False, color='r', hist_kwargs={"normed":True,
                                                                     "log":True}, label='tau')

        corner.corner(data_e, fig=figure, plot_contours=False, scale_hist=False, labels=labels,
                        plot_density = False, color='b', hist_kwargs={"normed":True,
                                                                     "log":True}, label='e&NC')

        red_line = mlines.Line2D([], [], color='red', label='tau')
        blue_line = mlines.Line2D([], [], color='blue', label='e&NC')

        plt.legend(handles=[red_line,blue_line], bbox_to_anchor=(0., 1.0, 1., .0), loc=4)
        plt.show()

def plot_condition(condition, vals_single, vals, recoPulse_timeList, recoPulse_chargeList, tDiff, file_num, frame_num, omkey):
    if condition:

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
        entries_in_bins = num[num_ampRatio > 0.2]
        bin_centers = bin_centers[num_ampRatio > 0.2]

        maxBinCenter = max(bin_centers)

        '''
        (x, y) values for the fit
        '''
        #x = bin_centers
        x = np.linspace(min(bin_centers)-30, max(bin_centers)+30, 1000)
        y_biGauss = biGauss(x, vals_single[1], vals_single[2], vals_single[3], vals_single[4])
        y_doublePeak = double_peak(x, vals[1], vals[2], vals[3], vals[4],
                                   vals[5], vals[6], vals[7], vals[8])

        plt.figure(figsize=(10,9))
        _ = plt.hist(timestamps, bins=bins, weights=max_charge, histtype='step', linewidth = 5)
        plt.plot(x, y_biGauss, '--', c = 'k', label = 'bigauss', linewidth=3)
        plt.plot(x, y_doublePeak, '--', c = 'r', label = 'double bigauss', linewidth=3)
        plt.axvline(x=vals[1], c='orange', label='Postion of First Gaussian')
        plt.axvline(x=vals[5], c='g', label='Postion of Second Gaussian')
        plt.axhline(y=vals[4], c='orange', label='Amplitude of First Gaussian')
        plt.axhline(y=vals[8], c='g', label='Amplitude of Second gaussian')
        plt.legend()
        plt.xlabel('Time(ns)', fontsize = 16)
        plt.title('Time Difference - ' +str(tDiff) + ' File num - ' + str(file_num) + ' Frame - ' + str(frame_num) + ' ' + str(omkey), fontsize=14)
