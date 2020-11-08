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

        '''
        Setting bin range
        '''

        bins_time  = self.bin_range(numBins_time, tDiff_t, tDiff_e)
        bins_width  = self.bin_range(numBins_wid, widDiff_t, widDiff_e)
        bins_skew  = self.bin_range(numBins_skew, skewDiff_t, skewDiff_e)
        bins_amp  = self.bin_range(numBins_amp, ampRatio_t, ampRatio_e)
        bins_LLH  = self.bin_range(numBins_LLH, LLHDiff_t, LLHDiff_e)

        print(bins_time.shape, bins_width.shape, bins_skew.shape, bins_amp.shape, bins_LLH.shape)
        print(tDiff_t.shape, LLHDiff_t.shape)

        '''
        plot
        '''

        #TimeDiff, LLH
        fig, axs = plt.subplots(1, 2)
        fig.set_figheight(8)
        fig.set_figwidth(20)

        a0 = axs[0].hist2d(tDiff_t, LLHDiff_t, cmap=plt.cm.Reds, norm=LogNorm(), normed=True, bins=(bins_time,bins_LLH))
        fig.colorbar(a0[3], ax=axs[0])
        axs[0].set_xlabel('Time Difference[ns]', fontsize = 18)
        axs[0].set_ylabel('2DelLLH', fontsize = 18)
        axs[0].set_xlim(range_time)
        axs[0].set_ylim(range_LLH)

        a1 = axs[1].hist2d(tDiff_e, LLHDiff_e, cmap=plt.cm.Blues, norm=LogNorm(), normed=True, bins=(bins_time,bins_LLH))
        fig.colorbar(a1[3], ax=axs[1])
        axs[1].set_xlabel('Time Difference[ns]', fontsize = 18)
        axs[1].set_ylabel('2DelLLH', fontsize = 18)
        axs[1].set_xlim(range_time)
        axs[1].set_ylim(range_LLH)
        plt.show()

        #TimeDiff, Skew
        fig, axs = plt.subplots(1, 2)
        fig.set_figheight(8)
        fig.set_figwidth(20)

        a0 = axs[0].hist2d(tDiff_t, skewDiff_t, cmap=plt.cm.Reds, norm=LogNorm(), normed=True, bins=(bins_time,bins_skew))
        fig.colorbar(a0[3], ax=axs[0])
        axs[0].set_xlabel('Time Difference[ns]', fontsize = 18)
        axs[0].set_ylabel('Skewness Difference', fontsize = 18)
        axs[0].set_xlim(range_time)
        axs[0].set_ylim(range_skew)

        a1 = axs[1].hist2d(tDiff_e, skewDiff_e, cmap=plt.cm.Blues, norm=LogNorm(), normed=True, bins=(bins_time,bins_skew))
        fig.colorbar(a1[3], ax=axs[1])
        axs[1].set_xlabel('Time Difference[ns]', fontsize = 18)
        axs[1].set_ylabel('Skewness Difference', fontsize = 18)
        axs[1].set_xlim(range_time)
        axs[1].set_ylim(range_skew)
        plt.show()

        #time, amp
        fig, axs = plt.subplots(1, 2)
        fig.set_figheight(8)
        fig.set_figwidth(20)

        a0 = axs[0].hist2d(tDiff_t, ampRatio_t, cmap=plt.cm.Reds, norm=LogNorm(), normed=True, bins=(bins_time,bins_amp))
        fig.colorbar(a0[3], ax=axs[0])
        axs[0].set_xlabel('Time Difference[ns]', fontsize = 18)
        axs[0].set_ylabel('Amplitude Ratio', fontsize = 18)
        axs[0].set_xlim(range_time)
        axs[0].set_ylim(range_amp)

        a1 = axs[1].hist2d(tDiff_e, ampRatio_e, cmap=plt.cm.Blues, norm=LogNorm(), normed=True, bins=(bins_time,bins_amp))
        fig.colorbar(a1[3], ax=axs[1])
        axs[1].set_xlabel('Time Difference[ns]', fontsize = 18)
        axs[1].set_ylabel('Amplitude Ratio', fontsize = 18)
        axs[1].set_xlim(range_time)
        axs[1].set_ylim(range_amp)
        plt.show()

        #time, wid
        fig, axs = plt.subplots(1, 2)
        fig.set_figheight(8)
        fig.set_figwidth(20)

        a0 = axs[0].hist2d(tDiff_t, widDiff_t, cmap=plt.cm.Reds, norm=LogNorm(), normed=True, bins=(bins_time,bins_width))
        fig.colorbar(a0[3], ax=axs[0])
        axs[0].set_xlabel('Time Difference[ns]', fontsize = 18)
        axs[0].set_ylabel('Width Difference[ns]', fontsize = 18)
        axs[0].set_xlim(range_time)
        axs[0].set_ylim(range_wid)

        a1 = axs[1].hist2d(tDiff_e, widDiff_e, cmap=plt.cm.Blues, norm=LogNorm(), normed=True, bins=(bins_time,bins_width))
        fig.colorbar(a1[3], ax=axs[1])
        axs[1].set_xlabel('Time Difference[ns]', fontsize = 18)
        axs[1].set_ylabel('Width Difference[ns]', fontsize = 18)
        axs[1].set_xlim(range_time)
        axs[1].set_ylim(range_wid)
        plt.show()

        #Wid, LLH
        fig, axs = plt.subplots(1, 2)
        fig.set_figheight(8)
        fig.set_figwidth(20)

        a0 = axs[0].hist2d(widDiff_t, LLHDiff_t, cmap=plt.cm.Reds, norm=LogNorm(), normed=True, bins=(bins_width,bins_LLH))
        fig.colorbar(a0[3], ax=axs[0])
        axs[0].set_xlabel('Width Difference[ns]', fontsize = 18)
        axs[0].set_ylabel('2DelLLH', fontsize = 18)
        axs[0].set_xlim(range_wid)
        axs[0].set_ylim(range_LLH)

        a1 = axs[1].hist2d(widDiff_e, LLHDiff_e, cmap=plt.cm.Blues, norm=LogNorm(), normed=True, bins=(bins_width,bins_LLH))
        fig.colorbar(a1[3], ax=axs[1])
        axs[1].set_xlabel('Width Difference[ns]', fontsize = 18)
        axs[1].set_ylabel('2DelLLH', fontsize = 18)
        axs[1].set_xlim(range_wid)
        axs[1].set_ylim(range_LLH)
        plt.show()

        #wid, skew
        fig, axs = plt.subplots(1, 2)
        fig.set_figheight(8)
        fig.set_figwidth(20)

        a0 = axs[0].hist2d(widDiff_t, skewDiff_t, cmap=plt.cm.Reds, norm=LogNorm(), normed=True, bins=(bins_width,bins_skew))
        fig.colorbar(a0[3], ax=axs[0])
        axs[0].set_xlabel('Width Difference[ns]', fontsize = 18)
        axs[0].set_ylabel('Skewness Difference', fontsize = 18)
        axs[0].set_xlim(range_wid)
        axs[0].set_ylim(range_skew)

        a1 = axs[1].hist2d(widDiff_e, skewDiff_e, cmap=plt.cm.Blues, norm=LogNorm(), normed=True, bins=(bins_width,bins_skew))
        fig.colorbar(a1[3], ax=axs[1])
        axs[1].set_xlabel('Width Difference[ns]', fontsize = 18)
        axs[1].set_ylabel('Skewness Difference', fontsize = 18)
        axs[1].set_xlim(range_wid)
        axs[1].set_ylim(range_skew)
        plt.show()

        #wid,amp
        fig, axs = plt.subplots(1, 2)
        fig.set_figheight(8)
        fig.set_figwidth(20)

        a0 = axs[0].hist2d(widDiff_t, ampRatio_t, cmap=plt.cm.Reds, norm=LogNorm(), normed=True, bins=(bins_width,bins_amp))
        fig.colorbar(a0[3], ax=axs[0])
        axs[0].set_xlabel('Width Difference[ns]', fontsize = 18)
        axs[0].set_ylabel('Amplitude Ratio', fontsize = 18)
        axs[0].set_xlim(range_wid)
        axs[0].set_ylim(range_amp)

        a1 = axs[1].hist2d(widDiff_e, ampRatio_e, cmap=plt.cm.Blues, norm=LogNorm(), normed=True, bins=(bins_width,bins_amp))
        fig.colorbar(a1[3], ax=axs[1])
        axs[1].set_xlabel('Width Difference[ns]', fontsize = 18)
        axs[1].set_ylabel('Amplitude Ratio', fontsize = 18)
        axs[1].set_xlim(range_wid)
        axs[1].set_ylim(range_amp)
        plt.show()

        #amp, LLH
        fig, axs = plt.subplots(1, 2)
        fig.set_figheight(8)
        fig.set_figwidth(20)

        a0 = axs[0].hist2d(ampRatio_t, LLHDiff_t, cmap=plt.cm.Reds, norm=LogNorm(), normed=True, bins=(bins_amp,bins_LLH))
        fig.colorbar(a0[3], ax=axs[0])
        axs[0].set_xlabel('Amplitude Ratio', fontsize = 18)
        axs[0].set_ylabel('2DelLLH', fontsize = 18)
        axs[0].set_xlim(range_amp)
        axs[0].set_ylim(range_LLH)

        a1 = axs[1].hist2d(ampRatio_e, LLHDiff_e, cmap=plt.cm.Blues, norm=LogNorm(), normed=True, bins=(bins_amp,bins_LLH))
        fig.colorbar(a1[3], ax=axs[1])
        axs[1].set_xlabel('Amplitude Ratio', fontsize = 18)
        axs[1].set_ylabel('2DelLLH', fontsize = 18)
        axs[1].set_xlim(range_amp)
        axs[1].set_ylim(range_LLH)
        plt.show()

        #amp, skew
        fig, axs = plt.subplots(1, 2)
        fig.set_figheight(8)
        fig.set_figwidth(20)

        a0 = axs[0].hist2d(ampRatio_t, skewDiff_t, cmap=plt.cm.Reds, norm=LogNorm(), normed=True, bins=(bins_amp,bins_skew))
        fig.colorbar(a0[3], ax=axs[0])
        axs[0].set_xlabel('Amplitude Ratio', fontsize = 18)
        axs[0].set_ylabel('Skewness Difference', fontsize = 18)
        axs[0].set_xlim(range_amp)
        axs[0].set_ylim(range_skew)

        a1 = axs[1].hist2d(ampRatio_e, skewDiff_e, cmap=plt.cm.Blues, norm=LogNorm(), normed=True, bins=(bins_amp,bins_skew))
        fig.colorbar(a1[3], ax=axs[1])
        axs[1].set_xlabel('Amplitude Ratio', fontsize = 18)
        axs[1].set_ylabel('Skewness Difference', fontsize = 18)
        axs[1].set_xlim(range_amp)
        axs[1].set_ylim(range_skew)
        plt.show()

        #skew, LLH
        fig, axs = plt.subplots(1, 2)
        fig.set_figheight(8)
        fig.set_figwidth(20)

        a0 = axs[0].hist2d(skewDiff_t, LLHDiff_t, cmap=plt.cm.Reds, norm=LogNorm(), normed=True, bins=(bins_skew,bins_LLH))
        fig.colorbar(a0[3], ax=axs[0])
        axs[0].set_xlabel('Skewness Difference', fontsize = 18)
        axs[0].set_ylabel('2DelLLH', fontsize = 18)
        axs[0].set_xlim(range_skew)
        axs[0].set_ylim(range_LLH)

        a1 = axs[1].hist2d(skewDiff_e, LLHDiff_e, cmap=plt.cm.Blues, norm=LogNorm(), normed=True, bins=(bins_skew,bins_LLH))
        fig.colorbar(a1[3], ax=axs[1])
        axs[1].set_xlabel('Skewness Difference', fontsize = 18)
        axs[1].set_ylabel('2DelLLH', fontsize = 18)
        axs[1].set_xlim(range_skew)
        axs[1].set_ylim(range_LLH)
        plt.show()
