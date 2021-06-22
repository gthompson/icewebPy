from obspy.core import read, Stream
import numpy as np

## imports for copying obspy spectrograms
from obspy.imaging.spectrogram import _nearest_pow_2
from matplotlib import mlab
import matplotlib.pyplot as plt

# colormaps
# https://docs.obspy.org/packages/autogen/obspy.imaging.cm.html
# obspy_sequential is same as viridis, which is default for matplotlib
from obspy.imaging.cm import viridis_white_r, obspy_divergent, pqlx, obspy_sequential

# adding colorbar
from mpl_toolkits.axes_grid1 import make_axes_locatable

class icewebSpectrogram:
    
    def __init__(self, stream=None, secsPerFFT=-1):
        """
        :type st: Stream
        :param st: an ObsPy Stream object. No detrending or filtering is done, so do that before calling this function.
        """
        self.stream = Stream()
        self.F = []
        self.T = []
        self.S = []
        self.precomputed = False
        if isinstance(stream, Stream):
            self.stream = stream
       
    def __str__(self):
        str = '\n\nicewebSpectrogram:\n'
        str += self.stream.__str__()
        str += '\nF: %d 1-D numpy arrays' % len(self.F)
        str += '\nT: %d 1-D numpy arrays' % len(self.T)
        str += '\nS: %d 2-D numpy arrays\n\n' % len(self.S)
        return str
    
    def precompute(self, secsPerFFT=None):
        """
    
        For each Trace in self.stream, call compute_one.
    
        :type secsPerFFT: int or float
        :param secsPerFFT: Window length for fft in seconds. If this parameter is too
            small, the calculation will take forever. If None, it defaults to
            ceil(sampling_rate/100.0).
        """

        # seconds to use for each FFT. 1 second if the event duration is <= 100 seconds, 6 seconds if it is 10-minutes
        if secsPerFFT is None:
            secsPerFFT = np.ceil((self.stream[0].stats.delta * self.stream[0].stats.npts)/100)
            #print('seconds per FFT = %.1f' % secsPerFFT)   
        
        for tr in self.stream:
            #[T, F, S] = icewebSpectrogram.compute_one(tr, wlen=secsPerFFT)
            [T, F, S] = compute_spectrogram(tr, wlen=secsPerFFT)
            self.T.append(T)
            self.F.append(F)
            self.S.append(S)
        self.precomputed = True
    
        return self    

    
    def plot(self, outfile=None, secsPerFFT=None, fmin=0.5, fmax=20.0, log=False, cmap=pqlx, clim=None, \
                      equal_scale=False, title=None, add_colorbar=False, precompute=False, dbscale=True ):   
        """
        For each Trace in a Stream, plot the seismogram and spectrogram. This results in 2*N subplots 
        on the figure, where N is the number of Trace objects.
        This is modelled after IceWeb spectrograms, which have been part of the real-time monitoring 
        system at the Alaska Volcano Observatory since March 1998
        MATLAB code for this exists in the GISMO toolbox at 
        https://github.com/geoscience-community-codes/GISMO/blob/master/applications/+iceweb/spectrogram_iceweb.m
    
        :type outfile: str
        :param outfile: String for the filename of output file     
        :type fmin: float
        :param fmin: frequency minimum to plot on spectrograms
        :type fmax: float
        :param fmax: frequency maximum to plot on spectrograms    
        :type log: bool
        :param log: Logarithmic frequency axis if True, linear frequency axis
            otherwise.
        :type cmap: :class:`matplotlib.colors.Colormap`
        :param cmap: Specify a custom colormap instance. If not specified, then the
            pqlx colormap is used.
        :type clim: [float, float]
        :param clim: colormap limits. adjust colormap to clip at lower and/or upper end.
            This overrides equal_scale parameter.
        :type equal_scale: bool
        :param equal_scale: Apply the same colormap limits to each spectrogram if True. 
            This requires more memory since all spectrograms have to be pre-computed 
            to determine overall min and max spectral amplitude within [fmin, fmax]. 
            If False (default), each spectrogram is individually scaled, which is best 
            for seeing details. equal_scale is overridden by clim if given.
        :type add_colorbar: bool
        :param add_colorbar: Add colorbars for each spectrogram (5% space will be created on RHS)        
        :type title: str
        :param title: String for figure super title
        :type dbscale: bool
        :param dbscale: If True 20 * log10 of color values is taken.         
        """
        
        if not self.precomputed and precompute:
            self = self.precompute(secsPerFFT=secsPerFFT)

        N = len(self.stream) # number of channels we are plotting        
        if N==0:
            print('Stream object is empty. Nothing to do')
            return
         
        fig, ax = plt.subplots(N*2, 1) # create fig and ax handles with approx positions for now
        fig.set_size_inches(5.76, 7.56) 
 
        if equal_scale: # calculate range of spectral amplitudes
            if self.precomputed:
                Smin, Smax = icewebSpectrogram.S_range(self, fmin=fmin, fmax=fmax)
            else:
                index_min = np.argmin(self.stream.max()) # find the index of largest Trace object
                #[T, F, S] = icewebSpectrogram.compute_one(self.stream[index_min], wlen=secsPerFFT)
                [T, F, S] = compute_spectrogram(self.stream[index_min], wlen=secsPerFFT)
                f_indexes = np.intersect1d(np.where(F>=fmin), np.where(F<fmax))
                S_filtered = S[f_indexes, :]
                Smax = np.nanmax(S_filtered)
                S_filtered[S_filtered == 0] = Smax
                Smin = np.nanmin(S_filtered)
                
            if Smin<Smax*1e-6: # impose a dynamic range limit of 1,000,000
                Smin=Smax*1e-6
                
            clim = (Smin, Smax)
        
        c = 0 # initialize channel number
        for c in range(N):
            tr = self.stream[c]
            if self.precomputed:
                T = self.T[c]
                F = self.F[c]
                S = self.S[c]                
            else:                
                #[T, F, S] = icewebSpectrogram.compute_one(tr, wlen=secsPerFFT)
                [T, F, S] = compute_spectrogram(tr, wlen=secsPerFFT)
        
            # fix the axes positions for this trace and spectrogram
            if add_colorbar and clim:
                spectrogramPosition, tracePosition = icewebSpectrogram.calculateSubplotPositions(N, c, 
                                                                       frameBottom = 0.13, totalHeight = 0.83)
                cax = fig.add_axes([spectrogramPosition[0], 0.04, spectrogramPosition[2], 0.025])
            else:
                spectrogramPosition, tracePosition = icewebSpectrogram.calculateSubplotPositions(N, c)
            ax[c*2].set_position(tracePosition)
            ax[c*2+1].set_position(spectrogramPosition)
        
            # plot the trace
            t = tr.times()
            ax[c*2].plot(t, tr.data, linewidth=0.5);
            ax[c*2].set_yticks(ticks=[]) # turn off yticks
        
            if log:
                # Log scaling for frequency values (y-axis)
                ax[c*2+1].set_yscale('log')
           
            # Plot spectrogram
            vmin = None
            vmax = None
            if clim:
                if dbscale:
                    vmin, vmax = amp2dB(clim)
                else:
                    vmin, vmax = clim
            if dbscale:
                S = amp2dB(S)
                                
            sgram_handle = ax[c*2+1].pcolormesh(T, F, S, vmin=vmin, vmax=vmax, cmap=cmap )
            ax[c*2+1].set_ylim(fmin, fmax)
        
            # Plot colorbar
            if add_colorbar:
                if clim: # Scaled. Add a colorbar at the bottom of the figure. Just do it once.
                    if c==0:
                        plt.colorbar(sgram_handle, cax=cax, orientation='horizontal') 
                        if dbscale:
                            cax.set_xlabel('dB relative to 1 m/s/Hz')
                            ''' 
                            # Failed attempt to add a second set of labels to top of colorbar
                            # It just wipes out the original xtick labels
                            cax2 = cax.twinx()
                            cax2.xaxis.set_ticks_position("top")
                            xticks = cax.get_xticks()
                            xticks = np.power(xticks/20, 10)
                            cax2.set_xticks(xticks)
                            cax2.set_xlabel('m/s/Hz')
                            '''
                        else:
                            cax.set_xlabel('m/s/Hz')
                else: # Unscaled, so each spectrogram has max resolution. Add a colorbar to the right of each spectrogram.
                    divider = make_axes_locatable(ax[c*2+1])
                    cax = divider.append_axes("right", size="5%", pad=0.05)
                    plt.colorbar(sgram_handle, cax=cax)
                    # also add space next to Trace
                    divider2 = make_axes_locatable(ax[c*2])
                    hide_ax = divider2.append_axes("right", size="5%", pad=0.05, visible=False)
                    
      
            # add a ylabel
            ax[c*2+1].set_ylabel('     ' + tr.stats.station + '.' + tr.stats.channel, rotation=80)
            
            # turn off xticklabels, except for bottom panel
            ax[c*2].set_xticklabels([])
            if c<N-1:
                ax[c*2+1].set_xticklabels([])
            
            # increment c and go to next trace (if any left)
            c += 1   
    
    
        ax[N*2-1].set_xlabel('Time [s]')
    
        if title:
            ax[0].set_title(title)
        
        # set the xlimits for each panel from the min and max time values we kept updating
        [min_t, max_t] = icewebSpectrogram.t_range(self)
        [min_T, max_T] = icewebSpectrogram.T_range(self) 
        #print(min_t, max_t, min_T, max_T) 
        for c in range(N): 
            ax[c*2].set_xlim(min_t, max_t)
            ax[c*2].grid(axis='x', linestyle = ':', linewidth=0.5)
            #ax[c*2+1].set_xlim(min_T, max_T)
            ax[c*2+1].set_xlim(min_t, max_t)
            ax[c*2+1].grid(True, linestyle = ':', linewidth=0.5)    
        
        # write plot to file
        if outfile:
            fig.savefig(outfile, dpi=100)

        return fig, ax

    def t_range(self):
        min_t = np.Inf
        max_t = -np.Inf    
        for tr in self.stream:
            t = tr.times()
            min_t = min([np.nanmin(t), min_t])
            max_t = max([np.nanmax(t), max_t])
        return min_t, max_t  
    
    def T_range(self):   
        min_T = np.Inf
        max_T = -np.Inf 
        for c in range(len(self.T)):
            T = self.T[c]
            min_T = min([np.nanmin(T), min_T])
            max_T = max([np.nanmax(T), max_T])  
        return min_T, max_T

    def S_range(self, fmin=0.5, fmax=20.0):
        Smin = np.Inf
        Smax = -np.Inf
        for c in range(len(self.T)):
            F = self.F[c]
            S = self.S[c]
                            
            # filter S between fmin and fmax and then update Smin and Smax
            f_indexes = np.intersect1d(np.where(F>=fmin), np.where(F<fmax))
            #print(F[f_indexes[0]], F[f_indexes[-1]])
            try:
                S_filtered = S[f_indexes, :]
            except:
                print('S_range failed. F is ',F.shape, ' S is ',S.shape)
            else:
                Smin = min([np.nanmin(S_filtered), Smin])
                Smax = max([np.nanmax(S_filtered), Smax])
        print('S ranges from %e to %e' % (Smin, Smax))
        return Smin, Smax    
    
   
        
    def calculateSubplotPositions(numchannels, channelNum, frameLeft=0.12, frameBottom=0.08, \
                              totalWidth = 0.8, totalHeight = 0.88, fractionalSpectrogramHeight = 0.8):
        """ Copied from the MATLAB/GISMO function """
    
        channelHeight = totalHeight/numchannels;
        spectrogramHeight = fractionalSpectrogramHeight * channelHeight;
        traceHeight = channelHeight - spectrogramHeight; 
        spectrogramBottom = frameBottom + (numchannels - channelNum - 1) * channelHeight; # only change was to subtract 1 here, since in Python channelnum goes from 0..N-1, not 1..N as in MATLAB
        traceBottom = spectrogramBottom + spectrogramHeight;
        spectrogramPosition = [frameLeft, spectrogramBottom, totalWidth, spectrogramHeight];
        tracePosition = [frameLeft, traceBottom, totalWidth, traceHeight];
    
        return spectrogramPosition, tracePosition  

def compute_spectrogram(tr, per_lap=0.9, wlen=None, mult=8.0):
    """
        Computes spectrogram of the input data.
        Modified from obspy.imaging.spectrogram because we want the plotting part
        of that method in a different function to this.

        :type tr: Trace
        :param tr: Trace object to compute spectrogram for
        :type per_lap: float
        :param per_lap: Percentage of overlap of sliding window, ranging from 0
        to 1. High overlaps take a long time to compute.
        :type wlen: int or float
        :param wlen: Window length for fft in seconds. If this parameter is too
            small, the calculation will take forever. If None, it defaults to
            (samp_rate/100.0).
        :type mult: float
        :param mult: Pad zeros to length mult * wlen. This will make the
            spectrogram smoother.

    """

    Fs = float(tr.stats.sampling_rate)
    y = tr.data
    npts = tr.stats.npts
        
    # set wlen from samp_rate if not specified otherwise
    if not wlen:
        wlen = Fs / 100.0
 
    # nfft needs to be an integer, otherwise a deprecation will be raised
    # XXX add condition for too many windows => calculation takes for ever
    nfft = int(_nearest_pow_2(wlen * Fs))
    if nfft > npts:
        nfft = int(_nearest_pow_2(npts / 8.0))

    if mult is not None:
        mult = int(_nearest_pow_2(mult))
        mult = mult * nfft
    nlap = int(nfft * float(per_lap))

    y = y - y.mean()

    # Here we call not plt.specgram as this already produces a plot
    # matplotlib.mlab.specgram should be faster as it computes only the
    # arrays
    S, F, T = mlab.specgram(y, Fs=Fs, NFFT=nfft, pad_to=mult, noverlap=nlap, mode='magnitude')
        
    return T, F, S

def amp2dB(X):
    return 20 * np.log10(X)