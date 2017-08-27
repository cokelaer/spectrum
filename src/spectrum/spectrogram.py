from spectrum import Periodogram, pmtm
import numpy as np


__all__ = ["Spectrogram"]


class Spectrogram(object):
    """Simple example of spectrogram

    .. plot::

        from spectrum import Spectrogram, dolphin_filename, readwav
        data, samplerate = readwav(dolphin_filename)

        p = Spectrogram(data, ws=128, W=4096, sampling=samplerate)
        p.periodogram()
        p.plot()

    .. warning:: this is a prototype and need careful checking about x/y axis


    """
    def __init__(self, signal, ws=128, W=4096, sampling=1, channel=1):
        if len(signal.shape) == 1:
            self.signal = signal
        else:
            self.signal = signal[:,channel-1]
        self.W = W
        self.ws = ws
        self._start_y = 10
        self.sampling = sampling
        self.duration = len(self.signal) / float(self.sampling)

    def plot(self, filename=None, vmin=None, vmax=None, cmap='jet_r'):
        import pylab
        pylab.clf()
        pylab.imshow(-np.log10(self.results[self._start_y:,:]), 
            origin="lower",
            aspect="auto", cmap=cmap, vmin=vmin, vmax=vmax)
        pylab.colorbar()

        # Fix xticks
        XMAX = float(self.results.shape[1])  # The max integer on xaxis
        xpos = list(range(0, int(XMAX), int(XMAX/5)))
        xx = [int(this*100)/100 for this in np.array(xpos) / XMAX * self.duration]
        pylab.xticks(xpos, xx, fontsize=16)

        # Fix yticks
        YMAX = float(self.results.shape[0])  # The max integer on xaxis
        ypos = list(range(0, int(YMAX), int(YMAX/5)))
        yy = [int(this) for this in np.array(ypos) / YMAX * self.sampling]
        pylab.yticks(ypos, yy, fontsize=16)

        #pylab.yticks([1000,2000,3000,4000], [5500,11000,16500,22000], fontsize=16)
        #pylab.title("%s echoes" %  filename.replace(".png", ""), fontsize=25)
        pylab.xlabel("Time (seconds)", fontsize=25)
        pylab.ylabel("Frequence (Hz)", fontsize=25)
        pylab.tight_layout()
        if filename:
            pylab.savefig(filename)

    def periodogram(self):
        W = self.W
        ws = self.ws
        N = int(len(self.signal)/ws)
        self.results = np.zeros((W*2+1, N-8))
        print("Duration: %s" % self.duration)
        print("W: %s" % W)
        print("ws: %s" % ws)
        print("Computing %s TFs" % N)
        for i in range(N-8):
            data = self.signal[i*ws:i*ws+W]
            p = Periodogram(data, sampling=self.sampling, NFFT=W*4)
            p()
            self.results[:,i] = p.psd
        print("done")

    def pmtm(self):
        W = self.W
        ws = self.ws
        N = int(len(self.signal) / ws)
        self.results = np.zeros((W+1, N-8))
        for i in range(N-8):
            data = self.signal[i*ws:i*ws+W]
            a = pmtm(data, 4, NFFT=W*4, show=False)
            Sk = np.mean(abs(a[0].transpose())**2 * a[1], axis=1)
            self.results[:, i] = Sk[0:self.W+1]
            print(i, N)
        print("done")


