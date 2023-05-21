import numpy as np
import math
from copy import deepcopy
from collections import deque
from collections import OrderedDict

class FLC():

    """ FLC filter class

    Attributes
    ----------
    n : int
        Number of harmonics
    X : ndarray
        Reference input vector
    W : ndarray
        Weights
    V : ndarray
        Angular  frequencies
    mu : float
        Adaptive filter gain
    f0 : float
        frequency
    """

    def __init__(self, n=5, mu=0.07, f0=8):
        """
        Parameters
        ----------
        n : int
            Number of harmonics
        mu : float
            Adaptive filter gain
        f0 : float
            Frequency of input signal
        """

        self.n = n
        self.mu = mu
        self.f0 = f0
        self.X = np.zeros(shape=(2, n))
        self.W = np.zeros(shape=(2, n))
        self.V = np.array(np.zeros([n]))

        for i in range(self.n):
            self.V[i] = i * 2 * math.pi * self.f0;


    def FLC(self, k, s):
        """ FLC filter

        Parameters
        ----------
        k : float
            Time instant
        s : float
            Reference signal

        Returns
        -------
        y : float
            Estimated signal
        """

        # Find reference input vector
        for i in range(self.n):
            self.X[0][i] = math.sin(self.V[i] * k)
            self.X[1][i] = math.cos(self.V[i] * k)
        # Find estimated signal
        y = np.dot(np.transpose(self.W[0]), self.X[0]) + np.dot(np.transpose(self.W[1]), self.X[1])

        err = s - y
        # Update weights

        self.W[0] += 2 * self.mu * self.X[0] * err
        self.W[1] += 2 * self.mu * self.X[1] * err

        return y

class WFLC():
    """ WFLC filter class

    Attributes
    ----------
    n : int
        Number of harmonics
    X : ndarray
        Reference input vector
    W : ndarray
        Weights
    V : ndarray
        Angular  frequencies
    mu : float
        Adaptive filter gain
    v0 : float
        ω0 fundamental angular frequency
    """

    def __init__(self, n=1, mu=0.001, mu0=0.0001, f0 = 6):
        """
        Parameters
        ----------
        n : int
            Number of harmonics
        mu : float
            Adaptive filter gain for reference input vector
        mu0 : float
            Adaptive filter gain for fundamental frequency
        """

        self.n = n
        self.mu = mu
        self.mu0 = mu0
        self.v0 = 2*math.pi*f0
        self.X = np.zeros(shape=(2, n))
        self.W = np.zeros(shape=(2, n))
        self.V = np.array(np.zeros([n]))
        self.estimatedFrequency = 0


    def WFLC(self,k, s):
        """ FLC filter

        Parameters
        ----------
        k : float
            Time instant
        s : float
            Reference signal

        Returns
        -------
        y : float
            Estimated signal
        """


        # Find reference input vector
        for i in range(self.n):
            self.X[0][i] = math.sin((i+1) * self.v0 * k)
            self.X[1][i] = math.cos((i+1) * self.v0 * k)
        # Find estimated signal
        y = np.dot(np.transpose(self.W[0]), self.X[0]) + np.dot(np.transpose(self.W[1]), self.X[1])

        err = s - y


        # Update fundamental angular frequency
        z = 0
        for i in range(self.n):
            z += (i+1)*(self.W[0][i]*self.X[1][i] - self.W[1][i]*self.X[0][i])
        self.v0 = self.v0 + 2*self.mu0*err*z

        # Update weights
        self.W[0] += 2 * self.mu * self.X[0] * err
        self.W[1] += 2 * self.mu * self.X[1] * err



        self.estimatedFrequency = self.v0 / (2 * math.pi)


        return y

class BMFLC():

    """ FLC filter class

    Attributes
    ----------
    n : int
        Number of harmonics
    X : ndarray
        Reference input vector
    W : ndarray
        Weights
    V : ndarray
        Angular  frequencies
    mu : float
        Adaptive filter gain
    f0 : float
        Starting frequency
    dF : float
        Frequrncey step
    """

    def __init__(self, mu=0.01, fmin=3, fmax=9, dF=0.1):
        """
        Parameters
        ----------
        n : int
            Number of harmonics
        mu : float
            Adaptive filter gain
        f0 : float
            Starting frequency
        """

        self.n = int((fmax - fmin) / dF) + 1
        self.mu = mu
        self.fmax = fmax
        self.fmin = fmin
        self.X = np.zeros(shape=(2, self.n))
        self.W = np.zeros(shape=(2, self.n))
        self.V = np.array(np.zeros([self.n]))
        self.estimatedFrequency = 0

        for i in range(self.n):
            self.V[i] = 2 * math.pi * (self.fmin + dF * i);


    def BMFLC(self, k, s):
        """ BMFLC filter

        Parameters
        ----------
        k : float
            Time instant
        s : float
            Reference signal

        Returns
        -------
        y : float
            Estimated signal
        """
        for i in range(self.n):
            self.X[0][i] = math.sin(self.V[i] * k)
            self.X[1][i] = math.cos(self.V[i] * k)

        y = np.dot(np.transpose(self.W[0]), self.X[0]) + np.dot(np.transpose(self.W[1]), self.X[1])

        err = s - y

        # Update weights
        for i in range(self.n):
            self.W[0][i] += 2 * self.mu * self.X[0][i] * err
            self.W[1][i] += 2 * self.mu * self.X[1][i] * err

        a=0
        b=0
        vest = 0

        for i in range(self.n):
            a += (self.W[0][i]**2+self.W[1][i]**2)*self.V[i]
            b += self.W[0][i] ** 2 + self.W[1][i] ** 2
        vest += a/b

        self.estimatedFrequency = vest/(2*math.pi)

        return y

class BMWFLC():
    """ FLC filter class

    Attributes
    ----------
    n : int
        Number of harmonics
    X : ndarray
        Reference input vector
    W : ndarray
        Weights
    V : ndarray
        Angular  frequencies
    mu : float
        Adaptive filter gain
    f_min : float
        Minimum frequency
    f_max : float
        Maximum frequency
    dF : float
        Frequency step
    """

    def __init__(self, mu=1.0, kappa = 0.009, g = 100, h = 0.00001,eta = 0.4 , f_min=6, f_max=7, dF=0.1, dT=0.01, Tp=2, alpha=0.67, beta=100, l=0.1, peaks_to_track = 2, adaptive_lr = True):
        """
        Parameters
        ----------
        n : int
            Number of harmonics
        mu : float
            Adaptive filter gain
        f_min : float
            Minimum frequency
        f_max : float
            Maximum frequency
        dT : float
            Sampling time in seconds
        Tp : float
            Width of memory window in seconds
        alpha : float
            Minimum amplification gain for memory window
        beta : float
            Multiplier for mu0 learning rate. The learning rate will start on the multiplied
            value of mu0, and decay to it reaches mu0.
        l : float
            Decay constant for mu0. l for lambda.
        d : float
            Scaling factor for mu. How much to scale mu with respect to max amplitude
            from the input signal
            mu = d/maxAmplitude
        g : int
            length of sliding window for max amplitude.
            If g = 100, the max amplitude from the last 100 samples
            from the input signal will be used to scale mu
        h : float
            factor to decrease mu0 with respect to mu.
            mu0 = mu*h
        eta: float
            Threshold when angular frequency should be reset for having
            to low magnitude
        """

        self.n = int((f_max - f_min) / dF) + 1
        self.f_max = f_max
        self.f_min = f_min
        self.X = np.zeros(shape=(2, self.n))
        self.W = np.zeros(shape=(2, self.n))
        self.V = np.array(np.zeros([self.n]))
        self.magnitudes = np.zeros(self.n)

        self.l = l
        for i in range(0,self.n):
            self.V[i] = 2 * math.pi * (self.f_min + dF * i)
        self.Vref = deepcopy(self.V)

        # Peak stuff
        self.allPeaksSorted = []
        self.trackedPeaksObj = {}
        self.peaks_to_track = peaks_to_track
        self.eta = eta

        # Memory window
        delta = (1 / dT) * Tp
        self.rho = (alpha) ** (1 / delta)

        # Dynamic learning rate mu
        self.adaptive_lr = adaptive_lr
        self.mu = mu
        self.kappa = kappa
        self._peakAmplitude = 0
        self._q_amplitudes = deque(g * [0], maxlen=g)
        self.beta = beta
        self.h = h



    def BMWFLC(self, k, s):
        """ BMFLC filter

        Parameters
        ----------
        k : float
            Time instant
        s : float
            Reference signal

        Returns
        -------
        y : float
            Estimated signal
        """

        # Adapt learning rates mu and mu0

        if  self.adaptive_lr:
            self.adapt_learningrate(k, s)

        for i in range(self.n):
            self.X[0][i] = math.sin(self.V[i] * k)
            self.X[1][i] = math.cos(self.V[i] * k)

        y = np.dot(np.transpose(self.W[0]), self.X[0]) + np.dot(np.transpose(self.W[1]), self.X[1])

        err = s - y

        # Update weights
        for i in range(self.n):
            self.W[0][i] = self.W[0][i] * self.rho + 2 * self.mu * self.X[0][i] * err
            self.W[1][i] = self.W[1][i] * self.rho + 2 * self.mu * self.X[1][i] * err

        # Find peaks
        self.find_peaks(err, k, s)
        # Create objects for the peaks being tracked
        self.create_tracked_peak_objects(k)
        # Update the angular frequencies for the tracked peaks
        self.update_V_for_peaks(err)


        return y

    def find_peaks(self, err, k, s):

        # Find all magnitudes
        magnitudes = []
        magnitudes_dict = {x: math.sqrt(self.W[0][x] ** 2 + self.W[1][x] ** 2) for x in range(self.n)}
        for key, m in magnitudes_dict.items():
            magnitudes.append(m)

        self.magnitudes = magnitudes

        # Reset angle when magnitudes gets small
        for n in range(self.n):
            if magnitudes[n] < self.eta and self.V[n] != self.Vref[n]:
                self.V[n] = self.Vref[n]

        # Key: position of peak, Value: Magnitude of peak
        peaksDict = {}

        magnitudePrev = 0
        magnitudeDiffPrev = 0
        # Find peaks
        for i in range(self.n):
            magnitudeDiff = magnitudes[i] - magnitudePrev
            if magnitudeDiff < 0 and magnitudeDiffPrev > 0:
                peaksDict[i - 1] = magnitudePrev
            # Check if last magnitude is a peak
            elif i == (self.n -1) and magnitudeDiff > 0:
                peaksDict[i] = magnitudes[i]

            magnitudeDiffPrev = magnitudeDiff
            magnitudePrev = magnitudes[i]

        # Sort peaks by magnitude
        peaksDict = OrderedDict(sorted(peaksDict.items(), key=lambda t: t[1], reverse = True))
        self.allPeaksSorted = list(peaksDict.items())


    def update_V_for_peaks(self, err):
        for i in range(self.peaks_to_track):
            if len(self.allPeaksSorted) > i:

                # Update V for peaks
                peakPos = self.allPeaksSorted[i][0]

                if self.allPeaksSorted[i][0] >= 0:
                    z = (peakPos + 1) * (
                    self.W[0][peakPos] * self.X[1][peakPos] - self.W[1][peakPos] * self.X[0][peakPos])

                    self.V[peakPos] = self.V[peakPos] + 2 * self.trackedPeaksObj[peakPos].mu0 * err * z

    def create_tracked_peak_objects(self, k):
        # Create peak objects for new peaks that are being tracked
        for i in range(self.peaks_to_track):
            if len(self.allPeaksSorted) > i:
                if self.allPeaksSorted[i][0] not in self.trackedPeaksObj:
                    self.trackedPeaksObj[self.allPeaksSorted[i][0]] = self.peak(mu0=self.mu*self.h, position=self.allPeaksSorted[i][0],
                                                                                magnitude=self.allPeaksSorted[i][1],
                                                                                V=self.V[self.allPeaksSorted[i][0]], l= self.l,
                                                                                beta = self.beta, k = k, h=self.h)

        # Remove the peak objects that are no longer being tracked
        keep_peak = False
        trackedPeaksCopy = deepcopy(self.trackedPeaksObj)
        for key, value in trackedPeaksCopy.items():
            for i in range(self.peaks_to_track):
                if len(self.allPeaksSorted) > i:
                    if self.allPeaksSorted[i][0] == key:
                        keep_peak = True
            if not keep_peak:
                if self.peaks_to_track == 1:
                    self.V[key] = self.Vref[key]
                self.trackedPeaksObj.pop(key, None)
            keep_peak = False

    def adapt_learningrate(self, k, s):
        self._q_amplitudes.append(s)
        self._peakAmplitude = max(self._q_amplitudes)
        if self._peakAmplitude > 0:
            self.mu = self.kappa / self._peakAmplitude

        if len(self.trackedPeaksObj) > 0:
            for key in self.trackedPeaksObj:
                self.trackedPeaksObj[key].adapt_m0(k, s, self.mu)

    class peak():

        def __init__(self, mu0, position, magnitude,V, l, beta, k, h):
            self.mu0 = mu0
            self.h = h
            self.position = position
            self.V = V
            self.influence = 1
            self.magnitude = magnitude
            self.l = l
            self.beta = beta
            self._k_created = k

        def adapt_m0(self, k, s, mu):
            self.mu0 = mu * self.h
            # learning rate decay
            self.mu0 += self.mu0 * math.exp(-1 * self.l * (k - self._k_created)) * self.beta
            #print(math.exp(-1 * self.l * (k - self._k_created)) * self.beta)

class EBMFLC():

    """ FLC filter class

    Attributes
    ----------
    n : int
        Number of harmonics
    X : ndarray
        Reference input vector
    W : ndarray
        Weights
    V : ndarray
        Angular  frequencies
    mu : float
        Adaptive filter gain
    f0 : float
        Starting frequency
    dF : float
        Frequrncey step
    """

    def __init__(self, mu=0.01, fmin=0,fmax=20,fa = 0,fb=2, fc=6, fd=8, dF=0.2, dT=0.01, Tp = 2, alpha = 0.05):
        """
        Parameters
        ----------
        n : int
            Number of harmonics
        mu : float
            Adaptive filter gain
        fmin : float
            Starting frequency of complete frequency range
        fmax : float
            Max frequency of the complete frequency range
        fa : float
            Starting frequency of voluntary motion
        fb : float
            End frequency for voluntary motion
        fc : float
            Start frequency for involuntary motion
        fd : float
            End frequency for involuntary motion
        dT : float
            Sampling time in seconds
        Tp : float
            Width of memory window in seconds
        """

        self.Na = int((fa - 0) / dF)
        self.Nb = int((fb - 0) / dF)
        self.Nc = int((fc - 0) / dF)
        self.Nd = int(round((fd - 0) / dF)) + 1
        self.n = int((fmax - fmin) / dF) + 1
        self.mu = mu
        self.fmax = fmax
        self.fmin = fmin
        self.fa = fa
        self.fb = fb
        self.fc = fc
        self.fd = fd
        self.X = np.zeros(shape=(2, self.n))
        self.Xi = np.zeros(shape=(2, self.Nd - self.Nc))
        self.W = np.zeros(shape=(2, self.n))
        self.Wi = np.zeros(shape=(2, self.Nd - self.Nc))
        self.V = np.array(np.zeros([self.n]))

        delta = (1/dT)*Tp
        self.rho = (alpha)**(1/delta)

        self.estimatedFrequency = 0

        for i in range(self.n):
            self.V[i] = 2 * math.pi * (self.fmin + dF * i);

        self.Vab = self.V[self.Na:self.Nb]
        self.Vcd = self.V[self.Nc:self.Nd]

    def EBMFLC(self,k, s):
        """ BMFLC filter

        Parameters
        ----------
        k : float
            Time instant
        s : float
            Reference signal

        Returns
        -------
        y : float
            Estimated signal
        """
        for i in range(self.n):
            self.X[0][i] = math.sin(self.V[i] * k)
            self.X[1][i] = math.cos(self.V[i] * k)

        m = np.dot(np.transpose(self.W[0]), self.X[0]) + np.dot(np.transpose(self.W[1]), self.X[1])

        err = s - m

        # Update weights
        for i in range(self.n):
            self.W[0][i] = self.W[0][i]*self.rho + 2 * self.mu * self.X[0][i] * err
            self.W[1][i] = self.W[1][i]*self.rho + 2 * self.mu * self.X[1][i] * err

        self.Wi[0] = self.W[0][self.Nc:self.Nd]
        self.Wi[1] = self.W[1][self.Nc:self.Nd]

        self.Xi[0] = self.X[0][self.Nc:self.Nd]
        self.Xi[1] = self.X[1][self.Nc:self.Nd]

        mi = np.dot(np.transpose(self.Wi[0]), self.Xi[0]) + np.dot(np.transpose(self.Wi[1]), self.Xi[1])


        # a=0
        # b=0
        # vest = 0
        # for i in range((self.Nd-self.Nc)):
        #     a += (self.Wi[0][i]**2+self.Wi[1][i]**2)*self.Vcd[i]
        #     for j in range((self.Nd-self.Nc)):
        #         b += self.Wi[0][i] ** 2 + self.Wi[1][j] ** 2
        #     vest += a/b
        #     a=0
        #     b=0
        # self.estimatedFrequency = vest/(2*math.pi)

        b = 0
        for j in range(self.Nd-self.Nc):
            b += self.Wi[0][j]**2 + self.Wi[1][j]**2

        if b == 0:
            raise ValueError("Denominator cannot be zero")

        a = 0
        vest = 0
        for i in range(self.Nd-self.Nc):
            a += (self.Wi[0][i]**2 + self.Wi[1][i]**2) * self.Vcd[i]
            vest += a / b
            a = 0

        self.estimatedFrequency = vest / (2 * math.pi)

        return mi
