#!/usr/bin/env python
"""
dedtmf.py

Python port of the dtmf-tone removal matlab code.

2014-01-23 Dan Ellis dpwe@ee.columbia.edu
"""

import numpy as np
# For filter
import scipy.signal
# For SRI's wavreading code
import scipy.io.wavfile as wav
# For command line
import sys

def frame(x, window, hop):
    """ Convert vector x into an array of successive window-point segments, 
        stepped by hop
        Done with stride_tricks, no copying, but have to lose the final 
        part-window 
    """

    nframes = 1 + int( np.floor( (len(x)-window)/hop ))
    shape = (nframes, window)
    strides = (x.strides[0] * hop, x.strides[0])
    return np.lib.stride_tricks.as_strided(x, shape=shape, strides=strides)

def ola(X, hop):
    """ Overlap-add rows of X into a single vector, stepping by hop """
    nw, W = np.shape(X)

    # How long should X be in the end?
    # one complete window, plus another hop's worth of points for each
    # additional window
    lenx = W + (nw-1)*hop
    Y = np.zeros(lenx)

    for i in range(nw):
        Y[i*hop : i*hop+W] += X[i,]

    return Y

# based on https://code.google.com/p/vocoder-lpc/source/browse/lpc.py

def autoco(x, p):
    """
    Calculate autocorrelation of x
    Return only the first p values, i.e. r_xx[0] to r_xx[p-1]
    """
    # Naive version
    #R = np.correlate(x, x, 'full')
    #R = R[len(x)-1:len(x)+p]
    # Better version (> twice as fast?)
    X = np.fft.rfft(np.r_[x, np.zeros(p-1)])
    rxx = np.fft.irfft(np.square(np.abs(X)))
    return rxx[:p]

def lpc_frame(x, p):
    """
    Fit p'th order LPC model to a single frame of data in x
    using Levinson-Durbin recursion
    """
    R = autoco(x, p+1)
    a = np.zeros((p+1, p+1))
    k = np.zeros(p+1)

    E = np.r_[R[0], np.empty(p)]
    for i in xrange(1, p + 1):
      c = 0
      for j in xrange(1, i):
         c += a[j, i-1] * R[i-j]
      k[i] = (R[i] - c) / E[i-1]

      a[i][i] = k[i]

      for j in xrange(1, i):
         a[j][i] = a[j][i-1] - k[i] * a[i-j][i-1]

      E[i] = (1 - k[i]**2) * E[i-1]

    fa = np.empty(p+1)
    fa[0] = 1
    for j in xrange(1, p+1):
      fa[j] = -a[j][p]
    return fa, E[p]

def lpc(X, order):
    """ Calculate LPC fit to each row of X, return as rows of A matrix, 
        plus vector of energies
    """
    rows, cols = np.shape(X)
    A = np.zeros( (rows, order+1) )
    E = np.zeros(rows)
    for row in range(rows):
        #print "frame ", row, " of ", rows
        A[row,], E[row] = lpc_frame(X[row,:], order)
    return A, E


def rootsbyrow(A):
    """ Calculate the roots for polynomials defined by each row of A
    """
    nr, nc = np.shape(A)
    R = np.zeros( (nr, nc-1), complex )
    for i in range(nr):
        R[i,] = np.roots(A[i,])
    return R

def polybyrow(R):
    """ Calculate the polynomial coefficients for rows of roots in R
    """
    nr, nc = np.shape(R)
    A = np.zeros((nr, nc+1))
    for i in range(nr):
        A[i,] = np.real(np.poly(R[i,]))
    return A

def sigmoid(X):
    """ simple sigmoid function """
    return 1./(1.+np.exp(-X))

########## main dedtmf ############

def dedtmf(X, lpc_order=40, win_pts=4096, hop_pts=256):
    """ Remove sustained, steady tones from an audio signal by finding strong 
        poles in LPC analysis on long-ish blocks, then inverse-filtering to 
        remove them.
    """
    # Default hop is half whatever window is
    if hop_pts == 0:
        hop_pts = win_pts/2

    # Break x into w-length blocks stepped by h
    Xf = frame(X, win_pts, hop_pts)

    print "lpc fitting frames..."
    # Fit LPC models to every frame
    Af, Ef = lpc(Xf, lpc_order)

    # Find all the poles of the per-frame LPC filters
    print "finding roots..."
    Rf = rootsbyrow(Af);
    # Magnitudes of each pole
    Mf = np.absolute(Rf)
    # Corresponding pure-phase parts
    Cf = Rf/Mf

    # Modify magnitude; ensure smaller than 1
    # Use a sigmoid to make poles close to 1 even closer
    poleradthresh = 0.98
    poleradtrans =  0.002
    #Mfm = (Mf > poleradthresh);
    Mfm = sigmoid( (Mf - poleradthresh)/poleradtrans )
    # New polynomial with exaggerated pole radii
    Bfe = polybyrow( Mfm * Cf )
    # Compensatory poles - just inside zeros
    polerad = 0.98
    Afe = polybyrow( polerad*Mfm * Cf);

    # Apply inverse filters to each block of X
    print "applying filters ..."
    nr, nc = np.shape(Xf)
    Xff = np.zeros( (nr, nc) )
    for i in range(nr):
        Xff[i,] = scipy.signal.lfilter(Bfe[i,], Afe[i,], Xf[i,])

    # Reconstruct with 50%-overlapped hanning window
    L = 2 * hop_pts
    pad = win_pts - L
    # Was an L point hanning, this is to make it flat
    win = np.r_[np.zeros(int(np.floor(pad/2))), 
                np.hanning(L+1),  
                np.zeros(int(np.ceil(pad/2)-1))]

    # Figure per-block energy ratio
    #ER = (np.sum(np.power(Xff*win, 2.0), axis=-1) /
    #      np.sum(np.power(Xf*win, 2.0), axis=-1))
    # Threshold
    #ER = np.minimum(ER, 1.0);
    # scale
    #Xff = repmat(ER, size(Xff,1), 1).* Xff;
    
    print "overlap-add ..."
    return ola(Xff*win, hop_pts)


############## Provide a command-line wrapper

def main(argv):
    """ Main routine to apply dedtmf from command line """
    if len(argv) != 3:
        raise NameError( ("Usage: ", argv[0], 
                          " inputsound.wav outsound.wav") )

    inwavfile = argv[1]
    outwavfile = argv[2]

    # Read in wav file
    srate, wavd = wav.read(inwavfile)
    # normalize short ints to floats of -1 / 1
    data = np.asfarray(wavd) / 32768.0  

    lendata = len(data)
    win_pts = 4096
    hop_pts = 256
    # Pad with zeros - half a window before, a whole window after 
    # just to make sure dropping the final part-window doesn't lose anything 
    prepad_pts = win_pts/2
    pdata = np.r_[np.zeros(prepad_pts),data,np.zeros(win_pts)]

    # Apply
    lpc_order = 40
    filtrd = dedtmf(pdata, lpc_order, win_pts, hop_pts)

    # Strip the padding
    udata = filtrd[prepad_pts:prepad_pts + lendata]

    # Write the wav out
    wav.write(outwavfile, srate, (32768.0*udata).astype(np.int16))
    print "Wrote ", outwavfile

# Actually run main
if __name__ == "__main__":
    main(sys.argv)



