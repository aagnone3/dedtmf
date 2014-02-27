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

    E = np.r_[R[0], np.zeros(p)]

    fa = np.zeros(p+1)
    fa[0] = 1

    if R[0] > 0:
      # Only do the recursion if this frame is nonempty
      for i in xrange(1, p + 1):
        c = 0
        for j in xrange(1, i):
           c += a[j, i-1] * R[i-j]
        k[i] = (R[i] - c) / E[i-1]

        a[i][i] = k[i]

        for j in xrange(1, i):
           a[j][i] = a[j][i-1] - k[i] * a[i-j][i-1]

        E[i] = (1 - k[i]**2) * E[i-1]

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
        if E[row] == 0:
            print "Empty frame: ", row, " of ", rows
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

def dedtmf(X, lpc_order=40, win_pts=4096, hop_pts=256, params={}):
    """ Remove sustained, steady tones from an audio signal by finding strong 
        poles in LPC analysis on long-ish blocks, then inverse-filtering to 
        remove them.
    """
    # Set up default params
    # poleradthresh is the threshold for "capturing" poles
    if 'poleradthresh' in params:
        poleradthresh = params['poleradthresh']
    else:
        poleradthresh = 0.98
    # poleradtrans is the transition width for pole capture
    if 'poleradtrans' in params:
        poleradtrans = params['poleradtrans']
    else:
        poleradtrans = 0.002
    # polerad is the radius of the new poles put behind the introduced 
    # circles to limit the notch width
    if 'polerad' in params:
        polerad = params['polerad']
    else:
        polerad = 0.98

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
    Cf = Rf/(Mf + (Mf==0))

    # Modify magnitude; ensure smaller than 1
    # Use a sigmoid to make poles close to 1 even closer
    # (defaults now set at top of function)
    #poleradthresh = 0.98
    #poleradtrans =  0.002
    #Mfm = (Mf > poleradthresh);
    Mfm = sigmoid( (Mf - poleradthresh)/poleradtrans )
    # New polynomial with exaggerated pole radii
    Bfe = polybyrow( Mfm * Cf )
    # Compensatory poles - just inside zeros
    #polerad = 0.98
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

def to_num(s):
    ''' Converts a string to an int or a float if that works, 
        else leaves it as a string'''
    try:
        return int(s)
    except ValueError:
        try:
            return float(s)
        except ValueError:
            return s

def parse_param_file(filename, defaults):
    """ Read a param file of "key value" lines.  Return a dictionary of 
        assignments.  Fill in defaults from default dict """
    result = defaults
    if len(filename):
        with open(filename, 'r') as f:
            for line in f:
                keyval = line.lstrip().split(None, 1)
                if len(keyval) > 0:
                    key = keyval[0]
                    if len(keyval) > 1:
                        val = keyval[1]
                    else:
                        val = []
                    if key != "" and key[0] != '#':
                        # Overwrite or create
                        result[key] = to_num(val)
    return result

def main(argv):
    """ Main routine to apply dedtmf from command line """
    if len(argv) != 3 and len(argv) != 4:
        raise NameError( ("Usage: ", argv[0], 
                          " inputsound.wav outsound.wav [paramfile]") )

    # Input and output files
    inwavfile = argv[1]
    outwavfile = argv[2]

    # Read in param file (if any) and/or setup default values
    paramfile = ''
    if len(argv) == 4:
        paramfile = argv[3]
    params = parse_param_file(paramfile, 
                              {'win_pts':4096, 
                               'hop_pts':256, 
                               'lpc_order':40, 
                               'poleradthresh':0.98, 
                               'poleradtrans':0.002, 
                               'polerad':0.98})
    #print params

    # Read in wav file
    srate, wavd = wav.read(inwavfile)
    # normalize short ints to floats of -1 / 1
    data = np.asfarray(wavd) / 32768.0  

    lendata = len(data)
    #win_pts = 4096
    #hop_pts = 256
    win_pts = params['win_pts']
    hop_pts = params['hop_pts']
    # Pad with zeros - half a window before, a whole window after 
    # just to make sure dropping the final part-window doesn't lose anything 
    prepad_pts = win_pts/2
    pdata = np.r_[np.zeros(prepad_pts), data, np.zeros(win_pts)]

    # Apply
    #lpc_order = 40
    lpc_order = params['lpc_order']
    filtrd = dedtmf(pdata, lpc_order, win_pts, hop_pts, params)

    # Strip the padding
    udata = filtrd[prepad_pts:prepad_pts + lendata]

    # Write the wav out
    wav.write(outwavfile, srate, (32768.0*udata).astype(np.int16))
    print "Wrote ", outwavfile

# Actually run main
if __name__ == "__main__":
    main(sys.argv)



