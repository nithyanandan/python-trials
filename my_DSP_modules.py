import sys as SYS
import numpy as NP
import my_operations as OPS

def FT1D(inp, ax=-1, use_real=False, shift=False, verbose=True):
    """
    ---------------------------------------------------------------------
    Compute FFT from Numpy. 
    Inputs:

    inp:    Input data (vector or array) to be Fourier transformed

    Keyword Inputs:

    ax:         Axis (scalar integer) over which FFT is performed. 
                Default = -1 (last axis)

    use_real:   [Boolean scalar] If True, compute only the positive
                frequency components using the real part of the data

    Fftoututs:    
    
    fftout: FFT of input data over the specified axes
    -------------------------------------------------------------------
    """

    inp = NP.asarray(inp)

    try:
        isinstance(inp, NP.ndarray)
        # type(inp) is numpy.ndarray
    except TypeError: 
        print 'Unable to convert to Numpy array data type'
        sys.exit(1) # Abort execution

    if use_real:
        inp = NP.real(inp)
        if verbose:
            print "Opted for FFT of real data. Hence performing numpy.rfft()."
            print "numpy.rfft() returns only positive frequencies."
        fftout = NP.fft.rfft(inp, axis=ax)
    else:
        fftout = NP.fft.fft(inp, axis=ax)

    if shift:
        fftout = NP.fft.fftshift(fftout, axes=ax)
    return fftout


def frequencies(length, delx=1.0, shift=False, use_real=False):
    """
    ----------------------------------------------------------------
    Compute frequency samples in the FFT

    Inputs:

    length:    Length of vector to be Fourier transformed

    Keyword Inputs:

    delx:        x-axis interval, used only in case of 1D inp.
                 Default = 1.0

    shift:       [Boolean scalar] True => Shift to center of frequencies

    use_real:    [Boolean scalar] True => Compute only positive 
                 frequencies using numpy.fft.rfftfreq() 

    Output:    
    
    freqs: Discrete frequencies in the output FFT
    ---------------------------------------------------------------
    """
    
    # try: 
    #     size(length) == 1 and isinstance(length, int)
    #     print type(length)
    # except: 
    #     print "length has to be a scalar positive integer."
    #     print "Aborted execution in my_DSP_modules.frequencies()"
    #     sys.exit(1) # Abort execution

    if use_real:
        freqs = NP.fft.rfftfreq(length, d=delx)
    else: 
        freqs = NP.fft.fftfreq(length, d=delx)
        if shift:
            freqs = NP.fft.fftshift(freqs)

    return freqs


def rfft_append(inp, axis=0):

    """
    ------------------------------------------------------------------
    Compute the negative frequency left out by numpy.rfft()
    and append in the right order to the output from numpy.rfft().

    Input:

    inp       Input data of any dimensions to which negative frequency 
              components have to be appended.

    Keyword Input: 

    axis      [scalar] Axis along which negative frequency components
              are to be appended. It has to be a scalar in the range
              0 to Ndim-1 where Ndim is the number of axes in the data.

    Output:

    Appended data along the axis specified. 
    -------------------------------------------------------------------
    """

    axis = axis[0]
    shp = NP.shape(inp)
    ndim = len(shp)

    if (axis < 0) or (axis >= ndim):
        print "Input data does not contain the axis specified."
        print "Aborted execution in my_operations.reverse()"
        sys.exit(1) 

    if shp[axis] == 1:
        return inp

    return NP.append(inp, NP.conj(OPS.reverse(inp, axis=axis, ind_range=[1,shp[axis]-2])), axis=axis)


def rfftfreq_append(rfft_freqs):

    """
    ------------------------------------------------------------
    Compute the negative frequencies for the output of 
    numpy.rfftfreq() and rearrange the frequencies in the correct
    order. 

    Input: 

    rfft_freqs      [Vector] Positive frequencies

    Output:

    Positive and negative frequencies computed from numpy.rfftfreq()
    made equal to the output of numpy.fftfreq()
    ------------------------------------------------------------
    """

    rfft_freqs = NP.asarray(rfft_freqs)
    return NP.append(rfft_freqs[:-1],-rfft_freqs[-1:0:-1],axis=0)
