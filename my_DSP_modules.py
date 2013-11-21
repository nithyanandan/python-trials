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

    oututs:    
    
    fftout: FFT of input data over the specified axes
    -------------------------------------------------------------------
    """

    try:
        inp
    except NameError:
        raise NameError('inp not defined. Aborting FT1D().')

    if not isinstance(inp, NP.ndarray):   # type(inp) is numpy.ndarray
        raise TypeError('Input array should be Numpy array data type')

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


def spectral_axis(length, delx=1.0, shift=False, use_real=False):
    """
    ----------------------------------------------------------------
    Compute spectral axis in the FFT

    Inputs:

    length:    Length of vector to be Fourier transformed

    Keyword Inputs:

    delx:        x-axis interval, used only in case of 1D inp.
                 Default = 1.0

    shift:       [Boolean scalar] True => Shift to center of frequencies

    use_real:    [Boolean scalar] True => Compute only positive 
                 frequencies using numpy.fft.rfftfreq() 

    Output:    
    
    spaxis: Discrete spectral axis in the output FFT
    ---------------------------------------------------------------
    """
    
    # try: 
    #     size(length) == 1 and isinstance(length, int)
    #     print type(length)
    # except: 
    #     print "length has to be a scalar positive integer."
    #     print "Aborted execution in my_DSP_modules.frequencies()"
    #     SYS.exit(1) # Abort execution

    if use_real:
        spaxis = NP.fft.rfftfreq(length, d=delx)
    else: 
        spaxis = NP.fft.fftfreq(length, d=delx)
        if shift:
            spaxis = NP.fft.fftshift(spaxis)

    return spaxis


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

    try:
        inp
    except NameError:
        raise NameError('inp undefined. Aborting rfft_append()')

    if not isinstance(inp, NP.ndarray):
        raise TypeError('inp should be Numpy array data type.')

    if isinstance(axis, (list, tuple, str)):
        raise TypeError('axis should be a scalar integer in the range 0 to Ndim-1')
    axis = int(axis)
    shp = NP.shape(inp)
    ndim = len(shp)

    if (axis < 0) or (axis >= ndim):
        raise ValueError("Input data does not contain the axis specified. Aborted execution in reverse()")

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

    try:
        rfft_freqs
    except NameError:
        raise NameError('Input rfft_freqs not specified. Aborting rfftfreq_append()')

    if not isinstance(rfft_freqs, (list, NP.ndarray)):
        raise TypeError('Input rfft_freqs should be a list or a 1D Numpy array')

    rfft_freqs = NP.asarray(rfft_freqs)

    return NP.append(rfft_freqs[:-1],-rfft_freqs[-1:0:-1],axis=0)
