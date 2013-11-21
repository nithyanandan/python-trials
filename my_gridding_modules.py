import numpy as NP
import scipy as SP

#############################################################################

def grid(rangelist, pad=None, spacing=None, pow2=True, verbose=True):
    """
    ------------------------------------------------------------------
    Produce a multi-dimensional grid.
    
    Inputs:
    rangelist [list of tuples] Each tuple is made of two
              elements, the min and max with min < max. One tuple for 
              each axis.
             
    pad       [Optional. Scalar or list] The padding (in same units as
              range) to be applied along the axes. Default=None implies
              no padding.
              
    spacing   [Optional. Scalar or list] The spacing for the grid
              along each of the axes. If not supplied, a default of 
              sqrt(max-min) is used. If a scalar is supplied, it applies
              for all axes. A list applies for each of the axes.
              
    pow2      [Optional, default=True] If set, the grid produced is a 
              power of 2 in all axes, which is useful for FFT.
              
    verbose   [Default=True]

    Outputs:

    tupleout  [List of tuples] A 4-element tuple for each axis. The 
              elements in each tuple are min, max, lengths, and spacing 
              (which could have been modified relative to input). 
              The precise grid points can be generated by using numpy's
              linspace using min, max and lengths.
    ------------------------------------------------------------------
    """

    for item in rangelist:
        if item[0] >= item[1]:
            raise ValueError('Data ranges provided not compatible with min < max. Exiting from grid().')

    if pad is None:
        pad = NP.zeros(len(rangelist))
    elif len(pad) == 1:
        pad = [pad] * len(rangelist)
    elif len(pad) > len(rangelist):
        pad = NP.asarray(pad[:len(rangelist)])
    elif (len(pad) > 1) and (len(pad) < len(rangelist)):
        if verbose is True:
            print 'Insufficient paddings provided compared to the number of data ranges.'
            print 'Assuming the remaining paddings to be zero.'
        pad += [0.0 for ranges in rangelist[len(pad):]]
    pad = NP.reshape(NP.asarray(pad),len(pad)) # Force it to be row vector
    pad.clip(min(NP.abs(pad))) # Remove any negative values for pad

    if spacing is None:
        if verbose is True:
            print 'No spacing provided. Setting defaults to sqrt(max-min)'
            print 'Final spacings could be different from defaults assumed.'
        spacing = [NP.sqrt(rangelist[i][1]-rangelist[i][0]+2*pad[i]) for i in range(len(rangelist))]
    elif len(spacing) == 1:
        if verbose is True:
            print 'Scalar value for spacing provided. Assuming spacing is identical along all axes.'
        spacing = [spacing] * len(rangelist)
    elif len(spacing) > len(rangelist):
        if verbose is True:
            print 'Too many values of spacing provided. Ignoring values for indices beyond the length of data ranges.'
        spacing = NP.asarray(spacing[:len(rangelist)])
    elif (len(spacing) > 1) and (len(spacing) < len(rangelist)):
        if verbose is True:
            print 'Insufficient spacings provided compared to the number of data ranges.'
            print 'Assuming the remaining spacings to be default spacing.'
            print 'Final spacings could be different from defaults assumed.'
        # spacing += [NP.sqrt(ranges[1]-ranges[0]) for ranges in rangelist[len(spacing):]]
        spacing += [NP.sqrt(rangelist[i][1]-rangelist[i][0]+2*pad[i]) for i in range(len(spacing),len(rangelist))]
    spacing = NP.asarray(spacing)
    spacing.clip(min(NP.abs(spacing)))

    rangelist = NP.asarray(rangelist)
    lengths = NP.ceil((rangelist[:,1]-rangelist[:,0]+2*pad)/spacing)+1
    lengths.astype(int)

    for i in xrange(len(lengths)): 
        if (lengths[i] % 2) == 0: lengths[i] += 1
        # make it odd number of bin edges enclsoing first
        # and last intervals so that the mid-point is one
        # of the bin edges.

    if pow2 is True:
        lengths = 2**NP.ceil(NP.log2(lengths))
        newspacing = (rangelist[:,1]-rangelist[:,0]+2*pad)/(lengths-2)
        offsets = rangelist[:,0]-pad+(lengths-2)*newspacing - (rangelist[:,1]+pad)
        tupleout = zip(rangelist[:,0]-pad-0.5*offsets-newspacing, rangelist[:,1]+pad+0.5*offsets, lengths, newspacing) # converts numpy arrays into a list of tuples
        # tupleout = tuple(map(tuple, NP.column_stack((rangelist[:,0]-pad-0.5*offsets-newspacing, rangelist[:,1]+pad+0.5*offsets, lengths, newspacing)))) # converts numpy arrays into a list of tuples
    else:
        offsets = rangelist[:,0]-pad+(lengths-1)*spacing - (rangelist[:,1]+pad)
        tupleout = zip(rangelist[:,0]-pad-0.5*offsets, rangelist[:,1]+pad+0.5*offsets, lengths, spacing) # converts numpy arrays into a list of tuples
        # tupleout = tuple(map(tuple, NP.column_stack((rangelist[:,0]-pad-0.5*offsets, rangelist[:,1]+pad+0.5*offsets, lengths, spacing)))) # converts numpy arrays into a list of tuples
    
    return tupleout

########################################################################

def grid_1d(range, pad=None, spacing=None, pow2=True, verbose=True):
    """
    ------------------------------------------------------------------
    1D wrapper for grid()
    
    Inputs:
    range    [2-element list or a tuple with 2 elements] min and max
             with min < max.
             
    pad      [Optional or Scalar] The padding (in same units as
             range) to be applied along the axis. Default=None 
             implies no padding.

    spacing  [Optional or Scalar] The spacing for the grid
             along the axis. If not supplied, a default of 
             sqrt(max-min) is used. 

    pow2     [Optional, default=True] If set, the grid produced is a 
             power of 2 in its axis, which is useful for FFT.

    verbose  [Default=True]

    Outputs:

    A numpy array with x-values on the grid.
    ------------------------------------------------------------------
    """

    if range is None:
        raise NameError('No range provided. Exiting from grid_1d()')
    if isinstance(range, (list,tuple)) is False:
        raise TypeError('A 2-element list or tuple containing range with min < max should be provided.')

    if isinstance(range, tuple):
        grid_info = grid([range], pad=[pad], spacing=[spacing], pow2=pow2, verbose=verbose)
    if isinstance(range, list):
        grid_info = grid([tuple(range)], pad=[pad], spacing=[spacing], pow2=pow2, verbose=verbose)
    return NP.linspace(grid_info[0][0], grid_info[0][1], int(grid_info[0][2]))

###########################################################################

def grid_2d(range, pad=None, spacing=None, pow2=True, verbose=True):
    """
    ------------------------------------------------------------------
    2D wrapper for grid()
    
    Inputs:
    range    [2-element list of tuples] Each tuple is made of two
             elements, the min and max with min < max.
             
    pad      [Optional. Scalar or list] The padding (in same units as
             range) to be applied along the two axes. Default=None 
             implies no padding.

    spacing  [Optional. Scalar or list] The spacing for the grid
             along each of the axes. If not supplied, a default of 
             sqrt(max-min) is used. If a scalar is supplied, it applies
             for all axes. A list applies for each of the axes.

    pow2     [Optional, default=True] If set, the grid produced is a 
             power of 2 in all axes, which is useful for FFT.

    verbose  [Default=True]

    Outputs:

    Two 2D numpy arrays. The first array with x-values on the grid, and 
    the second with y-values on the grid.
    ------------------------------------------------------------------
    """
    if range is None:
        raise NameError('No range provided. Exiting from grid_2d()')
    if isinstance(range, list) is False:
        raise TypeError('A 2-element list of tuples specifying range with min < max should be provided. Exiting from grid_2d()')
    elif isinstance(range, list) is True:
        if isinstance(range[0], tuple) is False:
            raise TypeError('A 2-element list of tuples specifying range with min < max should be provided. Exiting from grid_2d()')

    grid_info = grid(range, pad=pad, spacing=spacing, pow2=pow2, verbose=verbose)
    return NP.meshgrid(NP.linspace(grid_info[0][0], grid_info[0][1], int(grid_info[0][2])), NP.linspace(grid_info[1][0], grid_info[1][1], int(grid_info[1][2])))

###########################################################################

def conv_grid1d(xc, xkern, kernel, xgrid, method='NN'):
    """
    ----------------------------------------------------------------------
    Perform 1D gridding convolution.
    Inputs:

    xc:      [vector] x-coordinates of center of gridding function 

    xkern:   [vector] x-coordinates of gridding function

    kernel:  [vector] gridding function (or kernel). If kernel is a
             scalar, a nearest neighbour interpolation is used overriding
             the method requested for.

    xgrid:   [vector] x-coordinates of grid locations

    Keyword Inputs:

    method:  String indicating interpolation method. [Default = 'NN']
             'NN' => Nearest Neighbour 
             'SL' => Single linear
             'CS' => Cubic Spline

    Output:

    outdata: [vector] Gridded values at values of xgrid

    ----------------------------------------------------------------------
    """
    try:
        xc
    except NameError:
        raise NameError("Argument 'xc' not defined. Aborting conv_grid1d().")

    try:
        xkern
    except NameError:
        raise NameError("Argument 'xkern' not defined. Aborting conv_grid1d().")

    try:
        xgrid
    except NameError:
        raise NameError("Argument 'xgrid' not defined. Aborting conv_grid1d().")

    try:
        method
    except NameError:
        method='NN'

    if (method != 'SL') and (method != 'CS'):
        method = 'NN'

    try:
        kernel
    except NameError:
        print "Argument 'kernel' not defined. "
        if method == 'NN':
            print "Since method is Nearest Neighbor interpolation, "
            print "proceeding with kernel=1. "
            kernel = 1.0
        else:
            raise TypeError("Aborting conv_grid1d().")

    if not isinstance(kernel, (list, NP.ndarray)):
        print ' Kernel seems to be a scalar. Proceeding with Nearest Neighbour \n method of interpolation.'
        method = 'NN'

    if len(xkern) != len(kernel):
        raise IndexError(' Incompatible kernel coordinates. Verify their lengths are equal.\n Aborting conv_grid1d().')

    outdata = NP.zeros(len(xgrid))
    xckern = 0.5*(max(xkern)+min(xkern))
    xshift = xc - xckern
    for npoints in range(0,len(xc)):
        if method == 'SL':
            interp_function = SP.interpolate.interp1d(xkern+xshift[npoints], kernel, kind='slinear', fill_value=0.0)
        elif method == 'CS':
            interp_function = SP.interpolate.interp1d(xkern+xshift[npoints], kernel, kind='cubic', fill_value=0.0)
        else:
            interp_function = SP.interpolate.interp1d(xkern+xshift[npoints], kernel, kind='nearest', fill_value=0.0)
        outdata += interp_function(xgrid)
    return outdata
        
###########################################################################

def conv_grid2d(xc, yc, xkern, ykern, kernel, xgrid, ygrid, method='NN'):
    """
    ----------------------------------------------------------------------
    Perform gridding convolution.
    Inputs:

    xc:      [vector] x-coordinates of center of gridding function 

    yc:      [vector] y-coordinates of center of gridding function 

    xkern:   [vector] x-coordinates of gridding function

    ykern:   [vector] y-coordinates of gridding function

    kernel:  [vector] gridding function (or kernel). If kernel is a
             scalar, a nearest neighbour interpolation is used overriding
             the method requested for.

    xgrid:   [vector] x-coordinates of grid locations

    ygrid:   [vector] y-coordinates of grid locations
    
    Keyword Inputs:

    method:  String indicating interpolation method. [Default = 'NN']
             'NN' => Nearest Neighbour 
             'BL' => Bilinear
             'CS' => Cubic Spline

    Output:

    outdata: [vector] Gridded values at values of (xgrid, ygrid) 

    ----------------------------------------------------------------------
    """
    try:
        xc
    except NameError:
        raise NameError("Argument 'xc' not defined. Aborting conv_grid2d().")

    try:
        yc
    except NameError:
        raise NameError("Argument 'yc' not defined. Aborting conv_grid().")

    try:
        xkern
    except NameError:
        raise NameError("Argument 'xkern' not defined. Aborting conv_grid2d().")

    try:
        ykern
    except NameError:
        raise NameError("Argument 'ykern' not defined. Aborting conv_grid2d().")

    try:
        xgrid
    except NameError:
        raise NameError("Argument 'xgrid' not defined. Aborting conv_grid2d().")

    try:
        ygrid
    except NameError:
        raise NameError("Argument 'ygrid' not defined. Aborting conv_grid2d().")

    try:
        method
    except NameError:
        method='NN'

    if (method != 'BL') and (method != 'CS'):
        method = 'NN'

    try:
        kernel
    except NameError:
        print "Argument 'kernel' not defined. "
        if method == 'NN':
            print "Since method is Nearest Neighbor interpolation, "
            print "proceeding with kernel=1. "
            kernel = 1.0
        else:
            raise TypeError("Aborting conv_grid2d().")

    if not isinstance(kernel, (list, NP.ndarray)):
        print "Kernel is a scalar. Proceeding with Nearest Neighbour"
        print "method of interpolation."
        method = 'NN'

    if isinstance(xc, list) is False:
        xc = [xc]
    if isinstance(yc, list) is False:
        yc = [yc]

    if (len(xc) != len(yc)):
        raise IndexError(" Incompatible input location coordinates.\n Verify their lengths are equal. Aborting conv_grid2d().")

    if (len(xkern) != len(ykern)) or (len(xkern) != len(kernel)):
        raise IndexError(" Incompatible kernel coordinates. Verify their lengths are equal.\n Aborting conv_grid2d().")

    if (len(xgrid) != len(ygrid)):
        raise IndexError(" Incompatible grid lattice coordinates. Verify their lengths are equal.\n Aborting conv_grid2d().")

    outdata = NP.zeros(len(xgrid))
    xckern = 0.5*(max(xkern)+min(xkern))
    yckern = 0.5*(max(ykern)+min(ykern))
    xshift = xc - xckern
    yshift = yc - yckern
    for npoints in range(0,len(xc)):
        if method == 'BL':
            outdata += SP.interpolate.griddata(zip(xkern+xshift[npoints],ykern+yshift[npoints]), kernel, zip(xgrid,ygrid), method='linear', fill_value=0.0)
        elif method == 'CS':
            outdata += SP.interpolate.griddata(zip(xkern+xshift[npoints],ykern+yshift[npoints]), kernel, zip(xgrid,ygrid), method='cubic', fill_value=0.0)
        else:
            outdata += SP.interpolate.griddata(zip(xkern+xshift[npoints],ykern+yshift[npoints]), kernel, zip(xgrid,ygrid), method='nearest', fill_value=0.0)

    return outdata
        
###########################################################################
