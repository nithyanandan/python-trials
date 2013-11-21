from __future__ import division
import numpy as NP

try:
    from scipy.spatial import cKDTree as KDT
except ImportError:
    from scipy.spatial import KDTree as KDT

###############################################################

class Point:
    def __init__(self, xyz=None):
        if xyz is None:
            self.x, self.y, self.z = 0.0, 0.0, 0.0
        if isinstance(xyz, basestring): # input is str
            raise TypeError('Data type mismatch. Check input again.')
        elif isinstance(xyz, (tuple,list)): # input is list or tuple
            if len(xyz) == 0:
                self.x, self.y, self.z = 0.0, 0.0, 0.0
            elif len(xyz) == 1:
                self.x, self.y, self.z = float(xyz[0]), 0.0, 0.0
            elif len(xyz) == 2:
                self.x, self.y, self.z = float(xyz[0]), float(xyz[1]), 0.0
            else:
                self.x, self.y, self.z = float(xyz[0]), float(xyz[1]), float(xyz[2])
        else: # input is scalar
            self.x, self.y, self.z = float(xyz), 0.0, 0.0

    def __str__(self):
        return '({0}, {1}, {2})'.format(self.x, self.y, self.z)

    def __add__(self, other):
        try:
            other
        except NameError:
            print 'No object provided for addition.'
            return self
        
        if not isinstance(other, Point):
            raise TypeError('Object type is incompatible.')

        return Point((self.x+other.x, self.y+other.y, self.z+other.z))
 
    def __sub__(self, other):
        try:
            other
        except NameError:
            print 'No object provided for subtraction.'
            return self
        
        if not isinstance(other, Point):
            raise TypeError('Object type is incompatible.')

        return Point((self.x-other.x, self.y-other.y, self.z-other.z))

    # def __mul__(self, value):
    #     return Point(value*self.x, value*self.y, value*self.z)

    def __mul__(self, other): # Dot product
        try:
            other
        except NameError:
            print 'No object provided for product.'
            return self
        
        if not isinstance(other, Point):
            raise TypeError('Object type is incompatible.')

        return self.x*other.x + self.y*other.y + self.z*other.z

    # __rmul__ = __mul__

    def __rmul__(self, value): # scalar product
        try:
            value
        except NameError:
            print 'No scalar value provided for product.'
            return self
        
        if not isinstance(value, (int,float)):
            raise TypeError('Value is incompatible for product.')

        return Point((value*self.x, value*self.y, value*self.z))

#####################################################################

def altaz2dircos(altaz, units=None):
    """
    -----------------------------------------------------------------------
    Convert altitude and azimuth to direction cosines

    Inputs:
    
    altaz:       Altitude and Azimuth as a list of tuples or Nx2 Numpy array
    
    Keyword Inputs:

    units:       [Default = 'radians'] Units of altitude and azimuth. Could
                 be radians or degrees. 

    Output:
    
    dircos:      Direction cosines corresponding to altitude and azimuth.
                 The first axis corresponds to local East, second to local
                 North and the third corresponds to zenith.
    -----------------------------------------------------------------------
    """

    try:
        altaz
    except NameError:
        print 'No altitude or azimuth specified. Check inputs.'

    if not isinstance(altaz, NP.ndarray):
        if not isinstance(altaz, list):
            if not isinstance(altaz, tuple):
                raise TypeError('altaz should be a 2-element tuple, a list of 2-element tuples, a 2-element list of scalars, or a Nx2 numpy array. Check inputs.')
            elif len(altaz) > 2:
                altaz = (altaz[0], altaz[1])
            elif len(altaz) < 2:
                raise TypeError('altaz should be a 2-element tuple, a list of 2-element tuples, a 2-element list of scalars, or a Nx2 numpy array. Check inputs.')
        elif isinstance(altaz[0], tuple):
            for elem in altaz:
                if not isinstance(elem, tuple):
                    raise TypeError('altaz should be a 2-element tuple, a list of 2-element tuples, a 2-element list of scalars, or a Nx2 numpy array. Check inputs.')
                if len(elem) < 2:
                    raise TypeError('altaz should be a 2-element tuple, a list of 2-element tuples, a 2-element list of scalars, or a Nx2 numpy array. Check inputs.')
            altaz = [((elem[0],elem[1]) if len(elem) > 2 else elem) for elem in altaz]
        elif len(altaz) != 2:
            raise TypeError('altaz should be a 2-element tuple, a list of 2-element tuples, a 2-element list of scalars, or a Nx2 numpy array. Check inputs.')
    else:
        if altaz.shape[1] < 2:
            raise TypeError('altaz should be a 2-element tuple, a list of 2-element tuples, a 2-element list of scalars, or a Nx2 numpy array. Check inputs.')
        elif altaz.shape[1] > 2:
            altaz = altaz[:,:2]

    altaz = NP.asarray(altaz)

    if units is None: units = 'radians'
    if units == 'degrees':
        altaz = NP.radians(altaz)

    if (NP.any(altaz[:,0] > NP.pi/2)) or (NP.any(altaz[:,0] < 0.0)):
        raise ValueError('Altitude(s) should lie between 0 and 90 degrees. Check inputs and units.')

    altaz[:,1] = NP.pi/2 - altaz[:,1] # Convert azimuth (measured from north towards east) to angle measured from east towards north
    l = NP.cos(altaz[:,0])*NP.cos(altaz[:,1]) # towards east
    m = NP.cos(altaz[:,0])*NP.sin(altaz[:,1]) # towards north
    n = NP.sin(altaz[:,0])                    # towards zenith
    return NP.asarray(zip(l,m,n))

#######################################################################

def dircos2altaz(dircos, units=None):
    """
    -----------------------------------------------------------------------
    Convert direction cosines to altitude and azimuth

    Inputs:
    
    dircos:      directin cosines as a list of tuples or Nx3 Numpy array
    
    Keyword Inputs:

    units:       [Default = 'radians'] Units of altitude and azimuth. Could
                 be radians or degrees. 

    Output:
    
    altaz:       Altitude and azimuth corresponding to direction cosines.
                 The first axis corresponds to local East, second to local
                 North and the third corresponds to zenith.
    -----------------------------------------------------------------------
    """

    try:
        dircos
    except NameError:
        print 'No direction cosines specified. Check inputs.'

    if not isinstance(dircos, NP.ndarray):
        if not isinstance(dircos, list):
            if not isinstance(dircos, tuple):
                raise TypeError('dircos should be a 3-element tuple, a list of 3-element tuples, a 3-element list of scalars, or a Nx3 numpy array. Check inputs.')
            elif len(dircos) > 3:
                dircos = (dircos[0], dircos[1], dircos[2])
            elif len(dircos) < 3:
                raise TypeError('dircos should be a 3-element tuple, a list of 3-element tuples, a 3-element list of scalars, or a Nx3 numpy array. Check inputs.')
        elif isinstance(dircos[0], tuple):
            for elem in dircos:
                if not isinstance(elem, tuple):
                    raise TypeError('dircos should be a 3-element tuple, a list of 3-element tuples, a 3-element list of scalars, or a Nx3 numpy array. Check inputs.')
                if len(elem) < 3:
                    raise TypeError('dircos should be a 3-element tuple, a list of 3-element tuples, a 3-element list of scalars, or a Nx3 numpy array. Check inputs.')
            dircos = [((elem[0],elem[1],elem[2]) if len(elem) > 3 else elem) for elem in dircos]
        elif len(dircos) != 3:
            raise TypeError('dircos should be a 3-element tuple, a list of 3-element tuples, a 3-element list of scalars, or a Nx3 numpy array. Check inputs.')
    else:
        if dircos.shape[1] < 3:
            raise TypeError('dircos should be a 3-element tuple, a list of 3-element tuples, a 3-element list of scalars, or a Nx3 numpy array. Check inputs.')
        elif dircos.shape[1] > 3:
            dircos = dircos[:,:3]

    dircos = NP.asarray(dircos)

    eps = 1.0e-10
    if ((NP.any(NP.abs(dircos[:,0]) > 1.0)) or
        (NP.any(NP.abs(dircos[:,1]) > 1.0)) or
        (NP.any(NP.abs(dircos[:,2]) > 1.0)) or
        (NP.any(NP.absolute(NP.sqrt(NP.sum(dircos**2,axis=1))-1.0) > eps))):

        raise ValueError('Individual components of direction cosines should lie between 0 and 1 and the direction cosines should have unit magnitudes. Check inputs.')

    altaz = NP.empty((dircos.shape[0],3))
    altaz[:,0] = NP.pi/2 - NP.arccos(dircos[:,2]) # Altitude/elevation
    altaz[:,1] = NP.pi/2 - NP.arctan2(dircos[:,1],dircos[:,0]) # Azimuth (measured from North)

    if units is None: units = 'radians'
    if units == 'degrees':
        altaz = NP.degrees(altaz)

    return altaz

#######################################################################

def hadec2altaz(hadec, latitude, units=None):
    """
    -----------------------------------------------------------------------
    Convert HA and declination to altitude and azimuth

    Inputs:
    
    hadec:       HA and declination as a list of tuples or Nx2 Numpy array
    
    latitude:    Latitude of the observatory. 

    Keyword Inputs:

    units:       [Default = 'radians'] Units of HA, dec and latitude. Could
                 be radians or degrees. 

    Output:
    
    altaz:       Altitude and azimuth corresponding to HA and dec at the 
                 given latitude. Units are identical to those in input.
    -----------------------------------------------------------------------
    """

    try:
        hadec
    except NameError:
        print 'No HA or declination specified. Check inputs.'

    try:
        latitude
    except NameError:
        print 'No latitude specified. Check inputs.'

    if isinstance(latitude, (list, tuple, str)):
        raise TypeError('Latitude should be a scalar number.')

    if not isinstance(hadec, NP.ndarray):
        if not isinstance(hadec, list):
            if not isinstance(hadec, tuple):
                raise TypeError('hadec should be a 2-element tuple, a list of 2-element tuples, a 2-element list of scalars, or a Nx2 numpy array. Check inputs.')
            elif len(hadec) > 2:
                hadec = (hadec[0], hadec[1])
            elif len(hadec) < 2:
                raise TypeError('hadec should be a 2-element tuple, a list of 2-element tuples, a 2-element list of scalars, or a Nx2 numpy array. Check inputs.')
        elif isinstance(hadec[0], tuple):
            for elem in hadec:
                if not isinstance(elem, tuple):
                    raise TypeError('hadec should be a 2-element tuple, a list of 2-element tuples, a 2-element list of scalars, or a Nx2 numpy array. Check inputs.')
                if len(elem) < 2:
                    raise TypeError('hadec should be a 2-element tuple, a list of 2-element tuples, a 2-element list of scalars, or a Nx2 numpy array. Check inputs.')
            hadec = [((elem[0],elem[1]) if len(elem) >= 2 else elem) for elem in hadec]
        elif len(hadec) != 2:
            raise TypeError('hadec should be a 2-element tuple, a list of 2-element tuples, a 2-element list of scalars, or a Nx2 numpy array. Check inputs.')
    else:
        if len(hadec.shape) < 2:
            hadec = NP.asarray(hadec).reshape(-1,hadec.shape[0])

        if hadec.shape[1] < 2:
            raise TypeError('hadec should be a 2-element tuple, a list of 2-element tuples, a 2-element list of scalars, or a Nx2 numpy array. Check inputs.')
        elif hadec.shape[1] > 2:
            hadec = hadec[:,0:2]

    hadec = NP.asarray(hadec)

    if units is None: units = 'radians'
    if units == 'degrees':
        hadec = NP.radians(hadec)
        latitude = NP.radians(latitude)

    if NP.any(NP.absolute(hadec[:,1]) > NP.pi/2):
        raise ValueError('Declination(s) should lie between -90 and 90 degrees. Check inputs and units.')

    if NP.absolute(latitude) > NP.pi/2:
        raise ValueError('Latitude should lie between -90 and 90 degrees. Check inputs and units.')

    altitude = NP.arcsin( NP.sin(hadec[:,1])*NP.sin(latitude) + NP.cos(hadec[:,1])*NP.cos(latitude)*NP.cos(hadec[:,0]) )
    azimuth = NP.arccos( (NP.sin(hadec[:,1])-NP.sin(altitude)*NP.sin(latitude))/(NP.cos(altitude)*NP.cos(latitude)) )

    # Need to make sure the values are in the conventional range
    azimuth = NP.where(NP.sin(hadec[:,0])<0.0, azimuth, 2.0*NP.pi-azimuth)
    if units == 'degrees':
        altitude = NP.degrees(altitude)
        azimuth = NP.degrees(azimuth)
    return NP.asarray(zip(altitude, azimuth))

#########################################################################

def altaz2hadec(altaz, latitude, units=None):
    """
    -----------------------------------------------------------------------
    Convert altitude and azimuth to HA and declination
    Same transformation function as hadec2altaz, replace azimuth with HA
        and altitude with declination.

    Inputs:
    
    altaz:       Altitude and azimtuh as a list of tuples or Nx2 Numpy array
    
    latitude:    Latitude of the observatory. 

    Keyword Inputs:

    units:       [Default = 'radians'] Units of HA, dec and latitude. Could
                 be radians or degrees. 

    Output:
    
    hadec:       HA and declination corresponding to altitude and azimuth
                 given latitude. Units are identical to those in input.
    -----------------------------------------------------------------------
    """

    try:
        altaz
    except NameError:
        print 'No altitude or azimuth specified. Check inputs.'

    try:
        latitude
    except NameError:
        print 'No latitude specified. Check inputs.'

    if isinstance(latitude, (list, tuple, str)):
        raise TypeError('Latitude should be a scalar number.')

    if not isinstance(altaz, NP.ndarray):
        if not isinstance(altaz, list):
            if not isinstance(altaz, tuple):
                raise TypeError('altaz should be a 2-element tuple, a list of 2-element tuples, a 2-element list of scalars, or a Nx2 numpy array. Check inputs.')
            elif len(altaz) > 2:
                altaz = (altaz[0], altaz[1])
            elif len(altaz) < 2:
                raise TypeError('altaz should be a 2-element tuple, a list of 2-element tuples, a 2-element list of scalars, or a Nx2 numpy array. Check inputs.')
        elif isinstance(altaz[0], tuple):
            for elem in altaz:
                if not isinstance(elem, tuple):
                    raise TypeError('altaz should be a 2-element tuple, a list of 2-element tuples, a 2-element list of scalars, or a Nx2 numpy array. Check inputs.')
                if len(elem) < 2:
                    raise TypeError('altaz should be a 2-element tuple, a list of 2-element tuples, a 2-element list of scalars, or a Nx2 numpy array. Check inputs.')
            altaz = [((elem[0],elem[1]) if len(elem) >= 2 else elem) for elem in altaz]
        elif len(altaz) != 2:
            raise TypeError('altaz should be a 2-element tuple, a list of 2-element tuples, a 2-element list of scalars, or a Nx2 numpy array. Check inputs.')
    else:
        if altaz.shape[1] < 2:
            raise TypeError('altaz should be a 2-element tuple, a list of 2-element tuples, a 2-element list of scalars, or a Nx2 numpy array. Check inputs.')
        elif altaz.shape > 2:
            altaz = altaz[:,0:2]

    altaz = NP.asarray(altaz)

    if units is None: units = 'radians'
    if units == 'degrees':
        altaz = NP.radians(altaz)
        latitude = NP.radians(latitude)

    if (NP.any(altaz[:,0] > NP.pi/2)) or (NP.any(altaz[:,0] < 0.0)):
        raise ValueError('Altitude(s) should lie between 0 and 90 degrees. Check inputs and units.')

    if NP.absolute(latitude) > NP.pi/2:
        raise ValueError('Latitude should lie between -90 and 90 degrees. Check inputs and units.')

    dec = NP.arcsin( NP.sin(altaz[:,0])*NP.sin(latitude) + NP.cos(altaz[:,1])*NP.cos(latitude)*NP.cos(altaz[:,0]) )
    ha = NP.arccos( (NP.sin(altaz[:,0])-NP.sin(dec)*NP.sin(latitude))/(NP.cos(dec)*NP.cos(latitude)) )
 
    if units == 'degrees':
        ha *= 180.0/NP.pi
        dec *= 180.0/NP.pi
    return NP.asarray(zip(ha, dec))

#########################################################################

def enu2xyz(enu, latitude, units='radians'):
    """
    -----------------------------------------------------------------------
    Convert local ENU coordinates to equatorial XYZ coordinates.

    Inputs:
    
    enu:         local ENU coordinates as a list of tuples or Nx3 Numpy array
    
    latitude:    Latitude of the observatory. 

    Keyword Inputs:

    units:       [Default = 'radians'] Units of latitude. Could
                 be radians or degrees. 

    Output:
    
    xyz:         Equatorial XYZ coordinates corresponding to local ENU coordinates
                 given latitude. Units are identical to those in input.
    -----------------------------------------------------------------------
    """

    try:
        enu
    except NameError:
        print 'No baselines specified. Check inputs.'

    try:
        latitude
    except NameError:
        print 'No latitude specified. Check inputs.'

    if isinstance(latitude, (list, tuple, str)):
        raise TypeError('Latitude should be a scalar number.')

    if not isinstance(enu, NP.ndarray):
        if not isinstance(enu, list):
            if not isinstance(enu, tuple):
                raise TypeError('enu should be a 3-element tuple, a list of 3-element tuples, a 3-element list of scalars, or a Nx3 numpy array. Check inputs.')
            elif len(enu) > 3:
                enu = (enu[0], enu[1], enu[2])
            elif len(enu) < 3:
                raise TypeError('enu should be a 3-element tuple, a list of 3-element tuples, a 3-element list of scalars, or a Nx3 numpy array. Check inputs.')
        elif isinstance(enu[0], tuple):
            for elem in enu:
                if not isinstance(elem, tuple):
                    raise TypeError('enu should be a 3-element tuple, a list of 3-element tuples, a 3-element list of scalars, or a Nx3 numpy array. Check inputs.')
                if len(elem) < 3:
                    raise TypeError('enu should be a 3-element tuple, a list of 3-element tuples, a 3-element list of scalars, or a Nx3 numpy array. Check inputs.')
            enu = [((elem[0],elem[1]) if len(elem) >= 3 else elem) for elem in enu]
        elif len(enu) != 3:
            raise TypeError('enu should be a 3-element tuple, a list of 3-element tuples, a 3-element list of scalars, or a Nx3 numpy array. Check inputs.')
    else:
        if enu.shape[1] < 3:
            raise TypeError('enu should be a 3-element tuple, a list of 3-element tuples, a 3-element list of scalars, or a Nx3 numpy array. Check inputs.')
        elif enu.shape > 3:
            enu = enu[:,0:3]

    enu = NP.asarray(enu)

    if units == 'degrees':
        latitude = NP.radians(latitude)

    rotation_matrix = NP.asarray([[0.0, -NP.sin(latitude), NP.cos(latitude)],
                                  [1.0, 0.0,               0.0             ],
                                  [0.0, NP.cos(latitude),  NP.sin(latitude)]])

    xyz = NP.dot(enu, rotation_matrix.T)

    return xyz

#########################################################################

def xyz2enu(xyz, latitude, units='radians'):
    """
    -----------------------------------------------------------------------
    Convert equatorial XYZ coordinates to local ENU coordinates.

    Inputs:
    
    xyz:         equatorial XYZ coordinates as a list of tuples or Nx3 Numpy array
    
    latitude:    Latitude of the observatory. 

    Keyword Inputs:

    units:       [Default = 'radians'] Units of latitude. Could
                 be radians or degrees. 

    Output:
    
    enu:         local ENU coordinates corresponding to equatorial XYZ coordinates
                 given latitude. Units are identical to those in input.
    -----------------------------------------------------------------------
    """

    try:
        xyz
    except NameError:
        print 'No baselines specified. Check inputs.'

    try:
        latitude
    except NameError:
        print 'No latitude specified. Check inputs.'

    if isinstance(latitude, (list, tuple, str)):
        raise TypeError('Latitude should be a scalar number.')

    if not isinstance(xyz, NP.ndarray):
        if not isinstance(xyz, list):
            if not isinstance(xyz, tuple):
                raise TypeError('xyz should be a 3-element tuple, a list of 3-element tuples, a 3-element list of scalars, or a Nx3 numpy array. Check inputs.')
            elif len(xyz) > 3:
                xyz = (xyz[0], xyz[1], xyz[2])
            elif len(xyz) < 3:
                raise TypeError('xyz should be a 3-element tuple, a list of 3-element tuples, a 3-element list of scalars, or a Nx3 numpy array. Check inputs.')
        elif isinstance(xyz[0], tuple):
            for elem in xyz:
                if not isinstance(elem, tuple):
                    raise TypeError('xyz should be a 3-element tuple, a list of 3-element tuples, a 3-element list of scalars, or a Nx3 numpy array. Check inputs.')
                if len(elem) < 3:
                    raise TypeError('xyz should be a 3-element tuple, a list of 3-element tuples, a 3-element list of scalars, or a Nx3 numpy array. Check inputs.')
            xyz = [((elem[0],elem[1]) if len(elem) >= 3 else elem) for elem in xyz]
        elif len(xyz) != 3:
            raise TypeError('xyz should be a 3-element tuple, a list of 3-element tuples, a 3-element list of scalars, or a Nx3 numpy array. Check inputs.')
    else:
        if xyz.shape[1] < 3:
            raise TypeError('xyz should be a 3-element tuple, a list of 3-element tuples, a 3-element list of scalars, or a Nx3 numpy array. Check inputs.')
        elif xyz.shape > 3:
            xyz = xyz[:,0:3]

    xyz = NP.asarray(xyz)

    if units == 'degrees':
        latitude = NP.radians(latitude)

    rotation_matrix = NP.asarray([[0.0               , 1.0, 0.0],
                                  [-NP.sin(-latitude), 0.0, NP.cos(-latitude)],
                                  [NP.cos(-latitude) , 0.0, NP.sin(-latitude)]])

    enu = NP.dot(xyz, rotation_matrix.T)

    return enu

#########################################################################

# def angular_ring(skypos, angles, npoints=100, skyunits='radec', angleunits='degrees'):

#     skypos = NP.radians(skypos).reshape(-1,2)
#     angles = NP.asarray(angles).reshape(1,-1)
#     if angleunits == 'degrees':
#         angles = NP.radians(angles)

#     azangle = 2.0 * NP.pi * NP.arange(npoints)/npoints 

#     radec = NP.empty((npoints, 2, skypos.shape[0]))

#     for i in range(skypos.shape[0]):
#         local_enu = NP.hstack((NP.sin(angles[i])*NP.cos(azangle).reshape(-1,1), NP.sin(angles[i])*NP.sin(azangle).reshape(-1,1), NP.cos(angles[i])*NP.ones((npoints,1))))
#         rotation_matrix = NP.asarray([[0.0, -NP.sin(skypos[i,1]), NP.cos(skypos[i,1])],
#                                       [1.0, 0.0,                  0.0              ],
#                                       [0.0, NP.cos(skypos[i,1]),  NP.sin(skypos[i,1])]])
#         xyz = NP.dot(local_enu, rotation_matrix.T)
#         ring_dec = NP.degrees(NP.arcsin(xyz[:,2]))
#         ring_az = NP.degrees(NP.arctan(xyz[:,1], xyz[:,0]))
#         ring_ra = NP.degrees(skypos[:,0]) + ring_az

#         radec[:,0,i] = ring_ra
#         radec[:,1,i] = ring_dec

#         return radec

#     # local_east = NP.repeat(NP.cos(azangle).reshape(npoints,-1), len(angles), axis=1) * NP.repeat(NP.sin(angles).reshape(-1,len(angles)), npoints, axis=0)
#     # local_north = NP.repeat(NP.sin(azangle).reshape(npoints,-1), len(angles), axis=1) * NP.repeat(NP.sin(angles).reshape(-1,len(angles)), npoints, axis=0)
#     # local_up = NP.repeat(NP.cos(angles).reshape(-1,len(angles)), npoints, axis=0) 

#     # Rotate the "north" and "up" to equatorial coordinates by the declination provided



    

    

#########################################################################

def sph2xyz(lon, lat, rad=None):
    """
    --------------------------------------------------------------------
    Inputs:

    lon [scalar or vector] longitude in degrees.  

    lat [scalar or vector] latitude in degrees. Same size as lon.

    rad [Optional. scalar or vector] radius. Same size as lon and lat.
        Default = 1.0

    Outputs:

    x   [scalar or vector] x-coordinates. Same size as lon and lat

    y   [scalar or vector] y-coordinates. Same size as lon and lat

    z   [scalar or vector] z-coordinates. Same size as lon and lat
    --------------------------------------------------------------------
    """

    try:
        lon, lat
    except NameError:
        raise NameError('lon and/or lat not defined in sph2xyz().')

    lonr = NP.radians(lon)
    latr = NP.radians(lat)

    if lonr.size != latr.size:
        raise ValueError('lon and lat should have same size in sph2xyz().')

    if rad is None:
        rad = NP.ones(lonr.size)
    else:
        rad = NP.asarray(rad)

    if lonr.size != rad.size:
        raise ValueError('rad must have same size as lon and lat in sph2xyz().')

    x = rad * NP.cos(lonr) * NP.cos(latr)
    y = rad * NP.sin(lonr) * NP.cos(latr)
    z = rad * NP.sin(latr)

    return x, y, z

#########################################################################

def sphdist(lon1, lat1, lon2, lat2):
    """
    --------------------------------------------------------------------
    Returns great circle distance.  

    Uses vicenty distance formula - a bit slower than others, but
    numerically stable.

    Inputs: 

    lon1 [scalar or vector] Longtitude in first list in degrees

    lat1 [scalar or vector] Latitude in first list in degrees.
         Must be of same size as lon1 if vector

    lon2 [scalar or vector] Longtitude in second list in degrees.
         Must be of same size as lon1 if vector

    lat2 [scalar or vector] Latitude in second list in degrees.
         Must be of same size as lon1 if vector

    Outputs:

    Angular distance (in degrees) subtended on the great circle between
    the given set of points. Same size as the inputs.
    --------------------------------------------------------------------
    """

    try:
        lon1, lat1, lon2, lat2
    except NameError:
        raise NameError('At least one of lon1, lat1, lon2, lat2 undefined in sphdist().')

    lon1 = NP.asarray(lon1)
    lat1 = NP.asarray(lat1)

    if lon1.size != lat1.size:
        raise ValueError('lon1 and lat1 must be of same size in sphdist().')

    lon2 = NP.asarray(lon2)
    lat2 = NP.asarray(lat2)

    if lon2.size != lat2.size:
        raise ValueError('lon2 and lat2 must be of same size in sphdist().')

    if lon1.size != lon2.size:
        if lon1.size == 1:
            lon1 = lon1 * NP.ones(lon2.size)
            lat1 = lat1 * NP.ones(lat2.size)
        elif lon2.size == 1:
            lon2 = lon2 * NP.ones(lon1.size)
            lat2 = lat2 * NP.ones(lat1.size)
        else:
            raise ValueError('first and second sets of coordinates must have same size.')

    # terminology from the Vicenty formula - lambda and phi and
    # "standpoint" and "forepoint"

    lambs = NP.radians(lon1)
    phis = NP.radians(lat1)
    lambf = NP.radians(lon2)
    phif = NP.radians(lat2)

    dlamb = lambf - lambs

    numera = NP.cos(phif) * NP.sin(dlamb)
    numerb = NP.cos(phis) * NP.sin(phif) - NP.sin(phis) * NP.cos(phif) * NP.cos(dlamb)
    numer = NP.hypot(numera, numerb)
    denom = NP.sin(phis) * NP.sin(phif) + NP.cos(phis) * NP.cos(phif) * NP.cos(dlamb)

    return NP.degrees(NP.arctan2(numer, denom))

#########################################################################

def spherematch(lon1, lat1, lon2=None, lat2=None, matchrad=None,
                nnearest=0, maxmatches=-1):
    """
    ---------------------------------------------------------------------------
    Finds matches in one catalog to another. 

    Parameters
    lon1 : array-like
         Longitude-like (RA, etc.) coordinate in degrees of the first catalog
    lat1 : array-like
        Latitude-like (Dec, etc.) coordinate in degrees of the first catalog (shape of array must match `lon1`)
    lon2 : array-like
        Latitude-like (RA, etc.) coordinate in degrees of the second catalog
    lat2 : array-like
        Latitude-like (Dec, etc.) in degrees of the second catalog (shape of array must match `ra2`)
    matchrad : float or None, optional
        How close (in degrees) a match has to be to count as a match.  If None,
        all nearest neighbors for the first catalog will be returned. 
    nnearest : int, optional
        The nth neighbor to find.  
    maxmatches : int, optional
        Maximum number of matches to find. If maxmatches > 0, the code finds 
        all matches up to maxmatches satisfying matchrad and nnearest is
        ignored.

    Returns
    -------
    m1 : int array
        Indices into the first catalog of the matches. 
    m2 : int array
        Indices into the second catalog of the matches. 
    d12 : float array
        Distance (in degrees) between the matches 
    -------------------------------------------------------
    """

    try:
        lon1, lat1
    except NameError:
        raise NameError('lon1 and/or lat1 not defined. Aborting spherematch()')
    
    if not isinstance(lon1, (list, NP.ndarray)):
        lon1 = NP.asarray(lon1)
    if not isinstance(lat1, (list, NP.ndarray)):
        lat1 = NP.asarray(lat1)

    try:
        nnearest = int(nnearest)
    except TypeError:
        raise TypeError('nnearest should be a non-negative integer.')

    try:
        maxmatches = int(maxmatches)
    except TypeError:
        raise TypeError('maxmatches should be a non-negative integer.')

    if matchrad is None:
        if maxmatches > 0:
            nnearest = 0
            print 'No matchrad specified. Will determine all the {0} nearest neighbours.'.format(maxmatches)
        else:
            maxmatches = -1
            if nnearest <= 0:
                nnearest = 1
            print 'No matchrad specified. Will determine the nearest neighbour # {0}.'.format(nnearest)
    elif not isinstance(matchrad, (int,float)):
        raise TypeError('matchrad should be a scalar number.')
    elif matchrad > 0.0:
        matchrad_cartesian = 2.0*NP.sin(0.5*matchrad*NP.pi/180.0)
        if maxmatches >= 0:
            nnearest = 0
        else:
            if nnearest <= 0:
                nnearest = 1
            print 'maxmatches is negative. Will determine the nearest neighbour # {0}.'.format(nnearest)            
    else:
        raise ValueError('matchrad is not positive.')

    self_match = False
    if (lon2 is None) and (lat2 is None):
        self_match = True
    elif lon2 is None:
        lon2 = lon1
    elif lat2 is None:
        lat2 = lat2

    if lon1.shape != lat1.shape:
        raise ValueError('lon1 and lat1 should be of same length')
    if not self_match:
        if lon2.shape != lat2.shape:
            raise ValueError('lon2 and lat2 should be of same length')
    
    if lon1.size == 1:
        lon1 = NP.asarray(lon1).reshape(1)
        lat1 = NP.asarray(lat1).reshape(1)

    if lon2.size == 1:
        lon2 = NP.asarray(lon2).reshape(1)
        lat2 = NP.asarray(lat2).reshape(1)

    x1, y1, z1 = sph2xyz(lon1.ravel(), lat1.ravel())

    # this is equivalent to, but faster than just doing NP.array([x1, y1, z1])
    coords1 = NP.empty((x1.size, 3))
    coords1[:, 0] = x1
    coords1[:, 1] = y1
    coords1[:, 2] = z1

    if (lon2 is not None) and (lat2 is not None):
        x2, y2, z2 = sph2xyz(lon2.ravel(), lat2.ravel())

        # this is equivalent to, but faster than just doing NP.array([x1, y1, z1])
        coords2 = NP.empty((x2.size, 3))
        coords2[:, 0] = x2
        coords2[:, 1] = y2
        coords2[:, 2] = z2
    
    if maxmatches == 0:
        kdt1 = KDT(coords1)
        if not self_match:
            kdt2 = KDT(coords2)
            ngbr_of_first_in_second = kdt1.query_ball_tree(kdt2, matchrad_cartesian)
            m1 = [i for i in xrange(len(ngbr_of_first_in_second)) for j in ngbr_of_first_in_second[i] if ngbr_of_first_in_second[i] != []]
            m2 = [j for i in xrange(len(ngbr_of_first_in_second)) for j in ngbr_of_first_in_second[i] if ngbr_of_first_in_second[i] != []]
            d12 = sphdist(lon1[m1], lat1[m1], lon2[m2], lat2[m2])
        else:
            ngbr_of_first_in_itself = kdt1.query_ball_tree(kdt1, matchrad_cartesian)
            m1 = [i for i in xrange(len(ngbr_of_first_in_itself)) for j in ngbr_of_first_in_itself[i] if ngbr_of_first_in_itself[i] != [] and i != j]
            m2 = [j for i in xrange(len(ngbr_of_first_in_itself)) for j in ngbr_of_first_in_itself[i] if ngbr_of_first_in_itself[i] != [] and i != j]
            d12 = sphdist(lon1[m1], lat1[m1], lon1[m2], lat1[m2])
    else:
        if not self_match:
            kdt2 = KDT(coords2)
            if matchrad is None:
                m2 = kdt2.query(coords1, max(maxmatches,nnearest))[1]
            else:
                m2 = kdt2.query(coords1, max(maxmatches,nnearest), distance_upper_bound=matchrad_cartesian)[1]
        else:
            kdt1 = KDT(coords1)
            if matchrad is None:
                m2 = kdt1.query(coords1, max(maxmatches,nnearest)+1)[1]
            else:
                m2 = kdt2.query(coords1, max(maxmatches,nnearest)+1, distance_upper_bound=matchrad_cartesian)[1]

        if nnearest > 0:
            m1 = NP.arange(lon1.size)
            if nnearest > 1:
                m2 = m2[:,-1]
        else: 
            m1 = NP.repeat(NP.arange(lon1.size).reshape(lon1.size,1), maxmatches, axis=1).flatten()
            m2 = m2[:,-maxmatches:].flatten() # Extract the last maxmatches columns

        if not self_match:
            msk = m2 < lon2.size
        else:
            msk = m1 < lon1.size
        m1 = m1[msk]
        m2 = m2[msk]

        if not self_match:
            d12 = sphdist(lon1[m1], lat1[m1], lon2[m2], lat2[m2])
        else:
            d12 = sphdist(lon1[m1], lat1[m1], lon1[m2], lat1[m2])

        if matchrad is not None:
            msk = d12 <= matchrad
            m1 = m1[msk]
            m2 = m2[msk]
            d12 = d12[msk]

    return m1, m2, d12        

#########################################################################

def baseline_generator(antenna_locations, auto=False, conjugate=False):
    """
    -------------------------------------------------------------------
    Generate baseline from antenna locations.

    Inputs:

    antenna_locations: List of tuples containing antenna coordinates, 
                       or list of instances of class Point containing
                       antenna coordinates, or Numpy array (Nx3) array
                       with each row specifying an antenna location.

    Input keywords:

    auto:              [Default=False] If True, compute zero spacings of
                       antennas with themselves.

    conjugate:         [Default=False] If True, compute conjugate 
                       baselines.

    Output:

    baseline_locations: Baseline locations in the same data type as 
                        antenna locations (list of tuples, list of 
                        instances of class Point or Numpy array of size
                        Nb x 3 with each row specifying one baseline 
                        vector)

    -------------------------------------------------------------------
    """

    try:
        antenna_locations
    except NameError:
        print 'No antenna locations supplied. Returning from baseline_generator()'
        return None

    inp_type = 'tbd'

    if not isinstance(antenna_locations, NP.ndarray):
        if isinstance(antenna_locations, list):
            if isinstance(antenna_locations[0], Point):
                inp_type = 'loo' # list of objects
            elif isinstance(antenna_locations[0], tuple):
                inp_type = 'lot' # list of tuples
                antenna_locations = [(tuple(loc) if len(loc) == 3 else (tuple([loc[0],0.0,0.0]) if len(loc) == 1 else (tuple([loc[0],loc[1],0.0]) if len(loc) == 2 else (tuple([loc[0],loc[1],loc[2]]))))) for loc in antenna_locations if len(loc) != 0] # Remove empty tuples and validate the data range and data type for antenna locations. Force it to have three components for every antenna location.
        elif isinstance(antenna_locations, Point):
            if not auto:
                print 'No non-zero spacings found since auto=False.'
                return None
            else:
                return Point()
        elif isinstance(antenna_locations, tuple):
            if not auto:
                print 'No non-zero spacings found since auto=False.'
                return None
            else:
                return (0.0,0.0,0.0)
        else:
            if not auto:
                print 'No non-zero spacings found since auto=False.'
                return None
            else:
                return (0.0,0.0,0.0)
    else:
        inp_type = 'npa' # A numpy array
        if antenna_locations.shape[0] == 1:
            if not auto:
                print 'No non-zero spacings found since auto=False.'
                return None
            else:
                return NP.zeros(1,3)
        else:
            if antenna_locations.shape[1] > 3:
                antenna_locations = antenna_locations[:,:3]
            elif antenna_locations.shape[1] < 3:
                antenna_locations = NP.hstack((antenna_locations, NP.zeros(antenna_locations.shape[0],3-antenna_locations.shape[1])))

    if isinstance(antenna_locations, list):
        num_ants = len(antenna_locations)
    else:
        num_ants = antenna_locations.shape[0]

    if inp_type == 'loo':
        if auto:
            baseline_locations = [antenna_locations[j]-antenna_locations[i] for i in xrange(0,num_ants) for j in xrange(0,num_ants) if j >= i]
        else:
            baseline_locations = [antenna_locations[j]-antenna_locations[i] for i in range(0,num_ants) for j in range(0,num_ants) if j > i]                
        if conjugate:
            baseline_locations += [antenna_locations[j]-antenna_locations[i] for i in xrange(0,num_ants) for j in xrange(0,num_ants) if j < i]
    elif inp_type == 'lot':
        if auto:
            baseline_locations = [tuple((antenna_locations[j][0]-antenna_locations[i][0], antenna_locations[j][1]-antenna_locations[i][1], antenna_locations[j][2]-antenna_locations[i][2])) for i in xrange(0,num_ants) for j in xrange(0,num_ants) if j >= i]
        else:
            baseline_locations = [tuple((antenna_locations[j][0]-antenna_locations[i][0], antenna_locations[j][1]-antenna_locations[i][1], antenna_locations[j][2]-antenna_locations[i][2])) for i in xrange(0,num_ants) for j in xrange(0,num_ants) if j > i]
        if conjugate:
            baseline_locations += [tuple((antenna_locations[j][0]-antenna_locations[i][0], antenna_locations[j][1]-antenna_locations[i][1], antenna_locations[j][2]-antenna_locations[i][2])) for i in xrange(0,num_ants) for j in xrange(0,num_ants) if j < i]
    elif inp_type == 'npa':
        if auto:
            baseline_locations = [antenna_locations[i,:]-antenna_locations[j,:] for i in xrange(0,num_ants) for j in xrange(0,num_ants) if j >= i]
        else:
            baseline_locations = [antenna_locations[i,:]-antenna_locations[j,:] for i in xrange(0,num_ants) for j in xrange(0,num_ants) if j > i]        
        if conjugate:
            baseline_locations += [antenna_locations[i,:]-antenna_locations[j,:] for i in xrange(0,num_ants) for j in xrange(0,num_ants) if j < i]        
        baseline_locations = NP.asarray(baseline_locations)

    return baseline_locations

##########################################################################


