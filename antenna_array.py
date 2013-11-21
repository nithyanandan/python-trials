import numpy as NP
import my_DSP_modules as DSP
import geometry as GEOM
import my_gridding_modules as GRD

#####################################################################  

# def FT(inp, ax=0, use_real=False, shift=False):
    
#     try:
#         inp
#     except NameError:
#         raise NameError('inp not defined. Aborting FT().')

#     return DSP.FT1D(inp, ax=0, use_real=False, shift=False)

#####################################################################  

def XC(inp1, inp2=None):

    try:
        inp1
    except NameError:
        raise NameError('inp1 not defined. Aborting XC().')

    if not isinstance(inp1, (list, tuple, NP.ndarray, int, float, complex)):
        raise TypeError('inp1 is of the wrong data type. Check inputs again. Aborting XC().')

    inp1 = NP.asarray(inp1)

    if inp2 is None:
        inp2 = inp1
    elif not isinstance(inp2, (list, tuple, int, float, complex, NP.ndarray)):
        raise TypeError('inp2 has incompatible data type. Verify inputs. Aborting XC().')

    inp2 = NP.asarray(inp1)

    zero_pad_length = 2**NP.ceil(NP.log2(len(inp1)+len(inp2)-1))-(len(inp1)+len(inp2)-1)

    return NP.roll(NP.append(NP.correlate(inp1, inp2, mode='full'), NP.zeros(zero_pad_length)), -(len(inp2)-1))   # zero pad and shift to ensure identical results as FX() operation

#####################################################################  

def spaxis(length, resolution=1.0, use_real=False, shift=False):
    
    try:
        length
    except NameError:
        raise NameError('Input length not defined. Aborting spaxis().')
        
    if not isinstance(resolution, (int, float)):
        raise TypeError('Input resolution must be a positive scalar integer or floating point number. Aborting spaxis().')
    elif resolution < 0.0:
        raise ValueError('Input resolution must be positive. Aborting spaxis().')

    return DSP.spectral_axis(length, delx=resolution, use_real=False, shift=False)

#####################################################################  

class PolInfo:

    def __str__(self):
        return ' Instance of class "{0}" in module "{1}" \n type: {2} \n flag (X): {3} \n flag (Y): {4} '.format(self.__class__.__name__, self.__module__, self.flag_x, self.flag_y)

    def __init__(self, Ext, Eyt, flag_x=False, flag_y=False, pol_type='Linear'):
        self.Ext = NP.asarray(Ext)
        self.Eyt = NP.asarray(Eyt)
        self.flag_y = flag_y
        self.flag_y = flag_y
        self.pol_type = pol_type
        self.Exf, self.Eyf = self.F()

    def F(self, pol=None):
        if pol is None:
            self.Exf = DSP.FT1D(self.Ext, ax=0, use_real=False, shift=False)
            self.Eyf = DSP.FT1D(self.Eyt, ax=0, use_real=False, shift=False)
        elif pol == 'X' or pol == 'Y' or pol == 'x' or pol == 'y':
            if pol ==' X' or pol == 'x':
                self.Exf = DSP.FT1D(self.Ext, ax=0, use_real=False, shift=False)
            else:
                self.Eyf = DSP.FT1D(self.Eyt, ax=0, use_real=False, shift=False)
        else:
            raise ValueError('Polarization string unrecognized. Verify inputs. Aborting PolInfo.F()')

    def update(self, Ext=None, Eyt=None, flag_x=False, flag_y=False, pol_type='Linear'):
        
        if not Ext is None:
            self.Ext = NP.asarray(Ext)
            self.F(pol='X')

        if not Eyt is None:
            self.Eyt = NP.asarray(Eyt)
            self.F(pol='Y')

        if not flag_x is None: self.flag_x = flag_x
        if not flag_y is None: self.flag_y = flag_y
        if not pol_type is None: self.pol_type = pol_type


#####################################################################

# class Antenna_old:

#     def __str__(self):
#         return ' Instance of class "{0}" in module "{1}" \n label: {2} \n location: {3}'.format(self.__class__.__name__, self.__module__, self.label, self.pos.__str__())

#     def __init__(self, label, pos, Et, t, wtspos=None, gridwts=None):
#         self.label = label
#         # self.pos = GEOM.Point(pos.x, pos.y, pos.z)
#         self.pos = pos
#         self.Et = NP.asarray(Et)
#         self.t = NP.asarray(t)
#         self.Ef = FT(self.Et, ax=0, use_real=False, shift=False)
#         self.f = spaxis(len(self.Et), resolution=self.t[1]-self.t[0])
#         if (gridwts != None) and (wtspos != None):
#             if len(wtspos) == gridwts:
#                 self.gridwts = NP.asarray(gridwts)
#                 self.wtspos = wtspos
#                 # x_max, y_max, z_max = max(x_max, self.pos.x + max(wtspos.x)), max(y_max, self.pos.y + max(wtspos.y)), max(z_max, self.pos.z + max(wtspos.z))
#                 # x_min, y_min, z_min = min(x_min, self.pos.x - min(wtspos.x)), min(y_min, self.pos.y - min(wtspos.y)), min(z_min, self.pos.z - min(wtspos.z))
#         else:
#             print 'Dimension mismatch of gridding weights and location of the weights or gridding information not provided. Gridding information not initialized for antenna "{0}"'.format(self.label)
#             print 'Proceeding without grid information.'


#     def F(self):
#         return FT(self.Et, ax=0, use_real=False, shift=False)

#     def channels(self):
#         return spaxis(len(self.Et), resolution=self.t[1]-self.t[0])

# #####################################################################  

class Antenna:

    def __str__(self):
        return ' Instance of class "{0}" in module "{1}" \n label: {2} \n location: {3}'.format(self.__class__.__name__, self.__module__, self.label, self.pos.__str__())

    def __init__(self, label, pos, Ext, Eyt, t, timestamp, wtspos_x=None, wtspos_y=None, wts_x=None, wts_y=None, flag_x=False, flag_y=False, pol_type='Linear'):

        self.label = label

        if isinstance(pos, (tuple, GEOM.Point)):
            self.pos = pos
        else:
            raise TypeError('Antenna position must be a 3-element tuple or an instance of GEOM.Point')

        if (len(t) != len(Ext)) or (len(t) != len(Eyt)):
            raise IndexError('The size of voltage streams do not match that of time sequence.')

        self.pol = PolInfo(Ext, Eyt, flag_x, flag_y, pol_type)
        self.t = NP.asarray(t)
        self.timestamp = timestamp
        self.f = self.channels()

        self.gridpos_x = None
        self.gridpos_y = None
        self.gridwts_x = None
        self.gridwts_y = None

        if (wts_x != None) and (wtspos_x != None):
            if len(wtspos_x) == len(wts_x):
                self.wts_x = NP.asarray(wts_x)
                self.wtspos_x = NP.asarray(wtspos_x)
        else:
            print 'Dimension mismatch of gridding weights and location of the weights or gridding information not provided. Gridding information not initialized for polarization X in antenna {0}'.format(self.label)

        if (wts_y != None) and (wtspos_y != None):
            if len(wtspos_y) == len(wts_y):
                self.wts_y = NP.asarray(wts_y)
                self.wtspos_y = NP.asarray(wtspos_y)
        else:
            print 'Dimension mismatch of gridding weights and location of the weights or gridding information not provided. Gridding information not initialized for polarization Y in antenna {0}'.format(self.label)

    def channels(self):
        return spaxis(len(self.t), self.t[1]-self.t[0])

    def update(self, Ext=None, Eyt=None, t=None, timestamp=None, label=None, pos=None, wtspos_x=None, wts_x=None, wtspos_y=None, wts_y=None, flag_x=None, flag_y=None):

        if not label is None: self.label = label
        if not pos is None: self.pos = pos
        if not timestamp is None: self.timestamp = timestamp

        if not t is None:
            self.t = t
            self.f = self.channels()           

        if (not flag_x is None) or (not flag_y is None) or (not Ext is None) or (not Eyt is None):
            self.pol.update(Ext=Ext, Eyt=Eyt, wtspos_x=wtspos_x, wts_x=wts_x, wtspos_y=wtspos_y, wts_y=wts_y, flag_x=flag_x, flag_y=flag_y)

        if (not wts_x is None) and (not wtspos_x is None):
            if len(wtspos_x) == len(wts_x):
                self.wts_x = NP.asarray(wts_x)
                self.wtspos_x = NP.asarray(wtspos_x)
        else:
            print 'Dimension mismatch of gridding weights and location of the weights or gridding information not provided. Gridding information not updated for polarization X in antenna {0}'.format(self.label)

        if (not wts_y is None) and (not wtspos_y is None):
            if len(wtspos_y) == len(wts_y):
                self.wts_y = NP.asarray(wts_y)
                self.wtspos_y = NP.asarray(wtspos_y)
        else:
            print 'Dimension mismatch of gridding weights and location of the weights or gridding information not provided. Gridding information not updated for polarization Y in antenna {0}'.format(self.label)
       

#####################################################################  

class AntennaPair(Antenna):
    def __init__(self, A1, A2):
        self.A1, self.A2 = A1, A2
        self.label = A1.label+'-'+A2.label
        self.pos = A1.pos-A2.pos
        self.Et = XC(self.A1.Et, self.A2.Et)
        # self.t = (A1.t[1]-A1.t[0])*NP.asarray(range(0,len(self.Et)))
        self.Ef = DSP.FT1D(self.Et, ax=0, use_real=False, shift=False)
        self.f = spaxis(len(self.Et), resolution=A1.t[1]-A1.t[0])
        self.t = NP.fft.fftfreq(len(self.Ef),self.f[1]-self.f[0])        
        self.mode = 'XF'
    
    def __str__(self):
        return ' Instance of class "{0}" in module "{1}" \n label: {2} \n baseline: {3}'.format(self.__class__.__name__, self.__module__, self.label, self.pos.__str__())

    def XF(self):
        self.Ef = DSP.FT1D(self.Et, ax=0, use_real=False, shift=False)
        self.f = spaxis(len(self.Et), resolution=self.t[1]-self.t[0])
        self.mode = 'XF'

    def FX(self):
        # Ensure power of 2 and sufficient Nyquist length with zero padding to get identical results to XF-correlator output
        zero_pad_length1 = 2**NP.ceil(NP.log2(len(self.A1.Et)+len(self.A2.Et)-1))-len(self.A1.Et)   # zero pad length for first sequence
        zero_pad_length2 = 2**NP.ceil(NP.log2(len(self.A1.Et)+len(self.A2.Et)-1))-len(self.A2.Et)   # zero pad length for second sequence
        self.Ef = DSP.FT1D(NP.append(self.A1.Et, NP.zeros(zero_pad_length1)), ax=0, use_real=False, shift=False)*NP.conj(DSP.FT1D(NP.append(self.A2.Et, NP.zeros(zero_pad_length2)), ax=0, use_real=False, shift=False))  # product of FT
        self.f = spaxis(int(2**NP.ceil(NP.log2(len(self.A1.Et)+len(self.A2.Et)-1))), resolution=self.A1.t[1]-self.A1.t[0])
        self.mode = 'FX'
        self.t = NP.fft.fftfreq(len(self.Ef),self.f[1]-self.f[0])

#####################################################################  

# class AntennaArray_old:

#     # x_max, y_max, z_max = 0.0, 0.0, 0.0
#     # x_min, y_min, z_min = 0.0, 0.0, 0.0

#     # antlist = {}

#     def __init__(self):
#         self.antennas = {}
#         self.ants_xmin, self.ants_xmax, self.ants_ymin, self.ants_ymax = 0.0, 0.0, 0.0, 0.0
#         self.grid_xmin, self.grid_xmax, self.grid_ymin, self.grid_ymax = 0.0, 0.0, 0.0, 0.0
#         self.grid = 0.0
#         self.aperture = 0.0
#         self.E_fields_pol1 = 0.0
#         self.E_fields_pol2 = 0.0

#     def __str__(self):
#         printstr = ' Instance of class "{0}" in module "{1}".\n Holds the following "Antenna" class instances with labels:\n '.format(self.__class__.__name__, self.__module__)
#         printstr += '  '.join(self.antennas.keys())
#         printstr += '\n Antenna array bounds: tblc = [{0}, {1}], ttrc = [{2}, {3}]\n '.format(self.ants_xmin, self.ants_ymin, self.ants_xmax, self.ants_ymax)
#         printstr += '\n Grid bounds: tblc = [{0}, {1}], ttrc = [{2}, {3}]\n '.format(self.grid_xmin, self.grid_ymin, self.grid_xmax, self.grid_ymax)
#         return printstr

#     def add_antenna(self, A):
#         if A.label in self.antennas:
#             print "Antenna already included in the list of antennas."
#             print "For updating, use the update() method."
#         else:
#             self.antennas[A.label] = A
#             print 'Antenna "{0}" added to the list of antennas.'.format(A.label)

#     def remove_antenna(self, A=None, label=None):
#         if (A == None) and (label == None):
#             print 'No antenna specified for removal.'
#             pass
#         if (A.label in self.antennas) or (label in self.antennas):
#             if A.label in self.antennas:
#                 del self.antennas[A.label]
#                 print 'Antenna "{0}" removed from the list of antennas.'.format(A.label)
#             if label in self.antennas:
#                 del self.antennas[label]
#                 print 'Antenna "{0}" removed from the list of antennas.'.format(label)
#         else:
#             print 'Antenna(s) specified for removal not in the list of antennas.'

#     def antenna_grid(self, wavelength=1.0, spacing=0.5, xypad=None, pow2=True):


#         self.ants_xmin = min(ant.pos.x - min(ant.wtspos, key=lambda xyz: xyz[0])[0] for ant in self.antennas.values())
#         self.ants_xmax = max(ant.pos.x + max(ant.wtspos, key=lambda xyz: xyz[0])[0] for ant in self.antennas.values())
#         self.ants_ymin = min(ant.pos.y - min(ant.wtspos, key=lambda xyz: xyz[1])[1] for ant in self.antennas.values())
#         self.ants_ymax = max(ant.pos.y + max(ant.wtspos, key=lambda xyz: xyz[1])[1] for ant in self.antennas.values())

#         self.grid_xmin = self.ants_xmin/wavelength - xypad[0]
#         self.grid_xmax = self.ants_xmax/wavelength + xypad[0]
#         self.grid_ymin = self.ants_ymin/wavelength - xypad[1]
#         self.grid_ymax = self.ants_ymax/wavelength + xypad[1]

#     def grid_convolve():
#         pass
#         # conv_grid(self.pos.x, self.pos.y, self.wtspos.x, self.wtspos.y, self.gridwts, )

#     def ungrid(self):
#         pass

##############################################################################

class AntennaArray:

    # x_max, y_max, z_max = 0.0, 0.0, 0.0
    # x_min, y_min, z_min = 0.0, 0.0, 0.0

    # antlist = {}

    def __init__(self):
        self.antennas = {}
        self.ants_xmin, self.ants_xmax, self.ants_ymin, self.ants_ymax = 0.0, 0.0, 0.0, 0.0
        self.grid_xmin, self.grid_xmax, self.grid_ymin, self.grid_ymax = 0.0, 0.0, 0.0, 0.0
        self.gridx, self.gridy = 0.0, 0.0
        self.grid_illumination_polx = 0.0
        self.grid_illumination_poly = 0.0
        self.grid_Exf = 0.0
        self.grid_Eyf = 0.0
        

    def __str__(self):
        printstr = ' Instance of class "{0}" in module "{1}".\n Holds the following "Antenna" class instances with labels:\n '.format(self.__class__.__name__, self.__module__)
        printstr += '  '.join(self.antennas.keys())
        printstr += '\n Antenna array bounds: tblc = [{0}, {1}], ttrc = [{2}, {3}]\n '.format(self.ants_xmin, self.ants_ymin, self.ants_xmax, self.ants_ymax)
        printstr += '\n Grid bounds: tblc = [{0}, {1}], ttrc = [{2}, {3}]\n '.format(self.grid_xmin, self.grid_ymin, self.grid_xmax, self.grid_ymax)
        return printstr

    def __add__(self, others):
        retval = self.copy()
        if isinstance(others, AntennaArray):
            # for k,v in others.antennas.items():
            for k,v in others.antennas.iteritems():
                if k in retval.antennas:
                    print "Antenna {0} already included in the list of antennas.".format(k)
                    print "For updating, use the update() method. Ignoring antenna {0}".format(k)
                else:
                    retval.antennas[k] = v
                    print 'Antenna "{0}" added to the list of antennas.'.format(k)
        elif isinstance(others, dict):
            # for item in others.values():
            for item in others.itervalues():
                if isinstance(item, Antenna):
                    if item.label in retval.antennas:
                        print "Antenna {0} already included in the list of antennas.".format(item.label)
                        print "For updating, use the update() method. Ignoring antenna {0}".format(item.label)
                    else:
                        retval.antennas[item.label] = item
                        print 'Antenna "{0}" added to the list of antennas.'.format(item.label)
        elif isinstance(others, list):
            for i in range(0,len(others)):
                if isinstance(others[i], Antenna):
                    if others[i].label in retval.antennas:
                        print "Antenna {0} already included in the list of antennas.".format(others[i].label)
                        print "For updating, use the update() method. Ignoring antenna {0}".format(others[i].label)
                    else:
                        retval.antennas[others[i].label] = others[i]
                        print 'Antenna "{0}" added to the list of antennas.'.format(others[i].label)
                else:
                    print 'Element \# {0} is not an instance of class Antenna.'.format(i)
        elif isinstance(others, Antenna):
            if others.label in retval.antennas:
                print "Antenna {0} already included in the list of antennas.".format(others.label)
                print "For updating, use the update() method. Ignoring antenna {0}".format(others[i].label)
            else:
                retval.antennas[others.label] = others
                print 'Antenna "{0}" added to the list of antennas.'.format(others.label)
        else:
            print 'Input(s) is/are not instance(s) of class Antenna.'

        return retval

    def __radd__(self, others):
        return self.__add__(others)

    def __sub__(self, others):
        retval = self.copy()
        if isinstance(others, dict):
            for item in others.values():
                if isinstance(item, Antenna):
                    if item.label not in retval.antennas:
                        print "Antenna {0} does not exist in the list of antennas.".format(item.label)
                    else:
                        del retval.antennas[item.label]
                        print 'Antenna "{0}" removed from the list of antennas.'.format(item.label)
        elif isinstance(others, list):
            for i in range(0,len(others)):
                if others[i] in retval.antennas:
                    del retval.antennas[others[i]]
                    print 'Antenna {0} removed from the list of antennas.'.format(others[i])
                elif isinstance(others[i], Antenna):
                    if others[i].label in retval.antennas:
                        del retval.antennas[others[i].label]
                        print 'Antenna {0} removed from the list of antennas.'.format(others[i].label)
                    else:
                        print "Antenna {0} does not exist in the list of antennas.".format(others[i].label)
                else:
                    print 'Element \# {0} has no matches in the list of antennas.'.format(i)                        
        elif others in retval.antennas:
            del retval.antennas[others]
            print 'Antenna "{0}" removed from the list of antennas.'.format(others)
        elif isinstance(others, Antenna):
            if others.label in retval.antennas:
                del retval.antennas[others.label]
                print 'Antenna "{0}" removed from the list of antennas.'.format(others.label)
            else:
                print "Antenna {0} does not exist in the list of antennas.".format(others.label)
        else:
            print 'No matches found in existing list of antennas.'

        return retval


    def add_antennas(self, A=None):
        if A is None:
            print 'No antenna(s) supplied.'
        elif isinstance(A, (list, Antenna)):
            self = self.__add__(A)
        else:
            print 'Input(s) is/are not instance(s) of class Antenna.'


    def remove_antennas(self, A=None):
        if A is None:
            print 'No antenna specified for removal.'
        else:
            self = self.__sub__(A)


    def grid(self, wavelength=1.0, spacing=0.5, xypad=None, pow2=True):

        # Change itervalues() to values() when porting to Python 3.x

        ants_xmin_polx = min(ant.pos.x - min(ant.wtspos_x, key=lambda xyz: xyz[0])[0] for ant in self.antennas.itervalues())
        ants_xmax_polx = max(ant.pos.x + max(ant.wtspos_x, key=lambda xyz: xyz[0])[0] for ant in self.antennas.itervalues())
        ants_ymin_polx = min(ant.pos.y - min(ant.wtspos_x, key=lambda xyz: xyz[1])[1] for ant in self.antennas.itervalues())
        ants_ymax_polx = max(ant.pos.y + max(ant.wtspos_x, key=lambda xyz: xyz[1])[1] for ant in self.antennas.itervalues())

        ants_xmin_poly = min(ant.pos.x - min(ant.wtspos_y, key=lambda xyz: xyz[0])[0] for ant in self.antennas.itervalues())
        ants_xmax_poly = max(ant.pos.x + max(ant.wtspos_y, key=lambda xyz: xyz[0])[0] for ant in self.antennas.itervalues())
        ants_ymin_poly = min(ant.pos.y - min(ant.wtspos_y, key=lambda xyz: xyz[1])[1] for ant in self.antennas.itervalues())
        ants_ymax_poly = max(ant.pos.y + max(ant.wtspos_y, key=lambda xyz: xyz[1])[1] for ant in self.antennas.itervalues())

        self.ants_xmin = min(ants_xmin_polx, ants_xmin_poly)
        self.ants_xmax = max(ants_xmax_polx, ants_xmax_poly)
        self.ants_ymin = min(ants_ymin_polx, ants_ymin_poly)
        self.ants_ymax = max(ants_ymax_polx, ants_ymax_poly)

        # self.grid_xmin = (self.ants_xmin - xypad[0])/wavelength
        # self.grid_xmax = (self.ants_xmax + xypad[0])/wavelength
        # self.grid_ymin = (self.ants_ymin - xypad[1])/wavelength
        # self.grid_ymax = (self.ants_ymax + xypad[1])/wavelength

        self.gridx, self.gridy = GRD.grid_2d([(self.ants_xmin/wavelength, self.ants_xmax/wavelength),(self.ants_ymin/wavelength, self.ants_ymax/wavelength)], pad=xypad/wavelength, spacing=0.5, pow2=True)

    def grid_convolve(self):
        grid_illumination_polx = 0.0
        grid_illumination_poly = 0.0
        x_wtspos_polx = []
        y_wtspos_polx = []
        x_wtspos_poly = []
        y_wtspos_poly = []

        for ant in self.antennas.itervalues():
            x_wtspos_polx += [ant.wtspos_x.x]
            y_wtspos_polx += [ant.wtspos_x.y]
            x_wtspos_poly += [ant.wtspos_y.x]
            y_wtspos_poly += [ant.wtspos_y.y]

        for label, ant in self.antennas.iteritems():
            grid_illumination_polx = GRD.conv_grid2d(ant.pos.x, ant.pos.y, ant.wtspos_x, ant.wtspos_y, ant.wts_x, ant.wts_y, self.gridx, self.gridy, method='CS')
            grid_illumination_poly = GRD.conv_grid2d(ant.pos.x, ant.pos.y, ant.wtspos_x, ant.wtspos_y, ant.wts_x, ant.wts_y, self.gridx, self.gridy, method='CS')
        

    # def ungrid(self):
    #     pass

    def update(additions=None, removals=None):
        if not additions is None:
            self.add_antennas(additions)
        if not removals is None:
            self.remove)antennas(removals)

        self.grid()

#############################################################################

