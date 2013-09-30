import sys as SYS
import numpy as NP
import my_operations as OPS
import my_DSP_modules as DSP
import geometry as GEOM

def FT(inp, ax=0, use_real=False, shift=False):
    return DSP.FT1D(inp, ax=0, use_real=False, shift=False)

def XC(inp1, inp2=0.0):
    if (len(inp2) < 1):
        inp2 = inp1
    if (len(inp2) == 1) and (inp2 == 0.0):
        inp2 = inp1
    zero_pad_length = 2**NP.ceil(NP.log2(len(inp1)+len(inp2)-1))-(len(inp1)+len(inp2)-1)
    return NP.roll(NP.append(NP.correlate(inp1, inp2, mode='full'), NP.zeros(zero_pad_length)), -(len(inp2)-1))   # zero pad and shift to ensure identical results as FX() operation

def chans(length, resolution=1.0, use_real=False, shift=False):
    return DSP.frequencies(length, delx=resolution, use_real=False, shift=False)

class Antenna:
    def __init__(self, label, pos, Et, t):
        self.label = label
        # self.pos = GEOM.Point(pos.x, pos.y, pos.z)
        self.pos = pos
        self.Et = NP.asarray(Et)
        self.t = NP.asarray(t)
        self.Ef = FT(self.Et, ax=0, use_real=False, shift=False)
        self.f = chans(len(self.Et), resolution=self.t[1]-self.t[0])

    def __str__(self):
        return ' Instance of class "{0}" in module "{1}" \n label: {2} \n location: {3}'.format(self.__class__.__name__, self.__module__, self.label, self.pos.__str__())

    def F(self):
        return FT(self.Et, ax=0, use_real=False, shift=False)

    def channels(self):
        return chans(len(self.Et), resolution=self.t[1]-self.t[0])

class AntennaPair(Antenna):
    def __init__(self, A1, A2):
        self.A1, self.A2 = A1, A2
        self.label = A1.label+'-'+A2.label
        self.pos = A1.pos-A2.pos
        self.Et = XC(self.A1.Et, self.A2.Et)
        # self.t = (A1.t[1]-A1.t[0])*NP.asarray(range(0,len(self.Et)))
        self.Ef = FT(self.Et, ax=0, use_real=False, shift=False)
        self.f = chans(len(self.Et), resolution=A1.t[1]-A1.t[0])
        self.t = NP.fft.fftfreq(len(self.Ef),self.f[1]-self.f[0])        
        self.mode = 'XF'
    
    def __str__(self):
        return ' Instance of class "{0}" in module "{1}" \n label: {2} \n baseline: {3}'.format(self.__class__.__name__, self.__module__, self.label, self.pos.__str__())

    def XF(self):
        self.Ef = FT(self.Et, ax=0, use_real=False, shift=False)
        self.f = chans(len(self.Et), resolution=self.t[1]-self.t[0])
        self.mode = 'XF'

    def FX(self):
        # Ensure power of 2 and sufficient Nyquist length with zero padding to get identical results to XF-correlator output
        zero_pad_length1 = 2**NP.ceil(NP.log2(len(self.A1.Et)+len(self.A2.Et)-1))-len(self.A1.Et)   # zero pad length for first sequence
        zero_pad_length2 = 2**NP.ceil(NP.log2(len(self.A1.Et)+len(self.A2.Et)-1))-len(self.A2.Et)   # zero pad length for second sequence
        self.Ef = FT(NP.append(self.A1.Et, NP.zeros(zero_pad_length1)), ax=0, use_real=False, shift=False)*NP.conj(FT(NP.append(self.A2.Et, NP.zeros(zero_pad_length2)), ax=0, use_real=False, shift=False))  # product of FT
        self.f = chans(int(2**NP.ceil(NP.log2(len(self.A1.Et)+len(self.A2.Et)-1))), resolution=self.A1.t[1]-self.A1.t[0])
        self.mode = 'FX'
        self.t = NP.fft.fftfreq(len(self.Ef),self.f[1]-self.f[0])
