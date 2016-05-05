# -*- coding: utf-8 -*-
#
# Copyright 2016 University of Nevada, Reno
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
# Author: Mark Williams, Nevada Seismological Laboratory, UNR
#
"""
eqclustering

Python translation of MATLAB scripts/functions

- Mark C. Williams (2016) Nevada Seismological Lab

Port of a bunch of MATLAB code -- requires numpy

Defines one class -- BPTree, all functionality can be accessed via methods on
this instance.
"""
import os
import datetime
import numpy as np

def _ts(dt):
    """Float timestamp from python datetime"""
    return (dt-datetime.datetime(1970, 01, 01, 00, 00, 00)).total_seconds()

# I/O --- remove or do something ---------------------------------------------
def _input_line(line):
    """
    Return map from line
    """
    l = line.split()
    r = {
        "year": int(l[0]),
        "month": int(l[1]),
        "day": int(l[2]),
        "hour": int(l[3]),
        "minute": int(l[4]),
        "second": float(l[5]),
        "mag": float(l[9]),
        "lat": float(l[6]),
        "lon": float(l[7]),
        "dep": float(l[8]),
        #"evid": int(l[10]),
    }
    r['timestamp'] = _ts(datetime.datetime(r['year'], r['month'], r['day'],
        r['hour'], r['minute'], 0)) + r['second']
    return r


def read_input_file(f):
    return [_input_line(l) for l in f if l]
#------------------------------------------------------------------------------

def normalmix_1D(x, mu1i, mu2i, sigma1i, sigma2i):
    """
    Estimates parameters of a 2-component 1-D Gaussian mixture model
    (uses the EM algorithm)
    
    Inputs
    ------
    x : the data sample to be analysed
    mu#i,sigma#i : initial parameters of components
    
    Outputs
    -------
    c, mu1, mu2, sigma1, sigma2
        c : the maximum likelihood boundary between the hats
        mu#, sigma# : are the estimated paraeters of components
     
    """
    Ai = 0.5
    A2i = 0.5
    epsilon = 1
    p1a = np.zeros(x.size)
    p2a = np.zeros(x.size)
    p1 = np.zeros(x.size)
    p2 = np.zeros(x.size)
    epsilon0 = 0.001

    while epsilon >  epsilon0:
        p1a = Ai * (1/(np.sqrt(2*np.pi*sigma1i**2))*np.exp(-(((x-mu1i)/sigma1i)**2)/2))
        p2a = A2i * (1/(np.sqrt(2*np.pi*sigma2i**2))*np.exp(-(((x-mu2i)/sigma2i)**2)/2))

        p1 = p1a / (p1a+p2a)
        p2 = p2a / (p1a+p2a)

        A = np.sum(p1)/len(x)
        A2 = np.sum(p2)/len(x)

        mu1 = np.sum(p1*x)/np.sum(p1)
        mu2 = np.sum(p2*x)/np.sum(p2)

        sigma1 = np.sqrt(np.sum(p1*(x-mu1)**2)/np.sum(p1))
        sigma2 = np.sqrt(np.sum(p2*(x-mu2)**2)/np.sum(p2))

        epsilon1 = np.abs(Ai-A)
        epsilon2 = np.abs(mu1i-mu1)
        epsilon3 = np.abs(mu2i-mu2)
        epsilon4 = np.abs(sigma1i-sigma1)
        epsilon5 = np.abs(sigma2i-sigma2)
        epis = [epsilon1, epsilon2, epsilon3, epsilon4, epsilon5]
        epsilon = np.amax(epis)

        Ai = A
        A2i = A2
        mu1i = mu1
        mu2i = mu2
        sigma1i = sigma1
        sigma2i = sigma2

    aa = -2*sigma1**2 + 2*sigma2**2
    bb = (2*mu2 * 2*sigma1**2) + (-2*mu1 * 2*sigma2**2)
    cc = (-mu2**2 * 2*sigma1**2) + (mu1**2* 2*sigma2**2) - (np.log(sigma2/sigma1) * 2*sigma1**2 * 2*sigma2**2)
    c1 = (-bb + np.sqrt(bb**2 - 4*aa*cc)) / (2*aa)
    c2 = (-bb - np.sqrt(bb**2 - 4*aa*cc)) / (2*aa)

    if mu1 < c1 < mu2:
        c = c1
    else:
        c = c2
    return c, mu1, mu2, sigma1, sigma2, A, A2, p1, p2


def descend_all(parent, i):
    """
    Return all children from parent of index i
    """
    d = [i]
    for child in np.arange(parent.size)[parent==i]:
        d += descend_all(parent, child)
    return d


class BPFunction(object):
    """
    Baiesi-Paczuski metric function
    
    Inputs
    ------
    time  : the ocurrence time of events in years
    mag   : event magnitude
    lon   : event longitudes
    lat   : event latitudes
    depth : event depths
    b     : GR slope
    df    : fractal dimension of epicenters

    Outputs (all are the same size as the EQ catalog):
    -------
    P : resulting tree (parent-pointer format)
    n : nearest-neighbor BP-distance
    D : nearest-neighbor spatial distance
    T : nearest-neighbor temporal distance
    M : magniude of the parent event
    
    Notes
    -----
    Reference: PHYSICAL REVIEW E 69, 066106 (2004)

    """
    # Currently not used
    tmin = 0 # seconds
    dmin = 0 # meters
    
    @staticmethod
    def dist(pnt, lats, lons, deps=None):
        """
        Surface distance (in km) between points given as [lat,lon,depth]
        If depth is not specified, the surface distance is computed
        """
        R0 = 6.3673e6
        
        lat1, lon1 = pnt[:2]
        lat2 = lats
        lon2 = lons

        D = R0 * np.arccos(
            np.sin(np.radians(lat1)) * np.sin(np.radians(lat2)) +
            np.cos(np.radians(lat1)) * np.cos(np.radians(lat2)) * np.cos(np.radians(lon1-lon2))
        )

        if len(pnt) == 3 and deps:
            dep1 = pnt[2] * 1000
            dep2 = deps * 1000
            D = np.sqrt(D**2 + (dep1-dep2)**2)
        
        return D/1000.0

    def __init__(self, tmin=0, dmin=0):
        self.tmin = tmin
        self.dmin = dmin

    def __call__(self, time, mag, lon, lat, depth, b, df, lag=None, iswaitbar=False):
        """
        Run fxn
        """
        _n = len(time)

        P = np.empty(_n) * np.NaN # tree
        n = np.empty(_n) * np.Inf # eta*
        T = np.empty(_n) * np.NaN # T used in eta*
        D = np.empty(_n) * np.NaN # R used in eta*
        M = np.empty(_n) * np.NaN # Mag of parent event in tree
        
        if lag is None:
            lag = np.Inf
        
        # Function to calulate "Eta(i,j)" distance
        def eta(t, d, m):
            return t * d**df * 10**(-b*m)
        
        for i in np.arange(time.size):
            # TODO: put in logger or toggle verbose info
            #print "{0}".format(i),
            
            # This mask only calulates distances for previous events. This is
            # needed for the implementation of the tree which is strictly
            # time-ordered, events are only associated with previous events.
            I = (time < time[i]) & (time >= time[i]-lag)
            
            if np.any(I):
                d = 1000 * self.dist([lat[i], lon[i]], lat[I], lon[I])  # TODO: compare to dmin?
                t =  time[i] - time[I]  # TODO: compare to tmin?
                nc = eta(t, d, mag[I])

                imin = nc.argmin()
                n[i] = nc[imin]
                P[i] = np.arange(len(I))[imin]
                D[i] = d[imin]
                T[i] = t[imin]
                M[i] = mag[I][imin]

        T = T/365.35/24/60/60
        D /= 1000 
        return P, n, D, T, M


class BPTree(object):
    """
    Class to create/store/analyze B-P tree structure
    """
    df = 1.6  # Fractal dimension of epicenters
    b = 1.0   # GR slope
    p = 0.5   # Weight parameter to define Tn and Rn
    
    bp = BPFunction()  # setup bp function

    @classmethod
    def from_file(cls, fname):
        # Load the catalog
        with open(fname) as f:
            clg = read_input_file(f)
        # Init arrays
        tm, lat, lon, mag, dep = zip(*[(c['timestamp'], c['lat'], c['lon'], 
            c['mag'], c['dep']) for c in clg])
        return cls(times=tm, mags=mag, lons=lon, lats=lat, depths=dep)

    def __init__(self, times=[], mags=[], lons=[], lats=[], depths=[]):
        self._time = np.array(times)
        self._mag = np.array(mags) 
        self._lon = np.array(lons) 
        self._lat = np.array(lats)
        self._dep = np.array(depths)
    
    def grow(self):
        """
        Call bp function to construct the tree from event distances

        TODO: enable appending to an existing tree?
        """
        tr = self.bp(self._time, self._mag, self._lon, self._lat, [], self.b, 
                self.df, None, False)
        self.P = tr[0]
        self.n = tr[1]
        self.D = tr[2]
        self.T = tr[3]
        self.M = tr[4]

    def normalized_time(self):
        return self.T * 10**(-self.p*self.b*self.M)  # Normalized time

    def normalized_distance(self):
        return self.D**self.df * 10**(-(1-self.p)*self.b*self.M) # Normalized distance
    
    def prune(self, c=None):
        """
        Prune a B-P tree to find clusters for a given cutoff eta value

        If no c is given, calulate one.
        """ 
        TN = self.normalized_time()
        RN = self.normalized_distance()

        if c is None:
            I = (np.isfinite(np.log10(TN*RN))) & (TN > 0) & (RN > 0)
            c, mu1, mu2, sigma1, sigma2, A, A2, p1, p2 = normalmix_1D(np.log10(TN[I]*RN[I]), -6, -3, 1, 1)
            #print "Calculated c:{0}".format(c)  # TODO: log debug
        
        P1 = self.P.copy()  # maybe not needed
        K = (np.log10(TN) + np.log10(RN) > c)  # Mainshock condition
        P1[K] = np.nan

        # TODO: fixed in bp function, can remove...
        #P1[0] = np.nan  # first event won't pass cut, but is always a root

        self.Ifor = np.array([])
        self.Iaft = np.array([])
        self.Imain = np.array([])
        self.Pfin = np.zeros(P1.size)

        K = np.isnan(P1) # weird, K should already be this...
        Ki = np.arange(K.size)[K]

        for i in np.arange(Ki.size):
            I0 = descend_all(P1, Ki[i])
            I = np.array(I0)
            # TODO: if not mag: imax = I[0]
            imax = self._mag[I].argmax()
            IM = I[imax]
            IF0 = I[self._time[I] < self._time[IM]]
            IA0 = I[self._time[I] > self._time[IM]]
            #print "Subtree: event:{0}, main:{2}, cluster:{1}".format(Ki[i], I0, IM)
            
            # Keep track of fore/after shocks
            self.Ifor = np.append(self.Ifor, IF0)
            self.Iaft = np.append(self.Iaft, IA0)
            self.Imain = np.append(self.Imain, IM)
            
            # Change parent of all cluster members to main shock ID
            self.Pfin[IF0] = IM
            self.Pfin[IA0] = IM
            self.Pfin[IM] = IM

        self.TN = TN
        self.RN = RN

    def output2file(self, fname):
        with open(fname, 'w') as f:
            for n in np.arange(self.P.size):
                f.write("{0} {1} {2} {3} {4} {5} {6} {7}\n".format(
                    n,
                    self._time[n],
                    self._lat[n],
                    self._lon[n],
                    self._dep[n],
                    self._mag[n],
                    self.P[n],
                    self.Pfin[n],
                    )
                )


def main():
    """
    Command line version
    """
    pass  # TODO: parse args and use default CSV style file format???


