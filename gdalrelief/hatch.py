#!/usr/bin/env python
# a bar plot with errorbars
import numpy as np
import numpy.ma as ma
from osgeo import gdal
import matplotlib.pyplot as plt
import gdalrelief.diff as diff
from scipy.signal import argrelextrema
import collections
from osgeo.gdalconst import *
from common import *


TESTRASTER_GOLOUSTNOYE="../data/Goloustnoye/ALTITUDE 1Trim.grd"
TESTRASTER_OLKHON="../data/Olkhon/dem.gtiff"
SAFE_GDAL_FORMAT="GTiff"

class Hatch(object):
    """Defines section data consisting of
    x1,y1, x2,y2 are points on the section line and
    d is a distance from the line
    """

    def __init__(self, x1,y1,  x2,y2, d):
        self.p1=(x1,y1)
        self.p2=(x2,y2)
        self.d=d

class RasterSection(RasterProcessor):
    """
    """

    def __init__(self, raster, hatch):
        RasterProcessor.__init__(self, raster)
        self.hatch=hatch

    def __call__(self, layer):
        for x,y in self.line():
            sec.append(band[x,y])
        return np.array(sec)

    def scan_line(self, layer, minheight, p1=None, p2=None):
        band=RasterProcessor.__call__(self, layer)
        self.current_band=band
        if p1 == None:
            p1=self.hatch.p1
        if p2 == None:
            p2=self.hatch.p2
        sec=collections.deque()
        for x,y in self.line(p1,p2,extra=True):
            sec.append((x,y,band[y,x]))
            r=(x,y)
        for x,y in self.line(p1,p2,extra=True, forward=False, current=False):
            sec.appendleft((x,y,band[y,x]))
            l=(x,y)
        a=np.array(sec)
        z=a[:,2]
        zv=z>0

        x=a[zv,0]
        y=a[zv,1]
        z=a[zv,2]

        zu=argrelextrema(z, np.greater)
        zd=argrelextrema(z, np.less)
        zu=list(zu[0]) # Numers of heights values
        zd=list(zd[0])
        iu=id=0
        # Removing small extremums

        #import pudb; pu.db
        while True:
            # FIXME, I'm wrong
            try:
                d=abs(z[zu[iu]]-z[zd[id]])>=minheight
            except IndexError:
                break
            q=zu[iu]<zd[id] # x-coord of lower point is early
            if d:
                if (q): # Compare x coords
                    iu+=1
                else:
                    id+=1
            else:
                zu.pop(iu)
                zd.pop(id)
                if (q):
                    if id>0:
                        id-=1
                else:
                    if iu>0:
                        iu-=1

        azu=np.empty(shape=(len(zu),3), dtype=float)
        azu[:,0]=x[zu]
        azu[:,1]=y[zu]
        azu[:,2]=z[zu]
        azd=np.empty(shape=(len(zd),3), dtype=float)
        azd[:,0]=x[zd]
        azd[:,1]=y[zd]
        azd[:,2]=z[zd]
        yield z,zu,zd,azu,azd,r,l

    def scan(self, layer, minheight):
        for z,zu,zd,azu,azd,r,l in self.scan_line(layer, minheight):
            yield z,zu,zd,azu,azd

    def line(self, p1, p2, extra=False, forward=True, current=True):
        def sign(x):
            if x<0:
                return -1
            elif x>0:
                return 1
            else:
                return 0

        band=self.current_band
        my,mx=band.shape
        # print (my,mx)

        # Brasenham

        x,y=p1
        dx,dy=(p2[0]-p1[0]), (p2[1]-p1[1])
        ex=ey=0
        ix=sign(dx)
        iy=sign(dy)
        if not forward:
            ix=-ix
            iy=-iy
        dx,dy=abs(dx),abs(dy)
        d=max(dx,dy)
        if current:
            yield (x,y)

        while True:
            if not extra:
                if x==p2[0] and y==p2[1]:
                    return
            ex+=dx
            ey+=dy
            ax=ex-d
            ay=ey-d
            if ax>=0:
                x+=ix
                ex=ax
            if ay>=0:
                y+=iy
                ey=ay
            if x<0 or y<0:
                return
            if x>=mx or y>=my:
                return
            yield (x,y)


# TEST

def test_1():
    # h=Hatch(100,100, 732,684, 5)
    h=Hatch(510,520, 500,500, 5)
    rs=RasterSection(raster=TESTRASTER_GOLOUSTNOYE, hatch=h)
    rs.info()
    for z,zu,zd,u,d in rs.scan(4, 20):
        print ("===", u,d)

    # print (sec)
    height=z
    #valid=height[height>0]
    valid=height

    print (valid)

    x1=np.arange(len(valid))  #np.linspace(0,832)
    y1=valid
    #x2=np.linspace(0,784)
    #y2=valid
    #plt.subplot(2,1,1)
    plt.plot(x1, y1, 'r.-')
    plt.plot(zu, z[zu], 'g.-')
    plt.plot(zd, z[zd], 'g.-')


    #plt.plot(x2, y2, 'r.-')
    plt.show()

if __name__=="__main__":
    # register all of the GDAL drivers
    gdal.AllRegister()

    test_1()

    quit()
