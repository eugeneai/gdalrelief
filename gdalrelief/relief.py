#!/usr/bin/env python
# a bar plot with errorbars
import numpy as np
from osgeo import gdal
import matplotlib.pyplot as plt
import gdalrelief.diff as diff
#import scipy

TESTRASTER="../data/Goloustnoye/ALTITUDE 1Trim.grd"

class Hatch(object):
    """Defines section data consisting of
    x1,y1, x2,y2 are points on the section line and
    d is a distance from the line
    """

    def __init__(self, x1,y1,  x2,y2, d):
        self.p1=(x1,y1)
        self.p2=(x2,y2)
        self.d=d

class RasterProcessor(object):
    def __init__(self, raster):
        self.raster=self.load(raster)

    def load(self, raster):
        if type(raster)==str:
            rc=gdal.Open(raster)
            if rc == None:
                raise ValueError("cannot open raster file")
            return rc
        else:
            return raster

    def __call__(self, layer):
        sec=[]
        band=self.raster.GetRasterBand(layer).ReadAsArray()
        return band
        print (len(band[0]))
        for x,y in self.line():
            sec.append(band[x,y])
        return np.array(sec)

    def info(self):
        src_ds=self.raster
        print ("\n")
        print ("-----------------------------------------")
        print ("[ RASTER BAND COUNT ]: ", src_ds.RasterCount)
        for band in range( src_ds.RasterCount ):
            band += 1
            print ("[ GETTING BAND ]: ", band)
            srcband = src_ds.GetRasterBand(band)
            if srcband is None:
                continue

            stats = srcband.GetStatistics( True, True )
            if stats is None:
                continue

            print ("[ STATS ] =  Minimum=%.3f, Maximum=%.3f, Mean=%.3f, StdDev=%.3f" % ( \
                        stats[0], stats[1], stats[2], stats[3] ))

        print ("-----------------------------------------")

    def display(self, raster, interpolation="nearest", between=None):
        fig, ax = plt.subplots()
        print (raster)
        if between:
            raster=np.where(raster < between[0], between[0], raster)
            raster=np.where(raster > between[1], between[1], raster)
        print (raster)
        image = raster
        (mx,my) = image.shape
        ax.imshow(image, cmap=plt.cm.gray, interpolation=interpolation)
        ax.set_title('Display of a raster')

        # Move left and bottom spines outward by 10 points
        ax.spines['left'].set_position(('outward', my))
        ax.spines['bottom'].set_position(('outward', mx))
        # Hide the right and top spines
        ax.spines['right'].set_visible(False)
        ax.spines['top'].set_visible(False)
        # Only show ticks on the left and bottom spines
        ax.yaxis.set_ticks_position('left')
        ax.xaxis.set_ticks_position('bottom')

        plt.show()


class RasterSection(RasterProcessor):
    """
    """

    def __init__(self, raster, hatch):
        RasterProcessor.__init__(self, raster)
        self.hatch=hatch

    def __call__(self, layer):
        """
        """
        band=RasterProcessor.__call__(self, layer)
        sec=[]
        for x,y in self.line():
            sec.append(band[x,y])
        return np.array(sec)

    def line(self):
        def sign(x):
            if x<0:
                return -1
            elif x>0:
                return 1
            else:
                return 0

        # Brasenham

        p1=self.hatch.p1
        p2=self.hatch.p2

        x,y=p1
        dx,dy=p2[0]-p1[0], p2[1]-p1[1]
        d=max(dx,dy)
        ex=ey=0
        ix=sign(dx)
        iy=sign(dy)

        while True:
            yield (x,y)
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

class RasterPlastics(RasterProcessor):
    """Figures out plastic data of the relief.
    """

    def __init__(self, *args, **kwargs):
        RasterProcessor.__init__(self, *args, **kwargs)
        self.gradient=None
        self.plastic=None

    def __call__(self, layer, method="simple", r=1):
        """
        """
        band=RasterProcessor.__call__(self, layer)
        gx, gy = diff.gradient(band)
        self.gradient=(gx,gy)
        if method=="simple":
            plastic=diff.agrad(gx, gy)
        elif method=="circle":
            plastic=diff.roundagrad(gx, gy, r=r)
        self.plastic=plastic
        return plastic

# TEST

def test_1():
    h=Hatch(0,0, 832,784, 5)
    rs=RasterSection(raster=TESTRASTER, hatch=h)
    rs.info()
    sec=rs(4)
    valid=sec[sec>0]
    print (valid)
    x1=np.arange(len(valid))  #np.linspace(0,832)
    y1=valid
    #x2=np.linspace(0,784)
    #y2=valid
    #plt.subplot(2,1,1)
    plt.plot(x1, y1, 'r.-')

    #plt.plot(x2, y2, 'r.-')
    plt.show()

def test_plastics():
    rp=RasterPlastics(raster=TESTRASTER)
    plastic=rp(4)
    rp.display(plastic, between=(-10,10), interpolation="none")

if __name__=="__main__":
    #test_1()
    test_plastics()
    quit()
