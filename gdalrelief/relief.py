#!/usr/bin/env python
# a bar plot with errorbars
import numpy as np
from osgeo import gdal
import matplotlib.pyplot as plt
import gdalrelief.diff as diff
from scipy.signal import argrelextrema
import collections

#TESTRASTER="../data/Goloustnoye/ALTITUDE 1Trim.grd"
TESTRASTER="../data/Olkhon/dem.gtiff"

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
        self.alphas={}

    def load(self, raster):
        if type(raster)==str:
            rc=gdal.Open(raster)
            if rc == None:
                raise ValueError("cannot open raster file")
            return rc
        else:
            return raster

    def __call__(self, layer):
        return self.band(layer)

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

    def display(self, raster, between=None, cmap=plt.cm.gray, raster_alpha=None, **kwargs):
        fig, ax = plt.subplots()
        #print (raster)
        if between:
            raster=np.where(raster <= between[0], np.nan, raster)
            raster=np.where(raster >= between[1], np.nan, raster)
        if raster_alpha != None:
            print(raster.shape,raster_alpha.shape)
            raster[raster_alpha]=np.nan
        # raster[raster<=0]=np.nan
        #print (raster)
        image = raster
        (mx,my) = image.shape
        # ax.imshow(image, cmap=plt.cm.gray, interpolation=interpolation)
        ax.imshow(image, cmap=cmap, **kwargs)
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

    def band(self, layer):
        band=self.raster.GetRasterBand(layer).ReadAsArray().astype(np.float)
        alpha=band<=0
        self.alphas[layer]=alpha
        band[alpha]=np.nan
        return band


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
        zu=list(zu[0])
        zd=list(zd[0])
        iu=id=0
        while True:
            try:
                d=abs(z[zu[iu]]-z[zd[id]])>=minheight
            except IndexError:
                break
            if d:
                if (zu[iu]<zd[id]):
                    iu+=1
                else:
                    id+=1
            else:
                if (iu<id):
                    zu.pop(iu)
                else:
                    zd.pop(id)
        azu=np.empty(shape=(len(zu),3), dtype=float)
        azu[:,0]=x[zu]
        azu[:,1]=y[zu]
        azu[:,2]=z[zu]
        azd=np.empty(shape=(len(zd),3), dtype=float)
        azd[:,0]=x[zd]
        azd[:,1]=y[zd]
        azd[:,2]=z[zd]
        yield azu,azd,r,l

    def scan(self, layer, minheight):
        for azu,azd,r,l in self.scan_line(layer, minheight):
            yield azu,azd

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


class RasterPlastics(RasterProcessor):
    """Figures out plastic data of the relief.
    """

    def __init__(self, *args, **kwargs):
        RasterProcessor.__init__(self, *args, **kwargs)
        self.gradient=None
        self.plastic=None

    def __call__(self, layer, method="simple", r=1, bitonal=False):
        """
        """
        band=RasterProcessor.__call__(self, layer)
        gx, gy = diff.gradient(band)
        self.gradient=(gx,gy)
        if method=="simple":
            plastic=diff.agrad(gx, gy)
        elif method=="circle":
            plastic=diff.roundagrad(gx, gy, r=r)
        if bitonal:
            plastic=np.where(plastic > 1, 1, plastic)
            plastic=np.where(plastic < -1, -1, plastic)
        self.plastic=plastic
        return plastic

# TEST

def test_1():
    # h=Hatch(100,100, 732,684, 5)
    h=Hatch(510,520, 500,500, 5)
    rs=RasterSection(raster=TESTRASTER, hatch=h)
    rs.info()
    for u,d in rs.scan(4, 300):
        print ("===", u,d)
    """
    print (sec)
    height=sec[:,2]
    valid=height[height>0]
    print (valid)
    x1=np.arange(len(valid))  #np.linspace(0,832)
    y1=valid
    #x2=np.linspace(0,784)
    #y2=valid
    #plt.subplot(2,1,1)
    plt.plot(x1, y1, 'r.-')

    #plt.plot(x2, y2, 'r.-')
    plt.show()
    """

def test_plastics():
    rp=RasterPlastics(raster=TESTRASTER)
    rp.info()
    layer=1
    plastic=rp(layer, method="simple", r=1, bitonal=False)

    alpha=rp.alphas[layer][2:-2,2:-2]
    plastic[alpha]=np.nan
    rp.display(plastic,
               between=(-100,100),
               interpolation="none",
               raster_alpha=alpha)
    rp.display(rp.band(layer), interpolation="none", cmap='gist_earth')

if __name__=="__main__":
    #test_1()
    test_plastics()
    quit()
