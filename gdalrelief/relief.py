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

# register all of the GDAL drivers
gdal.AllRegister()

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

    def display(self, raster, between=None, cmap=plt.cm.gray, layer=None, **kwargs):
        fig, ax = plt.subplots()
        #print (raster)
        if between:
            raster=np.where(raster <= between[0], np.nan, raster)
            raster=np.where(raster >= between[1], np.nan, raster)
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
        bnd=self.raster.GetRasterBand(layer)
        band=bnd.ReadAsArray()
        nan_value=bnd.GetNoDataValue()
        if nan_value == None:
            nan_value=np.amin(band)

        return ma.masked_values(band, nan_value)

    def save(self, filename, data, sx=0, sy=0, driver=None, nan_value=None):
        if driver == None:
            driver = self.raster.GetDriver()
        rows,cols = data.shape
        outGRID = driver.Create(filename, cols, rows, 1, GDT_Float32)
        if outGRID is None:
            driver = gdal.GetDriverByName(SAFE_GDAL_FORMAT)
            filename+=".gtiff"
            outGRID = driver.Create(filename, cols, rows, 1, GDT_Int32)
            if outGRID is None:
                raise RuntimeError("Could not create {}.".format(filename))
        else:
            if nan_value == None:
                nan_value = self.raster.GetRasterBand(1).GetNoDataValue()
        if nan_value == None:
            nan_value = -7000000
        outBand = outGRID.GetRasterBand(1)
        d=np.array(data.data)
        d[data.mask]=nan_value
        outBand.WriteArray(d, sx, sy)

        # flush data to disk, set the NoData value and calculate stats
        outBand.FlushCache()
        outBand.SetNoDataValue(nan_value)

        # georeference the image and set the projection
        outGRID.SetGeoTransform(self.raster.GetGeoTransform())
        outGRID.SetProjection(self.raster.GetProjection())
        outBand.FlushCache()
        del outBand
        del outGRID
        print ("Saved {}".format(filename))




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
        yield z,azu,azd,r,l

    def scan(self, layer, minheight):
        for z,azu,azd,r,l in self.scan_line(layer, minheight):
            yield z,azu,azd

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
        gx, gy = self.get_diffs(layer=layer)
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

    def get_diffs(self, band=None, layer=None):
        if layer != None:
            band=RasterProcessor.band(self, layer)
        return diff.gradient(band)

# TEST

def test_1():
    # h=Hatch(100,100, 732,684, 5)
    h=Hatch(510,520, 500,500, 5)
    rs=RasterSection(raster=TESTRASTER_GOLOUSTNOYE, hatch=h)
    rs.info()
    for z,u,d in rs.scan(4, 30):
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

    #plt.plot(x2, y2, 'r.-')
    plt.show()

def script_plastics(raster, name, layer, display=False):
    rp=RasterPlastics(raster) # =TESTRASTER)
    rp.info()

    lmin,lmax=-200,200

    r=1
    plastic=rp(layer, method="simple", r=1, bitonal=False)
    print (plastic)
    if display:
        rp.display(plastic, interpolation="none")
    # rp.display(rp.band(layer), interpolation="none", cmap='gist_earth')

    gx,gy=rp.get_diffs(layer=layer)

    rp.save("plastic-{}-simple-r-1-n.gtiff".format(name), plastic)
    rp.save("plastic-{}-diff-x.gtiff".format(name), gx)
    rp.save("plastic-{}-diff-y.gtiff".format(name), gy)

    for r in [5,10]:
        plastic=rp(layer, method="circle", r=r, bitonal=False)
        #rp.display(plastic, between=(lmin,lmax), interpolation="none")
        rp.save("plastic-{}-circle-r-{}-n.gtiff".format(name,r), plastic)



if __name__=="__main__":
    test_1()
    #script_plastics(TESTRASTER_OLKHON, name="olkhon", layer=1)
    #script_plastics(TESTRASTER_GOLOUSTNOYE, name="goloustnoye", layer=4)
    quit()
