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

if __name__=="__main__":
    from hatch import test_1
    gdal.AllRegister()

    test_1()

    quit()
