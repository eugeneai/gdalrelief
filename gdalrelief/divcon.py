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

# register all of the GDAL drivers

TESTRASTER_GOLOUSTNOYE="../data/Goloustnoye/ALTITUDE 1Trim.grd"
TESTRASTER_OLKHON="../data/Olkhon/dem.gtiff"
SAFE_GDAL_FORMAT="GTiff"


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
    # register all of the GDAL drivers
    gdal.AllRegister()

    script_plastics(TESTRASTER_OLKHON, name="olkhon", layer=1)
    script_plastics(TESTRASTER_GOLOUSTNOYE, name="goloustnoye", layer=4)
    quit()
