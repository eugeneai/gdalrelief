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
from scipy import interpolate
from matplotlib import cm
from mpl_toolkits.mplot3d.axes3d import get_test_data
from itertools import product


TESTRASTER_GOLOUSTNOYE="../data/Goloustnoye/ALTITUDE 1Trim.grd"
TESTRASTER_OLKHON="../data/Olkhon/dem.gtiff"
SAFE_GDAL_FORMAT="GTiff"

class Hatch(object):
    """Defines section data consisting of
    x1,y1, x2,y2 are points on the section line and
    d is a distance from the line
    """

    def __init__(self, x1=None,y1=None,  x2=None,y2=None, d=5.0, p1=None, p2=None):
        if p1 != None:
            self.p1=p1
        else:
            self.p1=(x1,y1)
        if p2 != None:
            self.p2=p2
        else:
            self.p2=(x2,y2)
        self.d=d

    def copy(self):
        return self.__class__(p1=self.p1, p2=self.p2)

class RasterSection(RasterProcessor):
    """
    """

    def __init__(self, raster, hatch=None):
        RasterProcessor.__init__(self, raster)
        self.set_hatch(hatch)

    def set_hatch(self, hatch):
        self.hatch=hatch

    def __call__(self, layer):
        for x,y in self.line():
            sec.append(band[x,y])
        return np.array(sec)

    def scan_line(self, layer, minheight, p1=None, p2=None, need_bounds=False):
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
        if need_bounds:
            yield z,zu,zd,azu,azd,r,l
        else:
            yield z,zu,zd,azu,azd

    def scan(self, layer, minheight):
        save=self.hatch.copy()
        n = 5
        c = 0
        d=self.hatch.d

        for z,zu,zd,azu,azd,r,l in self.scan_line(layer, minheight, need_bounds=True):
            yield z,zu,zd,azu,azd

        if r[0]<l[0]: # wrong way
            p=r; r=l; l=p

        lu=ld=l
        ru=rd=r

        while True:

            lu=lu[0],lu[1] - d
            ld=ld[0],ld[1] + d
            ru=ru[0],ru[1] - d
            rd=rd[0],rd[1] + d

            h=Hatch(d=d,p1=lu,p2=ru)
            self.set_hatch(h)
            for z,zu,zd,azu,azd in self.scan_line(layer, minheight):
                yield z,zu,zd,azu,azd

            h=Hatch(d=d,p1=ld,p2=rd)
            self.set_hatch(h)
            for z,zu,zd,azu,azd in self.scan_line(layer, minheight):
                yield z,zu,zd,azu,azd

            c+=1
            if c>n:
                break

        self.set_hatch(save)

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
            if x<0 or y<0 or x>=mx or y>=my:
                return
            yield (x,y)

    def scan_sq(self, layer):
        band=RasterProcessor.__call__(self, layer)
        a=self.current_band=band
        h,w=a.shape
        print (h,w)
        yield from self.ss(0,0,w,h)

    def ss(self, r, c, w,h):
        """
        """
        if w == 1 and h==1:
            yield from self.check_ext(r,c)
            return
        if w>1:
            w2=w//2
        if h>1:
            h2=h//2
        if w>1 and h>1:
            yield from self.ss(r,c,w2,h2)
            yield from self.ss(r,c+w2,w-w2,h2)
            yield from self.ss(r+h2,c,w2,h-h2)
            yield from self.ss(r+h2,c+w2,w-w2,h-h2)
        elif h>1:
            yield from self.ss(r,c,w,h2)
            yield from self.ss(r+h2,c,w,h-h2)
        elif w>1:
            yield from self.ss(r,c,w2,h)
            yield from self.ss(r,c+w2,w-w2,h)
        else:
            raise RuntimeError("wrong combinations")

    DD=((0,1),(0,-1),(1,0),(-1,0),(1,1),(1,-1),(-1,1),(-1,-1))

    def mat_scan(self, layer):
        band=RasterProcessor.__call__(self, layer)
        a=self.current_band=band
        h,w=a.shape
        print (h,w)
        mi=np.ones((h,w), dtype=bool)
        ma=np.ones((h,w), dtype=bool)

        for dr,dc in self.__class__.DD:
            self.d(a,mi,True , dr, dc)
            self.d(a,ma,False, dr, dc)

        return mi, ma

    def d(self, a,m,miop, dr,dc):
        b=a[:,:]
        c=a[:,:]
        q=m[:,:]
        if dc>0:
            c=c[:-1,:]
            q=q[:-1,:]
            b=b[1:,:]
        elif dc<0:
            c=c[1:,:]
            q=q[1:,:]
            b=b[:-1,:]

        if dr>0:
            c=c[:,:-1]
            q=q[:,:-1]
            b=b[:,1:]
        else:
            c=c[:,1:]
            q=q[:,1:]
            b=b[:,:-1]
        if miop:
            d=c<b
        else:
            d=c>b
        q &= d

        # if dc>0:
        #     c=c[:-1,:]
        #     q=q[:-1,:]
        #     b=b[1:,:]
        # elif dc<0:
        #     c=c[1:,:]
        #     q=q[1:,:]
        #     b=b[:-1,:]

        # if dr>0:
        #     c=c[:,:-1]
        #     q=q[:,:-1]
        #     b=b[:,1:]
        # else:
        #     c=c[:,1:]
        #     q=q[:,1:]
        #     b=b[:,:-1]
        # print (q)
        #return q


    def check_ext(self, r,c):
        """
        """
        a=self.current_band
        h,w=a.shape
        h=hmin=hmax=a[r,c]
        if np.isnan(h):
            return
        ma=mi=True
        for dc,dr in self.__class__.DD:
            cr=r+dr
            cc=c+dc
            if cc>=w:
                continue
            if cr>=h:
                continue
            if cc<0:
                continue
            if cr<0:
                continue
            print (cr,cc)
            h=a[cr,cc]
            if np.isnan(h):
                continue

            if h>hmax:
                ma=False
            if h<hmin:
                mi=False

        if ma:
            yield (True, r,c,hmax) # True means max
        if mi:
            yield (False, r,c,hmin) # False means min

# TEST

def test_1():
    # h=Hatch(100,100, 732,684, 5)
    h=Hatch(3,400, 700, 300, 50)
    rs=RasterSection(raster=TESTRASTER_GOLOUSTNOYE, hatch=h)
    rs.info()

    au=None
    ad=None
    for z,zu,zd,u,d in rs.scan(4, 20):
        #print ("===", u,d)
        if au != None:
            au = np.append(au, u, axis=0)
        else:
            au=u
        if ad != None:
            ad = np.append(ad, d, axis=0)
        else:
            ad=d

    print (au, ad)

    fu = interpolate.interp2d(au[:,0], au[:,0], au[:,0], kind='cubic')
    fd = interpolate.interp2d(ad[:,0], ad[:,0], ad[:,0], kind='cubic')

    if False:
        height=z
        valid=height

        x1=np.arange(len(valid))  #np.linspace(0,832)
        y1=valid

        plt.plot(x1, y1, 'r.-')
        plt.plot(zu, z[zu], 'g.-')
        plt.plot(zd, z[zd], 'g.-')
        #plt.plot(x2, y2, 'r.-')
        plt.show()

    if True:
        w,h=rs.current_band.shape
        prec=10
        xnew = np.arange(0, w, prec)
        ynew = np.arange(0, h, prec)
        znew = fd(xnew, ynew) # FIXME, please, I'm very stupid copypaster.
        xnew, ynew = np.meshgrid(xnew, ynew)

        fig = plt.figure(figsize=plt.figaspect(1))

        #---- First subplot
        ax = fig.add_subplot(1, 1, 1, projection='3d')

        surf = ax.plot_surface(xnew, ynew, znew, rstride=1, cstride=1, cmap=cm.coolwarm,
                               linewidth=0, antialiased=True)
        ax.set_zlim3d(420, 1500)

        fig.colorbar(surf, shrink=0.5, aspect=10)

        #---- Second subplot
        #ax = fig.add_subplot(1, 2, 2, projection='3d')
        #X, Y, Z = get_test_data(0.05)
        #ax.plot_wireframe(X, Y, Z, rstride=10, cstride=10)

        plt.show()

def test_2():
    # h=Hatch(100,100, 732,684, 5)
    h=Hatch(3,400, 700, 300, 50)
    rs=RasterSection(raster=TESTRASTER_GOLOUSTNOYE)
    rs.info()

    mi,ma=rs.mat_scan(4)
    a=rs.current_band
    h,w=a.shape
    print(mi,ma)

    au=[]
    ad=[]

    for r,c in product(range(h),range(w)):
        z=a[r,c]
        if mi[r,c]==True:
            ad.append((c,r,z))
        if ma[r,c]==True:
            au.append((c,r,z))

    au=np.array(au)
    ad=np.array(ad)

    interpol='cubic'
    fu = interpolate.interp2d(au[:,0], au[:,0], au[:,0], kind=interpol)
    fd = interpolate.interp2d(ad[:,0], ad[:,0], ad[:,0], kind=interpol)

    f=fu
    if True:
        w,h=rs.current_band.shape
        prec=10
        xnew = np.arange(0, w, prec)
        ynew = np.arange(0, h, prec)
        znew = f(xnew, ynew) # FIXME, please, I'm very stupid copypaster.
        xnew, ynew = np.meshgrid(xnew, ynew)

        fig = plt.figure(figsize=plt.figaspect(1))

        #---- First subplot
        ax = fig.add_subplot(1, 1, 1, projection='3d')

        surf = ax.plot_surface(xnew, ynew, znew, rstride=1, cstride=1, cmap=cm.coolwarm,
                               linewidth=0, antialiased=True)
        ax.set_zlim3d(420, 1500)

        fig.colorbar(surf, shrink=0.5, aspect=10)

        #---- Second subplot
        #ax = fig.add_subplot(1, 2, 2, projection='3d')
        #X, Y, Z = get_test_data(0.05)
        #ax.plot_wireframe(X, Y, Z, rstride=10, cstride=10)

        plt.show()






if __name__=="__main__":
    # register all of the GDAL drivers
    gdal.AllRegister()

    test_2()

    quit()
