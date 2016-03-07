from numpy import *
INVALID = 0             # Data threshould to detect physiscal data in grids

def watch(x, mx, nint): # Auxiliary function to watch a computational process
    d = int(mx / nint)
    if x % d == 0:
        return 1
    return 0

def gradient(a, threshould=INVALID, default_height=0.0):
    (mx,my) = a.shape # Get grid extends
    print ("Adjusting...")
    a = where(greater(a, threshould), a, default_height) # Originally default_height=455.0

    print ("Find a stream line...")
    b  = a[:-2, 1:-1]-a[2:,1:-1]
    c  = a[1:-1,:-2]-a[1:-1,2:]

    return (0.5*b, 0.5*c)   # gx, gy are the gradient

def agrad(gx, gy): # Relief plastics accounting
    s =  gx[2:, 1:-1] - gx[:-2, 1:-1]  + gy[1:-1, 2:] - gy[1:-1, :-2]
    return s

def roundagrad(gx, gy, r): # Relief plastics accounting using a circle
    # gx[r:-r, r:-r]
    #gy[r:-r, r:-r]
    def sim(x,y, r):
        print (x,y,r)
        mx, my=gx.shape # shapes
        mx-=r
        my-=r
        #return
        s =-gx[r+x:mx+x, r:-r] * x - gy[r:-r, r+y:my+y] * y
        s+= gx[r-x:mx-x, r:-r] * x - gy[r:-r, r+y:my+y] * y
        s+= gx[r-x:mx-x, r:-r] * x + gy[r:-r, r-y:my-y] * y
        s+=-gx[r+x:mx+x, r:-r] * x + gy[r:-r, r-y:my-y] * y

        s+=-gx[r+y:mx+y, r:-r] * y - gy[r:-r, r+x:my+x] * x
        s+= gx[r-y:mx-y, r:-r] * y - gy[r:-r, r+x:my+x] * x
        s+= gx[r-y:mx-y, r:-r] * y + gy[r:-r, r-x:my-x] * x
        s+=-gx[r+y:mx+y, r:-r] * y + gy[r:-r, r-x:my-x] * x
        s =-s

        if x==0:
            s/=2
        if y==x:
            s/=2
        return s/r/8
    x=0
    y=r
    d=3-2*y
    s=None
    while x <= y :
        if s is None:
            s=sim(x,y, r)
        else:
            s+=sim(x,y, r)
        if d<0:
            d=d+4*x+6
        else:
            d=d+4*(x-y)+10;
            y-=1
        x+=1
    return s

if __name__=="__main__":
        print (roundagrad(0, 1,1, 10))
