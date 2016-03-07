# -*- coding: cp1251 -*-
from Numeric import *
INVALID = 0             # порог данных, считающихся неправильными

def watch(x, mx, nint): # вспомогательная функция для обображения работы процесса
    d = int(mx / nint)
    if x % d == 0:
        return 1
    return 0
    
def stream(a, adjust = 0): # вычисление поля стоков
    (mx,my) = a.shape # узнать параметры грида (размер по x и по y)
    print "Adjusting..."
    a = where(greater(a, 0), a, 455.0)
                    
    print "Find a stream line..."
    b  = a[:-2, 1:-1]-a[2:,1:-1]
    c  = a[1:-1,:-2]-a[1:-1,2:]

    return (0.5*b, 0.5*c) # возврат поля векторов стока

def gradient(a, adjust = 0): # вычисление поля градиентов
    (mx,my) = a.shape # узнать параметры грида (размер по x и по y)
    print "Adjusting..."
    a = where(greater(a, 0), a, 455.0)
                    
    print "Find a stream line..."
    b  = a[2:,1:-1]-a[:-2, 1:-1]
    c  = a[1:-1,2:]-a[1:-1,:-2]
    return (0.5*b, 0.5*c) # возврат поля градиентов
    
def agrad(dem, gx, gy):
    s =  gx[2:, 1:-1] - gx[:-2, 1:-1]  + gy[1:-1, 2:] - gy[1:-1, :-2]
    return s

def roundagrad(dem, gx, gy, r):
    # gx[r:-r, r:-r]
    #gy[r:-r, r:-r]
    def sim(x,y, r):
        print x,y,r
        mx, my=gx.shape # shapes 
        mx-=r
        my-=r
        #return
        s=- gx[r+x:mx+x, r:-r] * x - gy[r:-r, r+y:my+y] * y
        s+=gx[r-x:mx-x, r:-r] * x - gy[r:-r, r+y:my+y] * y
        s+=gx[r-x:mx-x, r:-r] * x + gy[r:-r, r-y:my-y] * y
        s+=-gx[r+x:mx+x, r:-r] * x + gy[r:-r, r-y:my-y] * y
            
        s+=- gx[r+y:mx+y, r:-r] * y - gy[r:-r, r+x:my+x] * x
        s+=gx[r-y:mx-y, r:-r] * y - gy[r:-r, r+x:my+x] * x
        s+=gx[r-y:mx-y, r:-r] * y + gy[r:-r, r-x:my-x] * x
        s+=-gx[r+y:mx+y, r:-r] * y + gy[r:-r, r-x:my-x] * x
        s=-s
            
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
        print roundagrad(0, 1,1, 10)
        