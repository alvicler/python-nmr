#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import nmrglue as ng
import sys
import numpy as np
import pandas as pd
get_ipython().run_line_magic('matplotlib', 'qt5')
import matplotlib.pyplot as plt

#import matplotlib
#print('Python version ' + sys.version)
#print('Pandas version ' + pd.__version__)
#print('Matplotlib version ' + matplotlib.__version__)
import os
samples = len(os.listdir('urine'))


this = sys.modules[__name__] # this is now your current namespace
def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx # or array[idx]


i = 1
df = pd.DataFrame([])
yarr=[]
xarr=[]
## set the number of 1r files to open, always open first processed spectra ##
while i < samples+1:
    name="urine/Calvin_Gab_ATU"+str(i)+"/1/pdata/1"
    dic, data= ng.bruker.read_pdata(name)
    
    sf=float(dic["procs"]["SF"])    
    sfo1=float(dic["acqus"]["SFO1"])  
    o1=float(dic["acqus"]["O1"])
    hzppt=float(dic["acqus"]["SW_h"])/len(data)
    swh=float(dic["acqus"]["SW_h"])
    sr=o1+(sf-sfo1)*1000000.
    
    pts=int(sr//hzppt) # Calc pts to Calibrate 0ppm to Xscale
    data = ng.proc_base.rev(data)    # reverse the data
    
    #setattr(this, 'data%s' % i, data)
    
    #### scale x from pts to ppm ###
    ## Bin size for PCA##
    
    si=len(data) 
    xs=[]
    
    for j in range(0-pts,si-pts):
        hz=float(((o1-swh/2)+(hzppt*(j)))/sf)
        xs+=[hz]
    xs = np.asarray(xs)
    
    #setattr(this, 'xs%s' % i, xs)

    
    #xmin=xs.min()
    xmin=-.25
    #xmax=xs.max()
    xmax=1.4
    ## Bin size for PCA##
    xbin=(xmax-xmin)/5
    #xbin=.25
    k=1
    f=0
    a={}
    for j in np.arange(xmin,xmax, xbin): 
        f=j+xbin
        fpos=find_nearest(xs, f)
        jpos=find_nearest(xs, j)
        #print(jpos,fpos)
        
        peak = data[jpos:fpos]
        #peak_scale=xs[j:f]
        #if peak.sum()<0:
         #   a['slice.'+str(k)]=0 
        #else:
        #a[k]=peak.max().cumsum()
        a[k]=peak.sum()
        k+=1
    #setattr(this, 'databin%s' % i, a)
    b=pd.Series(a, name=i)
    df = df.append(b)
    yarr.append(data)
    xarr.append(xs)
    i += 1 ## index for number of spectra
    
df

