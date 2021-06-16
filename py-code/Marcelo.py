#!/usr/bin/env python
# coding: utf-8

# # Open 1r
# 
# Open all 1r @ samples with **exp=1** and **process=1** using **NMRGLUE PACKAGE**
# Reference all 1r file from **procs** and **acqus** from sample directory
# 
# ## PCA
# 
# Divide data into range, for PCA 
# * change **xmin/xmax** in **PPM** for the wanted range
# 
# **xbin** was calc for **5 points** in this range, means no more then **5 components at PCA was allowed**
# 
# * Change it if you want to (xbin=(xmax-xmin)/5
# The final data was send to ***df*** Variable
# 
# 

# In[55]:


import nmrglue as ng
import sys
import numpy as np
import pandas as pd
get_ipython().run_line_magic('matplotlib', 'qt5')
import matplotlib.pyplot as plt
import os
samples = os.listdir('Marcelo')


this = sys.modules[__name__] # this is now your current namespace
def find_nearest(array, value):
    array = np.asarray(array)
    idx = (np.abs(array - value)).argmin()
    return idx # or array[idx]


i = 0
df = pd.DataFrame([])
df1 =  pd.DataFrame([])
yarr=[]
xarr=[]
names=[]
## set the number of 1r files to open, always open first processed spectra ##
while i < len(samples):
    name="Marcelo/"+samples[i]+"/1/pdata/1"

    dic, data= ng.bruker.read_pdata(name)
    
    sf=float(dic["procs"]["SF"])    
    sfo1=float(dic["acqus"]["SFO1"])  
    o1=float(dic["acqus"]["O1"])
    hzppt=float(dic["acqus"]["SW_h"])/len(data)
    swh=float(dic["acqus"]["SW_h"])
    sr=o1+(sf-sfo1)*1000000.
    
    pts=int(sr//hzppt) # Calc pts to Calibrate 0ppm to Xscale
    data = ng.proc_base.rev(data)    # reverse the data
    
    #### scale x from pts to ppm ###
    ## Bin size for PCA##
    
    si=len(data) 
    xs=[]
    
    for j in range(0-pts,si-pts):
        hz=float(((o1-swh/2)+(hzppt*(j)))/sf)
        xs+=[hz]
    xs = np.asarray(xs)
    
    yarr.append(data)
    xarr.append(xs)
    names.append(name)
    
    df = pd.read_table(name+'/integrals.txt', delim_whitespace=True, skiprows = (1,3), index_col=0, header=2)
  
    
    if i == 0:
        df1=df
    else:
        df1=df1.append(df)
    i += 1 ## index for number of spectra

df1


#df1.to_excel("output.xlsx")  


# In[2]:


### Calc Collors for nmr spectra
plt.rcParams["figure.figsize"]=15,10
from cycler import cycler
NUM_COLORS = 10
txt=[]
cm = plt.get_cmap('jet')
for i in range(NUM_COLORS):
    colr=cm(i//3*4/NUM_COLORS)
    txt.append(colr)
plt.rcParams['axes.prop_cycle'] = cycler(color=txt)

samples = len(os.listdir('Marcelo'))

xarr=np.array(xarr)
yarr=np.array(yarr)
arrx = xarr.reshape(samples,si).T
arry = yarr.reshape(samples,si).T
scalemim=data.min()-1000000
scalemax=data.max()/0.75
label= list(range(1, samples+1))

init=0


fig = plt.figure()
ax = fig.add_subplot(111)

ax.plot(arrx[init:si],arry[init:si])
label=names
# decorate axes
#ax.set_yticklabels([])
ax.set_xlabel("$^{1}H\ ppm$")
ax.set_title("Stack Plot of $"+str(samples)+"$ NMR Urine Samples")
ax.invert_xaxis()
ax.set_xlim(14, -0.5)
ax.set_ylim(scalemim,scalemax)
ax.legend(label,ncol=1,loc="upper left")
plt.grid(color='gray', linestyle='-.', linewidth=.25)
#fig.savefig('figure_nmrglue.svg')


# In[78]:


name="Marcelo/Gasoleo-100920-bp75/4/pdata/1"

dic, data= ng.bruker.read_pdata(name)
    
sf=float(dic["procs"]["SF"])    
sfo1=float(dic["acqus"]["SFO1"])  
o1=float(dic["acqus"]["O1"])
hzppt=float(dic["acqus"]["SW_h"])/len(data)
swh=float(dic["acqus"]["SW_h"])
sr=o1+(sf-sfo1)*1000000.
    
pts=int(sr//hzppt) # Calc pts to Calibrate 0ppm to Xscale
data = ng.proc_base.rev(data)    # reverse the data
    
#### scale x from pts to ppm ###
## Bin size for PCA##
si=len(data) 
xs=[]
for j in range(0-pts,si-pts):
        hz=float(((o1-swh/2)+(hzppt*(j)))/sf)
        xs+=[hz]
xs = np.asarray(xs)
    
arry=np.array(data).T
arrx=np.array(xs).T
names=name

init=0
factor=2
fig = plt.figure()
ax = fig.add_subplot(111)

ax.plot(arrx[init:si],arry[init:si]/factor,linewidth=.35)
label=names
# decorate axes
#ax.set_yticklabels([])
# Add labels to the plot
anote=[
    ['furane',[109.6094,143.0968]],
    ['2-Methylfurane',[152.6181, 152.6181, 111.555,106.3261,13.2976]],
    ['2-Ethylfurane',[158.2671, 141.6036,110.9273, 104.7653, 21.4574, 12.7028]],
    ['2,5-Dimethylfurane',[149.1741, 105.3859, 12.8116]],
    ['2,3,5-Trimethylfuran',[145.4946,142.8696,114.7653,109.1919,13.3763,11.1383,9.8918]],
    ['Benzofuran',[154.6271,144.6128,127.3870,124.0520,122.5424,120.8599,111.2968,106.6386]],
    ['3-Furaldehyde *****',[185.3,153.1,146.2,129.9,107.5]],
    ['2-Ethyl-5-methylfuran *****',[156.44,140.66,110.01,104.67,30.02,21.40, 13.72 ]],
    ['Furfural',[178.2550,154.1614,149.3534,122.4469,113.5815]]
    ]
#print(anote[0][1])    
i=0
for note,w in anote[:][:]:
    #print(z) 
    for x in w:
        #print(x,y)
        peakidx=find_nearest(arrx, x)
        y=arry[peakidx]/factor*1.2
        ax.annotate(" ------>   "+note, ha='center', xy=(x,y),
                    color='k', rotation='vertical');
    i+=1 


ax.axvspan(30.3796, 29.4559,color='r', alpha=0.2)# Aceton
ax.axvspan(205.967, 206.707,color='r', alpha=0.2)# Aceton

ax.set_xlabel("$^{1}H\ ppm$")
ax.set_title("Plot of $"+str(names)+"$ NMR")
ax.invert_xaxis()
ax.set_xlim(215, -5)
ax.set_ylim(scalemim,scalemax)
#ax.legend(label,loc="upper left")
plt.grid(color='gray', linestyle='-.', linewidth=.25)
#fig.savefig('figure_nmrglue.svg')


# In[ ]:





# In[ ]:




