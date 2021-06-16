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

# In[3]:


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
    
    setattr(this, 'data%s' % i, data)
    
    #### scale x from pts to ppm ###
    ## Bin size for PCA##
    
    si=len(data) 
    xs=[]
    
    for j in range(0-pts,si-pts):
        hz=float(((o1-swh/2)+(hzppt*(j)))/sf)
        xs+=[hz]
    xs = np.asarray(xs)
    
    setattr(this, 'xs%s' % i, xs)

    
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
    setattr(this, 'databin%s' % i, a)
    b=pd.Series(a, name=i)
    df = df.append(b)
    yarr.append(data)
    xarr.append(xs)
    i += 1 ## index for number of spectra
    
df


# # PCA
# 
# Runs PCA itself
# 
# To run ICA change principalComponents = **pca**.fit_transform(x) to principalComponents = **ica**.fit_transform(x)
# 
# to unscale data uncomment x=df and comment line x = StandardScaler().fit_transform(df)
# 
# if you desired to change scale, change x = **StandardScaler()**.fit_transform(df) funcition to **MinMaxScaler()** or **RobustScaler()**
# 
# 

# In[4]:


# Enabling the `widget` backend.
# This requires jupyter-matplotlib a.k.a. ipympl.
# ipympl can be install via pip or conda.
get_ipython().run_line_magic('matplotlib', 'qt5')
import matplotlib.pyplot as plt

plt.rcParams["figure.figsize"]=15,10
#plt.fig.canvas.header_visible = False # Hide the Figure name at the top of the figure
#from sklearn.decomposition import PCA
from sklearn import decomposition
from sklearn import datasets
from sklearn.decomposition import FastICA, PCA
from sklearn.preprocessing import MinMaxScaler
from sklearn.preprocessing import StandardScaler
from sklearn.preprocessing import RobustScaler

#df2=df
#df2.to_csv (r'dataframe.csv',  header=True)# Compute ICA
ica = FastICA(n_components=5)


x = StandardScaler().fit_transform(df)
#x=df
pca = decomposition.PCA(n_components=5)

principalComponents = pca.fit_transform(x)
principalDf = pd.DataFrame(data = principalComponents
             , columns = ['PC1', 'PC2','PC3', 'PC4','PC5'])
#print(principalDf['principal component 1'])
#fig = plt.figure(figsize = (8,8))
fig, ([ax1, ax2], [ax3, ax4],[ax5,ax6],[ax7,ax8],[ax9,ax10]) = plt.subplots(5, 2,sharex=True, figsize=(10,10))
#ax1.set_xlabel('PC1')
#ax2.set_xlabel('PC1')
#ax3.set_xlabel('PC1')
#ax4.set_xlabel('PC1')

#ax1.set_ylabel('PC2')
#ax2.set_ylabel('PC3')
#ax3.set_ylabel('PC4')
#ax4.set_ylabel('PC5')

ax1.set_title('PC1 x PC2')
ax2.set_title('PC1 x PC3')
ax3.set_title('PC1 x PC4')
ax4.set_title('PC1 x PC5')
ax5.set_title('PC2 x PC3')
ax6.set_title('PC2 x PC4')
ax7.set_title('PC2 x PC5')
ax8.set_title('PC3 x PC4')
ax8.set_title('PC3 x PC5')
ax8.set_title('PC4 x PC5')

i=1
for x,y,z,w,u in zip(principalDf['PC1'], principalDf['PC2'],principalDf['PC3'],principalDf['PC4'],principalDf['PC5']):
    label = f"{i}"
    ax1.scatter(x,y,s=30,marker="o",facecolors='none', edgecolors='darkorange')
    ax1.annotate(label, # this is the text
                 (x,y), # this is the point to label
                 textcoords="offset points", # how to position the text
                 xytext=(0,10), # distance from text to points (x,y)
                 ha='center', # horizontal alignment can be left, right or center
                 c='bisque') 
    ax2.scatter(x,z,s=30,marker="o",facecolors='none', edgecolors='r')
    ax2.annotate(label, # this is the text
                 (x,z), # this is the point to label
                 textcoords="offset points", # how to position the text
                 xytext=(0,10), # distance from text to points (x,y)
                 ha='center', # horizontal alignment can be left, right or center
                 c='r') 
    ax3.scatter(x,w,s=30,marker="o",facecolors='none', edgecolors='g')
    ax3.annotate(label, # this is the text
                 (x,w), # this is the point to label
                 textcoords="offset points", # how to position the text
                 xytext=(0,10), # distance from text to points (x,y)
                 ha='center', # horizontal alignment can be left, right or center
                 c='g') 
    ax4.scatter(x,u,s=50,marker="o",facecolors='none', edgecolors='k')
    ax4.annotate(label, # this is the text
                 (x,u), # this is the point to label
                 textcoords="offset points", # how to position the text
                 xytext=(0,10), # distance from text to points (x,y)
                 ha='center', # horizontal alignment can be left, right or center
                 c='k') 
    ax5.scatter(y,z,s=30,marker="o",facecolors='none', edgecolors='c')
    ax5.annotate(label, # this is the text
                 (y,z), # this is the point to label
                 textcoords="offset points", # how to position the text
                 xytext=(0,10), # distance from text to points (x,y)
                 ha='center', # horizontal alignment can be left, right or center
                 c='c') 
    ax6.scatter(y,w,s=30,marker="o",facecolors='none', edgecolors='m')
    ax6.annotate(label, # this is the text
                 (y,w), # this is the point to label
                 textcoords="offset points", # how to position the text
                 xytext=(0,10), # distance from text to points (x,y)
                 ha='center', # horizontal alignment can be left, right or center
                 c='m') 
    ax7.scatter(y,u,s=30,marker="o",facecolors='none', edgecolors='r')
    ax7.annotate(label, # this is the text
                 (y,u), # this is the point to label
                 textcoords="offset points", # how to position the text
                 xytext=(0,10), # distance from text to points (x,y)
                 ha='center', # horizontal alignment can be left, right or center
                 c='g') 
    ax8.scatter(z,w,s=30,marker="o",facecolors='none', edgecolors='g')
    ax8.annotate(label, # this is the text
                 (z,w), # this is the point to label
                 textcoords="offset points", # how to position the text
                 xytext=(0,10), # distance from text to points (x,y)
                 ha='center', # horizontal alignment can be left, right or center
                 c='r') 
    ax9.scatter(y,u,s=30,marker="o",facecolors='none', edgecolors='b')
    ax9.annotate(label, # this is the text
                 (y,u), # this is the point to label
                 textcoords="offset points", # how to position the text
                 xytext=(0,10), # distance from text to points (x,y)
                 ha='center', # horizontal alignment can be left, right or center
                 c='g') 
    color = [str(item/255.) for item in principalDf['PC4']]
    ax10.scatter(z,w, marker="o",facecolors='none', edgecolors='c')
    ax10.annotate(label, # this is the text
                 (z,w), # this is the point to label
                 textcoords="offset points", # how to position the text
                 xytext=(0,10), # distance from text to points (x,y)
                 ha='center', # horizontal alignment can be left, right or center
                 c='m') 
    #ax10.plt.colorbar()
    i+=1



plt.show()


# # Stackplot of all sample spectrum

# In[31]:


### Calc Collors for nmr spectra
plt.rcParams["figure.figsize"]=15,10
from cycler import cycler
NUM_COLORS = 25
txt=[]
cm = plt.get_cmap('jet')
for i in range(NUM_COLORS):
    colr=cm(i//3*4/NUM_COLORS)
    txt.append(colr)
plt.rcParams['axes.prop_cycle'] = cycler(color=txt)



xarr=np.array(xarr)
yarr=np.array(yarr)
arrx = xarr.reshape(samples,si).T
arry = yarr.reshape(samples,si).T
scalemim=data.min()/100
scalemax=data.max()/0.75
label= list(range(1, samples+1))

init=0


fig = plt.figure()
ax = fig.add_subplot(111)

ax.plot(arrx[init:si],arry[init:si])

# decorate axes
ax.set_yticklabels([])
ax.set_xlabel("$^{1}H\ ppm$")
ax.set_title("Stack Plot of $"+str(samples)+"$ NMR Urine Samples")
ax.invert_xaxis()
ax.set_xlim(10, -0.5)
ax.set_ylim(scalemim,scalemax)
ax.legend(label,ncol=10,loc="upper right")
plt.grid(color='gray', linestyle='-.', linewidth=.25)
#fig.savefig('figure_nmrglue.svg')


# # Scaled PCA with ellipsis of 95%
# 
# * Change **perc** from 95% (0.95) to your choice

# In[32]:


import scipy, random
data = df
pca = decomposition.PCA(n_components = 5)
scaler = StandardScaler()
scaler.fit(data)
data = scaler.transform(data)
pcaFit = pca.fit(data)
dataProject = pcaFit.transform(data)

perc=0.9500
#Calculate ellipse bounds and plot with scores
theta = np.concatenate((np.linspace(-np.pi, np.pi, 10), np.linspace(np.pi, -np.pi, 10)))
circle = np.array((np.cos(theta), np.sin(theta)))
sigma = np.cov(np.array((dataProject[:, 0], dataProject[:, 1])))
ed = np.sqrt(scipy.stats.chi2.ppf(perc, 2))
ell = np.transpose(circle).dot(np.linalg.cholesky(sigma) * ed)
a, b = np.max(ell[: ,0]), np.max(ell[: ,1]) #95% ellipse bounds
t = np.linspace(0, 2 * np.pi, 100)


fig = plt.figure(figsize = (8,8))
ax = fig.add_subplot() 
ax.set_title('PCA - '+str(perc*100)+' %', fontsize = 20)

ax.set_xlabel('PC1', fontsize = 15)
ax.set_ylabel('PC2', fontsize = 15)

ax.plot(a * np.cos(t), b * np.sin(t), color = 'red')
ax.grid(color = 'lightgray', linestyle = '--')

i=1
for x,y in zip(dataProject[:, 0], dataProject[:, 1]):
    ax.scatter(x, y,s=50,marker="o",facecolors='none', edgecolors='b')
    label = f"{i}"
    ax.annotate(label, # this is the text
                 (x,y), # this is the point to label
                 textcoords="offset points", # how to position the text
                 xytext=(0,10), # distance from text to points (x,y)
                 ha='center', # horizontal alignment can be left, right or center
                 c='k') 
    i+=1


plt.show()


# # 3D PCA (PC1,2,3) PLOT

# In[14]:


# libraries
import pandas as pd
import numpy as np
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
import seaborn as sns
 

# Run The PCA
pca = PCA(n_components=3)
pca.fit(df)
 
# Store results of PCA in a data frame
result=pd.DataFrame(pca.transform(df), columns=['PCA%i' % i for i in range(3)], index=df.index)
 
# Plot initialisation
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(result['PCA0'], result['PCA1'], result['PCA2'], cmap="Set2_r", s=60)
 
# make simple, bare axis lines through space:
xAxisLine = ((min(result['PCA0']), max(result['PCA0'])), (0, 0), (0,0))
ax.plot(xAxisLine[0], xAxisLine[1], xAxisLine[2], 'r')
yAxisLine = ((0, 0), (min(result['PCA1']), max(result['PCA1'])), (0,0))
ax.plot(yAxisLine[0], yAxisLine[1], yAxisLine[2], 'r')
zAxisLine = ((0, 0), (0,0), (min(result['PCA2']), max(result['PCA2'])))
ax.plot(zAxisLine[0], zAxisLine[1], zAxisLine[2], 'r')
 
# label the axes
ax.set_xlabel("PC1")
ax.set_ylabel("PC2")
ax.set_zlabel("PC3")
ax.set_title("PCA on the iris data set")
plt.show()


# # PLOT Gaussian, Lorentzian and Voigt lines

# In[33]:


import numpy as np
from scipy.special import wofz
import pylab

def G(x, alpha):
    """ Return Gaussian line shape at x with HWHM alpha """
    return np.sqrt(np.log(2) / np.pi) / alpha                             * np.exp(-(x / alpha)**2 * np.log(2))

def L(x, gamma):
    """ Return Lorentzian line shape at x with HWHM gamma """
    return gamma / np.pi / (x**2 + gamma**2)

def V(x, alpha, gamma):
    """
    Return the Voigt line shape at x with Lorentzian component HWHM gamma
    and Gaussian component HWHM alpha.

    """
    sigma = alpha / np.sqrt(2 * np.log(2))

    return np.real(wofz((x + 1j*gamma)/sigma/np.sqrt(2))) / sigma                                                           /np.sqrt(2*np.pi)

alpha, gamma = 0.1, 0.1
x = np.linspace(-0.8,0.8,1000)
pylab.plot(x, G(x, alpha), ls=':', label='Gaussian')
pylab.plot(x, L(x, gamma), ls='--', label='Lorentzian')
pylab.plot(x, V(x, alpha, gamma), label='Voigt')
pylab.xlim(-0.8,0.8)
pylab.legend()
pylab.show()


# # Align 2 spectra 

# In[16]:


### Align region of spec based on max of one peak
import pandas as pd
import numpy as np
from scipy import signal
from matplotlib import pylab as plt
from scipy.signal import chirp, find_peaks, peak_widths

# regiao espectral
f1=4.0
f2=4.1
# regiao em pts espectro 1
f1xpos=find_nearest(xs37, f1)
f2xpos=find_nearest(xs37, f2)

# maximos naquela regiao
f1h=data37[f1xpos:f2xpos].max()
f2h=data64[f1xpos:f2xpos].max()


# posicao em pts destes maximos
f1pos=find_nearest(data37, f1h)
f2pos=find_nearest(data64, f2h)

#calcula em ppm a distancia dos pontos
shft01=xs37[f2pos]-xs37[f1pos]
#print(shft01)
peaks, _ = find_peaks(data37[f1xpos:f2xpos])
p2=xs37[f1xpos:f2xpos]
p2=p2[peaks]
#print(p2)
results_half = peak_widths(data37[f1xpos:f2xpos], peaks, rel_height=0.5)

results_full = peak_widths(data37[f1xpos:f2xpos], peaks, rel_height=1)
alpha, gamma = 0.1, 0.002


p3=data37[f1xpos:f2xpos]
p3=p3[peaks]
maxval=p3.max()
print(maxval)
maxpts=p3.argmax()
print(p2[6])
x = np.linspace(-0.8*p2,0.8*p2,100000)
#muda o eixo x do segundo espectro
xss=xs37-shft01

plt.plot(xs37[f1xpos:f2xpos], data37[f1xpos:f2xpos], label='sample 37')
plt.plot(xs64[f1xpos:f2xpos], data64[f1xpos:f2xpos], label='sample 64')
plt.plot(xss[f1xpos:f2xpos], data64[f1xpos:f2xpos], ls='dashed', label='aligned')
plt.plot(p2, p3, "o", label='peak')
plt.plot(x+p2,L(x, gamma)*p3/160, ls="--", label='Sim.sample 37 Lorentzian')
plt.xlim(3.95, 4.15)
plt.legend()


# # Old Stack method, manual inserction

# In[17]:


data=np.stack((data1,data2,data3,data4,data5,data6,data7,data8,data9,data10,
               data11,data12,data13,data14,data15,data16,data17,data18,data19,data20,
               data21,data22,data23,data24,data25,data26,data27,data28,data29,data30,
               data31,data32,data33,data34,data35,data36,data37,data38,data39,data40,
               data41,data42,data43,data44,data45,data46,data47,data48,data49,data50,
               data51,data52,data53,data54,data55,data56,data57,data58,data59,data60,
               data61,data62,data63,data64,data65,data66,data67,data68,data69,data70,
               data71,data72,data73,data74,data75,data76,data77,data78,data79,data80,
               data81,data82,data83,data84,data85,data86,data87,data88,data89,data90,
               data91,data92,data93,data94,data95,data96,data97,data98,data99,data100), axis=1)

xs=np.stack((xs1,xs2,xs3,xs4,xs5,xs6,xs7,xs8,xs9,xs10,
               xs11,xs12,xs13,xs14,xs15,xs16,xs17,xs18,xs19,xs20,
               xs21,xs22,xs23,xs24,xs25,xs26,xs27,xs28,xs29,xs30,
               xs31,xs32,xs33,xs34,xs35,xs36,xs37,xs38,xs39,xs40,
               xs41,xs42,xs43,xs44,xs45,xs46,xs47,xs48,xs49,xs50,
               xs51,xs52,xs53,xs54,xs55,xs56,xs57,xs58,xs59,xs60,
               xs61,xs62,xs63,xs64,xs65,xs66,xs67,xs68,xs69,xs70,
               xs71,xs72,xs73,xs74,xs75,xs76,xs77,xs78,xs79,xs80,
               xs81,xs82,xs83,xs84,xs85,xs86,xs87,xs88,xs89,xs90,
               xs91,xs92,xs93,xs94,xs95,xs96,xs97,xs98,xs99,xs100),axis=1)




scalemim=data.min()/100
scalemax=data.max()/0.75

init=0


fig = plt.figure()
ax = fig.add_subplot(111)

ax.plot(xs[init:si],data[init:si])

# decorate axes
#ax.set_yticklabels([])
ax.set_xlabel("$^{1}H\ ppm$")
ax.invert_xaxis()
ax.set_xlim(10, -0.5)
ax.set_ylim(scalemim,scalemax)
#fig.savefig('figure_nmrglue.svg')


# In[18]:


import numpy as np
import scipy.signal as signal
import matplotlib.pyplot as plt

fs = si*2

amp = 2*np.sqrt(2)
freq = 1270
noise_power = 0.1 
time = xs37
x = amp*np.sin(2*np.pi*freq*time)
x=data37
x += np.random.normal(scale=np.sqrt(noise_power), size=time.shape)

f, Pper_spec = signal.periodogram(x, fs, 'flattop', scaling='spectrum')

plt.semilogy(f, Pper_spec)
plt.xlabel('frequency [Hz]')
plt.ylabel('PSD')
plt.grid()
plt.show()


# In[19]:


## Plot real component, imaginary component, and envelope for a 5 Hz pulse, sampled at 100 Hz for 2 seconds:

from scipy import signal
import matplotlib.pyplot as plt
t = np.linspace(-1, 1, 2 * 500, endpoint=False)
i, q, e = signal.gausspulse(t, fc=5, retquad=True, retenv=True)
plt.plot(t, i, t, q, t, e, '--')


# In[20]:


from scipy import signal
from scipy.fft import fft, fftshift
import matplotlib.pyplot as plt

M = si
tau = 3.0
window = data37
#plt.plot(window)
plt.title("Exponential Window (tau=3.0)")
plt.ylabel("Amplitude")
plt.xlabel("Sample")

A = fft(window, si) / (len(window)/2.0)
freq = np.linspace(-0.5, 0.5, len(A))
response = 20 * np.log10(np.abs(fftshift(A / abs(A).max())))
plt.plot(freq, response)


# In[ ]:


import matplotlib.pyplot as plt
import numpy as np

NUM_COLORS = 101

cm = plt.get_cmap('tab20c_r')
fig = plt.figure()
ax = fig.add_subplot(111)
for i in range(NUM_COLORS):
    lines = ax.plot(np.arange(10)*(i+1))
    lines[0].set_color(cm(i//3*3.0/NUM_COLORS))
    #lines[0].set_linewidth(i%3 + 1)
    print(cm(i//3*3.0/NUM_COLORS),',')
#fig.savefig('moreColors.png')
plt.show()


# In[ ]:


#from bokeh.plotting import figure, show
# create a new plot with a title and axis labels
p = figure(title="Simple line example", x_axis_label='ppm', y_axis_label='u.a',x_range=(10, -1),y_range=(-1e5, 1e6),)
# add a line renderer with legend and line thickness to the plot
p.line(arrx, arry, legend_label="Temp.", line_width=2)
# show the results
show(p)


# In[ ]:




