#!/usr/bin/env python
import nmrglue as ng
import matplotlib.pyplot as plt
plt.rcParams["figure.figsize"]=15,10
import numpy as np
import sys
this = sys.modules[__name__] # this is now your current namespace


i = 1
xs=[]

while i < 61:
    name="urine/Calvin_Gab_ATU"+str(i)+"/1/pdata/1"

    dic, data= ng.bruker.read_pdata(name)

    sf=float(dic["procs"]["SF"])
    sfo1=float(dic["acqus"]["SFO1"])
    o1=float(dic["acqus"]["O1"])
    hzppt=float(dic["acqus"]["SW_h"])/len(data)
    sr=o1+(sf-sfo1)*1000000.
    pts=int(sr//hzppt)
    data=np.roll(data, pts) # Calibrate 0ppm

    data = ng.proc_base.rev(data)               # reverse the data
    setattr(this, 'data%s' % i, data)
    i += 1

si=len(data)
swh=float(dic["acqus"]["SW_h"])
bf1=float(dic["acqus"]["BF1"])

j=0
xs=[]
for j in range(0,si):
    hz=float(((o1-swh/2)+(hzppt*j))/sfo1)
    xs+=[round(hz,5)]
xs = np.asarray(xs)
#print(type(my_array))
bin=0.05 #ppm
i=1.0

intreg=[]
xmax=xs.max()/1
xmin=xs.min()/1

print(xmin)
print(xmax)
i=1
for j in np.arange(xmin,xmax-bin, bin):
    #print(j)
    intreg+=[[i]+[j]+[j+bin]]
    i=i+1
peak_list=intreg
#for name, start, end in peak_list:
#    min = uc(start, "ppm")
#    max = uc(end, "ppm")
#    if min > max:
#        min, max = max, min

    # extract the peak
#    peak = data[min:max + 1]
#    peak_scale = ppm_scale[min:max + 1]
#    tup = (name, peak_scale[0], peak_scale[-1], peak.sum())
#    print(tup)



data=np.stack((data1,data2,data3,data4,data5,data6,data7,data8,data9,data10,
               data11,data12,data13,data14,data15,data16,data17,data18,data19,data20,
               data21,data22,data23,data24,data25,data26,data27,data28,data29,data30,
               data31,data32,data33,data34,data35,data36,data37,data38,data39,data40,
               data41,data42,data43,data44,data45,data46,data47,data48,data49,data50,
               data51,data52,data53,data54,data55,data56,data57,data58,data59,data60), axis=1)


scalemim=data.min()*0-10000
scalemax=data.max()/10

init=0
Real=data

fig = plt.figure()
ax = fig.add_subplot(111)
ax.plot(xs[init:si],Real[init:si])
#ax.plot(Real2[init:si])
# decorate axes
ax.set_yticklabels([])
ax.set_xlabel("1H ppm")
ax.invert_xaxis()
ax.set_xlim(10, -0.25)
ax.set_ylim(scalemim,scalemax)
plt.show()
#fig.savefig('figure_nmrglue.svg')
fig = plt.figure()


