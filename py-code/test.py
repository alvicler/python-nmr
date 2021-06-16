import numpy as np

import csv
pts=[]
aq=[]
with open('test.txt') as csvDataFile:
    csvReader = csv.reader(csvDataFile)
    for row in csvReader:
        pts.append(row[0])
        aq.append(row[1])


from scipy.fft import fft, ifft, fftfreq, fftshift


# size of fid
N=len(aq)
print('FID size = ' + str(N));

#convert string into integer!
aqs = list(map(int, aq))
x = list(map(int, pts))[:N//2]

#split odd and even points of FID file
recr=aqs[0::2]
reci=aqs[1::2]

#calc FFT of real and imaginary FID
Ri=fft(recr)[:N//2]
Ii=fft(reci)[:N//2]

#Check size of everything 
#print('Ri = ' + str(len(Ri)));
#print('Ii = ' + str(len(Ii)));
#print('recr = ' + str(len(recr)));
#print('reci = ' + str(len(reci)));



# Fase correction
# acording to Journal of Magnetic Resonance 158 (2002) 164–168


phc0=(np.pi/180)*90
phc1=(np.pi/180)*-6
i=0

for ROR in Ri.real:
    
    phase=phc0+phc1*i/(N/2)
    Ri.real[i]*np.cos(phase).append()
    #R0I[i]=-Ri.imag[i]*np.sin(phase)
    #I0R[i]=Ii.real[i]*np.cos(phase) 
    #I0I[i]=Ii.imag[i]*np.sin(phase)
    i=i+1
print(Ri.real[0]*np.cos(phase))
print(ROR)
#phase=phc0+phc1*i/N

#R0R=Ri.real[0:N]*cos(phase)
#I0R=-Ri.imag[0:N]*sin(phase)

#I0R=Ii.real[0:N]*cos(phase)
#I0I=Ii.imag[0:N]*sin(phase)

#
Mx=R0R
Mxi=R0I
My=I0R
Myi=I0I


T = 1/N

print('xf= '+str(len(xf)));
t1=(My+Myi-Mx+Mxi)
t2=(Mx+Mxi+My-Myi)

import matplotlib.pyplot as plt
#plt.plot(x[60:2000], recr[60:2000])
#plt.plot(x[60:2000], reci[60:2000])
#plt.grid()
#plt.show()
plt.plot(x, t1)
plt.grid()
plt.show()
plt.plot(x, t2)
plt.grid()
plt.show()
plt.plot(x, t1+t2)
plt.grid()
plt.show()


