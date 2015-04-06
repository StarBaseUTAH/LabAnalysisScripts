import numpy as np
from matplotlib import pyplot as plt 
from scipy.optimize import curve_fit

def line(x,m,b):
    return m*x+b

def sqrt(x,A):
    return A*np.sqrt(x)

e=1.619e-19
Gamp=5000.
f=1.e6 #Frequency of pulses used

#Enter in values from extractPhotonNumber.py
voltage=['900V']
avgInt=np.array([1.26,3.10,3.25,5.22,7.03,9.01,10.98,14.77,16.65])
sigInt=np.array([0.229,0.295,0.356,0.416,0.469,0.523,0.584,0.632,0.674])
Vdc=1e-3*np.array([4.91,12.6,12.6,20.7,28.1,36.0,43.9,59.7,67.0]) #(V)
Nph=np.array([30.5,83.1,110.6,156.9,224.8,296.7,353.6,546.2,610.4])*f #ph/sec

#First check to see if you have reasonable gain results
a=np.polyfit(Vdc,Nph,1)
Gpmt=1./(e*Gamp*a[0])
print("Calculated Gain of the PMT: ",Gpmt)

#Fit sqrt function to Vdc vs. Deviation in Integral Vals
sqfit,sqsuccess=curve_fit(sqrt,Vdc,sigInt)
xvalss=np.arange(0,np.max(Vdc),0.0001)
theorySqrt=sqrt(xvalss,sqfit[0])
plt.figure()
plt.plot(xvalss*1e3,theorySqrt)
plt.plot(Vdc*1e3,sigInt,'.')
plt.xlabel('DC Voltage (mV)')
plt.ylabel('Deviation in Trace Integrals (a.u.)')

xvals=np.arange(0,np.max(Vdc),0.00001)
theoryCalibration=line(xvals,a[0],a[1])
plt.figure()
plt.plot(xvals*1e3,theoryCalibration*1e-9)
plt.plot(Vdc*1e3,Nph*1e-9,'.') 
plt.xlabel('DC Voltage (mV)')
plt.ylabel('Photon Count Rate (ph/s) * $10^9$')
plt.title('Photon Rate Estimator for HV = '+voltage[0])
plt.show()
