import numpy as np
from matplotlib import pyplot as plt 
from scipy.optimize import curve_fit
from scipy.integrate import simps
from scipy.stats import mode

def resp(t,A,m,sigR,b):
    if m<0:
        m=1e6
    m_fact=np.math.factorial(np.int(m))
    return A*(np.sqrt(m+1)/(m_fact*sigR))*(np.sqrt(m+1)*t/sigR)**m * np.exp(-1*np.sqrt(m+1)*t/sigR)+b

filename="15Msamp_950mV_1100V" #Filename to be analyzed
print("Working on filename: ",filename)
searchRange=400. #Number of samples to look for first pulse
numberOfPulses=10000.

f_led=1.0e6 #Frequency of full pulse signal (Hz)
sampRate=250.e6 #Sampling Rate 
f_ledC=sampRate/f_led #Number of samples per full pulse 
pulseWidth=100. #Number of samples per pulse width. Follows duty cycle
extWidth=pulseWidth/0.9 #Extraction width around each pulse. Decided that 90% should contain actual pulse..
binNumbers=50. #Number of bins to use for histogram

#Load in time series data
timeSeries=np.loadtxt(filename)
timeSeries=np.reshape(timeSeries,np.size(timeSeries))
tvals=np.arange(0,np.size(timeSeries))

#Identify baseline. Use weighted mean of most frequent values
hist=np.histogram(timeSeries)
numbMostFreq=7 #Number of points to use to determine baseline
mcVals,numbOcc=np.zeros(numbMostFreq),np.zeros(numbMostFreq)
copyArr=timeSeries
for i in np.arange(0,numbMostFreq):
    mcVals[i],numbOcc[i]=mode(copyArr)[0][0],mode(copyArr)[1][0]
    copyArr=copyArr[copyArr!=mcVals[i]]

baseline=np.sum(mcVals*numbOcc)/np.sum(numbOcc)
timeSeries=timeSeries-baseline

#Identify first pulse by min. amplitude in first amount of samples where there should be a pulse
firstMinIndex=np.argmin(timeSeries[searchRange:2*searchRange])
firstMinIndex=firstMinIndex+searchRange

#Intial guesses for fit
Ainit=timeSeries[firstMinIndex]
sigRinit=1.0
minit=8.

#Extract region around first pulse and fit
lowCutoff=np.int(firstMinIndex-extWidth/2)
highCutoff=np.int(firstMinIndex+extWidth/2+1)
firstPulse=timeSeries[lowCutoff:highCutoff]
timeIndex=np.arange(lowCutoff,highCutoff)
distArray=np.zeros((numberOfPulses,np.size(firstPulse)-1))

#Display and check first pulse
plt.plot(firstPulse)
plt.title('Could this be a LED pulse/baseline correct(=0)?')
plt.show()

#Iterate through data finding pulse integrals.
integralVals=np.zeros(numberOfPulses)
for i in np.arange(0,numberOfPulses): 
    #Find region containing next pulse
    if i==0:
        currentPulse=firstPulse
    currentMinIndex=np.argmin(timeSeries[lowCutoff:highCutoff])+lowCutoff
    lowCutoff=np.int(currentMinIndex+f_ledC-extWidth/2)
    highCutoff=np.int(currentMinIndex+f_ledC+extWidth/2)
        #print np.size(timeSeries[lowCutoff:highCutoff])
        #        plt.plot(timeSeries[lowCutoff:highCutoff])
        #        plt.show()
    integralVals[i]=simps(timeSeries[lowCutoff:highCutoff])
    #print(np.size(timeSeries[lowCutoff:highCutoff]))
    distArray[i,:]=timeSeries[lowCutoff:highCutoff]

histSetup=np.histogram(integralVals,bins=binNumbers)
intXVALS=histSetup[1][0:np.size(histSetup[0])]
intHVALS=histSetup[0]
print("Average pulse intgral: ",np.mean(integralVals))
print("Deviation in pulse integral: ",np.std(integralVals))
print("Average of time series data: ",np.mean(timeSeries))
print("Number of photo-electrons/pulse: ",(np.mean(integralVals)/np.std(integralVals))**2)

#Local Variables
binNumbers2=200
minCts=-50
maxCts=5

distSetup=np.histogram(distArray,bins=binNumbers2,range=(minCts,maxCts))
distXVALS=distSetup[1][0:np.size(distSetup[0])]
distYVALS=distSetup[0]

plt.plot(distXVALS[distYVALS!=0],distYVALS[distYVALS!=0])
#plt.plot(intXVALS,intYVALS,'.')
plt.ylim(-1,10000)
plt.title('Photo-Electron Pulse Height Distribution')
plt.show()

plt.plot(intXVALS,intHVALS,'.')
plt.xlabel('Integral bin value')
plt.ylabel('# of Occurances')
plt.show()
