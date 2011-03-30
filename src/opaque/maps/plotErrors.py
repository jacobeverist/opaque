import csv
import pylab
import math


errorFile = open("errors.csv", 'r')
errorReader = csv.reader(open('errors.csv', 'rb'), delimiter=',')

header = errorReader.next()

print header

xyErrors1 = []
angErrors1 = []
xyErrors2 = []
angErrors2 = []
newPose= []
currPose = None
t = []

for row in errorReader:

    if currPose == None:
        t = [0]
    else:
        t.append(t[-1]+1)

    if currPose != int(row[0]):
        currPose = int(row[0])
        newPose.append(t[-1])
        
    xyErrors1.append(float(row[2]))
    angErrors1.append(float(row[3]))
    xyErrors2.append(float(row[4]))
    angErrors2.append(float(row[5]))

t = range(0,len(xyErrors1))    


#pylab.plot(t,xyErrors2, color=(0.5,0.5,1.0), linewidth='1')
#pylab.plot(t,xyErrors1, color='r', linewidth='1')
pylab.plot(t,angErrors2, color=(0.5,0.5,1.0), linewidth='1')
pylab.plot(t,angErrors1, color='r', linewidth='1')
pylab.axhline(y=0, xmin=0, xmax=1)

print newPose

for v in newPose:
    pylab.axvline(x = v, color='0.5')
    
pylab.ylim(-math.pi/4., math.pi/4.)  

pylab.show()