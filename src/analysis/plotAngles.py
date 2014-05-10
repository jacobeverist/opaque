import csv
import pylab
import math


errorReader = csv.reader(open('angles.txt', 'rb'), delimiter=',')

header = errorReader.next()

print header

angles = []
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
        
    angles.append(float(row[2]))
    xyErrors1.append(float(row[3]))
    angErrors1.append(float(row[4]))
    xyErrors2.append(float(row[5]))
    angErrors2.append(float(row[6]))

t = range(0,len(xyErrors1))    


#pylab.plot(t,xyErrors2, color=(0.5,0.5,1.0), linewidth='1')
#pylab.plot(t,xyErrors1, color='r', linewidth='1')
pylab.plot(t,angErrors2, color=(0.5,0.5,1.0), linewidth='1')
pylab.plot(t,angErrors1, color='r', linewidth='1')
pylab.plot(angles, color='g')
pylab.axhline(y=0, xmin=0, xmax=1)

print newPose

for v in newPose:
    pylab.axvline(x = v, color='0.5')
    
pylab.ylim(-math.pi/4., math.pi/4.)  

pylab.show()