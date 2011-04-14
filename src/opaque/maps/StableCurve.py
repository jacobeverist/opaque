
from math import sqrt
import scipy.interpolate
from random import gauss
from copy import copy
from math import floor, asin, acos
from time import time
from numpy import linspace, concatenate

from Pose import Pose
import pca_module
import scipy

"Gross Posture Approximation Curve"

class StableCurve:

    def __init__ (self, points, rotated = True, smooth = 0.5, kp = 2):
        self.pointSet = points
        self.smoothNess = smooth
        self.kp = kp
        
        self.rotated = rotated
        self.rotAngle = 0.0
        self.rotatedPoints = []
        self.rotateProfile = Pose([0.0,0.0,0.0])
        
        newP = []
        
        if self.rotated:
            self.rotatedPoints, self.rotAngle = self.horizontalizePosture(points)
            #print "rotAngle =", self.rotAngle
            for p in self.rotatedPoints:
                newP.append(p)
        else:
            for p in points:
                newP.append(p)
        
        # unzip the points
        intArray = [[],[]]

        " perturb the points to prevent degeneracies "
        for p in newP:
            intArray[0].append(p[0] + gauss(0.0,0.0001))
            intArray[1].append(p[1] + gauss(0.0,0.0001))
            
        foreArray = [[],[]]
        backArray = [[],[]]
        for i in range(20):
            foreArray[0].append(intArray[0][i])
            foreArray[1].append(intArray[1][i])
        for i in range(20,39):
            backArray[0].append(intArray[0][i])
            backArray[1].append(intArray[1][i])

        self.foreSpline = scipy.interpolate.UnivariateSpline(foreArray[0], foreArray[1], s = 2.0)    
        self.backSpline = scipy.interpolate.UnivariateSpline(backArray[0], backArray[1], s = 2.0)
        

    def getKnots(self):

        x1 = self.foreSpline.get_knots()
        y1 = self.foreSpline(x1)
        x2 = self.backSpline.get_knots()
        y2 = self.backSpline(x2)
        
        x = concatenate((x1, x2))
        y = concatenate((y1,y2))

        #print x1, x2
        #print x
        #x = self.spline.get_knots()
        #y = self.spline(x)

        xP = []
        yP = []
        if self.rotated:
            for i in range(len(x)):
                newPoint = self.rotateProfile.convertLocalToGlobal([x[i],y[i]])
                xP.append(newPoint[0])
                yP.append(newPoint[1])

            return xP, yP


        return x, y
    
    def getPlot(self):
        
        #knots = self.spline.get_knots()    
        #x = linspace(knots[0],knots[-1],100)
        #y = self.spline(x)


        knots1 = self.foreSpline.get_knots()
        knots2 = self.backSpline.get_knots()

        x1 = linspace(knots1[0],knots1[-1],100)
        y1 = self.foreSpline(x1)

        x2 = linspace(knots2[0],knots2[-1],100)
        y2 = self.backSpline(x2)


        #x = x1 + x2
        #y = y1 + y2
        x = concatenate((x1, x2))
        y = concatenate((y1,y2))
        xP = []
        yP = []
        if self.rotated:
            for i in range(len(x)):
                newPoint = self.rotateProfile.convertLocalToGlobal([x[i],y[i]])
                xP.append(newPoint[0])
                yP.append(newPoint[1])

            return xP, yP
        
        return x,y

    def horizontalizePosture(self, posture):
    
        x_list = []
        y_list = []
    
        for p in posture:
            x_list.append(p[0])
            y_list.append(p[1])
    
        cov_a = scipy.cov(x_list,y_list)
    
        loadings = []
    
        " NOTE:  seems to create opposing colinear vectors if data is colinear, not orthogonal vectors "
    
        try:
            scores, loadings, E = pca_module.nipals_mat(cov_a, 2, 0.000001, False)
    
        except:
            raise
    
        if len(loadings) < 2:
            raise
    
    

        highVarVec = loadings[0]
        lowVarVec = loadings[1]
        
        
        angle = acos(highVarVec[0])
        if asin(highVarVec[1]) < 0:
            angle = -angle
        
        self.rotateProfile = Pose([0.0,0.0,angle])
        
        rotatedPosture = []
        #print "len(posture) =", len(posture)
        for i in range(len(posture)):
            rotatedPosture.append(self.rotateProfile.convertGlobalPoseToLocal(posture[i]))
            #print "CHECK BOO"    
            
        " enforce condition that 0'th node has greater X than -1'th node so we do not have 180 degree errors "
        if rotatedPosture[0][0] < rotatedPosture[-1][0]:
            #print "ROTATING StableCurve by 180 degrees"

            " rotate by 180 degrees "
            highVarVec = [-highVarVec[0], -highVarVec[1]]
            lowVarVec = [-lowVarVec[0], -lowVarVec[1]]
            
            angle = acos(highVarVec[0])
            if asin(highVarVec[1]) < 0:
                angle = -angle
            
            self.rotateProfile = Pose([0.0,0.0,angle])
            
            rotatedPosture = []
            for i in range(len(posture)):
                rotatedPosture.append(self.rotateProfile.convertGlobalPoseToLocal(posture[i]))
                #print "CHECK FOO"
                
    
        #print "posture:", posture
        #print "rotatedPosture:", rotatedPosture
        "return the first vector returned from PCA because this has the highest variance"
        return rotatedPosture, angle
      
    def getPoses(self):
        
        knots1 = self.foreSpline.get_knots()
        
        xMid1 = (knots1[0] + knots1[1])/2.0
        yMid1 = self.foreSpline(xMid1)

        x1 = xMid1 + 0.3
        y1 = self.foreSpline(x1)
        
        newPoint1 = [xMid1, yMid1]
        newPoint2 = [x1, y1]

        if self.rotated:
            newPoint1 = self.rotateProfile.convertLocalToGlobal(newPoint1)
            newPoint2 = self.rotateProfile.convertLocalToGlobal(newPoint2)

        # tangent vector
        vec = [newPoint2[0] - newPoint1[0], newPoint2[1] - newPoint1[1]]

        # normalize
        mag = sqrt(vec[0]**2 + vec[1]**2)
        vec = [vec[0]/mag, vec[1]/mag]

        angle = acos(vec[0])
        if asin(vec[1]) < 0:
            angle = -angle
        
        pose1 = [newPoint1[0], newPoint1[1], angle]

        knots2 = self.backSpline.get_knots()
                 
        xMid2 = (knots2[0] + knots2[1])/2.0
        yMid2 = self.backSpline(xMid2)

        x2 = xMid2 + 0.3
        y2 = self.backSpline(x2)
        
        newPoint1 = [xMid2, yMid2]
        newPoint2 = [x2, y2]

        if self.rotated:
            newPoint1 = self.rotateProfile.convertLocalToGlobal(newPoint1)
            newPoint2 = self.rotateProfile.convertLocalToGlobal(newPoint2)

        # tangent vector
        vec = [newPoint2[0] - newPoint1[0], newPoint2[1] - newPoint1[1]]

        # normalize
        mag = sqrt(vec[0]**2 + vec[1]**2)
        vec = [vec[0]/mag, vec[1]/mag]

        angle = acos(vec[0])
        if asin(vec[1]) < 0:
            angle = -angle
        
        pose2 = [newPoint1[0], newPoint1[1], angle]
        
    
        return pose1, pose2

if __name__ == '__main__':
    
    dirName = "../../testData/poseTest"
    numPoses = 16
    
    localPostures = []
    centerCurves = []
    
    estPoses = []
    gndPoses = []

    for i in range(numPoses):

        f = open(dirName + "/posture%04u.txt" % i, 'r')
        localPostures.append(eval(f.read().rstrip()))
        f.close()
    print len(localPostures), numPoses
    print time()
    for i in range(numPoses):    
        " compute a fitted curve of a center point "
        centerCurves.append(SplineFit(localPostures[i], smooth = 0.5, kp = 2))
        vec = centerCurves[i].getUVector(0.5)
        pnt = centerCurves[i].getU(0.5)

    print time()
        

