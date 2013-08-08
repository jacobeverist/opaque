
from math import sqrt
import scipy.interpolate
import random
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
        
        random.seed(0)        

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
            try:
                intArray[0].append(p[0] + gauss(0.0,0.0001))
            except:
                intArray[0].append(p[0])

            try:
                intArray[1].append(p[1] + gauss(0.0,0.0001))
            except:
                intArray[1].append(p[1])
            
            
        foreArray = [[],[]]
        backArray = [[],[]]
        for i in range(20):
            foreArray[0].append(intArray[0][i])
            foreArray[1].append(intArray[1][i])
        for i in range(20,39):
            backArray[0].append(intArray[0][i])
            backArray[1].append(intArray[1][i])

        #print "foreArray:", foreArray
        #print "backArray:", backArray

        self.foreSpline = scipy.interpolate.UnivariateSpline(foreArray[0], foreArray[1], s = 2.5)    
        #print "foreSpline =", self.foreSpline
        self.backSpline = scipy.interpolate.UnivariateSpline(backArray[0], backArray[1], s = 2.5)
        
        #print "backSpline =", self.backSpline

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
    
    
    def getUniformSamples(self):
        pass
    
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


        newPoints1 = []
        newPoints2 = []
        if self.rotated:
            for i in range(len(x1)):
                newPoints1.append(self.rotateProfile.convertLocalToGlobal([x1[i],y1[i]]))
            for i in range(len(x2)):
                newPoints2.append(self.rotateProfile.convertLocalToGlobal([x2[i],y2[i]]))
        else:
            for i in range(len(x1)):
                newPoints1.append([x1[i],y1[i]])
            for i in range(len(x2)):
                newPoints2.append([x2[i],y2[i]])

        # tangent vector
        vec1 = [newPoints1[1][0] - newPoints1[2][0], newPoints1[1][1] - newPoints1[2][1]]
        vec2 = [newPoints2[-2][0] - newPoints2[-3][0], newPoints2[-2][1] - newPoints2[-3][1]]

        #print "vec1:", vec1
        #print "vec2:", vec2

        # normalize
        mag1 = sqrt(vec1[0]**2 + vec1[1]**2)
        vec1 = [vec1[0]/mag1, vec1[1]/mag1]
        mag2 = sqrt(vec2[0]**2 + vec2[1]**2)
        vec2 = [vec2[0]/mag2, vec2[1]/mag2]

        #print "mag1:", mag1
        #print "mag2:", mag2
        #print "vec1:", vec1
        #print "vec2:", vec2

        P1 = [newPoints1[0][0], newPoints1[0][1]]
        P2 = [newPoints2[-1][0], newPoints2[-1][1]]
        #print "P1:", P1
        #print "P2:", P2

        newP1 = [P1[0] + 3*vec1[0], P1[1] + 3*vec1[1]]
        newP2 = [P2[0] + 3*vec2[0], P2[1] + 3*vec2[1]]

        #print "newP1:", newP1
        #print "newP2:", newP2

        #xExtr1_1 = knots1[0] + 10.0
        #xExtr1_2 = knots1[-1] - 10.0

        #xExtr2_1 = knots2[0] + 10.0
        #xExtr2_2 = knots2[-1] - 10.0
        
        #print xExtr1_1, self.foreSpline(xExtr1_1)
        #print xExtr2_2, self.backSpline(xExtr2_2)
 
        newPoints1.insert(0,newP1)
        newPoints2.append(newP2)
        
        return newPoints1 + newPoints2

        #x = concatenate((x1, x2))
        #y = concatenate((y1,y2))
        #xP = []
        #yP = []
        #if self.rotated:
        #    for i in range(len(x)):
        #        newPoint = self.rotateProfile.convertLocalToGlobal([x[i],y[i]])
        #        xP.append(newPoint[0])
        #        yP.append(newPoint[1])

        #    return xP, yP
        
        #return x, y

    
    def getPoints(self):
        x1 = linspace(self.rotatedPoints[0][0],self.rotatedPoints[9][0],20)
        y1 = self.foreSpline(x1)

        x2 = linspace(self.rotatedPoints[30][0],self.rotatedPoints[38][0],20)
        y2 = self.backSpline(x2)

        #print "x1:", x1
        #print "y1:", y1

        newPoints1 = []
        newPoints2 = []
        if self.rotated:
            for i in range(len(x1)):
                newPoints1.append(self.rotateProfile.convertLocalToGlobal([x1[i],y1[i]]))
            for i in range(len(x2)):
                newPoints2.append(self.rotateProfile.convertLocalToGlobal([x2[i],y2[i]]))
        else:
            for i in range(len(x1)):
                newPoints1.append([x1[i],y1[i]])
            for i in range(len(x2)):
                newPoints2.append([x2[i],y2[i]])

        return newPoints1 + newPoints2
      
    def getTips(self):

        
        x1 = self.rotatedPoints[0][0]
        y1 = self.foreSpline(x1)
        
        newPoint1 = [x1, y1]

        if self.rotated:
            newPoint1 = self.rotateProfile.convertLocalToGlobal(newPoint1)

        x2 = self.rotatedPoints[-1][0]
        y2 = self.backSpline(x2)
        
        newPoint2 = [x2, y2]

        if self.rotated:
            newPoint2 = self.rotateProfile.convertLocalToGlobal(newPoint2)
            
        return newPoint1,  newPoint2
    
      
    def getPoses(self):
        
        knots1 = self.foreSpline.get_knots()
        #print "knots1:", knots1
        
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
        #print "knots2:", knots2
                 
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
        

