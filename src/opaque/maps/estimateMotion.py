from SplineFit import SplineFit
from GPACurve import GPACurve
from StableCurve import StableCurve
from Pose import Pose
import computeCovar
from math import *
#from copy import copy
from functions import normalizeAngle
import pylab
import numpy
import gen_icp
import scipy.optimize
from time import time

from scipy.interpolate import *

import scipy
import pca_module
import random

import os

estPlotCount = 0


def makeGuess2(centerCurve1, centerCurve2, stepDist, posture1, posture2, originPose =[0.0,0.0,0.0]):

    " get the curve root vectors "
    vec1 = centerCurve1.getUVector(0.3)
    vec2 = centerCurve2.getUVector(0.5)
    
    " negate vectors since they're pointing backwards instead of forwards "
    vec1 = [-vec1[0],-vec1[1]]
    vec2 = [-vec2[0],-vec2[1]]
    
    
    vecPoint1 = centerCurve1.getU(0.3)
    vecPoint2 = centerCurve2.getU(0.5)
    
    angle1 = acos(vec1[0])
    if asin(vec1[1])  < 0:
        angle1 = -angle1

    angle2 = acos(vec2[0])
    if asin(vec2[1]) < 0:
        angle2 = -angle2
    
    curve1Pose = [vecPoint1[0], vecPoint1[1], angle1]
    curve2Pose = [vecPoint2[0], vecPoint2[1], angle2]
    
    
    " root nodes, start at [0,0,0] "
    
    " we want to project stepDist in the direction of angle 1 "
    xCurveLocalDes2 = curve1Pose[0] + stepDist * cos(curve1Pose[2])
    yCurveLocalDes2 = curve1Pose[1] + stepDist * sin(curve1Pose[2])
    
    originProfile = Pose(originPose)
    
    #oP = originProfile.convertLocalToGlobal([curve1Pose[0],curve1Pose[1]])
    #xP = [oP[0]]
    #yP = [oP[1]]
    #pylab.scatter(xP,yP, color='k')

    #nP = originProfile.convertLocalToGlobal([xCurveLocalDes2,yCurveLocalDes2])
    #xP = [nP[0]]
    #yP = [nP[1]]
    #pylab.scatter(xP,yP, color='b')

    " convert the desired intersection point on curve 1 into global coordinates "
    poseProfile1 = Pose([0.0,0.0,0.0])
    curveDesOffset2 = [xCurveLocalDes2, yCurveLocalDes2, angle1]
    curveDesGlobal2 = poseProfile1.convertLocalOffsetToGlobal(curveDesOffset2)
    
    " now convert this point into a pose, and perform the inverse transform using curve2Pose "
    desGlobalPose2 = Pose(curveDesGlobal2)
    
    #print "des:", curveDesGlobal2
    
    " perform inverse offset from the destination pose "
    negCurve2Pose = desGlobalPose2.doInverse(curve2Pose)
    #print "curve2Pose:", curve2Pose
    #print "negCurve2Pose:", negCurve2Pose
    
    resPose2 = desGlobalPose2.convertLocalOffsetToGlobal(negCurve2Pose)

    resProfile = Pose(resPose2)
    #print "resProfile:", resProfile.convertLocalOffsetToGlobal(curve2Pose)

    #print "resPose2:", resPose2
    #print "rev:", resProfile.convertLocalOffsetToGlobal(curve2Pose)
    desGlobalPose2Rev = Pose(resProfile.convertLocalOffsetToGlobal(curve2Pose))
    negCurveDesOffset2 = desGlobalPose2Rev.doInverse(curveDesOffset2)
    #print desGlobalPose2Rev.convertLocalOffsetToGlobal(negCurveDesOffset2)

    localOffset = poseProfile1.convertGlobalPoseToLocal(resPose2)
    #print "localOffset:", localOffset

    
    #plotOffset(originPose, localOffset, posture1, posture2)
       
    return [localOffset[0], localOffset[1], localOffset[2]]

    " we want to put xStep, yStep on the point of vecPoint2 "
    " convert xStep, yStep from 1-local coordinates to 2-local coordinates "

    " now compute the offset terms "    
    #xT = cos(0)*(xStep) + sin(0)*(yStep)
    #yT = -sin(0)*(xStep) + cos(0)*(yStep)
    #pT = normalizeAngle(angle2-angle1)

    #return [-xT, -yT, -pT]
    #return [xT, yT, pT]
    
    #return globalCurve2

def horizontalizePosture(posture):

    x_list = []
    y_list = []

    for p in posture:
        x_list.append(p[0])
        y_list.append(p[1])

    cov_a = scipy.cov(x_list,y_list)

    loadings = []

    # NOTE:  seems to create opposing colinear vectors if data is colinear, not orthogonal vectors

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
    
    rotateProfile = Pose([0.0,0.0,angle])
    
    rotatedPosture = []
    for i in range(len(posture)):
        rotatedPosture.append(rotateProfile.convertGlobalPoseToLocal(posture[i]))

    # return the second vector returned from PCA because this has the least variance (orthogonal to plane)
    #return highVarVec
    
    return rotatedPosture, angle

def cost_func(angle, posture1, posture2):


    poseProfile1 = Pose([0.0,0.0,0.0])

    posture1_offset = []
    for p in posture1:
        posture1_offset.append(p)

    pose2 = poseProfile1.convertLocalOffsetToGlobal([0.0,0.0,angle])
    
    poseProfile2 = Pose(pose2)
    
    posture2_offset = []
    for p in posture2:
        nP = poseProfile2.convertLocalOffsetToGlobal(p)
        posture2_offset.append(nP)   

    indices1 = range(0,20)
    indices2 = range(20,39) 

    cost1 = 0.0
    cost2 = 0.0

    for j in indices1:
        p1 = posture1_offset[j]
        p2 = posture2_offset[j]
        
        cost1 += sqrt((p2[0]-p1[0])**2 + (p2[1]-p1[1])**2)

    for j in indices2:
        p1 = posture1_offset[j]
        p2 = posture2_offset[j]
        
        cost2 += sqrt((p2[0]-p1[0])**2 + (p2[1]-p1[1])**2)

    if cost1 > cost2:
        return cost2
    else:
        return cost1



def correctOrientation3(posture1, posture2):

    curve1 = StableCurve(posture1, rotated=True)
    curve2 = StableCurve(posture2, rotated=True)

    forePose1, backPose1 = curve1.getPoses()
    forePose2, backPose2 = curve2.getPoses()

    originPoint = [0.0,0.0]

    #handleLocs = [0.3, 0.7]
    handleLocs = [0.1, 0.2, 0.3, 0.4, 0.6, 0.7, 0.8, 0.9]
    handlePoints1 = curve1.getUSet(handleLocs)
    handlePoints2 = curve2.getUSet(handleLocs)

    hypotheses = []
    
    for i in range(len(handleLocs)):
        vec1 = [handlePoints1[i][0], handlePoints1[i][1]]
        mag = sqrt(vec1[0]*vec1[0] + vec1[1]*vec1[1])
        vec1[0] /= mag
        vec1[1] /= mag
        
        angle1 = acos(vec1[0])
        if asin(vec1[1])  < 0:
            angle1 = -angle1
        
        vec2 = [handlePoints2[i][0], handlePoints2[i][1]]
        mag = sqrt(vec2[0]*vec2[0] + vec2[1]*vec2[1])
        vec2[0] /= mag
        vec2[1] /= mag

        angle2 = acos(vec2[0])
        if asin(vec2[1])  < 0:
            angle2 = -angle2

        hypotheses.append(normalizeAngle(angle1 - angle2))

    poseProfile1 = Pose([0.0,0.0,0.0])

    posture1_offset = []
    for p in posture1:
        nP = poseProfile1.convertLocalOffsetToGlobal(p)
        posture1_offset.append(nP)

    indices = []
    for i in range(len(handleLocs)):
        if handleLocs[i] < 0.5:
            indices.append(range(0,20))
        else:
            indices.append(range(20,39))

    costs = []
    for i in range(len(handleLocs)):
        pose2 = poseProfile1.convertLocalOffsetToGlobal([0.0,0.0,hypotheses[i]])
        
        poseProfile2 = Pose(pose2)
        
        posture2_offset = []
        for p in posture2:
            nP = poseProfile2.convertLocalOffsetToGlobal(p)
            posture2_offset.append(nP)

        " measure the cost of side 1 and side 2 "
        
        cost = 0.0
        for j in indices[i]:
            p1 = posture1_offset[j]
            p2 = posture2_offset[j]
            
            cost += sqrt((p2[0]-p1[0])**2 + (p2[1]-p1[1])**2)

        costs.append(cost)

    minCost = 1e100
    minIndex = 0

    for i in range(len(costs)):
        if costs[i] < minCost:
            minCost = costs[i]
            minIndex = i

    return hypotheses[minIndex], minCost

    if costs[0] > costs[1]:
        return hypotheses[1], costs[1]
    else:
        return hypotheses[0], costs[0]

        
    #return hypotheses[0]
    
    
def correctOrientation2(posture1, posture2):

    angle = scipy.optimize.fmin(cost_func, 0.0, [posture1, posture2], disp = 0)

    return angle, 0.0


def correctOrientation(posture1, posture2):

    #poseProfile1 = Pose(originPose)

    #posture1_offset = []
    #for p in posture1:
    #    nP = poseProfile1.convertLocalToGlobal(p)
    #    posture1_offset.append(nP)
        
    #pose2 = poseProfile1.convertLocalOffsetToGlobal(offset)
    
    #poseProfile2 = Pose(pose2)
    
    #posture2_offset = []
    #for p in posture2:
    #    nP = poseProfile2.convertLocalToGlobal(p)
    #    posture2_offset.append(nP)
    
    #curve1 = SplineFit(posture1, smooth = 0.5, kp = 2)
    #curve2 = SplineFit(posture2, smooth = 0.5, kp = 2)
    curve1 = GPACurve(posture1, rotated=True)
    curve2 = GPACurve(posture2, rotated=True)

    originPoint = [0.0,0.0]

    #handleLocs = [0.3, 0.7]
    handleLocs = [0.1, 0.2, 0.3, 0.4, 0.6, 0.7, 0.8, 0.9]
    handlePoints1 = curve1.getUSet(handleLocs)
    handlePoints2 = curve2.getUSet(handleLocs)

    hypotheses = []
    
    for i in range(len(handleLocs)):
        vec1 = [handlePoints1[i][0], handlePoints1[i][1]]
        mag = sqrt(vec1[0]*vec1[0] + vec1[1]*vec1[1])
        vec1[0] /= mag
        vec1[1] /= mag
        
        angle1 = acos(vec1[0])
        if asin(vec1[1])  < 0:
            angle1 = -angle1
        
        vec2 = [handlePoints2[i][0], handlePoints2[i][1]]
        mag = sqrt(vec2[0]*vec2[0] + vec2[1]*vec2[1])
        vec2[0] /= mag
        vec2[1] /= mag

        angle2 = acos(vec2[0])
        if asin(vec2[1])  < 0:
            angle2 = -angle2

        hypotheses.append(normalizeAngle(angle1 - angle2))

    poseProfile1 = Pose([0.0,0.0,0.0])

    posture1_offset = []
    for p in posture1:
        nP = poseProfile1.convertLocalOffsetToGlobal(p)
        posture1_offset.append(nP)

    indices = []
    for i in range(len(handleLocs)):
        if handleLocs[i] < 0.5:
            indices.append(range(0,20))
        else:
            indices.append(range(20,39))

    costs = []
    for i in range(len(handleLocs)):
        pose2 = poseProfile1.convertLocalOffsetToGlobal([0.0,0.0,hypotheses[i]])
        
        poseProfile2 = Pose(pose2)
        
        posture2_offset = []
        for p in posture2:
            nP = poseProfile2.convertLocalOffsetToGlobal(p)
            posture2_offset.append(nP)

        " measure the cost of side 1 and side 2 "
        
        cost = 0.0
        for j in indices[i]:
            p1 = posture1_offset[j]
            p2 = posture2_offset[j]
            
            cost += sqrt((p2[0]-p1[0])**2 + (p2[1]-p1[1])**2)

        costs.append(cost)

    minCost = 1e100
    minIndex = 0

    for i in range(len(costs)):
        if costs[i] < minCost:
            minCost = costs[i]
            minIndex = i

    return hypotheses[minIndex], minCost

    if costs[0] > costs[1]:
        return hypotheses[1], costs[1]
    else:
        return hypotheses[0], costs[0]

        
    #return hypotheses[0]
    


def plotPosture(posture1, posture2):
    global estPlotCount
 
    xP = []
    yP = []
    for p in posture1:
        xP.append(p[0])
        yP.append(p[1])
        
    pylab.plot(xP,yP, color='b')

    #splPoints1 = curve1.getPlot()
       
    #xP = []
    #yP = []
    #for p in splPoints1:
    #    xP.append(p[0])
    #    yP.append(p[1])
        
    #pylab.plot(xP,yP, color='b')    

    xP = []
    yP = []
    for p in posture2:
        xP.append(p[0])
        yP.append(p[1])
        
    pylab.plot(xP,yP, color='r')

    #splPoints2 = curve2.getPlot()
       
    #xP = []
    #yP = []
    #for p in splPoints2:
    #    xP.append(p[0])
    #    yP.append(p[1])
        
    #pylab.plot(xP,yP, color='r') 
    pylab.xlim(-4,4)
    pylab.ylim(-4,4)
    pylab.savefig("plotPosture%04u.png" % estPlotCount)
    pylab.clf()      
    
    estPlotCount += 1


def plotPoses(pose1, pose2, posture1, posture2):
    global estPlotCount


            
    poseProfile1 = Pose(pose1)
    poseProfile2 = Pose(pose2)
 
    xP = []
    yP = []
    for p in posture1:
        nP = poseProfile1.convertLocalToGlobal(p)
        xP.append(nP[0])
        yP.append(nP[1])
        
    pylab.plot(xP,yP, color='b')

    xP = []
    yP = []
    for p in posture2:
        nP = poseProfile2.convertLocalToGlobal(p)
        xP.append(nP[0])
        yP.append(nP[1])
        
    pylab.plot(xP,yP, color='r') 

    curve1 = SplineFit(posture1, smooth = 0.5, kp = 2)
    curve2 = SplineFit(posture2, smooth = 0.5, kp = 2)

    upoints = numpy.arange(0,1.0,0.01)
    splPoints1 = curve1.getUSet(upoints)
    splPoints2 = curve2.getUSet(upoints)
       
    xP = []
    yP = []
    for p in splPoints1:
        nP = poseProfile1.convertLocalToGlobal(p)
        xP.append(nP[0])
        yP.append(nP[1])
        
    pylab.plot(xP,yP, color='b')    

    xP = []
    yP = []
    for p in splPoints2:
        nP = poseProfile2.convertLocalToGlobal(p)
        xP.append(nP[0])
        yP.append(nP[1])
        
    pylab.plot(xP,yP, color='r')    
    
    pylab.xlim(-4,4)
    pylab.ylim(-4,4)
    pylab.savefig("plotCenter%04u.png" % estPlotCount)
    pylab.clf()      
    
    #estPlotCount += 1


def plotOffsetAndGnd(originPose, offset, newPosture, gndPose, gndPosture):
    
    global estPlotCount
  
    gndProfile = Pose(gndPose)

    gndPosture_offset = []
    for p in gndPosture:
        nP = gndProfile.convertLocalOffsetToGlobal(p)
        gndPosture_offset.append(nP)
        
    xP = []
    yP = []
    for p in gndPosture_offset:
        xP.append(p[0])
        yP.append(p[1])
    pylab.plot(xP,yP, color='b')
    
    
    originProfile = Pose(originPose)    
    correctedPose = originProfile.convertLocalOffsetToGlobal(offset)
    correctedProfile = Pose(correctedPose)
    posture_offset = []
    for p in newPosture:
        nP = correctedProfile.convertLocalOffsetToGlobal(p)
        posture_offset.append(nP)
    
    xP = []
    yP = []
    for p in posture_offset:
        xP.append(p[0])
        yP.append(p[1])
        
    pylab.plot(xP,yP, color='r')    


    distError = sqrt((correctedPose[0]-gndPose[0])**2 + (correctedPose[1]-gndPose[1])**2)
    angleError = normalizeAngle(correctedPose[2]-gndPose[2])
    #print distError, angleError

    pylab.title("XY Error = %3.2f, Ang Error = %3.2f" % (distError, angleError))
    
    pylab.xlim(-4,4)
    pylab.ylim(-4,4)

    pylab.savefig("plotPosture%04u.png" % estPlotCount)
    pylab.clf()        

    estPlotCount += 1


def plotOffset(originPose, offset, posture1, posture2, cost = None):
    
    global estPlotCount
  
    poseProfile1 = Pose(originPose)

    posture1_offset = []
    for p in posture1:
        nP = poseProfile1.convertLocalOffsetToGlobal(p)
        posture1_offset.append(nP)
        
    xP = []
    yP = []
    for p in posture1_offset:
        #nP = poseProfile1.convertLocalToGlobal(p)
        xP.append(p[0])
        yP.append(p[1])
        
    pylab.plot(xP,yP, color='b')
    
    pose2 = poseProfile1.convertLocalOffsetToGlobal(offset)
    
    poseProfile2 = Pose(pose2)
    
    posture2_offset = []
    for p in posture2:
        nP = poseProfile2.convertLocalOffsetToGlobal(p)
        posture2_offset.append(nP)
    
    xP = []
    yP = []
    for p in posture2_offset:
        #nP = poseProfile2.convertLocalToGlobal(p)
        xP.append(p[0])
        yP.append(p[1])
        
    pylab.plot(xP,yP, color='r')    


    #FIXME:  not the same curve from spline fit at different configurations


    #curve1 = SplineFit(posture1_offset, smooth = 0.5, kp = 2)
    #curve2 = SplineFit(posture2_offset, smooth = 0.5, kp = 2)
    #curve1 = SplineFit(posture1_offset, smooth = 0.5, kp = 5)
    #curve2 = SplineFit(posture2_offset, smooth = 0.5, kp = 5)
    curve1 = GPACurve(posture1_offset, rotated=True)
    curve2 = GPACurve(posture2_offset, rotated=True)

    upoints = numpy.arange(0,1.0,0.01)
    splPoints1 = curve1.getUSet(upoints)
    splPoints2 = curve2.getUSet(upoints)

    handleLocs = [0.3,0.7]
    handlePoints1 = curve1.getUSet(handleLocs)
    handlePoints2 = curve2.getUSet(handleLocs)

    xP = []
    yP = []
    for p in splPoints1:
        #nP = poseProfile1.convertLocalToGlobal(p)
        xP.append(p[0])
        yP.append(p[1])
        
    pylab.plot(xP,yP, color='k')    

    xP = []
    yP = []
    for p in splPoints2:
        #nP = poseProfile2.convertLocalToGlobal(p)
        xP.append(p[0])
        yP.append(p[1])
        
    pylab.plot(xP,yP, color='k')    
    
    pylab.scatter([originPose[0]],[originPose[1]], color = 'k')

    xP = []
    yP = []
    for i in range(2):
        xP.append(handlePoints1[i][0])
        yP.append(handlePoints1[i][1])

        xP.append(handlePoints2[i][0])
        yP.append(handlePoints2[i][1])

    pylab.scatter(xP,yP, color = 'k')

    if cost:
        pylab.title("Cost = %3.2f" % cost)

    pylab.xlim(-4,4)
    pylab.ylim(-4,4)

    pylab.savefig("plotPosture%04u.png" % estPlotCount)
    pylab.clf()        

    estPlotCount += 1


def splineTest(posture):

    global estPlotCount

    " perturb the points to prevent degeneracies "
    #xP = []
    #yP = []
    #for p in posture:
    #    xP.append(p[0] + random.gauss(0.0,0.0001))
    #    yP.append(p[1] + random.gauss(0.0,0.0001))
    
    #spline1 = UnivariateSpline(xP, yP, s = 2.0)
    #knots = spline1.get_knots()
    #xMid = (knots[0] + knots[1])/2.0
    #spline3 = LSQUnivariateSpline(xP, yP, [0.0, 0.2, 0.4, 0.6, 0.8, 1.0])
    #spline3 = LSQUnivariateSpline(xP, yP, [xMid-0.05, xMid+0.05])


    #print spline1.get_knots()

    xP = []
    yP = []
    for p in posture:
        xP.append(p[0])
        yP.append(p[1])
    pylab.plot(xP,yP, color='b')

    filledPosture = []
    for i in range(len(posture)-1):
        posture[i]
        posture[i+1]

    curve1 = GPACurve(posture, rotated = True)

    x, y = curve1.getPlot()
        
    #x = numpy.linspace(-3,3,100)
    #y = spline1(x)    
    pylab.plot(x,y, color='r')
    
    #xP = curve1.spline.get_knots()
    #yP = curve1.spline(xP)
    
    xP, yP = curve1.getKnots()
    
    pylab.scatter(xP,yP, color='k')

    #x = numpy.linspace(-3,3,100)
    #y = spline2(x)    
    #pylab.plot(x,y, color='g')
    
    pylab.xlim(-4,4)
    pylab.ylim(-4,4)
    
    pylab.savefig("plotPosture%04u.png" % estPlotCount)
    pylab.clf()
    
    estPlotCount += 1
    
    
    #print spline1.get_knots()

if __name__ == '__main__':

    dirName = "../../testData/fixedPoseTest3/"
    #dirName = "../../"
    #numPoses = 12
    numPoses = 35
    #numPoses = 1
    
    poseData = []

    for i in range(numPoses):
        
        files = os.listdir(dirName)
        estNames = []
        gndNames = []
        postureNames = []
        
        name1 = "estpose_%03u" % i
        name2 = "gndpose_%03u" % i
        name3 = "posture_snap_%03u" % i
        
        for nameF in files:
            if  name1 == nameF[0:11]:
                estNames.append(nameF)

            if  name2 == nameF[0:11]:
                gndNames.append(nameF)

            if  name3 == nameF[0:16]:
                postureNames.append(nameF)

        estNames.sort()
        gndNames.sort()
        postureNames.sort()
        
        #print len(estNames), len(gndNames), len(postureNames)

        localPostures= []
        localCurves = []        
        estPoses = []
        gndPoses = []

        for j in range(len(postureNames)):
                
            f = open(dirName + postureNames[j], 'r')
            localPostures.append(eval(f.read().rstrip()))
            f.close()
            
            " compute a fitted curve of a center point "
            #localCurves.append(GPACurve(localPostures[j], rotated=True))

            #splineTest(localPostures[j])
    
            f = open(dirName + estNames[j], 'r')
            estPoses.append(eval(f.read().rstrip()))
            f.close()
    
            f = open(dirName + gndNames[j], 'r')
            gndPoses.append(eval(f.read().rstrip()))
            f.close()

        poseData.append([estPoses, gndPoses, localPostures, localCurves])

    #for poseIndex in range(numPoses):
    #    for j in range(len(poseData[poseIndex][2])-1):
    #        posture = poseData[poseIndex][2][j]
            
    #        splineTest(posture)

    #exit()

    for poseIndex in range(numPoses):


        originPosture = poseData[poseIndex][2][0]
        originCurve = StableCurve(originPosture)
        originForePose, originBackPose = originCurve.getPoses()
         
        originForeProfile = Pose(originForePose)
        originBackProfile = Pose(originBackPose)
        
        originForePosture = []
        originBackPosture = []
        for i in range(len(originPosture)):
            originForePosture.append(originForeProfile.convertGlobalPoseToLocal(originPosture[i]))
        for i in range(len(originPosture)):
            originBackPosture.append(originBackProfile.convertGlobalPoseToLocal(originPosture[i]))

        for j in range(len(poseData[poseIndex][2])-1):


            #plotOffset([0.0,0.0,0.0], result, localPostures[poseIndex], localPostures[poseIndex+1])
 
            localPosture = poseData[poseIndex][2][j+1]
            localCurve = StableCurve(localPosture)

            localForePose, localBackPose = localCurve.getPoses()
             
            localForeProfile = Pose(localForePose)
            localBackProfile = Pose(localBackPose)
            
            localForePosture = []
            localBackPosture = []
            for i in range(len(localPosture)):
                localForePosture.append(localForeProfile.convertGlobalPoseToLocal(localPosture[i]))
            for i in range(len(localPosture)):
                localBackPosture.append(localBackProfile.convertGlobalPoseToLocal(localPosture[i]))
        
            foreCost = 0.0
            for j in range(0,20):
                p1 = originForePosture[j]
                p2 = localForePosture[j]
                
                foreCost += sqrt((p2[0]-p1[0])**2 + (p2[1]-p1[1])**2)           
            
            backCost = 0.0
            for j in range(20,39):
                p1 = originBackPosture[j]
                p2 = localBackPosture[j]
                
                backCost += sqrt((p2[0]-p1[0])**2 + (p2[1]-p1[1])**2)           

            if foreCost > backCost:
                plotPosture(originBackPosture, localBackPosture)
                
            else:
                plotPosture(originForePosture, localForePosture)
            
            
            #t1 = time()
            #angle, cost = correctOrientation3(originAIRPosture, airPosture2)
            #t2 = time()                    
            
            #gndPoseCurr = poseData[poseIndex][1][j+1]
            #globalGndProfile = Pose(gndPoseCurr)
            #globalGndAIRPose = globalGndProfile.convertLocalOffsetToGlobal(localAIRPose2)
            
            #tempProfile = Pose(globalOriginAIRPose)    
            #correctedPose = tempProfile.convertLocalOffsetToGlobal([0.0,0.0,angle])
            #uncorrectedPose = tempProfile.convertLocalOffsetToGlobal([0.0,0.0,0.0])
               
            #distError = sqrt((correctedPose[0]-globalGndAIRPose[0])**2 + (correctedPose[1]-globalGndAIRPose[1])**2)
            #angleError = normalizeAngle(correctedPose[2]-globalGndAIRPose[2])
            #distError2 = sqrt((uncorrectedPose[0]-globalGndAIRPose[0])**2 + (uncorrectedPose[1]-globalGndAIRPose[1])**2)
            #angleError2 = normalizeAngle(uncorrectedPose[2]-globalGndAIRPose[2])
            #print "%d, %d, %1.5f, %1.5f, %1.5f, %1.5f, %1.5f" % (poseIndex, j+1, angle, distError, angleError, distError2, angleError2)

    exit()            
 
    
    #for j in range(numPoses):
    #    posture_example = poseData[j][2][0]
    #    curve_example = poseData[j][3][0]

    #    rotatedPosture, rotAngle = horizontalizePosture(posture_example)
        
    #    regularCurve = GPACurve(posture_example, rotated=False)
    #    unrotatedCurve = GPACurve(posture_example, rotated=True)
    #    rotatedCurve = GPACurve(rotatedPosture, rotated=True)
        
    #    plotPosture(posture_example, curve_example)
    #    plotPosture(posture_example, regularCurve)
    #    plotPosture(posture_example, unrotatedCurve)
    #    plotPosture(rotatedPosture, rotatedCurve)
    
    #    print regularCurve.rotAngle, unrotatedCurve.rotAngle, rotatedCurve.rotAngle
    
    #print poseData[0]
    #exit()

    naives = []
    for i in range(numPoses):
        naives.append([0.0,0.0,0.0])

    " plot estimated positions "
    for i in range(numPoses-1):
        pass
        #plotOffset(gndPoses[i], naives[i], localPostures[i], localPostures[i+1])
        #plotPoses(gndPoses[i], gndPoses[i+1], localPostures[i], localPostures[i+1])

    #samples = []
    #for i in range(len(centerCurves)):
    #    samples.append(centerCurves[i].getUniformSamples())

    print "Pose No, Snapshot No, Correction Angle, Corrected XY Error, Corrected Ang Error, Uncorrected XY Error, Uncorrected Ang Error"
    
    for poseIndex in range(len(naives)):
        #result = gen_icp.motionICP(samples[poseIndex], samples[poseIndex+1], naives[poseIndex], plotIter = True, n1 = poseIndex, n2 = poseIndex+1)
        
        originCurve = poseData[poseIndex][3][0]
        originAIRPose = originCurve.getPose()
         
        originAIRProfile = Pose(originAIRPose)
 
        originPosture = poseData[poseIndex][2][0]
        
        originAIRPosture = []
        for i in range(len(originPosture)):
            originAIRPosture.append(originAIRProfile.convertGlobalPoseToLocal(originPosture[i]))

        gndPose = poseData[poseIndex][1][0]
        globalOriginProfile = Pose(gndPose)
        globalOriginAIRPose = globalOriginProfile.convertLocalOffsetToGlobal(originAIRPose)

        #plotOffset([0.0,0.0,0.0], result, localPostures[poseIndex], localPostures[poseIndex+1])
        for j in range(len(poseData[poseIndex][2])-1):

            curve1 = poseData[poseIndex][3][j]
            curve2 = poseData[poseIndex][3][j+1]

            localAIRPose1 = curve1.getPose() 
            localAIRPose2 = curve2.getPose()
    
            localAIRProfile1 = Pose(localAIRPose1)
            localAIRProfile2 = Pose(localAIRPose2)
            
            posture1 = poseData[poseIndex][2][j]
            posture2 = poseData[poseIndex][2][j+1]
            
            
            airPosture1 = []
            airPosture2 = []
            for i in range(len(posture1)):
                airPosture1.append(localAIRProfile1.convertGlobalPoseToLocal(posture1[i]))
                airPosture2.append(localAIRProfile2.convertGlobalPoseToLocal(posture2[i]))
            
            #plotOffset([0.0,0.0,0.0], [0.0,0.0,0.0], airPosture1, airPosture2)
            t1 = time()
            angle, cost = correctOrientation(originAIRPosture, airPosture2)
            t2 = time()        
            #print poseIndex, "cost:", cost, globalOriginAIRPose
            #plotOffset([0.0,0.0,0.0], [0.0,0.0,0.0], originAIRPosture, airPosture2)
            #if cost > 1.0:
            #    plotOffset([0.0,0.0,0.0], [0.0,0.0,angle], originAIRPosture, airPosture2, cost)
            #plotOffset(globalOriginAIRPose, [0.0,0.0,angle], originAIRPosture, airPosture2, cost)
            #plotOffset([0.0,0.0,0.0], [0.0,0.0,angle], originAIRPosture, airPosture2, cost)
            #plotOffset(globalOriginAIRPose, [0.0,0.0,0.0], originAIRPosture, airPosture2, localAIRPose2[2])
            
            #print angle
            
            
            gndPoseCurr = poseData[poseIndex][1][j+1]
            globalGndProfile = Pose(gndPoseCurr)
            globalGndAIRPose = globalGndProfile.convertLocalOffsetToGlobal(localAIRPose2)
            
            #plotOffsetAndGnd(globalOriginAIRPose, [0.0,0.0,angle], airPosture2, globalGndAIRPose, airPosture2)

               
            tempProfile = Pose(globalOriginAIRPose)    
            correctedPose = tempProfile.convertLocalOffsetToGlobal([0.0,0.0,angle])
            uncorrectedPose = tempProfile.convertLocalOffsetToGlobal([0.0,0.0,0.0])
               
            distError = sqrt((correctedPose[0]-globalGndAIRPose[0])**2 + (correctedPose[1]-globalGndAIRPose[1])**2)
            angleError = normalizeAngle(correctedPose[2]-globalGndAIRPose[2])
            distError2 = sqrt((uncorrectedPose[0]-globalGndAIRPose[0])**2 + (uncorrectedPose[1]-globalGndAIRPose[1])**2)
            angleError2 = normalizeAngle(uncorrectedPose[2]-globalGndAIRPose[2])
            #print "Pose %d, Snapshot %d, XY Error = %3.2f, Ang Error = %3.2f" % (poseIndex, j+1, distError, angleError)
            print "%d, %d, %1.5f, %1.5f, %1.5f, %1.5f, %1.5f" % (poseIndex, j+1, angle, distError, angleError, distError2, angleError2)
            

    
    
