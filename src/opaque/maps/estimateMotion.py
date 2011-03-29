from SplineFit import SplineFit
from Pose import Pose
import computeCovar
from math import *
#from copy import copy
from functions import normalizeAngle
import pylab
import numpy
import gen_icp
from time import time

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
    curve1 = SplineFit(posture1, smooth = 0.5, kp = 5)
    curve2 = SplineFit(posture2, smooth = 0.5, kp = 5)

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
        nP = poseProfile1.convertLocalToGlobal(p)
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
            nP = poseProfile2.convertLocalToGlobal(p)
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

def plotOffset(originPose, offset, posture1, posture2, cost = None):
    
    global estPlotCount
  
    poseProfile1 = Pose(originPose)

    posture1_offset = []
    for p in posture1:
        nP = poseProfile1.convertLocalToGlobal(p)
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
        nP = poseProfile2.convertLocalToGlobal(p)
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
    curve1 = SplineFit(posture1_offset, smooth = 0.5, kp = 5)
    curve2 = SplineFit(posture2_offset, smooth = 0.5, kp = 5)

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
    
    pylab.scatter([0.0],[0.0], color = 'k')

    xP = []
    yP = []
    for i in range(2):
        xP.append(handlePoints1[i][0])
        yP.append(handlePoints1[i][1])

        xP.append(handlePoints2[i][0])
        yP.append(handlePoints2[i][1])

    pylab.scatter(xP,yP, color = 'k')

    #for i in range(0,20):

    #    p1 = posture1_offset[i]
    #    p2 = posture2_offset[i]

    #    xP = [p1[0],p2[0]]
    #    yP = [p1[1],p2[1]]

    #    pylab.plot(xP,yP,color='k')

        
    """
    xP = []
    yP = []
    for p in splPoints:
        #nP = p
        nP = poseProfile1.convertLocalToGlobal(p)

        xP.append(nP[0])
        yP.append(nP[1])
    pylab.plot(xP,yP, color='b')


    localPosture = currNode2.localPosture
    splPoints = currNode2.splPoints

    estPose2 = poseProfile1.convertLocalOffsetToGlobal(gndOffset) 
    poseProfile2 = Pose(estPose2)
    poseProfile2 = Pose(pose2)
    #print pose2, estPose2


    
    xP = []
    yP = []
    for p in splPoints:
        #nP = gen_icp.dispOffset(p,firstGuess)
        #nP = gen_icp.dispOffset(p,gndOffset)
        nP = poseProfile2.convertLocalToGlobal(p)
        xP.append(nP[0])
        yP.append(nP[1])
    pylab.plot(xP,yP, color='r')

    """
    if cost:
        pylab.title("Cost = %3.2f" % cost)

    pylab.xlim(-4,4)
    pylab.ylim(-4,4)

    pylab.savefig("plotPosture%04u.png" % estPlotCount)
    pylab.clf()        

    estPlotCount += 1

if __name__ == '__main__':

    dirName = "../../testData/fixedPoseTest2/"
    numPoses = 8
    
    poseData = []
 
    for i in range(numPoses):
        
        files = os.listdir(dirName)
        estNames = []
        gndNames = []
        postureNames = []
        
        name1 = "estpose_%02u" % i
        name2 = "gndpose_%02u" % i
        name3 = "posture_snap_%02u" % i
        
        for nameF in files:
            if  name1 == nameF[0:10]:
                estNames.append(nameF)

            if  name2 == nameF[0:10]:
                gndNames.append(nameF)

            if  name3 == nameF[0:15]:
                postureNames.append(nameF)

        estNames.sort()
        gndNames.sort()
        postureNames.sort()
        
        print len(estNames), len(gndNames), len(postureNames)

        localPostures= []
        localCurves = []        

        for j in range(len(postureNames)):
                
            f = open(dirName + postureNames[j], 'r')
            localPostures.append(eval(f.read().rstrip()))
            f.close()
            
            " compute a fitted curve of a center point "
            localCurves.append(SplineFit(localPostures[j], smooth = 0.5, kp = 2))

        f = open(dirName + estNames[0], 'r')
        estpose = eval(f.read().rstrip())
        f.close()

        f = open(dirName + gndNames[0], 'r')
        gndpose = eval(f.read().rstrip())
        f.close()

        poseData.append([estpose, gndpose, localPostures, localCurves])

    #print poseData[0]
    #exit()

    naives = []
    for i in range(numPoses-1):
        naives.append([0.0,0.0,0.0])

    " plot estimated positions "
    for i in range(numPoses-1):
        pass
        #plotOffset(gndPoses[i], naives[i], localPostures[i], localPostures[i+1])
        #plotPoses(gndPoses[i], gndPoses[i+1], localPostures[i], localPostures[i+1])

    #samples = []
    #for i in range(len(centerCurves)):
    #    samples.append(centerCurves[i].getUniformSamples())
    
    for poseIndex in range(len(naives)):
        #result = gen_icp.motionICP(samples[poseIndex], samples[poseIndex+1], naives[poseIndex], plotIter = True, n1 = poseIndex, n2 = poseIndex+1)
        
        originCurve = poseData[poseIndex][3][0]
        vec = originCurve.getUVector(0.5)
        vecPoint = originCurve.getU(0.5)
        
        angle = acos(vec[0])
        if asin(vec[1]) < 0:
            angle = -angle
        
        originAIRPose = [vecPoint[0], vecPoint[1], angle]   
        originAIRProfile = Pose(originAIRPose)
 
        originPosture = poseData[poseIndex][2][0]
        
        originAIRPosture = []
        for i in range(len(originPosture)):
            originAIRPosture.append(originAIRProfile.convertGlobalPoseToLocal(originPosture[i]))


        #plotOffset([0.0,0.0,0.0], result, localPostures[poseIndex], localPostures[poseIndex+1])
        for j in range(len(poseData[poseIndex][2])-1):

            curve1 = poseData[poseIndex][3][j]
            curve2 = poseData[poseIndex][3][j+1]

            vec = curve1.getUVector(0.5)
            vecPoint = curve1.getU(0.5)
            
            angle = acos(vec[0])
            if asin(vec[1]) < 0:
                angle = -angle
            
            localAIRPose1 = [vecPoint[0], vecPoint[1], angle]    
        
            vec = curve2.getUVector(0.5)
            vecPoint = curve2.getU(0.5)
            
            angle = acos(vec[0])
            if asin(vec[1]) < 0:
                angle = -angle
        
            localAIRPose2 = [vecPoint[0], vecPoint[1], angle]    
    
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
            print "cost:", cost
            #plotOffset([0.0,0.0,0.0], [0.0,0.0,0.0], originAIRPosture, airPosture2)
            #if cost > 1.0:
            #    plotOffset([0.0,0.0,0.0], [0.0,0.0,angle], originAIRPosture, airPosture2, cost)
            plotOffset([0.0,0.0,0.0], [0.0,0.0,angle], originAIRPosture, airPosture2, cost)
               
    
    
    
