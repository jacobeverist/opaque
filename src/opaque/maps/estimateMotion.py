from SplineFit import SplineFit
from Pose import Pose
import computeCovar
from math import *
#from copy import copy
from functions import normalizeAngle
import pylab
import numpy
import gen_icp

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
    
    estPlotCount += 1

def plotOffset(originPose, offset, posture1, posture2):
    
    global estPlotCount
  
    poseProfile1 = Pose(originPose)

    xP = []
    yP = []
    for p in posture1:
        nP = poseProfile1.convertLocalToGlobal(p)
        xP.append(nP[0])
        yP.append(nP[1])
        
    pylab.plot(xP,yP, color='b')
    
    pose2 = poseProfile1.convertLocalOffsetToGlobal(offset)
    
    poseProfile2 = Pose(pose2)
    
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
    pylab.xlim(-4,4)
    pylab.ylim(-4,4)

    pylab.savefig("plotCenter%04u.png" % estPlotCount)
    pylab.clf()        

    estPlotCount += 1

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
        
        " compute a fitted curve of a center point "
        centerCurves.append(SplineFit(localPostures[i], smooth = 0.5, kp = 2))
        
        
        f = open(dirName + "/estpose%04u.txt" % i, 'r')
        estPose = eval(f.read().rstrip())
        f.close()
        
        estPoses.append(estPose)
    
        f = open(dirName + "/gndpose%04u.txt" % i, 'r')
        gndPose = eval(f.read().rstrip())
        f.close()
    
        gndPoses.append(gndPose)


    gnds = []
    ests = []
    naives = []

    for i in range(numPoses-1):
        gndOffset, estOffset = computeCovar.getMotionConstraint(gndPoses[i], gndPoses[i+1], estPoses[i], estPoses[i+1])
        gnds.append(gndOffset)
        ests.append(estOffset)
        
        naiveOffset = makeGuess2(centerCurves[i], centerCurves[i+1], 0.2, localPostures[i], localPostures[i+1], originPose=gndPoses[i])
        naives.append(naiveOffset)
        
    #print gnds[0], naives[0]

    #print computeCovar.computeCovar(gnds, ests)
    print computeCovar.computeCovar(gnds, gnds)
    print computeCovar.computeCovar(gnds, naives)

    for i in range(numPoses-1):
        pass
        #print "x:", gnds[i][0], naives[i][0]    
        #print "y:", gnds[i][1], naives[i][1]    
        #print "p:", gnds[i][2], naives[i][2]
        print abs(gnds[i][0] - naives[i][0]), abs(gnds[i][1] - naives[i][1]), abs(gnds[i][2] - naives[i][2])
    
    " plot estimated positions "
    #plotOffset(gndPoses[0], gnds[0], localPostures[0], localPostures[1])
    #plotPoses(gndPoses[0], gndPoses[1], localPostures[0], localPostures[1])
    for i in range(numPoses-1):
        pass
        #plotOffset(gndPoses[i], naives[i], localPostures[i], localPostures[i+1])
        #plotPoses(gndPoses[i], gndPoses[i+1], localPostures[i], localPostures[i+1])

    #offset = [sqrt(0.5),sqrt(0.5),0.0]

    " starting position "
    initPose = [sqrt(0.5),sqrt(0.5),1.0]
    testPose = Pose(initPose)

    " offset and its inverse "
    offset = [sqrt(0.5),sqrt(0.2),pi/2+pi/4.5]
    invOffset = testPose.doInverse(offset)

    " new pose from inverse offset "
    globalPose = testPose.convertLocalOffsetToGlobal(invOffset)
    
    " return to original pose using forward offset "
    testPose2 = Pose(globalPose)
    globalPose2 = testPose2.convertLocalOffsetToGlobal(offset)
    
    " compare that the two are equal "
    #print initPose, globalPose2
    
    samples = []
    for i in range(len(centerCurves)):
        samples.append(centerCurves[i].getUniformSamples())
    
    """   
    pylab.clf()
    for i in range(len(samples)):
        
        xP = []
        yP = []
        for p in samples[i]:
            xP.append(p[0])
            yP.append(p[1])
            
        pylab.plot(xP,yP)
        pylab.scatter(xP,yP)
        
    pylab.show()
    """

    #fixed_offset = gen_icp.motionICP(samples[7], samples[8], naives[7], plotIter = True, n1 = 7, n2 = 8)
    #plotOffset(gndPoses[7], fixed_offset, localPostures[7], localPostures[8])
    
    results = []
    for i in range(len(samples)-1):
        results.append(gen_icp.motionICP(samples[i], samples[i+1], naives[i], plotIter = True, n1 = i, n2 = i+1))
        
    for i in range(len(results)):
        plotOffset(gndPoses[i], results[i], localPostures[i], localPostures[i+1])
        plotOffset(gndPoses[i], gnds[i], localPostures[i], localPostures[i+1])

               
    
    
    