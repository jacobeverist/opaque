import random
from copy import deepcopy
from SplineFit import SplineFit
import pylab
from math import sqrt
from Pose import Pose
import gen_icp
from math import acos, asin

from scipy.spatial import KDTree, cKDTree

"""
def mysum(lvals):
    s2 = 0
    s = 0
    for e in l:
        s += e
        s2 += e * e
    return (s, s2)

N = len(lvals)

s, s2 = mysum(lvals)

var = (s2 - (s*s) / N ) / N
mean = s / N

"""



class ParticleFilter:

    def __init__(self, initPath, nodeHash):
    
        self.numParticles = 20

        " pathID, distMean "
        self.initPose = [0, 0.1]
        
        self.rootPose = [-2.0,0.0]
 
        self.nodeHash = nodeHash
        #self.distMove = 0.628
        #self.distVar = 0.0964       
        self.distMove = 1.0
        self.distVar = 0.5
        
        self.updateCount = 0
        
        self.paths = {}
        self.paths[self.updateCount] = deepcopy(initPath)


        p0 = self.paths[self.updateCount][0]
        pN = self.paths[self.updateCount][-1]
        
        dist0 = sqrt((self.rootPose[0]-p0[0])*(self.rootPose[0]-p0[0]) + (self.rootPose[1]-p0[1])*(self.rootPose[1]-p0[1]))
        distN = sqrt((self.rootPose[0]-pN[0])*(self.rootPose[0]-pN[0]) + (self.rootPose[1]-pN[1])*(self.rootPose[1]-pN[1]))

        if dist0 > distN:
            self.paths[self.updateCount].reverse()
        
        " each particle:  [pathID, distMean, distVar] "
        pathID = self.initPose[0]
        distMean = self.initPose[1]
        
        particles = []
        for i in range(self.numParticles):
            particles.append([pathID, distMean])

        self.particleSnapshots = {}
        self.particleSnapshots[self.updateCount] = particles

        self.drawState()
       
        
    def update(self, newPath0, nodeID, isForward = True):

        " old path + old positions => x,y "
        " x,y => new path + old positions "
        " new path + old positions + move => new path + new positions "
        
        oldPath = self.paths[self.updateCount]
        oldSpline = SplineFit(oldPath, smooth=0.1)
        

        " convert particles to x,y positions "
        xyParticles = []
        partPoints = []
        oldParticles = self.particleSnapshots[self.updateCount]
        for part in oldParticles:
            pathID = part[0]
            distMean = part[1]

            point = oldSpline.getPointOfDist(distMean)
            #print "dist:", distMean, "point:", point

            newPart = [pathID, point, distMean]
            partPoints.append(point[:2])
            xyParticles.append(newPart)

        newPath = deepcopy(newPath0)
        p0 = newPath[0]
        pN = newPath[-1]
        
        dist0 = sqrt((self.rootPose[0]-p0[0])*(self.rootPose[0]-p0[0]) + (self.rootPose[1]-p0[1])*(self.rootPose[1]-p0[1]))
        distN = sqrt((self.rootPose[0]-pN[0])*(self.rootPose[0]-pN[0]) + (self.rootPose[1]-pN[1])*(self.rootPose[1]-pN[1]))

        if dist0 > distN:
            newPath.reverse()
        
        
        " new path "
        newSpline = SplineFit(newPath, smooth = 0.1)
        
        newSpline.precompute()
        splinePoints = newSpline.densePoints
        cartPoints = []
        for p in splinePoints:
            cartPoints.append(p[:2])
        tree = cKDTree(cartPoints)
        
        qResults = tree.query(partPoints)
        distances = qResults[0]
        indices = qResults[1]

        " fit old particles onto new path and compute their distance "
        newParticles = []        
        for k in range(len(xyParticles)):
            
            #distances[k]
            #indices[k]
            
            arcDist = newSpline.distPoints[indices[k]]

            part = xyParticles[k]
            newPart = [part[0], arcDist]

            " move forw0ard or backward "        
            if isForward:
                newPart[1] += self.distMove
            else:
                newPart[1] += -self.distMove

            newPart[1] = random.gauss(newPart[1],self.distVar)

            if newPart[1] <= newSpline.dist_u(1.0) and newPart[1] >= newSpline.dist_u(0.0):  
                newParticles.append(newPart)
            else:
                pass
                " add in random sample "


        " compute the medial axis for each pose "
        node1 = self.nodeHash[nodeID]

        medial1 = node1.getBestMedialAxis()
        medial1 = node1.medialLongPaths[0]
        
        #medial1 = computeHullAxis(nodeID, node1, tailCutOff = False)
        globalPath = newPath
        globalPathReverse = deepcopy(globalPath)
        globalPathReverse.reverse()
        

        
        estPose1 = node1.getGlobalGPACPose()        
        poseOrigin = Pose(estPose1)
        
        globalMedial = []
        for p in medial1:
            globalMedial.append(poseOrigin.convertLocalToGlobal(p))
        
        medialSpline1 = SplineFit(globalMedial, smooth=0.1)
        globalSpline = SplineFit(globalPath, smooth=0.1)
        globalSplineReverse = SplineFit(globalPathReverse, smooth=0.1)
        medialReverse = deepcopy(medial1)
        medialReverse.reverse()

        overlapMatch = []
        angleSum1 = 0.0
        angleSum2 = 0.0
        for i in range(0,len(globalMedial)):
            p_1 = globalMedial[i]
            p_2, j, minDist = gen_icp.findClosestPointInA(globalPath, p_1)
            if minDist < 0.5:
                overlapMatch.append((i,j,minDist))

                pathU1 = medialSpline1.findU(p_1)    
                pathU2 = globalSpline.findU(p_2)    
                pathU2_R = globalSplineReverse.findU(p_2)    

                pathVec1 = medialSpline1.getUVector(pathU1)
                pathVec2 = globalSpline.getUVector(pathU2)
                pathVec2_R = globalSplineReverse.getUVector(pathU2_R)

                val1 = pathVec1[0]*pathVec2[0] + pathVec1[1]*pathVec2[1]
                if val1 > 1.0:
                    val1 = 1.0
                elif val1 < -1.0:
                    val1 = -1.0
                ang1 = acos(val1)
                
                val2 = pathVec1[0]*pathVec2_R[0] + pathVec1[1]*pathVec2_R[1]
                if val2 > 1.0:
                    val2 = 1.0
                elif val2 < -1.0:
                    val2 = -1.0
                ang2 = acos(val2)

                angleSum1 += ang1
                angleSum2 += ang2
        
        " select global path orientation based on which has the smallest angle between tangent vectors "
        print i, "angleSum1 =", angleSum1, "angleSum2 =", angleSum2
        if angleSum1 > angleSum2:
            orientedGlobalPath = globalPathReverse
            isReverse = True
        else:
            orientedGlobalPath = globalPath
            isReverse = False
        

        #globalSplineReverse = SplineFit(globalPathReverse, smooth=0.1)

        orientedSpline = SplineFit(orientedGlobalPath, smooth=0.1)
        orientedSpline.precompute()
        cartPoints1 = []
        for p in orientedSpline.densePoints:
            cartPoints1.append(p[:2])
        treeO = cKDTree(cartPoints1)
        
        
        " sample the entire path for fitting likelihood "
        
        totalDist = orientedSpline.dist_u(1.0)
        
        """
        #totalDist = newSpline.dist_u(1.0)
        distInc = 0.1
        numInc = int(totalDist / distInc) + 1
        #uSet = [newSpline.dist_search(i*distInc)*0.001 for i in range(numInc)]
        uSet = [orientedSpline.dist_search(i*distInc)*0.001 for i in range(numInc)]

        args = []
        for i in range(len(uSet)):
    
            if isReverse:
                u1 = uSet[i]
            else:
                u1 = uSet[i]
                

            args.append([u1, 0.5, 0.0])
            
            # args.append([u1,u2,angGuess])

        print "ICP", len(args), "times"
        print args
        #exit()
        fitResults = gen_icp.batchGlobalICP(orientedGlobalPath, medial1, args)
        #fitResults = gen_icp.batchGlobalICP(newPath, medial1, args)
        print "RESULTS_SAMPLED:", len(fitResults)
 
        measResults1 = []
        queryPoints1 = []
        for k in range(len(fitResults)):
            resultPose, lastCost = fitResults[k]
            
            measResults1.append([resultPose,lastCost])
            queryPoints1.append(resultPose[:2])
        
        qResults1 = treeO.query(queryPoints1)
        distances1 = qResults1[0]
        indices1 = qResults1[1]
            
        measParticles1 = []
        for k in range(len(measResults1)):
            pos = measResults1[k][1]
            cost = measResults1[k][1]
            #arcDist = newSpline.distPoints[indices1[k]]
            arcDist = orientedSpline.distPoints[indices1[k]]
            
            if isReverse:
                arcDist = totalDist - arcDist

            print k*distInc, uSet[k], arcDist, cost

            #if measResults1[k][1] < 500:
            #   measParticles1.append([newPart[0],arcDist])
        """




                
        
        " measurement adjustment with ICP "       
        measResults = []
        queryPoints = []

        args = []
        for k in range(len(newParticles)):
            newPart = newParticles[k]
            
            dist = newPart[1]
            #uVal = newSpline.dist_search(dist) / 1000.

            if isReverse:
                dist = totalDist - dist

            
            uVal = orientedSpline.dist_search(dist) / 1000.
            #print dist, uVal2
            
            u2 = 0.5
            u1 = uVal
            angGuess = 0.0
            
            args.append([u1,u2,angGuess])
        
                
        results = gen_icp.batchGlobalICP(orientedGlobalPath, medial1, args)
        print "len(results) =", len(results), "len(newParticles) =", len(newParticles)
        
        
        for k in range(len(newParticles)):

            #resultPose, lastCost = gen_icp.globalOverlapICP_GPU2(args[k], orientedGlobalPath, medial1)
            #print "estimating branch pose", nodeID, "at",  resultPose[0], resultPose[1], resultPose[2]
            
            resultPose, lastCost = results[k]
            
            measResults.append([resultPose,lastCost])
            queryPoints.append(resultPose[:2])
        

        qResults2 = treeO.query(queryPoints)
        distances2 = qResults2[0]
        indices2 = qResults2[1]
            
        measParticles = []
        for k in range(len(newParticles)):
            newPart = newParticles[k]
            #arcDist = newSpline.distPoints[indices2[k]]
            arcDist = orientedSpline.distPoints[indices2[k]]
            
            if isReverse:
                arcDist = totalDist - arcDist

            print queryPoints[k], arcDist, newPart[1], measResults[k][1]

            if measResults[k][1] < 500:
                measParticles.append([newPart[0],arcDist])
        
        
        " generate random particles with weighted probabilities around existing particles "
        numGenParticles = self.numParticles - len(measParticles)
        numExistingParticles = len(measParticles)
        genParticles = []
        for k in range(numGenParticles):
            
            " randomly select particle out of remaining "
            selectParticleInd = random.randrange(0, numExistingParticles)
            " select its position with a variance "
            #newDist = random.gauss(measParticles[selectParticleInd][1], self.distVar)
            newDist = measParticles[selectParticleInd][1]
            
            " add particle "
            genParticles.append([measParticles[selectParticleInd][0], newDist])
        

        totalParticles = measParticles + genParticles
        
        #newSpline.findClosestPoint(resultPose)
            
            #self.makeGlobalMedialOverlapConstraint(nodeID, newPath, uVal, originU2 = 0.5)
        
        
        
        """
        " fit old particles onto new path and compute their distance "
        newParticles = []        
        for part in xyParticles:
            #print "part:", part
            pos = part[1]
            minDist, minU, closestPoint = newSpline.findClosestPoint(pos)
            arcDist = newSpline.dist(0.0, minU)

            #print "oldDist,newDist,pos,minU,closestPoint:", part[2], arcDist, pos, minU, closestPoint

            newPart = [part[0], arcDist]

            " move forward or backward "        
            if isForward:
                newPart[1] += self.distMove
            else:
                newPart[1] += -self.distMove

            newPart[1] = random.gauss(newPart[1],self.distVar)

            if newPart[1] <= newSpline.dist_u(1.0) and newPart[1] >= newSpline.dist_u(0.0):  
                newParticles.append(newPart)
        """
        
        self.updateCount += 1
        #self.particleSnapshots[self.updateCount] = newParticles
        self.particleSnapshots[self.updateCount] = totalParticles
        self.paths[self.updateCount]  = deepcopy(newPath)    


        self.drawState(nodeID)


    def drawState(self, nodeID = 0):
    
        pylab.clf()
    
        path = self.paths[self.updateCount]
        currSpline = SplineFit(path, smooth=0.1)

        splinePoints = currSpline.getUniformSamples()

        xP = []
        yP = []
        for p in splinePoints:
            xP.append(p[0])
            yP.append(p[1])
        
        pylab.plot(xP,yP, color='r')
        
        particles = self.particleSnapshots[self.updateCount]
        
        
        xP = []
        yP = []
        for p in particles:
            distVal = p[1]
            
            splinePoint = currSpline.getPointOfDist(distVal)
            xP.append(splinePoint[0])
            yP.append(splinePoint[1])
            
        pylab.scatter(xP,yP, color='b')
        
        
        originPoint = currSpline.getPointOfDist(0.0)
        pylab.scatter([originPoint[0]], [originPoint[1]], color='g')
        
        
        node = self.nodeHash[nodeID]
        
        gndPose1 = node.getGndGlobalGPACPose()
        medial1 = node.getBestMedialAxis()
        
        " set the origin of pose 1 "
        poseGnd = Pose(gndPose1)

        xP = []
        yP = []
        for p in medial1:
            p1 = poseGnd.convertLocalToGlobal(p)
            xP.append(p1[0])
            yP.append(p1[1])
        
        pylab.plot(xP,yP, color='k')
            

        for k in range(nodeID+1):
            print self.nodeHash[k].getGndGlobalGPACPose()
        
        
        pylab.xlim(-5,10)
        pylab.ylim(-8,8)
        
        pylab.savefig("particles_%04u.png" % self.updateCount)
    
