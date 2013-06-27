
from Pose import Pose
import math
import sys
from copy import copy
from functions import *
import random
from subprocess import Popen, PIPE
import Image
from medialaxis import computeMedialAxis
import graph
from LocalNode import getLongestPath
import gen_icp
from SplineFit import SplineFit
import pylab
import numpy
from operator import itemgetter
import hashlib
from Splices import batchGlobalMultiFit, getMultiDeparturePoint, orientPath

import alphamod

" 1) probability that a path is the same as another path "
" 2) new ID for each new path "
" 3) test fitness of being its own path and being the same as other path "
" 4) collapse path to old path if pass threshold "

"""
Path Data Structure

parentID: int, the parent path ID
branchNodeID: int, the node from which the branching point is computed
localJunctionPose: [float*3], the pose of the branching point in local coordinates of branchNodeID
sameProb: dict of floats for each pathID indicating probability this path is the same as older path 

"""

" Likeness Test "
" 1) % overlap of two paths via ICP "
" 2) successful overlap of member nodes translated to old path "
" 2a) maintain separate set of intra-path constraints for each path "
" 2b) measure X_theta error for each path hypothesis "

" New Path Condition "
" 1) any new departure automatically creates a new path "
" 2) a new departure that does not overlap correctly creates a new path "


" Overlap State "
" 1) perform best overlap of each splice "
" 2) perform ICP of each splice "
" 3) maintain path-state machine to keep track of overlaps and directions and likely configurations "
" 4) extract path-state from particle filter motion estimation and subsequent ICP fit on splice "


" FIXME:  Branch direction does not seem to be used properly yet.  Only used in getOverlapDeparture() "

"""

Things that are happening with the data set with branches 
-----------

Second branch detection fails because is deemed "similar" in departure angle to path 1.
The difference is about 1.01 radians but is just under the difference threshold.
Paradoxically, this node is not added to the similar path and is instead added to path 0.
This is because the similar path 1 is not part of the orderedOverlap set.

When the next instance of branch detection occurs, it is determined to be successful because of both dist and ang diff.
However, this new branch eventually is constrained to the junction where it actually belongs and on top of the node 
with the previous branch detection fail.  The mechanism of how this happens is?   NOTE:  This happened because of a bug in my code
that allowed cross-path constraints between member nodes (node path membership always returned True).  Remarkably, there were no map errors from this permissiveness.

In the event that the first node is successfully branch detected, subsequent branch nodes are added to existing path,
but not properly merged or constrained.

"""

" FIXME: branch from set of ordered overlap pathIDs (0,2) selects 2 as the parent, but should really be 0 "
" How did it overlap with 2 in the first place?  It departed from 0, but selected 2 as parent since it's the terminating pathID "
" 2 was marked positive by overlapCondition since it was parallel to 2 within the minMatchDist of 1.0 "

def computeBareHull(node1, sweep = False, static = False):
    
    if static:
        node1.computeStaticAlphaBoundary()

        a_data = node1.getAlphaBoundary(static=True)
        m = hashlib.md5()
        m.update(repr(a_data))
        #print "Paths: computeBareHull():", int(m.digest().encode('hex'),16)        
        
        
        a_data = decimatePoints(a_data)

        " convert hull points to GPAC coordinates before adding covariances "
        localGPACPose = node1.getLocalGPACPose()
        localGPACProfile = Pose(localGPACPose)
        
        a_data_GPAC = []
        for pnt in a_data:
            a_data_GPAC.append(localGPACProfile.convertGlobalToLocal(pnt))
        
    else:
                
        " Read in data of Alpha-Shapes without their associated covariances "
        node1.computeAlphaBoundary(sweep = sweep)
        a_data = node1.getAlphaBoundary(sweep = sweep)
        m = hashlib.md5()
        m.update(repr(a_data))
        #print "Paths: computeBareHull():", int(m.digest().encode('hex'),16)     
        
        a_data = decimatePoints(a_data)
        
        " convert hull points to GPAC coordinates "
        localGPACPose = node1.getLocalGPACPose()
        localGPACProfile = Pose(localGPACPose)
        
        a_data_GPAC = []
        for pnt in a_data:
            a_data_GPAC.append(localGPACProfile.convertGlobalToLocal(pnt))
    
    return a_data_GPAC


def computeHullAxis(nodeID, node2, tailCutOff = False):

    medial2 = node2.getBestMedialAxis()

    if tailCutOff:
        medial2 = node2.medialTailCuts[0]
    else:
        medial2 = node2.medialLongPaths[0]


    if node2.isBowtie:            
        hull2 = computeBareHull(node2, sweep = False, static = True)
        hull2.append(hull2[0])

    else:
        hull2 = computeBareHull(node2, sweep = False)
        hull2.append(hull2[0])
    
    
    return hull2, medial2




class Paths:
    
    
    def __init__(self, nodeHash):
            
        self.nodeHash = nodeHash
        
        self.pathIDs = 0

        self.paths = {0 : []}
        self.hulls = {0 : []}
        self.pathTermsVisited = {0: False}
        self.trimmedPaths  = {}
        self.pathClasses = {}
        self.pathClasses[0] = {"parentID" : None, "branchNodeID" : None, "localJunctionPose" : None, 
                            "sameProb" : {}, "nodeSet" : [], "globalJunctionPose" : None}
        
        self.mergeMapping = {}
        
        self.consistency = {}
        
        self.rootPoint = [-3.2, 0.0]
        
        self.pathGraph = graph.graph()
        self.joins = []
        
        self.graphCount = 0
        self.medialCount = 1000
        self.spliceCount = 0
        self.pathPlotCount = 0
        self.pathPlotCount2 = 0
        self.overlapPlotCount = 0
        self.overlapPlotCount2 = 0
        self.topCount = 0
        self.medialSoupCount = 0
        self.pathIDs += 1
        self.multiDepCount = 0
        self.colors = []
        for i in range(1000):
            self.colors.append((random.random(),random.random(),random.random()))
    
    def saveState(self, id):
        
        saveFile = ""
        
        saveFile += "self.pathClasses = " + repr(self.pathClasses) + "\n"
        saveFile += "self.consistency = " + repr(self.consistency) + "\n"
        
        saveFile += "self.pathTermsVisited = " + repr(self.pathTermsVisited) + "\n"
        saveFile += "self.pathIDs = " + repr(self.pathIDs) + "\n"        
    
        saveFile += "self.paths = " + repr(self.paths) + "\n"        
        saveFile += "self.hulls = " + repr(self.hulls) + "\n"        
        saveFile += "self.trimmedPaths = " + repr(self.trimmedPaths) + "\n"        

        f = open("pathStateSave_%04u.txt" % id, 'w')
        f.write(saveFile)
        f.close()
        
    def restoreState(self, dirName, numNodes):
        
        print "loading" + dirName + "/pathStateSave_%04u.txt" % (numNodes-1)
        f = open(dirName + "/pathStateSave_%04u.txt" % (numNodes-1), 'r')        
        saveStr = f.read()
        print saveStr
        f.close()
        
        saveStr = saveStr.replace('\r\n','\n')
        
        exec(saveStr)
        
        #print self.numNodes
        #print self.edgePriorityHash
            
    def generatePaths(self):
        
        " COMPUTE MEDIAL AXIS FROM UNION OF PATH-CLASSIFIED NODES "
        self.paths = {}
        self.hulls = {}
        
        pathIDs = self.getPathIDs()
        for k in pathIDs:
            print "computing path for node set", k, ":", self.getNodes(k)
            self.paths[k], self.hulls[k] = self.getTopology(k)

        #self.getUnionTopology()


        " compute the junction points between parent paths and child branches "
        for pathID in pathIDs:
            
            if self.pathClasses[pathID]["parentID"] != None:

                parentPathID = self.pathClasses[pathID]["parentID"]
                childPathID = pathID

                path1 = self.paths[parentPathID]
                path2 = self.paths[childPathID]
                

                globJuncPose = self.getBranchPoint(parentPathID, childPathID, self.paths)
                print "generated globJuncPose:",  globJuncPose
                
                self.pathClasses[childPathID]["globalJunctionPose"] = globJuncPose



        self.trimmedPaths = self.trimPaths(self.paths)

        " for each path, attempt to join with its parent path "
        self.joins = []
        self.junctions = {}
        for pathID in pathIDs:

            cPath = self.getPath(pathID)
            
            parentPathID = cPath["parentID"]
            
            " parent does not concern us "
            if parentPathID == None:
                continue
            
            junctionNodeID = cPath["branchNodeID"]
            localJunctionPoint = cPath["localJunctionPose"]
            poseOrigin = Pose(self.nodeHash[junctionNodeID].getEstPose())
            junctionPoint = poseOrigin.convertLocalToGlobal(localJunctionPoint)
            junctionPose = poseOrigin.convertLocalOffsetToGlobal(localJunctionPoint)

            junctionPose = self.getGlobalJunctionPose(pathID)
            #junctionPose = cPath["globalJunctionPose"]
            junctionPoint = [junctionPose[0],junctionPose[1]]

            print "node globJuncPose:",  junctionPose

            path1 = self.trimmedPaths[pathID]
            path2 = self.trimmedPaths[parentPathID]
            
            minDist1 = 1e100
            minI1 = 0        
            for i in range(len(path1)):
                pnt = path1[i]
                dist = sqrt((pnt[0]-junctionPoint[0])**2 + (pnt[1]-junctionPoint[1])**2)
            
                if dist < minDist1:
                    minDist1 = dist
                    minI1 = i
    
            minDist2 = 1e100
            minI2 = 0        
            for i in range(len(path2)):
                pnt = path2[i]
                dist = sqrt((pnt[0]-junctionPoint[0])**2 + (pnt[1]-junctionPoint[1])**2)
            
                if dist < minDist2:
                    minDist2 = dist
                    minI2 = i
            
            self.joins.append([(pathID, minI1),(parentPathID, minI2), junctionPoint])


            " get junctions " 
            #localPose = cPath["localJunctionPose"]
            branchNodeID = cPath["branchNodeID"]
            
            #poseOrigin = Pose(self.nodeHash[branchNodeID].getEstPose())
            #junctionPoint = poseOrigin.convertLocalToGlobal(localPose)

            self.junctions[pathID] = [branchNodeID, junctionPose, (parentPathID,minI2), path2[minI2], minI1]
            
            #self.junctions[pathID] = [branchNodeID, junctionPoint, (parentPathID,minI2), path2[minI2], minI1]

            #globalPose = self.nodeHash[branchNodeID].getGlobalGPACPose()
            #junctions[pathID] = [cPath["branchNodeID"], cPath["localJunctionPose"]]
            


        " create a tree with the node IDs and then stitch them together with joins "
        #pylab.clf()
        
        
        self.pathGraph = graph.graph()
        for pathID in pathIDs:
            path = self.trimmedPaths[pathID]
            for k in range(len(path)):
                self.pathGraph.add_node((pathID, k), path[k])

            xP = []
            yP = []            
            if len(path) > 0:
                pass
                #xP = [path[0][0]]
                #yP = [path[0][1]]
            for k in range(len(path)-1):
                self.pathGraph.add_edge((pathID, k), (pathID, k+1))
                #xP.append(path[k+1][0])
                #yP.append(path[k+1][1])
            
            #pylab.plot(xP,yP, color='k')
            
            cPath = self.getPath(pathID)            
            parentPathID = cPath["parentID"]
            
            " parent does not concern us "
            """
            if parentPathID != None:
                
                junctionNodeID = cPath["branchNodeID"]
                localJunctionPoint = cPath["localJunctionPose"]
                poseOrigin = Pose(self.nodeHash[junctionNodeID].getEstPose())
                junctionPoint = poseOrigin.convertLocalToGlobal(localJunctionPoint)

                #xP = [junctionPoint[0]]
                #yP = [junctionPoint[1]]
                #pylab.scatter(xP,yP, color='r')
            """
            
        " join with the junction in between the join points "
        for k in range(len(self.joins)):
            join = self.joins[k]
            #pathGraph.add_node(k, join[2])
            #pathGraph.add_edge(join[0], k)
            #pathGraph.add_edge(k, join[1])
            self.pathGraph.add_edge(join[0], join[1])

            pID1, k1 = join[0]
            pID2, k2 = join[1]

            #xP = [self.trimmedPaths[pID1][k1][0], self.trimmedPaths[pID2][k2][0]]
            #yP = [self.trimmedPaths[pID1][k1][1], self.trimmedPaths[pID2][k2][1]]
            #pylab.plot(xP,yP, color='k')

        #pylab.title("%s" % repr(pathIDs))
        #pylab.savefig("splices_%04d.png" % self.graphCount)
        self.graphCount += 1

        self.topDict = {}
        
        
        " get terminals "
        self.terminals = {}
        for pathID in pathIDs:
            path = self.trimmedPaths[pathID]
            
            if len(path) > 0:
                if pathID == 0:
                    
                    dist1 = sqrt((self.trimmedPaths[pathID][0][0] - self.rootPoint[0])**2 + (self.trimmedPaths[pathID][0][1] - self.rootPoint[1])**2)
                    dist2 = sqrt((self.trimmedPaths[pathID][len(path)-1][0] - self.rootPoint[0])**2 + (self.trimmedPaths[pathID][len(path)-1][1] - self.rootPoint[1])**2)

                    if dist1 < dist2:
                        self.terminals[0] = [(pathID,0), path[0]]
                        self.terminals[1] = [(pathID,len(path)-1), path[len(path)-1]]
                        self.topDict["t%u" % 0] = (pathID,0)
                        self.topDict["t%u" % 1] = (pathID,len(path)-1)
                        self.rootPoint = self.trimmedPaths[pathID][0]

                    else:
                        self.terminals[1] = [(pathID,0), path[0]]
                        self.terminals[0] = [(pathID,len(path)-1), path[len(path)-1]]
                        self.topDict["t%u" % 1] = (pathID,0)
                        self.topDict["t%u" % 0] = (pathID,len(path)-1)
                        self.rootPoint = self.trimmedPaths[pathID][len(path)-1]
                    
                else:
                    
                    self.topDict["j%u" % pathID] = self.junctions[pathID][2]
                    minI1 = self.junctions[pathID][4]
                    
                    " determine which side is junction and which side is terminal"
                    if minI1 > len(path)-1 - minI1:    
                        self.terminals[pathID+1] = [(pathID,0), path[0]]
                        self.topDict["t%u" % (pathID+1)] = (pathID,0)
                    else:
                        self.terminals[pathID+1] = [(pathID,len(path)-1), path[len(path)-1]]
                        self.topDict["t%u" % (pathID+1)] = (pathID,len(path)-1)
                    #self.terminals[pathID+1] = [(pathID,len(path)-1), path[len(path)-1]]
                    #self.topDict["t%u" % (pathID+1)] = (pathID,len(path)-1)


    def getGlobalJunctionPose(self, pathID):

        try:
            globJuncPose = self.pathClasses[pathID]["globalJunctionPose"]
        except:
            
            branchNodeID = self.pathClasses[pathID]["branchNodeID"]
            
            if branchNodeID == None:
                return None
            
            localJunctionPose = self.pathClasses[pathID]["localJunctionPose"]
            
            poseOrigin = Pose(self.nodeHash[branchNodeID].getEstPose())
            globJuncPose = poseOrigin.convertLocalOffsetToGlobal(localJunctionPose)

        return globJuncPose
            
    def getNodes(self, pathID):
        
        return self.pathClasses[pathID]["nodeSet"]
    
    def getPathIDs(self):
        #print "keys:", self.pathClasses.keys()
        return self.pathClasses.keys()
    
    def getParentPathID(self, pathID):
        
        parentID = self.pathClasses[pathID]["parentID"]
    
        return parentID
    
    def getPath(self, pathID):
        return self.pathClasses[pathID]
    
    def getAllSplices(self, plotIter = False):

        print "junctions:", self.junctions
        print "terminals:", self.terminals

        " determine which paths are leaves "
        
        pathIDs = self.getPathIDs()
        isAParent = {}
        parents = {}
        pathIDGraph = graph.graph()

        for pathID in pathIDs:
            isAParent[pathID] = False
            pathIDGraph.add_node(pathID, [])
        
        for pathID in pathIDs:    
            parentPathID = self.getPath(pathID)["parentID"]
            parents[pathID] = parentPathID
            if parentPathID != None:
                isAParent[parentPathID] = True
                pathIDGraph.add_edge(pathID, parentPathID)
        
        leafCount = 0
        for pathID in pathIDs:
            if not isAParent[pathID]:
                leafCount += 1

        #self.pathClasses[0] = {"parentID" : None, "branchNodeID" : None, "localJunctionPose" : None, 
        

        " build topological graph "
        topGraph = graph.graph()
        
        topNodes = []
                       
        for k, term in self.terminals.iteritems():
            topGraph.add_node(term[0], term[1])
            print "add_node", term[0], term[1]
            topNodes.append(term[0])
        for k, junc in self.junctions.iteritems():
            topGraph.add_node(junc[2], junc[3])
            print "add_node", junc[2], junc[3]
            topNodes.append(junc[2])
            
        
        for pathID in pathIDs:
            pathIndex = []
            if pathID == 0:
                term = self.terminals[0]
                graphNodeKey = term[0]
                pathIndex.append(graphNodeKey[1])
                
                term = self.terminals[1]
                graphNodeKey = term[0]
                pathIndex.append(graphNodeKey[1])
            else:
                term = self.terminals[pathID+1]
                graphNodeKey = term[0]
                pathIndex.append(graphNodeKey[1])
                
            for k, junc in self.junctions.iteritems():
                graphNodeKey = junc[2]
                if graphNodeKey[0] == pathID:
                    pathIndex.append(graphNodeKey[1])
                        
            pathIndex.sort()
            
            print "pathIndex:", pathIndex
            
            for k in range(len(pathIndex)-1):
                print "add_edge", (pathID,pathIndex[k]),(pathID,pathIndex[k+1])
                topGraph.add_edge((pathID,pathIndex[k]),(pathID,pathIndex[k+1]))

            if pathID != 0:
                graphNodeKey = self.junctions[pathID][2]
                print "add_edge", graphNodeKey,(pathID,pathIndex[0])
                topGraph.add_edge(graphNodeKey,(pathID,pathIndex[0]))
                

        results = {}
        for i in range(len(pathIDs)):
            pathID1 = pathIDs[i]
            for j in range(i,len(pathIDs)):
                
                pathID2 = pathIDs[j]
                
                if pathID1 == pathID2:
                    "term1 to term 2"

                    termPaths = []
                    if pathID1 == 0:
                        termIDs = [self.topDict["t%u" % (pathID1)], self.topDict["t%u" % (pathID1+1)]]
                    else:
                        termIDs = [self.topDict["j%u" % (pathID1)], self.topDict["t%u" % (pathID1+1)]]
                    termPaths.append(termIDs)  

                    orderedPathIDs = [pathID1]
                else:
                    
                    startPathID = pathID1
                    endPathID = pathID2
        
                    shortestPathSpanTree, shortestDist = pathIDGraph.shortest_path(endPathID)
                    nextPathID = shortestPathSpanTree[startPathID]
            
                    orderedPathIDs = [startPathID, nextPathID]
                    while nextPathID != endPathID:
                        nextPathID = shortestPathSpanTree[nextPathID]
                        orderedPathIDs.append(nextPathID)            

                    
                    " start at first path ID "
                    termPaths = []
                    if parents[orderedPathIDs[0]] == orderedPathIDs[1]:
                        " going to parent "
                        termIDs = [self.topDict["t%u" % (orderedPathIDs[0]+1)]]
                        termPaths.append(termIDs)
                    
                    elif parents[orderedPathIDs[1]] == orderedPathIDs[0]:
                        " going to child "
                        termIDs = [self.topDict["t%u" % (orderedPathIDs[0]+1)]]
                        termPaths.append(termIDs)

                        if orderedPathIDs[0] == 0:
                            termIDs = [self.topDict["t%u" % 0]]
                            termPaths.append(termIDs)
                        else:
                            termIDs = [self.topDict["j%u" % orderedPathIDs[0]]]
                            termPaths.append(termIDs)
                    else:
                        print parents
                        print orderedPathIDs
                        print termPaths
                        raise
                        
                    " transitions between middle path IDs "
                    for k in range(len(orderedPathIDs)-1):
                        if parents[orderedPathIDs[k]] == orderedPathIDs[k+1]:
                            " child -> parent:  prev -> junction "
                            for ids in termPaths:
                                ids.append(self.topDict["j%u" % orderedPathIDs[k]])

                        elif parents[orderedPathIDs[k+1]] == orderedPathIDs[k]:
                            " parent -> child: prev -> junction "
                            for ids in termPaths:
                                ids.append(self.topDict["j%u" % orderedPathIDs[k+1]])

                        else:
                            print parents
                            print orderedPathIDs
                            print termPaths
                            print k
                            raise
                    
                    " last path ID "
                    if parents[orderedPathIDs[-2]] == orderedPathIDs[-1]:
                        " child -> parent:  prev -> junction "
                        oldTermPaths = termPaths
                        termPaths = []
                        " split into 2 possibilities into parent path ID "
                        for ids in oldTermPaths:
                            ids1 = deepcopy(ids)
                            ids2 = deepcopy(ids)
                            if orderedPathIDs[-1] == 0:
                                ids1.append(self.topDict["t%u" % orderedPathIDs[-1]])
                                ids2.append(self.topDict["t%u" % (orderedPathIDs[-1]+1) ])                                
                            else:
                                ids1.append(self.topDict["j%u" % orderedPathIDs[-1]])
                                ids2.append(self.topDict["t%u" % (orderedPathIDs[-1]+1) ])
                            
                            termPaths.append(ids1)
                            termPaths.append(ids2)

                    elif parents[orderedPathIDs[-1]] == orderedPathIDs[-2]:
                        " parent -> child: prev -> junction "
                        for ids in termPaths:
                            ids.append(self.topDict["t%u" % (orderedPathIDs[-1]+1) ])
                                        
                    else:
                        print parents
                        print orderedPathIDs
                        print termPaths
                        raise
                
                finalResults = []
                
                for termPath in termPaths:

                    # self.joins.append([(pathID, minI1),(parentPathID, minI2), junctionPoint])

                    joinPairs = []

                    startKey = termPath[0]
                    endKey = termPath[-1]

                    shortestPathSpanTree, shortestDist = self.pathGraph.shortest_path(endKey)
                    currNode = shortestPathSpanTree[startKey]                    
                    splicedPath = []
                    while currNode != endKey:
                        splicedPath.append(self.pathGraph.get_node_attributes(currNode))
                        nextNode = shortestPathSpanTree[currNode]
                        
                        for join in self.joins:
                            if join[0] == currNode and join[1] == nextNode:
                                joinPairs.append((len(splicedPath)-1,len(splicedPath)))
                            elif join[1] == currNode and join[0] == nextNode:
                                joinPairs.append((len(splicedPath)-1,len(splicedPath)))
                            
                        currNode = nextNode
                        
                    splicedPath.append(self.pathGraph.get_node_attributes(currNode))


                    " if there are any joins, we must interpolate a smooth transition "
                    lastIndex = 0
                    newPath = []
                    for pair in joinPairs:
                        index1 = pair[0]
                        index2 = pair[1]
                        cPoints = [splicedPath[index1-1], splicedPath[index1], splicedPath[index2], splicedPath[index2+1]]               
                        spline = SplineFit(cPoints, smooth=0.0, kp=3)
                        points = spline.getUniformSamples(spacing = 0.1)
                        
                        newPath += splicedPath[lastIndex:index1]
                        newPath += points
                        
                        lastIndex = index2+2

                    newPath += splicedPath[lastIndex:]


                    sPath = {}
                    sPath['orderedPathIDs'] = orderedPathIDs
                    sPath['path'] = newPath
                    sPath['termPath'] = termPath
                    
                    finalResults.append(sPath)
                    #finalResults.append([termPath,splicedPath])
                    
                results[(pathID1,pathID2)] = finalResults
                #results.append(termPaths)
        
        print "results:"
        for k, result in results.iteritems():
            print k, ":"
            for sPath in result:
                print sPath['orderedPathIDs'], sPath['termPath'], len(sPath['path'])


        if plotIter:
            pylab.clf()
            
            for k,path in  self.trimmedPaths.iteritems():
                print "path has", len(path), "points"
                xP = []
                yP = []
                for p in path:
                    xP.append(p[0])
                    yP.append(p[1])
    
                pylab.plot(xP,yP, color = self.colors[k])
    
    
            for k, result in results.iteritems():
                for sPath in result:
                    path = sPath['path']
    
                    xP = []
                    yP = []
                    for p in path:
                        xP.append(p[0])
                        yP.append(p[1])
                    pylab.plot(xP,yP, color='b')


            for k in pathIDs:
                xP = []
                yP = []
                
                for p in self.paths[k]:
                    xP.append(p[0])
                    yP.append(p[1])
                pylab.plot(xP,yP, color=self.colors[k], linewidth=1)


            #pylab.xlim(-4,4)
            #pylab.ylim(-4,4)
    
            #self.drawWalls()
    
            pylab.title("Spliced Paths, pathIDs = %s" %  self.getPathIDs())
            pylab.savefig("splicedPath_%04u.png" % self.spliceCount)
            self.spliceCount += 1
        


        return results, self.terminals, self.junctions

    
    def isNodeExist(self, nodeID, pathID):
        return nodeID in self.pathClasses[pathID]["nodeSet"]

    def findSplice2(self, nodeID):

        splicePaths, spliceTerms, splicePathIDs = self.getSplicesByNearJunction(nodeID)
        
        node2 = self.nodeHash[nodeID]
        hull2, medial2 = computeHullAxis(nodeID, node2, tailCutOff = False)
        estPose2 = node2.getGlobalGPACPose()        

        resultSet = []

        for k in range(len(splicePaths)):
            
            currPath = splicePaths[k]
            pathIDs = splicePathIDs[k]
            results = getMultiDeparturePoint(currPath, medial2, estPose2, estPose2, pathIDs, nodeID)

            resultSet.append(results+(k,))

            "departurePoint1, angle1, isInterior1, isExist1, dist1, departurePoint2, angle2, isInterior2, isExist2, dist2, contigFrac, overlapSum, angDiff2"
        
        resultSet = sorted(resultSet, key=itemgetter(10), reverse=True)
        
        print "findSplice2 resultSet:"
        for result in resultSet:
            print result


        spliceIndex = resultSet[0][13]
        overlappedPathIDs = splicePathIDs[spliceIndex]
        orderedPathIDs = self.getPathOrdering(nodeID, overlappedPathIDs)

        print "returned ordered", orderedPathIDs, "from unordered", overlappedPathIDs
        
        return orderedPathIDs
        

    def findSplice(self, nodeID):

        node2 = self.nodeHash[nodeID]
  
        hull2, medial2 = computeHullAxis(nodeID, node2, tailCutOff = False)

        
        estPose2 = node2.getGlobalGPACPose()        
                    
        " set the initial guess "
        poseOrigin = Pose(estPose2)
        
        medialOffset = []
        for pnt in medial2:
            medialOffset.append(poseOrigin.convertLocalToGlobal(pnt))
     

        medialSpline2 = SplineFit(medialOffset, smooth=0.1)
        globalMedial = medialSpline2.getUniformSamples(spacing = 0.04)
     
     
        #allSplices = self.getAllSplices()
        allSplices, terminals, junctions = self.getAllSplices()
        
        selects = []
        
        for pathPair, sPaths in allSplices.iteritems():
            
            for sPath in sPaths:
                termPath1 = sPath['termPath']
                
                print "termPath1:", termPath1
        
                startKey = termPath1[0]
                endKey = termPath1[-1]
        
                shortestPathSpanTree, shortestDist = self.pathGraph.shortest_path(endKey)
                currNode = shortestPathSpanTree[startKey]                    
                path1 = []
                while currNode != endKey:
                    path1.append(self.pathGraph.get_node_attributes(currNode))
                    currNode = shortestPathSpanTree[currNode]
                path1.append(self.pathGraph.get_node_attributes(currNode))

                matchCount = self.getPathOverlapCondition2(globalMedial, path1, nodeID, -1)
                print "matchCount =", matchCount
                pathIDs = sPath['orderedPathIDs']
                selects.append((matchCount, pathIDs))


                if False:
                    departurePoint1, angle1, isInterior1, isExist1, dist1, departurePoint2, angle2, isInterior2, isExist2, dist2 = self.getPathDeparture2(globalMedial, path1, 0, 0, termPathIDs1 = termPath1)
                    
                    print 'departure results:'
                    print "interiors:", isInterior1, isInterior2
                    print "exists:", isExist1, isExist2

                    if not isInterior1 and not isInterior2:
                        selects.append(termPath1)

        selects.sort(reverse = True)
        overlappedPaths = selects[0][1]

        print "unordered:", overlappedPaths

        orderedPathIDs = self.getPathOrdering(nodeID, overlappedPaths)

        " already ordered! "
        #orderedPathIDs = overlappedPaths

        print "ordered:", orderedPathIDs

        #print "findSplices(", nodeID,") = "
        #print selects
        
        if len(orderedPathIDs) == 0:
            raise
        
        return orderedPathIDs


    def getSplicesByNearJunction(self, nodeID):
        
        allSplices, terminals, junctions = self.getAllSplices(plotIter = True)

        initPose2 = self.nodeHash[nodeID].getGlobalGPACPose()
        
        print "junctions:", junctions
        print "initPose2:", initPose2
        
        junctions2 = []
        for pathID, params in junctions.iteritems():
            junctionPose = params[1]
            
            print "checking", pathID, params
            
            dist2 = sqrt((junctionPose[0]-initPose2[0])**2 + (junctionPose[1]-initPose2[1])**2)
            
            print "dist2 =", dist2
            
            if dist2 < 3.0:
                junctions2.append((pathID,params[2][0]))
            
        "self.junctions[pathID] = [branchNodeID, junctionPoint, (parentPathID,minI2), path2[minI2], minI1]"
        
        closePathID2 = self.getClosestPath(initPose2)

        print "closePathID2:", closePathID2
        
        pathSet2 = [closePathID2]

        junctionSet = junctions2
        junctionSet = list(set(junctionSet))

        print "junctions2:", junctions2

        for junc in junctionSet:
            pathSet2.append(junc[0])
            pathSet2.append(junc[1])

        print "pathSet2:", pathSet2

        pathSet2 = list(set(pathSet2))
        
        pathSet2.sort()

        print "pathSet2:", pathSet2
        
        termSet2 = []
        for pathID in pathSet2:
            if pathID == 0:
                termSet2.append(terminals[0][0])
                termSet2.append(terminals[1][0])
            else:
                termSet2.append(terminals[pathID+1][0])
                
        for pathID in pathSet2:
            parentID = self.getParentPathID(pathID)

            if parentID != None and not parentID in pathSet2:
                termSet2.append(junctions[pathID][2])
            

        

        splicePaths = []
        spliceTerms = []
        splicePathIDs = []

        print "termSet2:", termSet2
        print "allSplices:"

        for k, result in allSplices.iteritems():
            print "id:", k
            if k[0] in pathSet2 and k[1] in pathSet2:
                
                for sPath in result:
                    termPath = sPath['termPath']
                    print "termPath:", termPath
                                                
                    if termPath[0] in termSet2 and termPath[-1] in termSet2:
                        splicePaths.append(sPath['path'])
                        spliceTerms.append(termPath)
                        splicePathIDs.append(sPath['orderedPathIDs'])
                        


        return splicePaths, spliceTerms, splicePathIDs

    def checkNodeConsistency(self, nodeID1, nodeID2):

        nodeSet = [nodeID1, nodeID2]

        resultSet = []
        
        for nodeID in nodeSet:

            pose0 = self.nodeHash[nodeID].getGlobalGPACPose()
            hull0, medial0 = computeHullAxis(nodeID, self.nodeHash[nodeID], tailCutOff = False)

            splicePaths, termSet2, pathSet2 = self.getSplicesByNearJunction(nodeID)
   
            for path in splicePaths:        

                orientedSplicePath = orientPath(path, medial0, pose0)    
                
                " (departurePoint1, angle1, isInterior1, isExist1, dist1, departurePoint2, angle2, isInterior2, isExist2, dist2, contigFrac, overlapSum, angDiff2 )"
                result0 = getMultiDeparturePoint(orientedSplicePath, medial0, pose0, pose0, [], nodeID, pathPlotCount = self.multiDepCount, plotIter = True)
                self.multiDepCount += 1

                
                isInterior1 = result0[2]
                isExist1 = result0[3]
                isInterior2 = result0[7]
                isExist2 = result0[8]
                
                contigFrac = result0[10]
                angDiff = result0[12]
    
                resultSet.append((nodeID, contigFrac, isInterior1, isInterior2, isExist1, isExist2))

        
        print "consistency for node", nodeID1, "and", nodeID2
        for result in resultSet:
            print "node", result[0], ":", result[1], result[2], result[3], result[4], result[5]




    def checkPathConsistency(self, targetPathID):

        #path = self.trimmedPaths[pathID]
        nodeSet = self.getNodes(targetPathID)

        resultSet = []

        for nodeID in nodeSet:
            pose0 = self.nodeHash[nodeID].getGlobalGPACPose()
            hull0, medial0 = computeHullAxis(nodeID, self.nodeHash[nodeID], tailCutOff = False)

            allSplices, terminals, junctions = self.getAllSplices(plotIter = True)

            initPose2 = self.nodeHash[nodeID].getGlobalGPACPose()
            
            print "junctions:", junctions
            print "initPose2:", initPose2
            
            junctions2 = []
            for pathID, params in junctions.iteritems():
                junctionPose = params[1]
                
                print "checking", pathID, params
                
                dist2 = sqrt((junctionPose[0]-initPose2[0])**2 + (junctionPose[1]-initPose2[1])**2)
                
                print "dist2 =", dist2
                
                if dist2 < 3.0:
                    junctions2.append((pathID,params[2][0]))
                
            "self.junctions[pathID] = [branchNodeID, junctionPoint, (parentPathID,minI2), path2[minI2], minI1]"
            
            closePathID2 = self.getClosestPath(initPose2)

            print "closePathID2:", closePathID2
            
            pathSet2 = [closePathID2]

            junctionSet = junctions2
            junctionSet = list(set(junctionSet))

            print "junctions2:", junctions2

            for junc in junctionSet:
                pathSet2.append(junc[0])
                pathSet2.append(junc[1])

            print "pathSet2:", pathSet2

            pathSet2 = list(set(pathSet2))
            
            pathSet2.sort()

            print "pathSet2:", pathSet2
            
            termSet2 = []
            for pathID in pathSet2:
                if pathID == 0:
                    termSet2.append(terminals[0][0])
                    termSet2.append(terminals[1][0])
                else:
                    termSet2.append(terminals[pathID+1][0])
                    
            for pathID in pathSet2:
                parentID = self.getParentPathID(pathID)

                if parentID != None and not parentID in pathSet2:
                    termSet2.append(junctions[pathID][2])
                

            

            splicePaths = []

            print "termSet2:", termSet2
            print "allSplices:"

            for k, result in allSplices.iteritems():
                print "id:", k
                if k[0] in pathSet2 and k[1] in pathSet2:
                    
                    for sPath in result:
                        termPath = sPath['termPath']
                        print "termPath:", termPath
                                                    
                        if termPath[0] in termSet2 and termPath[-1] in termSet2:
                            splicePaths.append(sPath['path'])
                            
            for path in splicePaths:        

                orientedSplicePath = orientPath(path, medial0, pose0)    
                
                " (departurePoint1, angle1, isInterior1, isExist1, dist1, departurePoint2, angle2, isInterior2, isExist2, dist2, contigFrac, overlapSum, angDiff2 )"
                result0 = getMultiDeparturePoint(orientedSplicePath, medial0, pose0, pose0, [], nodeID, pathPlotCount = self.multiDepCount, plotIter = True)
                self.multiDepCount += 1

                
                isInterior1 = result0[2]
                isExist1 = result0[3]
                isInterior2 = result0[7]
                isExist2 = result0[8]
                
                contigFrac = result0[10]
                angDiff = result0[12]
    
                resultSet.append((nodeID, contigFrac, isInterior1, isInterior2, isExist1, isExist2))

        
        print "consistency for path", targetPathID
        for result in resultSet:
            print "node", result[0], ":", result[1], result[2], result[3], result[4], result[5]

    
    def comparePaths(self):
 
        #allSplices = self.getAllSplices()
        allSplices, terminals, junctions = self.getAllSplices(plotIter = True)
        
        toBeMerged = []
        print "consistency:", self.consistency
        
        for pathPair,value in self.consistency.iteritems():
            pathID1 = pathPair[0]
            pathID2 = pathPair[1]

            path1 = self.getPath(pathID1)
            path2 = self.getPath(pathID2)
            
            parentID1 = path1["parentID"]
            parentID2 = path2["parentID"]
            
            
            if parentID1 != None:
                pathParent1 = self.getPath(parentID1)
                grandParentID1 = pathParent1["parentID"]
            else:
                pathParent1 = None
                grandParentID1 = None

            if parentID2 != None:
                pathParent2 = self.getPath(parentID2)
                grandParentID2 = pathParent2["parentID"]
            else:
                pathParent2 = None
                grandParentID2 = None
            
            
            self.consistency[pathPair] = 0.0     

            
            " Only consider case when two paths share the same parent "
            pathPairs = []
            
            rootPaths = []

            #if True:
            #if pathID1 == 0 and pathID2 == 3 or pathID1 == 3 and pathID2 == 0:

            print pathID1, pathID2, parentID1, parentID2
            
            if False and grandParentID1 == pathID2 and parentID1 != None:
                
                if parentID2 != None:
                    sPaths1 = allSplices[(parentID2, pathID1)]
                    sPaths2 = allSplices[(parentID2, pathID2)]
                else:
                    sPaths1 = allSplices[(pathID2, pathID1)]
                    sPaths2 = allSplices[(pathID2, pathID2)]

                for sPath1 in sPaths1:
                    termPath1 = sPath1['termPath']
                    print "termPath1:", termPath1
                    startKey = termPath1[0]
                    endKey = termPath1[-1]

                    shortestPathSpanTree, shortestDist = self.pathGraph.shortest_path(endKey)
                    currNode = shortestPathSpanTree[startKey]                    
                    path1 = []
                    while currNode != endKey:
                        path1.append(self.pathGraph.get_node_attributes(currNode))
                        currNode = shortestPathSpanTree[currNode]
                    path1.append(self.pathGraph.get_node_attributes(currNode))

                    rootPaths.append(path1)

                    for sPath2 in sPaths2:
                        termPath2 = sPath2['termPath']
                        
                        print "termPath2:", termPath2


                        startKey = termPath2[0]
                        endKey = termPath2[-1]
    
                        shortestPathSpanTree, shortestDist = self.pathGraph.shortest_path(endKey)
                        currNode = shortestPathSpanTree[startKey]                    
                        path2 = []
                        while currNode != endKey:
                            path2.append(self.pathGraph.get_node_attributes(currNode))
                            currNode = shortestPathSpanTree[currNode]
                        path2.append(self.pathGraph.get_node_attributes(currNode))

                        pathPairs.append((path1,path2,pathID1,pathID2))

            elif False and grandParentID2 == pathID1 and parentID2 != None:
                
                if parentID1 != None:
                    sPaths1 = allSplices[(parentID1, pathID1)]
                    sPaths2 = allSplices[(parentID1, pathID2)]
                else:
                    sPaths1 = allSplices[(pathID1, pathID1)]
                    sPaths2 = allSplices[(pathID1, pathID2)]

                for sPath1 in sPaths1:
                    termPath1 = sPath1['termPath']
                    print "termPath1:", termPath1
                    startKey = termPath1[0]
                    endKey = termPath1[-1]

                    shortestPathSpanTree, shortestDist = self.pathGraph.shortest_path(endKey)
                    currNode = shortestPathSpanTree[startKey]                    
                    path1 = []
                    while currNode != endKey:
                        path1.append(self.pathGraph.get_node_attributes(currNode))
                        currNode = shortestPathSpanTree[currNode]
                    path1.append(self.pathGraph.get_node_attributes(currNode))

                    rootPaths.append(path1)
                    
                    for sPath2 in sPaths2:
                        termPath2 = sPath2['termPath']
                        
                        print "termPath2:", termPath2


                        startKey = termPath2[0]
                        endKey = termPath2[-1]
    
                        shortestPathSpanTree, shortestDist = self.pathGraph.shortest_path(endKey)
                        currNode = shortestPathSpanTree[startKey]                    
                        path2 = []
                        while currNode != endKey:
                            path2.append(self.pathGraph.get_node_attributes(currNode))
                            currNode = shortestPathSpanTree[currNode]
                        path2.append(self.pathGraph.get_node_attributes(currNode))

                        pathPairs.append((path1,path2,pathID1,pathID2))

            elif parentID1 == pathID2:
                
                if parentID2 != None:
                    sPaths1 = allSplices[(parentID2, pathID1)]
                    sPaths2 = allSplices[(parentID2, pathID2)]
                else:
                    sPaths1 = allSplices[(pathID2, pathID1)]
                    sPaths2 = allSplices[(pathID2, pathID2)]

                for sPath1 in sPaths1:
                    termPath1 = sPath1['termPath']
                    print "termPath1:", termPath1

                    startKey = termPath1[0]
                    endKey = termPath1[-1] 
                    if len(termPath1) > 2:      
                        midKey = termPath1[1]
                    else:
                        midKey = None

                    shortestPathSpanTree, shortestDist = self.pathGraph.shortest_path(endKey)
                    currNode = shortestPathSpanTree[startKey]                    
                    path1 = []
                    midIndex = 0
                    while currNode != endKey:

                        if currNode == midKey:
                            midIndex = len(path1)

                        path1.append(self.pathGraph.get_node_attributes(currNode))
                        currNode = shortestPathSpanTree[currNode]
                        
                        
                    path1.append(self.pathGraph.get_node_attributes(currNode))
                    
                    print "path1:", len(path1), midIndex, midKey
                    
                    if midIndex > 150:
                        path1 = path1[midIndex-150:]

                    rootPaths.append(path1)
                    
                    for sPath2 in sPaths2:
                        termPath2 = sPath2['termPath']
                        
                        print "termPath2:", termPath2


                        startKey = termPath2[0]
                        endKey = termPath2[-1]
                        if len(termPath2) > 2:      
                            midKey = termPath2[1]
                        else:
                            midKey = None
                            
                        shortestPathSpanTree, shortestDist = self.pathGraph.shortest_path(endKey)
                        currNode = shortestPathSpanTree[startKey]                    
                        path2 = []
                        midIndex = 0
                        while currNode != endKey:
                            
                            if currNode == midKey:
                                midIndex = len(path2)
                                                            
                            path2.append(self.pathGraph.get_node_attributes(currNode))
                            currNode = shortestPathSpanTree[currNode]

                        path2.append(self.pathGraph.get_node_attributes(currNode))

                        print "path2:", len(path2), midIndex, midKey

                        if midIndex > 150:
                            path2 = path2[midIndex-150:]
                            
                        pathPairs.append((path1,path2,pathID1,pathID2))

            elif parentID2 == pathID1:
                
                if parentID1 != None:
                    sPaths1 = allSplices[(parentID1, pathID1)]
                    sPaths2 = allSplices[(parentID1, pathID2)]
                else:
                    sPaths1 = allSplices[(pathID1, pathID1)]
                    sPaths2 = allSplices[(pathID1, pathID2)]

                for sPath1 in sPaths1:

                    termPath1 = sPath1['termPath']
                    print "termPath1:", termPath1

                    startKey = termPath1[0]
                    endKey = termPath1[-1]
                    if len(termPath1) > 2:      
                        midKey = termPath1[1]
                    else:
                        midKey = None

                    shortestPathSpanTree, shortestDist = self.pathGraph.shortest_path(endKey)
                    currNode = shortestPathSpanTree[startKey]                    
                    path1 = []
                    midIndex = 0
                    
                    while currNode != endKey:
                        if currNode == midKey:
                            midIndex = len(path1)

                        path1.append(self.pathGraph.get_node_attributes(currNode))
                        currNode = shortestPathSpanTree[currNode]
                    path1.append(self.pathGraph.get_node_attributes(currNode))

                    print "path1:", len(path1), midIndex, midKey

                    if midIndex > 150:
                        path1 = path1[midIndex-150:]
                                        
                    rootPaths.append(path1)            
                    
                    for sPath2 in sPaths2:
                        termPath2 = sPath2['termPath']
                        
                        print "termPath2:", termPath2

                            
                        startKey = termPath2[0]
                        endKey = termPath2[-1]
                        if len(termPath2) > 2:      
                            midKey = termPath2[1]
                        else:
                            midKey = None
    
                        shortestPathSpanTree, shortestDist = self.pathGraph.shortest_path(endKey)
                        currNode = shortestPathSpanTree[startKey]                    
                        path2 = []
                        midIndex = 0
                        while currNode != endKey:
                            if currNode == midKey:
                                midIndex = len(path2)
                                
                            path2.append(self.pathGraph.get_node_attributes(currNode))
                            currNode = shortestPathSpanTree[currNode]
                                                                                        
                        path2.append(self.pathGraph.get_node_attributes(currNode))

                        print "path2:", len(path2), midIndex, midKey

                        if midIndex > 150:
                            path2 = path2[midIndex-150:]
                            
                        pathPairs.append((path1,path2,pathID1,pathID2))

            elif parentID1 == parentID2:
                
                values1 = allSplices[(parentID1,pathID1)]
                values2 = allSplices[(parentID2,pathID2)]
                
                print "sib:", pathID1, pathID2

                
                terms1 = []
                for sPath in values1:
                    orderedIDs = sPath['orderedPathIDs']
                    if len(orderedIDs) == 2:
                        termPath1 = sPath['termPath']
                        terms1.append(termPath1)
                        print termPath1
                
                terms2 = []
                for sPath in values2:
                    orderedIDs = sPath['orderedPathIDs']
                    if len(orderedIDs) == 2:
                        termPath2 = sPath['termPath']
                        terms2.append(termPath2)
                        print termPath2
                
                for k in range(len(terms1)):
                    
                    termPath1 = terms1[k]
                    termPath2 = terms2[k]                        
                    
                    startKey = termPath1[0]
                    endKey = termPath1[-1]

                    shortestPathSpanTree, shortestDist = self.pathGraph.shortest_path(endKey)
                    currNode = shortestPathSpanTree[startKey]                    
                    path1 = []
                    while currNode != endKey:
                        path1.append(self.pathGraph.get_node_attributes(currNode))
                        currNode = shortestPathSpanTree[currNode]
                    path1.append(self.pathGraph.get_node_attributes(currNode))

                    startKey = termPath2[0]
                    endKey = termPath2[-1]

                    shortestPathSpanTree, shortestDist = self.pathGraph.shortest_path(endKey)
                    currNode = shortestPathSpanTree[startKey]                    
                    path2 = []
                    while currNode != endKey:
                        path2.append(self.pathGraph.get_node_attributes(currNode))
                        currNode = shortestPathSpanTree[currNode]
                    path2.append(self.pathGraph.get_node_attributes(currNode))

                    pathPairs.append((path1,path2,pathID1,pathID2))
                    
                    rootPaths.append(path1)            


                pathPairs.append((None,None,pathID1,pathID2))
                pathPairs.append((None,None,pathID2,pathID1))
            
            print
            
            for pair in pathPairs:
                pathID1 = pair[2]
                pathID2 = pair[3]

                #if pathID1 == 3 or pathID2 == 3:
                if True:
                    
                    if pair[0] == None:
                        
                        print "comparing sibling paths", pathID1, pathID2
                        resultPose0, resultPose1, lastCost0, matchCount0, juncAngDiff = self.makeSiblingPathCompare(pathID1, pathID2, parentID1, plotIter = True)
                        
                        print "sibling compare result:", resultPose0, resultPose1, lastCost0, matchCount0, juncAngDiff
                        
                        if fabs(resultPose0[2]) < 0.5 and fabs(resultPose1[2]) < 0.5 and fabs(juncAngDiff) < 0.5:

                            print "queuing paths", pathID1, pathID2, "to be merged"

                            toBeMerged.append((pathID1, pathID2, resultPose0, rootPaths, resultPose1, lastCost0, matchCount0))
                    
                    else:
                        print "comparing paths", pathID1, pathID2
                        print pair[1][0], pair[0][0]
                        print pair[1][-1], pair[0][-1]
                        resultPose0, lastCost0, matchCount0 = self.makePathCompare(pair[1], pair[0], pathID2, pathID1, plotIter = True)
                    
                        if fabs(resultPose0[2]) < 0.5:
        
                            poseOrigin = Pose(resultPose0)
                            path1 = pair[1]
                            path1_offset = []
                            for p in path1:
                                result = poseOrigin.convertLocalToGlobal(p)
                                path1_offset.append(result)                    
                            
                            result0 = self.getPathDeparture(path1_offset, pair[0], pathID2, pathID1, plotIter = True)
                            print "getPathDeparture() = ", result0
        
                            isInterior1 = result0[2]
                            isInterior2 = result0[7]
        
                            if not isInterior1 and not isInterior2:
            
                                vals = allSplices[(pathID1,pathID1)]
                                sPath = vals[0]
                                termPath0 = sPath['termPath']
            
                                startKey = termPath0[0]
                                endKey = termPath0[-1]
                
                                shortestPathSpanTree, shortestDist = self.pathGraph.shortest_path(endKey)
                                currNode = shortestPathSpanTree[startKey]
                                path0 = []
                                while currNode != endKey:
                                    path0.append(self.pathGraph.get_node_attributes(currNode))
                                    currNode = shortestPathSpanTree[currNode]
                                path0.append(self.pathGraph.get_node_attributes(currNode))
                                
                                #cost = self.getPathOverlapCondition(self.trimmedPaths[pathID1], path1_offset, pathID1, pathID2)
                                cost = self.getPathOverlapCondition(path0, path1_offset, pathID1, pathID2, plotIter = True)                
                                #cost = self.getPathOverlapCondition(path2_offset, self.trimmedPaths[pathID1])
                                print "getPathOverlapCondition() =", cost
                                
                                if cost < 1e10:
                                    print pathID1, "and", pathID2, "are similar!"
                                    toBeMerged.append((pathID1, pathID2, resultPose0, rootPaths))
                                    if pathID1 < pathID2:
                                        self.consistency[(pathID1,pathID2)] = 1.0
                                    else:
                                        self.consistency[(pathID2,pathID1)] = 1.0
                                    print self.consistency
        
        return toBeMerged
                    
    
    def makeSiblingPathCompare(self, pathID1, pathID2, parentPathID, plotIter = False):


        globalPath1 = self.paths[pathID1]
        globalPath2 = self.paths[pathID2]
        parentPath = self.paths[parentPathID]

        #globalPath1_R = []
        #globalPath2_R = []
        #poseOrigin = Pose([0.0,0.0,pi/4.0])
        #for p in globalPath1:
        #    p1 = poseOrigin.convertLocalToGlobal(p)
        #    globalPath1_R.apppend(p1)

        #for p in globalPath2:
        #    p1 = poseOrigin.convertLocalToGlobal(p)
        #    globalPath2_R.apppend(p1)

        junctionPose1 = self.pathClasses[pathID1]["globalJunctionPose"]
        junctionPose2 = self.pathClasses[pathID2]["globalJunctionPose"]
        
        print "sibling pathIDs:", pathID1, pathID2
        print "junctionPose1:", junctionPose1
        print "junctionPose2:", junctionPose2
        
        " get closest point index onto the parent path to the child junctions "
        
        minDist_1 = 1e100
        minK_1 = 0
        minDist_2 = 1e100
        minK_2 = 0
        for k in range(len(parentPath)):
            p = parentPath[k]
        
            dist1 = sqrt((p[0]-junctionPose1[0])**2 + (p[1]-junctionPose1[1])**2)
            dist2 = sqrt((p[0]-junctionPose2[0])**2 + (p[1]-junctionPose2[1])**2)
        
            if dist1 < minDist_1:
                minDist_1 = dist1
                minK_1 = k

            if dist2 < minDist_2:
                minDist_2 = dist2
                minK_2 = k
        
        
        " select which path segments to use by checking if the other child path's junction point is contained within "
        
        if minK_1 < minK_2:
            pathSeg1 = parentPath[:minK_1+1]
            pathSeg2 = parentPath[minK_2:]

        else:        
            pathSeg2 = parentPath[:minK_2+1]
            pathSeg1 = parentPath[minK_1:]
        
        juncOrigin1 = Pose(junctionPose1)
        juncOrigin2 = Pose(junctionPose2)
        
        " align the two junctions "
        alignedJunctionPose2 = [junctionPose2[0], junctionPose2[1], junctionPose1[2]]
        localJunctionPose2 = juncOrigin1.convertGlobalPoseToLocal(alignedJunctionPose2)
        #localJunctionPose2 = juncOrigin1.convertGlobalPoseToLocal(junctionPose2)
        #deOffsetPose2 = juncOrigin1.doInverse(localJunctionPose2)
        #resultOffset = newJuncOrigin1.convertGlobalPoseToLocal(alignedJunctionPose2)
        
        print "alignedJunctionPose2:", alignedJunctionPose2
        print "localJunctionPose2:", localJunctionPose2
        
        alignJuncOrigin2 = Pose(alignedJunctionPose2)
        
        newPathSeg2 = []
        for k in range(len(pathSeg2)):
            p = pathSeg2[k]
            p2 = alignJuncOrigin2.convertGlobalToLocal(p)            
            #p2 = juncOrigin2.convertGlobalToLocal(p)            
            p1 = juncOrigin1.convertLocalToGlobal(p2)
            newPathSeg2.append(p1)
        
        
        if minK_1 < minK_2:
            stitchedPath = pathSeg1 + newPathSeg2
        else:
            stitchedPath = newPathSeg2 + pathSeg1
        
        
        if plotIter:
            
            pylab.clf()
            
            xP = []
            yP = []
            for p in stitchedPath:
                xP.append(p[0])
                yP.append(p[1])
                
            pylab.plot(xP,yP, color='k', zorder=10)
            
            xP = [junctionPose1[0], junctionPose2[0]]
            yP = [junctionPose1[1], junctionPose2[1]]
        
            pylab.scatter(xP,yP, color='k', zorder=10)
        
            xP = []
            yP = []
            for p in globalPath1:
                xP.append(p[0])
                yP.append(p[1])
            pylab.plot(xP,yP,color='r')

            xP = []
            yP = []
            for p in globalPath2:
                xP.append(p[0])
                yP.append(p[1])
            pylab.plot(xP,yP,color='b')


            xP = []
            yP = []
            for p in parentPath:
                xP.append(p[0])
                yP.append(p[1])
            pylab.plot(xP,yP,color='g')
        
        
            pylab.savefig("siblingCompare_%08u.png" % self.pathPlotCount2)
        
            self.pathPlotCount2 += 1
        
        #return


        " compute the medial axis for each pose "
        globalPath1 = stitchedPath
        globalPath2 = parentPath


        
        globalPath1Reverse = deepcopy(globalPath1)
        globalPath1Reverse.reverse()
        
        globalPath2Reverse = deepcopy(globalPath2)
        globalPath2Reverse.reverse()

        #globalMedial = globalPath1
        
        #medialSpline1 = SplineFit(globalMedial, smooth=0.1)
        globalSpline1 = SplineFit(globalPath1, smooth=0.1)
        globalSpline1Reverse = SplineFit(globalPath1Reverse, smooth=0.1)
        globalSpline2 = SplineFit(globalPath2, smooth=0.1)
        globalSpline2Reverse = SplineFit(globalPath2Reverse, smooth=0.1)


        overlapMatch = []
        angleSum1 = 0.0
        angleSum2 = 0.0
        for i in range(0,len(globalPath1)):
            p_1 = globalPath1[i]
            p_2, j, minDist = gen_icp.findClosestPointInA(globalPath2, p_1)
            if minDist < 0.5:
                overlapMatch.append((i,j,minDist))

                pathU1 = globalSpline1.findU(p_1)    
                pathU2 = globalSpline2.findU(p_2)    
                pathU2_R = globalSpline2Reverse.findU(p_2)    

                pathVec1 = globalSpline1.getUVector(pathU1)
                pathVec2 = globalSpline2.getUVector(pathU2)
                pathVec2_R = globalSpline2Reverse.getUVector(pathU2_R)

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
        #print i, "angleSum1 =", angleSum1, "angleSum2 =", angleSum2
        if angleSum1 > angleSum2:
            orientedGlobalPath = globalPath2Reverse
        else:
            orientedGlobalPath = globalPath2
            
        orientedGlobalSpline2 = SplineFit(orientedGlobalPath, smooth=0.1)


        globalSamples2 = orientedGlobalSpline2.getUniformSamples(spacing = 0.04)
        globalSamples1 = globalSpline1.getUniformSamples(spacing = 0.04)

        
        globalVar = []
        medialVar = []
        
        " compute the local variance of the angle "
        VAR_WIDTH = 40
        for i in range(len(globalSamples2)):
            
            lowK = i - VAR_WIDTH/2
            if lowK < 0:
                lowK = 0
                
            highK = i + VAR_WIDTH/2
            if highK >= len(globalSamples2):
                highK = len(globalSamples2)-1
            
            localSamp = []
            for k in range(lowK, highK+1):
                localSamp.append(globalSamples2[k][2])
            
            sum = 0.0
            for val in localSamp:
                sum += val
                
            meanSamp = sum / float(len(localSamp))
            

            sum = 0
            for k in range(len(localSamp)):
                sum += (localSamp[k] - meanSamp)*(localSamp[k] - meanSamp)
        
            varSamp = sum / float(len(localSamp))
            
            globalVar.append((meanSamp, varSamp))        

        for i in range(len(globalSamples1)):
            
            lowK = i - VAR_WIDTH/2
            if lowK < 0:
                lowK = 0
                
            highK = i + VAR_WIDTH/2
            if highK >= len(globalSamples1):
                highK = len(globalSamples1)-1
            
            localSamp = []
            for k in range(lowK, highK+1):
                localSamp.append(globalSamples1[k][2])
            
            sum = 0.0
            for val in localSamp:
                sum += val
                
            meanSamp = sum / float(len(localSamp))
            

            sum = 0
            for k in range(len(localSamp)):
                sum += (localSamp[k] - meanSamp)*(localSamp[k] - meanSamp)
        
            varSamp = sum / float(len(localSamp))
            
            medialVar.append((meanSamp, varSamp))        


        " now lets find closest points and save their local variances "            
        closestPairs = []
        #TERM_DIST = 0
        TERM_DIST = 20
        
        TERM_DIST1 = len(globalSamples1)/4
        TERM_DIST2 = len(globalSamples2)/4
        
        
        for i in range(TERM_DIST2, len(globalSamples2)-TERM_DIST2):
            pG = globalSamples2[i]
            minDist = 1e100
            minJ = -1
            for j in range(TERM_DIST1, len(globalSamples1)-TERM_DIST1):
                pM = globalSamples1[j]
                dist = math.sqrt((pG[0]-pM[0])**2 + (pG[1]-pM[1])**2)
                
                if dist < minDist:
                    minDist = dist
                    minJ = j
            
            juncDist = math.sqrt((junctionPose1[0]-pM[0])**2 + (junctionPose1[1]-pM[1])**2)
            #if minDist < 0.1:
            if True:
                closestPairs.append((i,minJ,minDist,globalVar[i][0],medialVar[minJ][0],globalVar[i][1],medialVar[minJ][1], juncDist))

        for j in range(TERM_DIST1, len(globalSamples1)-TERM_DIST1):
            pM = globalSamples1[j]
            minDist = 1e100
            minI = -1
            for i in range(TERM_DIST2, len(globalSamples2)-TERM_DIST2):
                pG = globalSamples2[i]
                dist = math.sqrt((pG[0]-pM[0])**2 + (pG[1]-pM[1])**2)
                
                if dist < minDist:
                    minDist = dist
                    minI = i
                    
            juncDist = math.sqrt((junctionPose1[0]-pM[0])**2 + (junctionPose1[1]-pM[1])**2)
            #if minDist < 0.1:
            if True:
                closestPairs.append((minI,j,minDist,globalVar[minI][0],medialVar[j][0],globalVar[minI][1],medialVar[j][1],juncDist))

        " remove duplicates "
        closestPairs = list(set(closestPairs))
        #closestPairs = list(closestPairs)
        
        " sort by lowest angular variance"
        closestPairs = sorted(closestPairs, key=itemgetter(2,5,6,7))
        #s = sorted(student_objects, key=attrgetter('age'))     # sort on secondary key
        #sorted(s, key=attrgetter('grade'), reverse=True)  
        #closestPairs.sort()
        
        #print len(closestPairs), "closest pairs"
        #for val in closestPairs:
        #    print val

        pathDict1 = self.getPath(pathID1)
        pathDict2 = self.getPath(pathID2)
        
        parentID1 = pathDict1["parentID"]
        parentID2 = pathDict2["parentID"]
        

        if False:
            globalJunctionPose1 = pathDict1["globalJunctionPose"]
            globalJunctionPose2 = pathDict2["globalJunctionPose"]


            originU2 = globalSpline1.findU(globalJunctionPose1)    
            originU1 = orientedGlobalSpline2.findU(globalJunctionPose1)
        
        elif len(closestPairs) > 0:

            globalJunctionPose1 = pathDict1["globalJunctionPose"]
            globalJunctionPose2 = pathDict2["globalJunctionPose"]            

            originU2 = globalSpline1.findU(globalJunctionPose1)    
            originU1 = orientedGlobalSpline2.findU(globalJunctionPose1)
            isFound = False
                        
            for pair in closestPairs:
                juncDist = pair[7]
                
                if juncDist < 1.0:
            
                    originU2 = globalSpline1.findU(globalSamples1[pair[1]])    
                    originU1 = orientedGlobalSpline2.findU(globalSamples2[pair[0]])
                    isFound = True
                    break
                
            if not isFound:
                print closestPairs
                raise
        else:
            raise
        
        u2 = originU2
        u1 = originU1
        angGuess = 0.0
        

        resultPose1, lastCost1, matchCount1 = gen_icp.pathOverlapICP([u1,u2,angGuess], orientedGlobalPath, globalPath1, plotIter = plotIter, n1 = pathID1, n2 = pathID2)

        print "matchCount,cost,result:", matchCount1, lastCost1, resultPose1

        globalPath1 = self.paths[pathID1]
        globalPath2 = self.paths[pathID2]
        parentPath = self.paths[parentPathID]

        if plotIter:
            
            pylab.clf()
            
            #xP = []
            #yP = []
            #for p in stitchedPath:
            #    xP.append(p[0])
            #    yP.append(p[1])
                
            #pylab.plot(xP,yP, color='k', zorder=10)
            
            xP = [junctionPose1[0], junctionPose2[0]]
            yP = [junctionPose1[1], junctionPose2[1]]
        
            pylab.scatter(xP,yP, color='k', zorder=10)

            poseOrigin = Pose(resultPose1)
            
            xP = []
            yP = []
            for p in stitchedPath:
                result = poseOrigin.convertLocalToGlobal(p)
                xP.append(result[0])
                yP.append(result[1])

            pylab.plot(xP,yP, color='k', zorder=10)
               
            xP = []
            yP = []
            for p in globalPath1:
                xP.append(p[0])
                yP.append(p[1])
            pylab.plot(xP,yP,color='r')

            xP = []
            yP = []
            for p in globalPath2:
                xP.append(p[0])
                yP.append(p[1])
            pylab.plot(xP,yP,color='b')


            xP = []
            yP = []
            for p in parentPath:
                xP.append(p[0])
                yP.append(p[1])
            pylab.plot(xP,yP,color='g')
        
        
            pylab.savefig("siblingCompare_%08u.png" % self.pathPlotCount2)
        
            self.pathPlotCount2 += 1
        

        
        #resultPose1 + 
        juncOrigin1 = Pose(junctionPose1)
        juncOrigin2 = Pose(junctionPose2)
        
        print "junctionPose1:", junctionPose1
        print "junctionPose2:", junctionPose2

        #resultPose1 = [resultPose1[0],resultPose1[1],0.0]
        print "resultPose1:", resultPose1
        
        
        " align the two junctions "
        alignedJunctionPose2 = [junctionPose2[0], junctionPose2[1], junctionPose1[2]]
        #localJunctionPose2 = juncOrigin1.convertGlobalPoseToLocal(alignedJunctionPose2)
        #deOffsetPose2 = juncOrigin1.doInverse(localJunctionPose2)
        
        print "alignedJunctionPose2:", alignedJunctionPose2
        
        offsetOrigin1 = Pose(resultPose1)
        newJuncPose1 = offsetOrigin1.convertLocalOffsetToGlobal(junctionPose1)

        print "newJuncPose1:", newJuncPose1

        #newJuncPose1 = juncOrigin1.convertLocalOffsetToGlobal(resultPose1)
        newJuncOrigin1 = Pose([newJuncPose1[0], newJuncPose1[1], 0.0])
        #resultOffset = newJuncOrigin1.convertGlobalPoseToLocal(alignedJunctionPose2)
        resultOffset = newJuncOrigin1.convertGlobalPoseToLocal([junctionPose2[0],junctionPose2[1], 0.0])
        
        print "resultOffset:", resultOffset
        
        resultPose2 = newJuncOrigin1.doInverse(resultOffset)

        #resultPose2[2] = resultPose1[2]
        


        def computeOffset(pose1, pose2):
        
            " corner points and orientations "
            corner1Pose = pose1
            corner2Pose = pose2
            
            " convert the desired intersection point on curve 1 into global coordinates "
            poseProfile1 = Pose([0.0,0.0,0.0])
                
            " now convert this point into a pose, and perform the inverse transform using corner2Pose "
            desGlobalPose2 = Pose(corner1Pose)
            
            " perform inverse offset from the destination pose "
            negCurve2Pose = desGlobalPose2.doInverse(corner2Pose)
            
            " relative pose between pose 1 and pose 2 to make corners coincide and same angle "
            resPose2 = desGlobalPose2.convertLocalOffsetToGlobal(negCurve2Pose)
            localOffset = poseProfile1.convertGlobalPoseToLocal(resPose2)
            
            return [localOffset[0], localOffset[1], localOffset[2]]
        
        
        resultPose2 = computeOffset(newJuncPose1, alignedJunctionPose2)
        #resultPose2 = computeOffset(newJuncPose1, junctionPose2)
        #fooPose2 = computeOffset(junctionPose2, newJuncPose1)

        print "resultPose2:", resultPose2

        #fooPose1 = [0.0,0.0,0.0]
        " array([[cos(self.estPose[2]), sin(self.estPose[2])],[-sin(self.estPose[2]),cos(self.estPose[2])]])  "
        #fooPose1[0] = resultPose2[0]*cos(resultPose1[2]) - resultPose2[1]*sin(resultPose1[2])
        #fooPose1[1] = resultPose2[0]*sin(resultPose1[2]) + resultPose2[1]*cos(resultPose1[2])
        #fooPose1[2] = resultPose1[2]
        
        #print "fooPose1:", fooPose1
        
        #fooPose2 = [0.0,0.0,0.0]
        " array([[cos(self.estPose[2]), sin(self.estPose[2])],[-sin(self.estPose[2]),cos(self.estPose[2])]])  "
        #fooPose2[0] = resultPose2[0]*cos(resultPose1[2]) + resultPose2[1]*sin(resultPose1[2])
        #fooPose2[1] = -resultPose2[0]*sin(resultPose1[2]) + resultPose2[1]*cos(resultPose1[2])
        #fooPose2[2] = resultPose1[2]
        
        #print "fooPose2:", fooPose2

        """
        alignJuncOrigin2 = Pose(alignedJunctionPose2)
        
        newPathSeg2 = []
        for k in range(len(pathSeg2)):
            p = pathSeg2[k]
            p2 = alignJuncOrigin2.convertGlobalToLocal(p)            
            #p2 = juncOrigin2.convertGlobalToLocal(p)            
            p1 = juncOrigin1.convertLocalToGlobal(p2)
            newPathSeg2.append(p1)


        poseOrigin = Pose(resultPose1)
        path1 = pair[1]
        path1_offset = []
        for p in path1:
            result = poseOrigin.convertLocalToGlobal(p)
            path1_offset.append(result)    

        """

        globalPath1 = self.paths[pathID1]
        globalPath2 = self.paths[pathID2]
        parentPath = self.paths[parentPathID]

        offsetOrigin1 = Pose(resultPose1)
        offsetOrigin2 = Pose(resultPose2)

        offJuncPose1 = offsetOrigin1.convertLocalOffsetToGlobal(junctionPose1)
        offJuncPose2 = offsetOrigin2.convertLocalOffsetToGlobal(junctionPose2)

        juncAngDiff = diffAngle(offJuncPose1[2],offJuncPose2[2])

        if plotIter:
            
            pylab.clf()
            
            xP = [junctionPose1[0], junctionPose2[0]]
            yP = [junctionPose1[1], junctionPose2[1]]
        
            pylab.scatter(xP,yP, color='k', zorder=9)
        
            xP = []
            yP = []

            for p in globalPath1:
                p1 = offsetOrigin1.convertLocalToGlobal(p)
                xP.append(p1[0])
                yP.append(p1[1])
            pylab.plot(xP,yP,color='r', zorder=8)

            offJuncPoint1 = offsetOrigin1.convertLocalToGlobal([junctionPose1[0],junctionPose1[1]])
            pylab.scatter([offJuncPoint1[0]],[offJuncPoint1[1]], color='r', zorder=10)
            print "offJuncPoint1:", offJuncPoint1

            xP = []
            yP = []
            #offsetOrigin2 = Pose(fooPose1)
            #offsetOrigin2 = Pose(resultOffset)
            for p in globalPath2:
                p1 = offsetOrigin2.convertLocalToGlobal(p)
                xP.append(p1[0])
                yP.append(p1[1])
            pylab.plot(xP,yP,color='b', zorder = 8)
            
            offJuncPoint2 = offsetOrigin2.convertLocalToGlobal([junctionPose2[0],junctionPose2[1]])
            pylab.scatter([offJuncPoint2[0]],[offJuncPoint2[1]], color='b', zorder=10)
            
            
            #xP = []
            #yP = []
            #offsetOrigin2 = Pose(fooPose2)
            #for p in globalPath2:
            #    p1 = offsetOrigin2.convertLocalToGlobal(p)
            #    xP.append(p1[0])
            #    yP.append(p1[1])
            #pylab.plot(xP,yP,color='m', zorder = 8)

            #offJuncPoint2 = offsetOrigin2.convertLocalToGlobal([junctionPose2[0],junctionPose2[1]])
            #pylab.scatter([offJuncPoint2[0]],[offJuncPoint2[1]], color='m', zorder=10)

            xP = []
            yP = []
            for p in parentPath:
                xP.append(p[0])
                yP.append(p[1])
            pylab.plot(xP,yP,color='g')
        
            pylab.title("resultPose = [%1.2f %1.2f %1.2f] [%1.2f %1.2f %1.2f] %1.2f\n" % (resultPose1[0], resultPose1[1], resultPose1[2], resultPose2[0], resultPose2[1], resultPose2[2], juncAngDiff))
        
            pylab.savefig("siblingCompare_%08u.png" % self.pathPlotCount2)
        
            self.pathPlotCount2 += 1
        

        


        return resultPose1, resultPose2, lastCost1, matchCount1, juncAngDiff
     



    def makePathCompare(self, globalPath1, globalPath2, pathID1, pathID2, plotIter = False):

        " compute the medial axis for each pose "
        
        #node1 = self.nodeHash[nodeID]
        #posture1 = node1.getStableGPACPosture()
        #hull1, medial1 = computeHullAxis(nodeID, node1, tailCutOff = True)
        #hull1, medial1 = computeHullAxis(nodeID, node1, tailCutOff = False)
        
        globalPath1Reverse = deepcopy(globalPath1)
        globalPath1Reverse.reverse()
        
        globalPath2Reverse = deepcopy(globalPath2)
        globalPath2Reverse.reverse()

        #globalMedial = globalPath1
        
        #medialSpline1 = SplineFit(globalMedial, smooth=0.1)
        globalSpline1 = SplineFit(globalPath1, smooth=0.1)
        globalSpline1Reverse = SplineFit(globalPath1Reverse, smooth=0.1)
        globalSpline2 = SplineFit(globalPath2, smooth=0.1)
        globalSpline2Reverse = SplineFit(globalPath2Reverse, smooth=0.1)


        overlapMatch = []
        angleSum1 = 0.0
        angleSum2 = 0.0
        for i in range(0,len(globalPath1)):
            p_1 = globalPath1[i]
            p_2, j, minDist = gen_icp.findClosestPointInA(globalPath2, p_1)
            if minDist < 0.5:
                overlapMatch.append((i,j,minDist))

                pathU1 = globalSpline1.findU(p_1)    
                pathU2 = globalSpline2.findU(p_2)    
                pathU2_R = globalSpline2Reverse.findU(p_2)    

                pathVec1 = globalSpline1.getUVector(pathU1)
                pathVec2 = globalSpline2.getUVector(pathU2)
                pathVec2_R = globalSpline2Reverse.getUVector(pathU2_R)

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
        #print i, "angleSum1 =", angleSum1, "angleSum2 =", angleSum2
        if angleSum1 > angleSum2:
            orientedGlobalPath = globalPath2Reverse
        else:
            orientedGlobalPath = globalPath2
            
        orientedGlobalSpline2 = SplineFit(orientedGlobalPath, smooth=0.1)


        globalSamples2 = orientedGlobalSpline2.getUniformSamples(spacing = 0.04)
        globalSamples1 = globalSpline1.getUniformSamples(spacing = 0.04)

        
        globalVar = []
        medialVar = []
        
        " compute the local variance of the angle "
        VAR_WIDTH = 40
        for i in range(len(globalSamples2)):
            
            lowK = i - VAR_WIDTH/2
            if lowK < 0:
                lowK = 0
                
            highK = i + VAR_WIDTH/2
            if highK >= len(globalSamples2):
                highK = len(globalSamples2)-1
            
            localSamp = []
            for k in range(lowK, highK+1):
                localSamp.append(globalSamples2[k][2])
            
            sum = 0.0
            for val in localSamp:
                sum += val
                
            meanSamp = sum / float(len(localSamp))
            

            sum = 0
            for k in range(len(localSamp)):
                sum += (localSamp[k] - meanSamp)*(localSamp[k] - meanSamp)
        
            varSamp = sum / float(len(localSamp))
            
            globalVar.append((meanSamp, varSamp))        

        for i in range(len(globalSamples1)):
            
            lowK = i - VAR_WIDTH/2
            if lowK < 0:
                lowK = 0
                
            highK = i + VAR_WIDTH/2
            if highK >= len(globalSamples1):
                highK = len(globalSamples1)-1
            
            localSamp = []
            for k in range(lowK, highK+1):
                localSamp.append(globalSamples1[k][2])
            
            sum = 0.0
            for val in localSamp:
                sum += val
                
            meanSamp = sum / float(len(localSamp))
            

            sum = 0
            for k in range(len(localSamp)):
                sum += (localSamp[k] - meanSamp)*(localSamp[k] - meanSamp)
        
            varSamp = sum / float(len(localSamp))
            
            medialVar.append((meanSamp, varSamp))        


        " now lets find closest points and save their local variances "            
        closestPairs = []
        #TERM_DIST = 0
        TERM_DIST = 20
        
        TERM_DIST1 = len(globalSamples1)/4
        TERM_DIST2 = len(globalSamples2)/4
        
        
        for i in range(TERM_DIST2, len(globalSamples2)-TERM_DIST2):
            pG = globalSamples2[i]
            minDist = 1e100
            minJ = -1
            for j in range(TERM_DIST1, len(globalSamples1)-TERM_DIST1):
                pM = globalSamples1[j]
                dist = math.sqrt((pG[0]-pM[0])**2 + (pG[1]-pM[1])**2)
                
                if dist < minDist:
                    minDist = dist
                    minJ = j
                    
            #if minDist < 0.1:
            if True:
                closestPairs.append((i,minJ,minDist,globalVar[i][0],medialVar[minJ][0],globalVar[i][1],medialVar[minJ][1]))

        for j in range(TERM_DIST1, len(globalSamples1)-TERM_DIST1):
            pM = globalSamples1[j]
            minDist = 1e100
            minI = -1
            for i in range(TERM_DIST2, len(globalSamples2)-TERM_DIST2):
                pG = globalSamples2[i]
                dist = math.sqrt((pG[0]-pM[0])**2 + (pG[1]-pM[1])**2)
                
                if dist < minDist:
                    minDist = dist
                    minI = i
                    
            #if minDist < 0.1:
            if True:
                closestPairs.append((minI,j,minDist,globalVar[minI][0],medialVar[j][0],globalVar[minI][1],medialVar[j][1]))

        " remove duplicates "
        closestPairs = list(set(closestPairs))
        #closestPairs = list(closestPairs)
        
        " sort by lowest angular variance"
        closestPairs = sorted(closestPairs, key=itemgetter(2,5,6))
        #s = sorted(student_objects, key=attrgetter('age'))     # sort on secondary key
        #sorted(s, key=attrgetter('grade'), reverse=True)  
        #closestPairs.sort()
        
        #print len(closestPairs), "closest pairs"
        #for val in closestPairs:
        #    print val

        pathDict1 = self.getPath(pathID1)
        pathDict2 = self.getPath(pathID2)
        
        parentID1 = pathDict1["parentID"]
        parentID2 = pathDict2["parentID"]
        
        if parentID1 == parentID2:
            
            """
            branchNodeID1 = pathDict1["branchNodeID"]
            branchNodeID2 = pathDict2["branchNodeID"]

            localJunctionPose1 = pathDict1["localJunctionPose"]
            localJunctionPose2 = pathDict2["localJunctionPose"]
            
            estPose1 = self.nodeHash[branchNodeID1].getEstPose()
            estPose2 = self.nodeHash[branchNodeID2].getEstPose()
            
            poseOrigin1 = Pose(estPose1)
            poseOrigin2 = Pose(estPose2)
            
            globalJunctionPose1 = poseOrigin1.convertLocalOffsetToGlobal(localJunctionPose1)
            globalJunctionPose2 = poseOrigin2.convertLocalOffsetToGlobal(localJunctionPose2)
            """


            globalJunctionPose1 = pathDict1["globalJunctionPose"]
            globalJunctionPose2 = pathDict2["globalJunctionPose"]


            originU2 = globalSpline1.findU(globalJunctionPose1)    
            originU1 = orientedGlobalSpline2.findU(globalJunctionPose2)
        
        elif len(closestPairs) > 0:
            originU2 = globalSpline1.findU(globalSamples1[closestPairs[0][1]])    
            originU1 = orientedGlobalSpline2.findU(globalSamples2[closestPairs[0][0]])

        else:
            raise
        
        u2 = originU2
        u1 = originU1
        angGuess = 0.0
        

        #resultPose2, lastCost2, matchCount2 = gen_icp.globalPathToNodeOverlapICP2([u2, uMedialOrigin2, 0.0], orientedSplicePath, medial2, plotIter = True, n1 = nodeID-1, n2 = -1, arcLimit = 0.1)
        #resultPose2, lastCost2, matchCount2 = gen_icp.globalOverlapICP_GPU2([u2, uMedialOrigin2, 0.0], orientedSplicePath, medial2, plotIter = True, n1 = nodeID-1, n2 = -1, arcLimit = 0.1)

        
        resultPose1, lastCost1, matchCount1 = gen_icp.pathOverlapICP([u1,u2,angGuess], orientedGlobalPath, globalPath1, plotIter = plotIter, n1 = pathID1, n2 = pathID2)
        #resultPose2, lastCost2, matchCount2 = gen_icp.pathOverlapICP([u1,u2+0.1,angGuess], orientedGlobalPath, globalPath1, plotIter = True, n1 = pathID1, n2 = pathID2)
        #resultPose3, lastCost3, matchCount3 = gen_icp.pathOverlapICP([u1,u2-0.1,angGuess], orientedGlobalPath, globalPath1, plotIter = True, n1 = pathID1, n2 = pathID2)
        #resultPose1, lastCost1, matchCount1 = gen_icp.pathOverlapICP([u1,u2,angGuess], orientedGlobalPath, globalPath1, plotIter = False, n1 = pathID1, n2 = pathID2)
        #resultPose2, lastCost2, matchCount2 = gen_icp.pathOverlapICP([u1,u2+0.1,angGuess], orientedGlobalPath, globalPath1, plotIter = False, n1 = pathID1, n2 = pathID2)
        #resultPose3, lastCost3, matchCount3 = gen_icp.pathOverlapICP([u1,u2-0.1,angGuess], orientedGlobalPath, globalPath1, plotIter = False, n1 = pathID1, n2 = pathID2)

        print "matchCount,cost,result:", matchCount1, lastCost1, resultPose1

        #print "matchCounts:", matchCount1, matchCount2, matchCount3
        #print "costs:", lastCost1, lastCost2, lastCost3
        #print "results:", resultPose1, resultPose2, resultPose3


        return resultPose1, lastCost1, matchCount1
     
     
    def delNode(self, nodeID, pathID):
        print "deleting node", nodeID, "from path", pathID
        self.pathClasses[pathID]["nodeSet"].remove(nodeID)
        
        
     
    def addNode(self, nodeID, pathID):
        
        print "adding node", nodeID, "to path", pathID
        self.pathClasses[pathID]["nodeSet"].append(nodeID)

        #self.generatePaths()
        #self.trimPaths(self.paths)        
        #self.comparePaths()
       
        """
        probDist = {}        
        for k in range(self.numPaths):
            probDist[k] = 0.0
        
        " distribute remainder evenly "
        if self.numPaths > 1:
            remainder = 1.0 - probA
            val = remainder/(self.numPaths-1)
            for k in range(self.numPaths):
                probDist[k] = val

        probDist[pathID] = probA
        
        
        self.pathProbs[nodeID] = probDist
        """
        
    def delPath(self, pathID, mergeTargetID):
        
        print "deleting path", pathID, "merging to path", mergeTargetID
        try: 
            del self.pathClasses[pathID]
            del self.pathTermsVisited[pathID]
            
            self.mergeMapping[pathID] = mergeTargetID

            
            keys = []
            
            for key, val in self.consistency.iteritems():
                
                if key[0] == pathID or key[1] == pathID:
                    keys.append(key)
            
            for key in keys:        
                del self.consistency[key]
                
            
            for k, pathClass in self.pathClasses.iteritems():
                
                if pathClass["parentID"] == pathID:
                    pathClass["parentID"] = mergeTargetID
                    
        except:
            pass
    
    def addPath(self, parentID, branchNodeID, localJunctionPose):

        print "addPath(", parentID, branchNodeID, localJunctionPose
        
        nodePose = self.nodeHash[branchNodeID].getEstPose()
        nodeOrigin = Pose(nodePose)
        globalJunctionPose = nodeOrigin.convertLocalOffsetToGlobal(localJunctionPose)
        
        
        oldPaths = self.getPathIDs()
        
        newPathID = self.pathIDs
        self.pathClasses[newPathID] = {"parentID" : parentID, "branchNodeID" : branchNodeID, "localJunctionPose" : localJunctionPose, 
                            "sameProb" : {}, "nodeSet" : [], "globalJunctionPose" : globalJunctionPose }        
                
        self.pathTermsVisited[newPathID] = False
        
        self.pathIDs += 1


        self.mergeMapping[newPathID] = None
        
        for oldID in oldPaths:
            
            path1 = self.getPath(oldID)
            path2 = self.getPath(newPathID)
            
            parentID1 = path1["parentID"]
            parentID2 = path2["parentID"]
            
            if parentID1 != parentID2:
                self.consistency[(oldID,newPathID)] = 0.0
            else:
                self.consistency[(oldID,newPathID)] = 0.5
        
        #self.generatePaths()
        #self.trimPaths(self.paths)        
        #self.comparePaths()
        print "consistency:", self.consistency
        
        return newPathID
 
    def getTopology(self, pathID):
        
        random.seed(0)        
 
        nodes = self.getNodes(pathID)

        junctionNodeID = self.pathClasses[pathID]["branchNodeID"]
        globalJunctionPoint = None

        if junctionNodeID != None:
        
            #localJunctionPoint = self.pathClasses[pathID]["localJunctionPose"]
            
            print "junctionNodeID:", junctionNodeID
            print "nodes:", self.nodeHash.keys()
            
            #poseOrigin = Pose(self.nodeHash[junctionNodeID].getEstPose())
            #globalJunctionPoint = poseOrigin.convertLocalOffsetToGlobal(localJunctionPoint)
            globalJunctionPoint = self.getGlobalJunctionPose(pathID) 
        
        def convertAlphaUniform(a_vert, max_spacing = 0.04):
            
            " make the vertices uniformly distributed "
            
            new_vert = []
        
            for i in range(len(a_vert)):
                p0 = a_vert[i]
                p1 = a_vert[(i+1) % len(a_vert)]
                dist = math.sqrt((p0[0]-p1[0])**2 + (p0[1]-p1[1])**2)
    
                vec = [p1[0]-p0[0], p1[1]-p0[1]]
                
                if dist == 0:
                    pass
                else:
                                        
                    vec[0] /= dist
                    vec[1] /= dist
                    
                    new_vert.append(copy(p0))
                    
                    if dist > max_spacing:
                        " cut into pieces max_spacing length or less "
                        numCount = int(math.floor(dist / max_spacing))
                        
                        for j in range(1, numCount+1):
                            newP = [j*max_spacing*vec[0] + p0[0], j*max_spacing*vec[1] + p0[1]]
                            new_vert.append(newP)
            
            return new_vert            

        medialPointSoup = []


        if nodes == []:
            return [], []


        if True:    
            for nodeID in nodes:
                estPose1 = self.nodeHash[nodeID].getGlobalGPACPose()        
        
                if self.nodeHash[nodeID].isBowtie:            
                    hull1 = computeBareHull(self.nodeHash[nodeID], sweep = False, static = True)
                else:
                    hull1 = computeBareHull(self.nodeHash[nodeID], sweep = False)
        
        
                m = hashlib.md5()
                m.update(repr(hull1))
                print nodeID, "hull1 =", int(m.digest().encode('hex'),16)
        
                " set the origin of pose 1 "
                poseOrigin = Pose(estPose1)
        
                xP = []
                yP = []    
                for k in range(len(hull1)):
     
                    p = hull1[k]
                    
                    m = hashlib.md5()
                    m.update(repr(p))
                    #print nodeID, "p =", int(m.digest().encode('hex'),16)               
                    
                    p1 = poseOrigin.convertLocalToGlobal(p)
                    
                    m = hashlib.md5()
                    m.update(repr(p1))
                    #print nodeID, "p1 =", int(m.digest().encode('hex'),16)
                    
                    medialPointSoup.append(p1)

        if False:
            pylab.clf()    
    
            for nodeID in nodes:
                estPose1 = self.nodeHash[nodeID].getGlobalGPACPose()        
        
                if self.nodeHash[nodeID].isBowtie:            
                    hull1 = computeBareHull(self.nodeHash[nodeID], sweep = False, static = True)
                else:
                    hull1 = computeBareHull(self.nodeHash[nodeID], sweep = False)
        
        
                m = hashlib.md5()
                m.update(repr(hull1))
                print nodeID, "hull1 =", int(m.digest().encode('hex'),16)
        
                " set the origin of pose 1 "
                poseOrigin = Pose(estPose1)
        
                xP = []
                yP = []    
                for k in range(len(hull1)):
     
                    p = hull1[k]
                    
                    m = hashlib.md5()
                    m.update(repr(p))
                    #print nodeID, "p =", int(m.digest().encode('hex'),16)               
                    
                    p1 = poseOrigin.convertLocalToGlobal(p)
                    
                    m = hashlib.md5()
                    m.update(repr(p1))
                    #print nodeID, "p1 =", int(m.digest().encode('hex'),16)
                    
                    medialPointSoup.append(p1)
                    xP.append(p1[0])
                    yP.append(p1[1])
    
                pylab.scatter(xP,yP)
    
            #pylab.xlim(-4,4)
            #pylab.ylim(-4,4)
            pylab.title("%d: %s" % (pathID, repr(nodes)))
            pylab.savefig("plotMedialSoup_%04u.png" % (self.medialSoupCount ))

            saveFile = "medialPointSoup = ["
            for k in range(len(medialPointSoup)):
                saveFile += repr(medialPointSoup[k]) + "\n"
            saveFile += "]\n"
    
            f = open("medialSoupSave_%04u.txt" % self.medialSoupCount, 'w')
            f.write(saveFile)
            f.close()
            self.medialSoupCount  += 1
        
        #nodeID = self.numNodes - 1
            
        #radius = 0.2
        radius = 0.2

        numPoints = len(medialPointSoup)

        isDone = False
        
        while not isDone:
    
            perturbPoints = []
            
            for p in medialPointSoup:
                p2 = copy(p)
                " add a little bit of noise to avoid degenerate conditions in CGAL "
                p2[0] += random.gauss(0.0,0.000001)
                p2[1] += random.gauss(0.0,0.000001)
    
                perturbPoints.append(p2)
        
            try:            
    
                vertices = alphamod.doAlpha(radius,perturbPoints)
                numVert = len(vertices)
                
                if numVert <= 2:
                    print "Failed, hull had only", numVert, "vertices"
                    raise
                
                isDone = True
            except:
                print "hull has holes!  retrying..."
                #print sArr     
        
        """       
        
        inputStr = str(numPoints) + " "

        " alpha shape circle radius "
        inputStr += str(radius) + " "

        #for p in medialPointSoup:
        #    p2 = copy(p)
        #    " add a little bit of noise to avoid degenerate conditions in CGAL "
        #    inputStr += str(p2[0]) + " " + str(p2[1]) + " "
        #
        #inputStr += "\n"


        isDone = False
        
        while not isDone:

            inputStr = str(numPoints) + " "
    
            " alpha shape circle radius "
            inputStr += str(radius) + " "
            
            for p in medialPointSoup:
                p2 = copy(p)
                " add a little bit of noise to avoid degenerate conditions in CGAL "
                p2[0] += random.gauss(0.0,0.000001)
                p2[1] += random.gauss(0.0,0.000001)
    
                inputStr += str(p2[0]) + " " + str(p2[1]) + " "
            
            inputStr += "\n"
            
            try:            
                " start the subprocess "
                if sys.platform == "win32":
                    subProc = Popen(["alpha2.exe"], stdin=PIPE, stdout=PIPE)
                else:
                    subProc = Popen(["./alpha2"], stdin=PIPE, stdout=PIPE)
                    
                
                " send input and receive output "
                sout, serr = subProc.communicate(inputStr)
        
                " convert string output to typed data "
                sArr = sout.split(" ")
                
                numVert = int(sArr[0])
                
                sArr = sArr[1:]
                print "hull has", numVert, "vertices"
                if numVert <= 2:
                    print "Failed, hull had only", numVert, "vertices"
                    raise
                
                vertices = []
                for i in range(len(sArr)/2):
                    vertices.append([float(sArr[2*i]), float(sArr[2*i + 1])])
                isDone = True
            except:
                print "hull has holes!  retrying..."
                #print sArr    

        """
        
        " cut out the repeat vertex "
        vertices = vertices[:-1]
        
        vertices = convertAlphaUniform(vertices)

        vertices.append(vertices[0])
        
        minX = 1e100
        maxX = -1e100
        minY = 1e100
        maxY = -1e100
        for p in vertices:
            if p[0] > maxX:
                maxX = p[0]
            if p[0] < minX:
                minX = p[0]
            if p[1] > maxY:
                maxY = p[1]
            if p[1] < minY:
                minY = p[1]
    
    
        " SPECIFY SIZE OF GRID AND CONVERSION PARAMETERS "
        PIXELSIZE = 0.05
        mapSize = 2*max(max(maxX,math.fabs(minX)),max(maxY,math.fabs(minY))) + 1
        pixelSize = PIXELSIZE
        numPixel = int(2.0*mapSize / pixelSize + 1.0)
        divPix = math.floor((2.0*mapSize/pixelSize)/mapSize)
        
        def realToGrid(point):
            indexX = int(math.floor(point[0]*divPix)) + numPixel/2 + 1
            indexY = int(math.floor(point[1]*divPix)) + numPixel/2 + 1
            return indexX, indexY
    
        def gridToReal(indices):
            i = indices[0]
            j = indices[1]
            point = [(i - numPixel/2 - 1)*(pixelSize/2.0) + pixelSize/2.0, (j - numPixel/2 - 1)*(pixelSize/2.0) + pixelSize/2.0]
            return point
    
        " CONVERT HULL TO GRID COORDINATES "
        gridHull = []
        for i in range(len(vertices)):
            p = vertices[i]
            gridHull.append(realToGrid(p))

        " BOUNDING BOX OF GRID HULL "    
        minX = 1e100
        maxX = -1e100
        minY = 1e100
        maxY = -1e100
        for p in gridHull:
            if p[0] > maxX:
                maxX = p[0]
            if p[0] < minX:
                minX = p[0]
            if p[1] > maxY:
                maxY = p[1]
            if p[1] < minY:
                minY = p[1]

        " COMPUTE MEDIAL AXIS OF HULL "
        resultImg = Image.new('L', (numPixel,numPixel))
        resultImg = computeMedialAxis(self.medialCount, numPixel,numPixel, 5, resultImg, len(gridHull[:-2]), gridHull[:-2])
        self.medialCount += 1

        " EXTRACT POINTS FROM GRID TO LIST "
        imgA = resultImg.load()    
        points = []
        for i in range(1, numPixel-1):
            for j in range(1, numPixel-1):
                if imgA[i,j] == 0:
                    points.append((i,j))
    
        " CREATE GRAPH NODES FOR EACH POINT "
        medialGraph = graph.graph()
        for p in points:
            medialGraph.add_node(p, [])
    
        " ADD EDGES BETWEEN NEIGHBORS "
        builtGraph = {}
        for i in range(2, numPixel-2):
            for j in range(2, numPixel-2):
                if imgA[i,j] == 0:
                    builtGraph[(i,j)] = []
                    for k in range(i-1, i+2):
                        for l in range(j-1, j+2):
                            if imgA[k,l] == 0 and not (k == i and l == j):
                                builtGraph[(i,j)].append((k,l))
                                medialGraph.add_edge((i,j), (k,l))
                                
        " COMPUTE MINIMUM SPANNING TREE "
        mst = medialGraph.minimal_spanning_tree()
        
        " INITIALIZE DATA DICTS FOR UNIDIRECTIONAL MST"
        uni_mst = {}
        for k, v in mst.items():
            uni_mst[k] = []

        
        " ADD EDGES TO DICT TREE REPRESENTATION "
        for k, v in mst.items():
            if v != None:
                uni_mst[k].append(v)
                uni_mst[v].append(k)

        " LOCATE ALL LEAVES "
        leaves = []
        for k, v in uni_mst.items():
            if len(v) == 1:
                leaves.append(k)
        
        " DELETE ALL NODES THAT ARE LEAVES, TO REMOVE SINGLE NODE BRANCHES "
        for leaf in leaves:
            medialGraph.del_node(leaf)    
            
        " RECOMPUTE MST "
        mst = medialGraph.minimal_spanning_tree()

        " AGAIN, CREATE OUR DATA STRUCTURES AND IDENTIFY LEAVES "
        uni_mst = {}
        isVisited = {}
        nodeSum = {}
        for k, v in mst.items():
            uni_mst[k] = []
            isVisited[k] = 0
            nodeSum[k] = 0
        
        for k, v in mst.items():
            if v != None:
                uni_mst[k].append(v)
                uni_mst[v].append(k)
                        
        " RECORD THE LEAVES AND JUNCTIONS "
        leaves = []
        junctions = []
        for k, v in uni_mst.items():
            if len(v) == 1:
                leaves.append(k)

            if len(v) > 2:
                junctions.append(k)

        print "junctionNodeID:", junctionNodeID
        print "junctions:", junctions
        
        " SAVE FOR LATER, IDENTIFIED BY THEIR INDEX NOW "
        numLeafs = len(leaves)
        allLeafs = []
        allJunctions = []
        for leaf in leaves:
            allLeafs.append(gridToReal(leaf))

        if junctionNodeID != None:
        
            print "finding theoretical junction:"
            minKey = None
            minCand = None
            minJuncDist = 1e100
            for k, v in uni_mst.items():
                gCand = gridToReal(k)
                juncDist = sqrt((globalJunctionPoint[0]-gCand[0])**2 + (globalJunctionPoint[1]-gCand[1])**2)
                
                if juncDist < minJuncDist:
                    minKey = k
                    minJuncDist = juncDist
                    minCand = gCand
        
            theoryVec = [cos(globalJunctionPoint[2]), sin(globalJunctionPoint[2])]

            theoryJunc = (minKey, minCand, minJuncDist, theoryVec)
            print "theoryJunc:", theoryJunc
            
            theoryJuncPoint = theoryJunc[1]
            theoryLeaf = []
            leafMag = 0.02
            for i in range(5):
                newPoint = [theoryVec[0]*leafMag + theoryJuncPoint[0], theoryVec[1]*leafMag + theoryJuncPoint[1]]
                theoryLeaf.append(newPoint)
                leafMag += 0.02


            print "theoryLeaf:", theoryLeaf
            
            if len(junctions) > 0:
                
                for junc in junctions:
                    gJunc = gridToReal(junc)
                    juncDist = sqrt((globalJunctionPoint[0]-gJunc[0])**2 + (globalJunctionPoint[1]-gJunc[1])**2)
                    allJunctions.append((junc, gJunc,juncDist))

            else:
                
                print "adding theoretical junction:", (minKey, minCand, minJuncDist)
                allJunctions.append(theoryJunc)                       
        
        print "allJunctions:", allJunctions
        
        " FIND ALL THE PATHS BETWEEN LEAVES "        
        nodePaths = {}
        for leaf in leaves:
            
            " INITIALIZE DATA DICTS FOR VISITED, PATH SO FAR, AND NODE HOP COUNT"
            isVisited = {}
            nodeSum = {}
            nodePath = {}
            for k, v in uni_mst.items():
                isVisited[k] = 0
                nodeSum[k] = 0
                nodePath[k] = []
    
            " PERFORM DFS FOR THIS LEAF "
            getLongestPath(leaf, 0, [], uni_mst, isVisited, nodePath, nodeSum)
    
            " SAVE THE RESULTING PATH DATA STRUCTURE "
            nodePaths[leaf] = nodePath
            

        " FOR EVERY PAIR OF LEAVES, SAVE ITS PATH IF ITS LONG ENOUGH "
        " should have X choose 2 combinations"
        longPaths = []
        
        MAX_LEN = 2
        for leaf1 in leaves:
            for leaf2 in leaves:
                if leaf1 < leaf2:
                    if len(nodePaths[leaf1][leaf2]) > MAX_LEN:
                        nPath = deepcopy(nodePaths[leaf1][leaf2])
                        juncIndices = []
                        for junc in allJunctions:
                            try:
                                index = nPath.index(junc[0])
                            except:
                                juncIndices.append(None)
                            else:
                                juncIndices.append(index)
                                
                        longPaths.append((len(nPath),nPath, juncIndices))


        if junctionNodeID != None:
            theoryMedialLongPaths = []
            theoryPaths = []
            print "theoryPaths:"
            for leaf1 in leaves:
                nPath = deepcopy(nodePaths[leaf1][theoryJunc[0]])                    
                theoryPaths.append((len(nPath), nPath))
                print leaf1, theoryJunc, len(nPath)
    
            theoryPaths.sort(reverse=True)
    
            for k in range(len(theoryPaths)):
                path = theoryPaths[k][1]
                realPath = []
                for p in path:
                    realPath.append(gridToReal(p))
                #realPath.reverse()
                
                
                theoryPaths[k] = realPath + theoryLeaf
                
                print "theoryPath(" , k , ")", len(realPath), realPath[-1], theoryLeaf
                
                theoryPath = theoryPaths[k]
                print "theoryPath:", len(theoryPath)
                
                leafPath = deepcopy(theoryPath)
                
                frontVec = [0.,0.]
                backVec = [0.,0.]
                indic = range(3)
                indic.reverse()
                
                for i in indic:
                    if i+2 < len(leafPath):
                        p1 = leafPath[i+2]
                        p2 = leafPath[i]
                        vec = [p2[0]-p1[0], p2[1]-p1[1]]
                        frontVec[0] += vec[0]
                        frontVec[1] += vec[1]
            
                        p1 = leafPath[-i-3]
                        p2 = leafPath[-i-1]
                        vec = [p2[0]-p1[0], p2[1]-p1[1]]
                        backVec[0] += vec[0]
                        backVec[1] += vec[1]
            
                frontMag = math.sqrt(frontVec[0]*frontVec[0] + frontVec[1]*frontVec[1])
                backMag = math.sqrt(backVec[0]*backVec[0] + backVec[1]*backVec[1])
            
                frontVec[0] /= frontMag
                frontVec[1] /= frontMag
                backVec[0] /= backMag
                backVec[1] /= backMag
            
                newP1 = (leafPath[0][0] + frontVec[0]*10, leafPath[0][1] + frontVec[1]*10)
                newP2 = (leafPath[-1][0] + backVec[0]*10, leafPath[-1][1] + backVec[1]*10)
            
                leafPath.insert(0,newP1)
                leafPath.append(newP2)
                
                medial2 = deepcopy(leafPath)
        
                " take the long length segments at tips of medial axis"
                edge1 = medial2[0:2]
                edge2 = medial2[-2:]
                
                frontVec = [edge1[0][0]-edge1[1][0], edge1[0][1]-edge1[1][1]]
                backVec = [edge2[1][0]-edge2[0][0], edge2[1][1]-edge2[0][1]]
                frontMag = math.sqrt(frontVec[0]*frontVec[0] + frontVec[1]*frontVec[1])
                backMag = math.sqrt(backVec[0]*backVec[0] + backVec[1]*backVec[1])
                
                frontVec[0] /= frontMag
                frontVec[1] /= frontMag
                backVec[0] /= backMag
                backVec[1] /= backMag
                
                " make a smaller version of these edges "
                newP1 = (edge1[1][0] + frontVec[0]*2, edge1[1][1] + frontVec[1]*2)
                newP2 = (edge2[0][0] + backVec[0]*2, edge2[0][1] + backVec[1]*2)
        
                edge1 = [newP1, edge1[1]]
                edge2 = [edge2[0], newP2]

                " find the intersection points with the hull "
                hull = vertices
                interPoints = []
                for k in range(len(hull)-1):
                    hullEdge = [hull[k],hull[k+1]]
                    isIntersect1, point1 = Intersect(edge1, hullEdge)
                    if isIntersect1:
                        interPoints.append(point1)
                        break
        
                for k in range(len(hull)-1):
                    hullEdge = [hull[k],hull[k+1]]
                    isIntersect2, point2 = Intersect(edge2, hullEdge)
                    if isIntersect2:
                        interPoints.append(point2)
                        break
                
                " replace the extended edges with a termination point at the hull edge "            
                medial2 = medial2[1:-2]
                
                if isIntersect1:
                    medial2.insert(0, point1)
                if isIntersect2:
                    medial2.append(point2)
    
                theoryMedialLongPaths.append(deepcopy(medial2))




    


        
        print "longPaths:"
        for longPath in longPaths:
            print longPath[2]
        
        " SORT FOR THE LONGEST TO SHORTEST "
        longPaths.sort(reverse=True)
        
        " REMOVE SIZE FROM TUPLE  "
        juncIndices = []
        for k in range(len(longPaths)):
            juncIndices.append(longPaths[k][2])
            longPaths[k] = longPaths[k][1]
        
        print "juncIndices:", juncIndices
        
        " GET THE LEAF INDEXES TO EACH PATH "
        leafPairs = []
        for path in longPaths:
            leaf1 = path[0]
            leaf2 = path[-1]    
            
            leafID1 = leaves.index(leaf1)
            leafID2 = leaves.index(leaf2)

            leafPairs.append((leafID1,leafID2))
        
        print "leafPairs:", leafPairs

        for k in range(len(longPaths)):
            path = longPaths[k]
            realPath = []
            for p in path:
                realPath.append(gridToReal(p))
            
            longPaths[k] = realPath
                    

        #longMedialWidths = []
        medialLongPaths = []
        juncAngSet = []
        #medialTailCuts = []
        #for longPath in longPaths:
            
        print len(longPaths), "long paths"
        for n in range(len(longPaths)):
            
            longPath = longPaths[n]
            print "longPath:", len(longPath)
            
            leafPath = deepcopy(longPath)
            
            frontVec = [0.,0.]
            backVec = [0.,0.]
            indic = range(3)
            indic.reverse()
            
            for i in indic:
                if i+2 < len(leafPath):
                    p1 = leafPath[i+2]
                    p2 = leafPath[i]
                    vec = [p2[0]-p1[0], p2[1]-p1[1]]
                    frontVec[0] += vec[0]
                    frontVec[1] += vec[1]
        
                    p1 = leafPath[-i-3]
                    p2 = leafPath[-i-1]
                    vec = [p2[0]-p1[0], p2[1]-p1[1]]
                    backVec[0] += vec[0]
                    backVec[1] += vec[1]
        
            frontMag = math.sqrt(frontVec[0]*frontVec[0] + frontVec[1]*frontVec[1])
            backMag = math.sqrt(backVec[0]*backVec[0] + backVec[1]*backVec[1])
        
            frontVec[0] /= frontMag
            frontVec[1] /= frontMag
            backVec[0] /= backMag
            backVec[1] /= backMag
        
            newP1 = (leafPath[0][0] + frontVec[0]*10, leafPath[0][1] + frontVec[1]*10)
            newP2 = (leafPath[-1][0] + backVec[0]*10, leafPath[-1][1] + backVec[1]*10)
        
            leafPath.insert(0,newP1)
            leafPath.append(newP2)
            
            medial2 = deepcopy(leafPath)
    
            " take the long length segments at tips of medial axis"
            edge1 = medial2[0:2]
            edge2 = medial2[-2:]
            
            frontVec = [edge1[0][0]-edge1[1][0], edge1[0][1]-edge1[1][1]]
            backVec = [edge2[1][0]-edge2[0][0], edge2[1][1]-edge2[0][1]]
            frontMag = math.sqrt(frontVec[0]*frontVec[0] + frontVec[1]*frontVec[1])
            backMag = math.sqrt(backVec[0]*backVec[0] + backVec[1]*backVec[1])
            
            frontVec[0] /= frontMag
            frontVec[1] /= frontMag
            backVec[0] /= backMag
            backVec[1] /= backMag
            
            " make a smaller version of these edges "
            newP1 = (edge1[1][0] + frontVec[0]*2, edge1[1][1] + frontVec[1]*2)
            newP2 = (edge2[0][0] + backVec[0]*2, edge2[0][1] + backVec[1]*2)
    
            edge1 = [newP1, edge1[1]]
            edge2 = [edge2[0], newP2]

            " find the intersection points with the hull "
            hull = vertices
            interPoints = []
            for k in range(len(hull)-1):
                hullEdge = [hull[k],hull[k+1]]
                isIntersect1, point1 = Intersect(edge1, hullEdge)
                if isIntersect1:
                    interPoints.append(point1)
                    break
    
            for k in range(len(hull)-1):
                hullEdge = [hull[k],hull[k+1]]
                isIntersect2, point2 = Intersect(edge2, hullEdge)
                if isIntersect2:
                    interPoints.append(point2)
                    break
            
            " replace the extended edges with a termination point at the hull edge "            
            medial2 = medial2[1:-2]
            
            if isIntersect1:
                medial2.insert(0, point1)
            if isIntersect2:
                medial2.append(point2)

            medialLongPaths.append(deepcopy(medial2))

            print "globalJunctionPoint:", globalJunctionPoint
            print "juncIndices:", juncIndices[n]

            juncAngs = []
            if globalJunctionPoint != None:
                for juncInd in juncIndices[n]:
                    if juncInd != None:
                        frontVec = [0.,0.]
                        backVec = [0.,0.]
                        indic = range(3)
                        indic.reverse()
                        
                        print "len(longPath):", len(longPath)
                        print "juncInd:", juncInd

                        #highIndex = 2+juncInd+2
                        #highMod = 0
                        #if highIndex >= len(longPath):
                        #    highMod = len(longPath) - highIndex - 1

                        highIndex = juncInd+4
                        highMod = highIndex
                        if highIndex+4 >= len(longPath):
                            highMod = len(longPath) - 5
                        
                        # FOOBAR
                         
                        #lowIndex = -3-juncInd-3
                        #lowMod = 0
                        #if -lowIndex > len(longPath):
                        #    lowMod = -len(longPath) - lowIndex

                        lowIndex = juncInd-4
                        lowMod = lowIndex
                        if lowIndex-4 < 0:
                            lowMod = 4
                        
                        print "highIndex, highMod:", highIndex, highMod
                        print "lowIndex, lowMod:", lowIndex, lowMod
                        
                        
                        for i in indic:
                            p1 = longPath[i+highMod]
                            p2 = longPath[i+2+highMod]
                            #p1 = longPath[i+juncInd+highMod]
                            #p2 = longPath[i+juncInd+2+highMod]
                            vec = [p2[0]-p1[0], p2[1]-p1[1]]
                            frontVec[0] += vec[0]
                            frontVec[1] += vec[1]
                    
                            p1 = longPath[-i+lowMod]
                            p2 = longPath[-i+lowMod-2]
                            #p1 = longPath[-i-juncInd-1+lowMod]
                            #p2 = longPath[-i-juncInd-3+lowMod]
                            vec = [p2[0]-p1[0], p2[1]-p1[1]]
                            backVec[0] += vec[0]
                            backVec[1] += vec[1]
                        
                        frontMag = math.sqrt(frontVec[0]*frontVec[0] + frontVec[1]*frontVec[1])
                        backMag = math.sqrt(backVec[0]*backVec[0] + backVec[1]*backVec[1])
                    
                        frontVec[0] /= frontMag
                        frontVec[1] /= frontMag
                        backVec[0] /= backMag
                        backVec[1] /= backMag
                        
                        print "frontVec:", frontVec
                        print "backVec:", backVec
                        
                        foreAng = acos(frontVec[0])
                        if frontVec[1] < 0.0:
                            foreAng = -foreAng
        
                        backAng = acos(backVec[0])
                        if backVec[1] < 0.0:
                            backAng = -backAng
        
                        print "foreAng:", foreAng
                        print "backAng:", backAng
        
                        frontError = normalizeAngle(globalJunctionPoint[2]-foreAng)
                        backError = normalizeAngle(globalJunctionPoint[2]-backAng)
                        juncAngs.append((frontError,backError))
                    else:
                        juncAngs.append(None)
            juncAngSet.append(juncAngs)
        
        juncDists = []
        for junc in allJunctions:
            juncDists.append(junc[2])
        
        if True:
            pylab.clf()
    
            for path in medialLongPaths:
                xP = []
                yP = []
                for p in path:
                    #p2 = gridToReal(p)
                    xP.append(p[0])
                    yP.append(p[1])
    
                pylab.plot(xP,yP)
            #pylab.scatter(xP,yP, color='b')

    
            if junctionNodeID != None:
    
                for path in theoryMedialLongPaths:
                    xP = []
                    yP = []
                    for p in path:
                        xP.append(p[0])
                        yP.append(p[1])
    
                    pylab.plot(xP,yP, color='k')

            globJuncPose = self.getGlobalJunctionPose(pathID)
            if globJuncPose != None:
                pylab.scatter([globJuncPose[0],], [globJuncPose[1],], color='k')

            
            xP = []
            yP = []
            for p in vertices:
                xP.append(p[0])
                yP.append(p[1])
            pylab.plot(xP,yP, color='r')
            
            sizes = []
            for path in longPaths:
                sizes.append(len(path))
            
            bufStr1 = ""
            for dist in juncDists:
                if dist != None:
                    bufStr1 += "%1.2f " % dist


            bufStr2 = ""
            for juncAngs in juncAngSet:
                for angs in juncAngs:
                    if angs != None:
                        bufStr2 += "%1.2f %1.2f " % (angs[0],angs[1])
            
            pylab.axis("equal")
            pylab.title("Path %d %s %s %s" % (pathID, sizes,bufStr1,bufStr2))
            pylab.savefig("medialOut2_%04u.png" % self.topCount)
            self.topCount += 1


        
        print "juncAngSet:", juncAngSet

        
        " searches for the branch from the parent junction "
        " selects splice of topology that has aligning junction "
        " does not select for distance or best medial axis representation of path "
        

        " sort by longest first "
        pathCands = []
        for k in range(len(longPaths)):
            pathCands.append((len(longPaths[k]), k))

        pathCands.sort(reverse=True)
        
        print "pathCands:", pathCands
        
        
        " FIXME:  Verify that leaf path exists for the junction point "
        
        " select longest path that has the best angle fit "
        if globalJunctionPoint != None:
            bestFit = -1
            minDiff = 1e100
            for cand in pathCands:
                #for k in range(len(juncAngSet)):
                k = cand[1]
                juncAngs = juncAngSet[k]
                
                for angs in juncAngs:

                    if angs != None:
                        angDiff1 = angs[0]
                        angDiff2 = angs[1]
                        
                        print "minDiff:", minDiff, angDiff1, angDiff2 
                        
                        if fabs(angDiff1) < minDiff:
                            minDiff = fabs(angDiff1)
                            bestFit = k
    
                        if fabs(angDiff2) < minDiff:
                            minDiff = fabs(angDiff2)
                            bestFit = k

                pass


            if bestFit != -1:
                
                if fabs(minDiff) < 1.047:

                    print "returning bestFit:", bestFit, minDiff
                    return medialLongPaths[bestFit], vertices

                else:
                    print "returning bestFit theory:", theoryJunc
                    return theoryMedialLongPaths[0], vertices
            
            else:
                print "not returning bestFit"

        print "returning longest fit"
                

        " find the longest path from each leaf"
        maxPairDist = 0
        maxPair = None
        maxPath = []
        for leaf in leaves:
    
            isVisited = {}
            nodeSum = {}
            nodePath = {}
            for k, v in uni_mst.items():
                isVisited[k] = 0
                nodeSum[k] = 0
                nodePath[k] = []
    
            getLongestPath(leaf, 0, [], uni_mst, isVisited, nodePath, nodeSum)
    
            maxDist = 0
            maxNode = None
            for k, v in nodeSum.items():
                if v > maxDist:
                    maxNode = k
                    maxDist = v
    
            if maxDist > maxPairDist:
                maxPairDist = maxDist
                maxPair = [leaf, maxNode]
                maxPath = nodePath[maxNode]

        " convert path points to real "
        realPath = []
        for p in maxPath:
            realPath.append(gridToReal(p))


        leafPath = deepcopy(realPath)
        
        frontVec = [0.,0.]
        backVec = [0.,0.]
        indic = range(3)
        #indic = range(16)
        indic.reverse()
    
        for i in indic:
            p1 = leafPath[i+2]
            p2 = leafPath[i]
            vec = [p2[0]-p1[0], p2[1]-p1[1]]
            frontVec[0] += vec[0]
            frontVec[1] += vec[1]
    
            p1 = leafPath[-i-3]
            p2 = leafPath[-i-1]
            vec = [p2[0]-p1[0], p2[1]-p1[1]]
            backVec[0] += vec[0]
            backVec[1] += vec[1]
    
        frontMag = math.sqrt(frontVec[0]*frontVec[0] + frontVec[1]*frontVec[1])
        backMag = math.sqrt(backVec[0]*backVec[0] + backVec[1]*backVec[1])
    
        frontVec[0] /= frontMag
        frontVec[1] /= frontMag
        backVec[0] /= backMag
        backVec[1] /= backMag
    
        newP1 = (leafPath[0][0] + frontVec[0]*10, leafPath[0][1] + frontVec[1]*10)
        newP2 = (leafPath[-1][0] + backVec[0]*10, leafPath[-1][1] + backVec[1]*10)
    
        leafPath.insert(0,newP1)
        leafPath.append(newP2)

        medial2 = deepcopy(leafPath)

        " take the long length segments at tips of medial axis"
        edge1 = medial2[0:2]
        edge2 = medial2[-2:]
        
        frontVec = [edge1[0][0]-edge1[1][0], edge1[0][1]-edge1[1][1]]
        backVec = [edge2[1][0]-edge2[0][0], edge2[1][1]-edge2[0][1]]
        frontMag = math.sqrt(frontVec[0]*frontVec[0] + frontVec[1]*frontVec[1])
        backMag = math.sqrt(backVec[0]*backVec[0] + backVec[1]*backVec[1])
        
        frontVec[0] /= frontMag
        frontVec[1] /= frontMag
        backVec[0] /= backMag
        backVec[1] /= backMag
        
        " make a smaller version of these edges "
        newP1 = (edge1[1][0] + frontVec[0]*2, edge1[1][1] + frontVec[1]*2)
        newP2 = (edge2[0][0] + backVec[0]*2, edge2[0][1] + backVec[1]*2)

        edge1 = [newP1, edge1[1]]
        edge2 = [edge2[0], newP2]
    
        " find the intersection points with the hull "
        hull = vertices
        interPoints = []
        for k in range(len(hull)-1):
            hullEdge = [hull[k],hull[k+1]]
            isIntersect1, point1 = Intersect(edge1, hullEdge)
            if isIntersect1:
                interPoints.append(point1)
                break

        for k in range(len(hull)-1):
            hullEdge = [hull[k],hull[k+1]]
            isIntersect2, point2 = Intersect(edge2, hullEdge)
            if isIntersect2:
                interPoints.append(point2)
                break
        
        " replace the extended edges with a termination point at the hull edge "            
        medial2 = medial2[1:-2]


        
        if isIntersect1:
            medial2.insert(0, point1)
        if isIntersect2:
            medial2.append(point2)


       
        return medial2, vertices
        #return realPath, vertices



    def getUnionTopology(self):
        
        random.seed(0)        
 
        allPathIDs = self.getPathIDs()

        nodes = []
        for pID in allPathIDs:
            nodes += self.getNodes(pID)
                                    
        
        def convertAlphaUniform(a_vert, max_spacing = 0.04):
            
            " make the vertices uniformly distributed "
            
            new_vert = []
        
            for i in range(len(a_vert)):
                p0 = a_vert[i]
                p1 = a_vert[(i+1) % len(a_vert)]
                dist = math.sqrt((p0[0]-p1[0])**2 + (p0[1]-p1[1])**2)
    
                vec = [p1[0]-p0[0], p1[1]-p0[1]]
                
                if dist == 0:
                    pass
                else:
                                        
                    vec[0] /= dist
                    vec[1] /= dist
                    
                    new_vert.append(copy(p0))
                    
                    if dist > max_spacing:
                        " cut into pieces max_spacing length or less "
                        numCount = int(math.floor(dist / max_spacing))
                        
                        for j in range(1, numCount+1):
                            newP = [j*max_spacing*vec[0] + p0[0], j*max_spacing*vec[1] + p0[1]]
                            new_vert.append(newP)
            
            return new_vert            

        medialPointSoup = []


        if nodes == []:
            return [], []


        if False:    
            for nodeID in nodes:
                estPose1 = self.nodeHash[nodeID].getGlobalGPACPose()        
        
                if self.nodeHash[nodeID].isBowtie:            
                    hull1 = computeBareHull(self.nodeHash[nodeID], sweep = False, static = True)
                else:
                    hull1 = computeBareHull(self.nodeHash[nodeID], sweep = False)
        
        
                m = hashlib.md5()
                m.update(repr(hull1))
                print nodeID, "hull1 =", int(m.digest().encode('hex'),16)
        
                " set the origin of pose 1 "
                poseOrigin = Pose(estPose1)
        
                xP = []
                yP = []    
                for k in range(len(hull1)):
     
                    p = hull1[k]
                    
                    m = hashlib.md5()
                    m.update(repr(p))
                    #print nodeID, "p =", int(m.digest().encode('hex'),16)               
                    
                    p1 = poseOrigin.convertLocalToGlobal(p)
                    
                    m = hashlib.md5()
                    m.update(repr(p1))
                    #print nodeID, "p1 =", int(m.digest().encode('hex'),16)
                    
                    medialPointSoup.append(p1)

        if True:
    
            pylab.clf()
            xP = []
            yP = []
            for nodeID in nodes:
                estPose1 = self.nodeHash[nodeID].getEstPose()        
                
                poseOrigin = Pose(estPose1)
        
                numPixel = self.nodeHash[nodeID].occMap.numPixel
                mapImage = self.nodeHash[nodeID].occMap.getMap()
                image = mapImage.load()
                
                for j in range(numPixel):
                    for k in range(numPixel):
                        if image[j,k] == 255:
                            pnt = self.nodeHash[nodeID].occMap.gridToReal([j,k])
                            p1 = poseOrigin.convertLocalToGlobal(pnt)
                            medialPointSoup.append(p1)
                            xP.append(p1[0])
                            yP.append(p1[1])
            
            pylab.scatter(xP,yP)
            pylab.title("%d: %s" % (-1, repr(nodes)))
            pylab.savefig("plotMedialSoup_%04u.png" % (self.topCount ))
        
        saveFile = "medialPointSoup = ["
        for k in range(len(medialPointSoup)):
            saveFile += repr(medialPointSoup[k]) + "\n"
        saveFile += "]\n"

        f = open("medialSoupSave_%04u.txt" % self.medialSoupCount, 'w')
        f.write(saveFile)
        f.close()
        self.medialSoupCount  += 1
        
        #nodeID = self.numNodes - 1
            
        #radius = 0.2
        radius = 0.05

        numPoints = len(medialPointSoup)

        isDone = False
        
        while not isDone:
    
            perturbPoints = []
            
            for p in medialPointSoup:
                p2 = copy(p)
                " add a little bit of noise to avoid degenerate conditions in CGAL "
                p2[0] += random.gauss(0.0,0.000001)
                p2[1] += random.gauss(0.0,0.000001)
    
                perturbPoints.append(p2)
        
            try:            
    
                vertices = alphamod.doAlpha(radius,perturbPoints)
                numVert = len(vertices)
                
                if numVert <= 2:
                    print "Failed, hull had only", numVert, "vertices"
                    raise
                
                isDone = True
            except:
                print "hull has holes!  retrying..."
                #print sArr     
       
        " cut out the repeat vertex "
        vertices = vertices[:-1]
        
        vertices = convertAlphaUniform(vertices)

        vertices.append(vertices[0])
        
        minX = 1e100
        maxX = -1e100
        minY = 1e100
        maxY = -1e100
        for p in vertices:
            if p[0] > maxX:
                maxX = p[0]
            if p[0] < minX:
                minX = p[0]
            if p[1] > maxY:
                maxY = p[1]
            if p[1] < minY:
                minY = p[1]
    
    
        " SPECIFY SIZE OF GRID AND CONVERSION PARAMETERS "
        PIXELSIZE = 0.05
        mapSize = 2*max(max(maxX,math.fabs(minX)),max(maxY,math.fabs(minY))) + 1
        pixelSize = PIXELSIZE
        numPixel = int(2.0*mapSize / pixelSize + 1.0)
        divPix = math.floor((2.0*mapSize/pixelSize)/mapSize)
        
        def realToGrid(point):
            indexX = int(math.floor(point[0]*divPix)) + numPixel/2 + 1
            indexY = int(math.floor(point[1]*divPix)) + numPixel/2 + 1
            return indexX, indexY
    
        def gridToReal(indices):
            i = indices[0]
            j = indices[1]
            point = [(i - numPixel/2 - 1)*(pixelSize/2.0) + pixelSize/2.0, (j - numPixel/2 - 1)*(pixelSize/2.0) + pixelSize/2.0]
            return point
    
        " CONVERT HULL TO GRID COORDINATES "
        gridHull = []
        for i in range(len(vertices)):
            p = vertices[i]
            gridHull.append(realToGrid(p))

        " BOUNDING BOX OF GRID HULL "    
        minX = 1e100
        maxX = -1e100
        minY = 1e100
        maxY = -1e100
        for p in gridHull:
            if p[0] > maxX:
                maxX = p[0]
            if p[0] < minX:
                minX = p[0]
            if p[1] > maxY:
                maxY = p[1]
            if p[1] < minY:
                minY = p[1]

        " COMPUTE MEDIAL AXIS OF HULL "
        resultImg = Image.new('L', (numPixel,numPixel))
        resultImg = computeMedialAxis(self.medialCount, numPixel,numPixel, 5, resultImg, len(gridHull[:-2]), gridHull[:-2])
        self.medialCount += 1

        " EXTRACT POINTS FROM GRID TO LIST "
        imgA = resultImg.load()    
        points = []
        for i in range(1, numPixel-1):
            for j in range(1, numPixel-1):
                if imgA[i,j] == 0:
                    points.append((i,j))
    
        " CREATE GRAPH NODES FOR EACH POINT "
        medialGraph = graph.graph()
        for p in points:
            medialGraph.add_node(p, [])
    
        " ADD EDGES BETWEEN NEIGHBORS "
        builtGraph = {}
        for i in range(2, numPixel-2):
            for j in range(2, numPixel-2):
                if imgA[i,j] == 0:
                    builtGraph[(i,j)] = []
                    for k in range(i-1, i+2):
                        for l in range(j-1, j+2):
                            if imgA[k,l] == 0 and not (k == i and l == j):
                                builtGraph[(i,j)].append((k,l))
                                medialGraph.add_edge((i,j), (k,l))
                                
        " COMPUTE MINIMUM SPANNING TREE "
        mst = medialGraph.minimal_spanning_tree()
        
        " INITIALIZE DATA DICTS FOR UNIDIRECTIONAL MST"
        uni_mst = {}
        for k, v in mst.items():
            uni_mst[k] = []

        
        " ADD EDGES TO DICT TREE REPRESENTATION "
        for k, v in mst.items():
            if v != None:
                uni_mst[k].append(v)
                uni_mst[v].append(k)

        " LOCATE ALL LEAVES "
        leaves = []
        for k, v in uni_mst.items():
            if len(v) == 1:
                leaves.append(k)
        
        " DELETE ALL NODES THAT ARE LEAVES, TO REMOVE SINGLE NODE BRANCHES "
        for leaf in leaves:
            medialGraph.del_node(leaf)    
            
        " RECOMPUTE MST "
        mst = medialGraph.minimal_spanning_tree()

        " AGAIN, CREATE OUR DATA STRUCTURES AND IDENTIFY LEAVES "
        uni_mst = {}
        isVisited = {}
        nodeSum = {}
        for k, v in mst.items():
            uni_mst[k] = []
            isVisited[k] = 0
            nodeSum[k] = 0
        
        for k, v in mst.items():
            if v != None:
                uni_mst[k].append(v)
                uni_mst[v].append(k)
                        
        " RECORD THE LEAVES AND JUNCTIONS "
        leaves = []
        junctions = []
        for k, v in uni_mst.items():
            if len(v) == 1:
                leaves.append(k)

            if len(v) > 2:
                junctions.append(k)

        
        " SAVE FOR LATER, IDENTIFIED BY THEIR INDEX NOW "
        numLeafs = len(leaves)
        allLeafs = []
        allJunctions = []
        for leaf in leaves:
            allLeafs.append(gridToReal(leaf))

         
        print "allJunctions:", allJunctions
        
        " FIND ALL THE PATHS BETWEEN LEAVES "        
        nodePaths = {}
        for leaf in leaves:
            
            " INITIALIZE DATA DICTS FOR VISITED, PATH SO FAR, AND NODE HOP COUNT"
            isVisited = {}
            nodeSum = {}
            nodePath = {}
            for k, v in uni_mst.items():
                isVisited[k] = 0
                nodeSum[k] = 0
                nodePath[k] = []
    
            " PERFORM DFS FOR THIS LEAF "
            getLongestPath(leaf, 0, [], uni_mst, isVisited, nodePath, nodeSum)
    
            " SAVE THE RESULTING PATH DATA STRUCTURE "
            nodePaths[leaf] = nodePath
            

        " FOR EVERY PAIR OF LEAVES, SAVE ITS PATH IF ITS LONG ENOUGH "
        " should have X choose 2 combinations"
        longPaths = []
        
        MAX_LEN = 2
        for leaf1 in leaves:
            for leaf2 in leaves:
                if leaf1 < leaf2:
                    if len(nodePaths[leaf1][leaf2]) > MAX_LEN:
                        nPath = deepcopy(nodePaths[leaf1][leaf2])
                        juncIndices = []
                        for junc in allJunctions:
                            try:
                                index = nPath.index(junc[0])
                            except:
                                juncIndices.append(None)
                            else:
                                juncIndices.append(index)
                                
                        longPaths.append((len(nPath),nPath, juncIndices))

            
        print "longPaths:"
        for longPath in longPaths:
            print longPath[2]
        
        " SORT FOR THE LONGEST TO SHORTEST "
        longPaths.sort(reverse=True)
        
        " REMOVE SIZE FROM TUPLE  "
        juncIndices = []
        for k in range(len(longPaths)):
            juncIndices.append(longPaths[k][2])
            longPaths[k] = longPaths[k][1]
        
        print "juncIndices:", juncIndices
        
        " GET THE LEAF INDEXES TO EACH PATH "
        leafPairs = []
        for path in longPaths:
            leaf1 = path[0]
            leaf2 = path[-1]    
            
            leafID1 = leaves.index(leaf1)
            leafID2 = leaves.index(leaf2)

            leafPairs.append((leafID1,leafID2))
        
        print "leafPairs:", leafPairs

        for k in range(len(longPaths)):
            path = longPaths[k]
            realPath = []
            for p in path:
                realPath.append(gridToReal(p))
            
            longPaths[k] = realPath
                    

        #longMedialWidths = []
        medialLongPaths = []
        juncAngSet = []
        #medialTailCuts = []
        #for longPath in longPaths:
            
        print len(longPaths), "long paths"
        for n in range(len(longPaths)):
            
            longPath = longPaths[n]
            print "longPath:", len(longPath)
            
            leafPath = deepcopy(longPath)
            
            frontVec = [0.,0.]
            backVec = [0.,0.]
            indic = range(3)
            indic.reverse()
            
            for i in indic:
                if i+2 < len(leafPath):
                    p1 = leafPath[i+2]
                    p2 = leafPath[i]
                    vec = [p2[0]-p1[0], p2[1]-p1[1]]
                    frontVec[0] += vec[0]
                    frontVec[1] += vec[1]
        
                    p1 = leafPath[-i-3]
                    p2 = leafPath[-i-1]
                    vec = [p2[0]-p1[0], p2[1]-p1[1]]
                    backVec[0] += vec[0]
                    backVec[1] += vec[1]
        
            frontMag = math.sqrt(frontVec[0]*frontVec[0] + frontVec[1]*frontVec[1])
            backMag = math.sqrt(backVec[0]*backVec[0] + backVec[1]*backVec[1])
        
            frontVec[0] /= frontMag
            frontVec[1] /= frontMag
            backVec[0] /= backMag
            backVec[1] /= backMag
        
            newP1 = (leafPath[0][0] + frontVec[0]*10, leafPath[0][1] + frontVec[1]*10)
            newP2 = (leafPath[-1][0] + backVec[0]*10, leafPath[-1][1] + backVec[1]*10)
        
            leafPath.insert(0,newP1)
            leafPath.append(newP2)
            
            medial2 = deepcopy(leafPath)
    
            " take the long length segments at tips of medial axis"
            edge1 = medial2[0:2]
            edge2 = medial2[-2:]
            
            frontVec = [edge1[0][0]-edge1[1][0], edge1[0][1]-edge1[1][1]]
            backVec = [edge2[1][0]-edge2[0][0], edge2[1][1]-edge2[0][1]]
            frontMag = math.sqrt(frontVec[0]*frontVec[0] + frontVec[1]*frontVec[1])
            backMag = math.sqrt(backVec[0]*backVec[0] + backVec[1]*backVec[1])
            
            frontVec[0] /= frontMag
            frontVec[1] /= frontMag
            backVec[0] /= backMag
            backVec[1] /= backMag
            
            " make a smaller version of these edges "
            newP1 = (edge1[1][0] + frontVec[0]*2, edge1[1][1] + frontVec[1]*2)
            newP2 = (edge2[0][0] + backVec[0]*2, edge2[0][1] + backVec[1]*2)
    
            edge1 = [newP1, edge1[1]]
            edge2 = [edge2[0], newP2]

            " find the intersection points with the hull "
            hull = vertices
            interPoints = []
            for k in range(len(hull)-1):
                hullEdge = [hull[k],hull[k+1]]
                isIntersect1, point1 = Intersect(edge1, hullEdge)
                if isIntersect1:
                    interPoints.append(point1)
                    break
    
            for k in range(len(hull)-1):
                hullEdge = [hull[k],hull[k+1]]
                isIntersect2, point2 = Intersect(edge2, hullEdge)
                if isIntersect2:
                    interPoints.append(point2)
                    break
            
            " replace the extended edges with a termination point at the hull edge "            
            medial2 = medial2[1:-2]
            
            if isIntersect1:
                medial2.insert(0, point1)
            if isIntersect2:
                medial2.append(point2)

            medialLongPaths.append(deepcopy(medial2))

            print "juncIndices:", juncIndices[n]

            juncAngs = []

        juncDists = []
        for junc in allJunctions:
            juncDists.append(junc[2])
        
        if True:
            pylab.clf()
    
            for path in medialLongPaths:
                xP = []
                yP = []
                for p in path:
                    #p2 = gridToReal(p)
                    xP.append(p[0])
                    yP.append(p[1])
    
                pylab.plot(xP,yP)
            #pylab.scatter(xP,yP, color='b')
            
            xP = []
            yP = []
            for p in vertices:
                xP.append(p[0])
                yP.append(p[1])
            pylab.plot(xP,yP, color='r')
            
            sizes = []
            for path in longPaths:
                sizes.append(len(path))
            
            bufStr1 = ""
            for dist in juncDists:
                if dist != None:
                    bufStr1 += "%1.2f " % dist


            bufStr2 = ""
            for juncAngs in juncAngSet:
                for angs in juncAngs:
                    if angs != None:
                        bufStr2 += "%1.2f %1.2f " % (angs[0],angs[1])
            
            pylab.axis("equal")
            pylab.title("Path %d %s %s %s" % (-1, sizes,bufStr1,bufStr2))
            pylab.savefig("medialOut2_%04u.png" % self.topCount)
            self.topCount += 1




    def trimPaths(self, paths):

        trimmedPaths = {}
        
        print "path lengths:"
        for k,v in paths.iteritems():
            print k, len(v)
        
        
        if len(paths) <= 1:
            for k, v in paths.iteritems():
                trimmedPaths[k] = v
                
            self.trimmedPaths = trimmedPaths
            return trimmedPaths

        " path0 is root, and path1 is a child of path0 "
        #parents = [None, [0, self.junctionNodeID, self.localJunctionPoint]]
        #parents = self.pathParents

        " information on parentage is required to determine which path belongs to which "

        " global junction point for parent 0 and child 1 "

        #print "pathParents:", parents
        pathIDs = self.getPathIDs()

        for pathID in pathIDs:
            
            #junctionPose = self.pathClasses[pathID]["localJunctionPose"]
            #poseOrigin2 = Pose(self.nodeHash[branchNodeID].getEstPose())
            #globalJunctionPoint = poseOrigin2.convertLocalToGlobal(localJunctionPoint)            

           
            
            if self.pathClasses[pathID]["parentID"] == None:
                trimmedPaths[pathID] = deepcopy(paths[pathID])

            else:
                parentPathID = self.pathClasses[pathID]["parentID"]
                childPathID = pathID

                path1 = paths[parentPathID]
                path2 = paths[childPathID]
                


                #secP1, secP2 = self.getOverlapDeparture(path1, path2)                
                secP1, secP2 = self.getOverlapDeparture(parentPathID, childPathID, paths)                

                minI_1 = 0        
                minI_2 = 0
                minDist_1 = 1e100
                minDist_2 = 1e100        
                for i in range(len(path2)):
                    pnt = path2[i]
                    dist1 = sqrt((pnt[0]-secP1[0])**2 + (pnt[1]-secP1[1])**2)
                    dist2 = sqrt((pnt[0]-secP2[0])**2 + (pnt[1]-secP2[1])**2)
                
                    if dist1 < minDist_1:
                        minDist_1 = dist1
                        minI_1 = i

                    if dist2 < minDist_2:
                        minDist_2 = dist2
                        minI_2 = i

                term0 = path2[0]
                termN = path2[-1]

                " smallest distance is the terminal point "
                dist0_1 = sqrt((term0[0]-secP1[0])**2 + (term0[1]-secP1[1])**2)
                distN_1 = sqrt((termN[0]-secP1[0])**2 + (termN[1]-secP1[1])**2)
                dist0_2 = sqrt((term0[0]-secP2[0])**2 + (term0[1]-secP2[1])**2)
                distN_2 = sqrt((termN[0]-secP2[0])**2 + (termN[1]-secP2[1])**2)
                print "terminal distances:", dist0_1, distN_1, dist0_2, distN_2
                
                distList = [dist0_1, distN_1, dist0_2, distN_2]
                
                minI = -1
                minDist = 1e100
                for i in range(len(distList)):
                    if distList[i] < minDist:
                        minDist = distList[i]
                        minI = i
                
                #junctionNodeID = self.pathClasses[pathID]["branchNodeID"]
                #localJunctionPoint = self.pathClasses[pathID]["localJunctionPose"]
                
                #poseOrigin = Pose(self.nodeHash[junctionNodeID].getEstPose())
                #junctionPoint = poseOrigin.convertLocalToGlobal(localJunctionPoint)
                
                junctionPoint = self.getGlobalJunctionPose(pathID)
                
                if minDist_1 < minDist_2:
                    junctionPoint_K = secP1
                    juncI_K = minI_1
                else:
                    junctionPoint_K = secP2
                    juncI_K = minI_2
                
                #junctionPoint = junctionPoint_K
                #juncI = juncI_K
                
                
                minDist2 = 1e100
                juncI = 0        
                for i in range(len(path2)):
                    pnt = path2[i]
                    dist = sqrt((pnt[0]-junctionPoint[0])**2 + (pnt[1]-junctionPoint[1])**2)
                
                    if dist < minDist2:
                        minDist2 = dist
                        juncI = i
                 
                 
                
                print "len(path1):", len(path1)
                print "len(path2):", len(path2)
                print "juncI:", juncI
                print "minDist:", minDist_1, minDist_2
                
                " now we have closest point to departure point. "
                " Which side is the departing side? "    
                #sum1 = self.getOverlapSum(path1, path2[0:minI2+1] + [junctionPoint])
                #sum2 = self.getOverlapSum(path1, [junctionPoint] + path2[minI2:])
                

                        
                if minI == 0:
                    "secP1 is terminal 0"
                    index = juncI-10
                    #index = juncI
                    if index < 1:
                        index = 1
                    #newPath2 = path2[:index+1] + [junctionPoint]
                    
                    #self.pathClasses[pathID]["globalJunctionPose"] = path2[index]

                    
                    newPath2 = path2[:index+1]
                    newPath2.reverse()
                
                elif minI == 1:
                    "secP1 is terminal N"
                    index = juncI+10
                    #index = juncI
                    if index >= len(path2)-1:
                        
                        " ensure at least 2 elements in path "
                        index = len(path2)-2

                    #self.pathClasses[pathID]["globalJunctionPose"] = path2[index]

                    #newPath2 = [junctionPoint] + path2[index:]
                    newPath2 = path2[index:]

                    
                elif minI == 2:
                    "secP2 is terminal 0"
                    index = juncI-10
                    #index = juncI
                    if index < 1:
                        index = 1

                    #self.pathClasses[pathID]["globalJunctionPose"] = path2[index]

                    #newPath2 = path2[:index+1] + [junctionPoint]
                    newPath2 = path2[:index+1]
                    newPath2.reverse()
                    
                elif minI == 3:
                    "secP2 is terminal N"
                    index = juncI+10
                    #index = juncI
                    if index >= len(path2)-1:
                        " ensure at least 2 elements in path "
                        index = len(path2)-2
                    
                    #self.pathClasses[pathID]["globalJunctionPose"] = path2[index]
                    
                    #newPath2 = [junctionPoint] + path2[index:]
                    newPath2 = path2[index:]
                
                else:
                    print "no terminal found"
                    raise
                                
                #if minI_1 > minI_2:
                #    path2[minI_2:minI_1+1]
                #else:
                #    path2[minI_1:minI_2+1]
                
                
                #if startI == 0:
                #    newPath2 = path2[startI:endI+1] + [junctionPoint]
                #    newPath2.reverse()
                #else:
                #    newPath2 = [junctionPoint] + path2[startI:endI+1]
                    

                #print "Sums:", sum1, sum2
                
                #if sum1 > sum2:
                #    index = minI2+1-10
                #    if index < 1:
                #        index = 1
                #    newPath2 = path2[0:index] + [junctionPoint]
                #    newPath2.reverse()
                #else:
                #    index = minI2+10
                #    if index >= len(path2):
                #        index = len(path2)-1
                #    newPath2 = [junctionPoint] + path2[index:]

                " convert path so that the points are uniformly distributed "
                
                
                max_spacing = 0.08
                #max_spacing = 0.04
                #max_spacing = 0.1
                newPath3 = []

                " make sure path is greater than 5 points "
                while len(newPath3) <= 5:
    
                    #print "len(newPath3) =", len(newPath3), "<= 5"
    
                    max_spacing /= 2
                    print "max_spacing =", max_spacing
                    
                    newPath3 = [copy(newPath2[0])]
                                        
                    for i in range(len(newPath2)-1):
                        p0 = newPath2[i]
                        p1 = newPath2[(i+1)]
                        dist = sqrt((p0[0]-p1[0])**2 + (p0[1]-p1[1])**2)
            
                        vec = [p1[0]-p0[0], p1[1]-p0[1]]
                        vec[0] /= dist
                        vec[1] /= dist
                        
                        
                        if dist > max_spacing:
                            " cut into pieces max_spacing length or less "
                            numCount = int(floor(dist / max_spacing))
                            
                            for j in range(1, numCount+1):
                                newP = [j*max_spacing*vec[0] + p0[0], j*max_spacing*vec[1] + p0[1]]
                                newPath3.append(newP)
    
                        newPath3.append(copy(p1))            
                
                #trimmedPaths.append(deepcopy(newPath3))
                trimmedPaths[pathID] = deepcopy(newPath3)

        """
        self.drawTrimmedPaths(self.trimmedPaths)
        """
        
        self.trimmedPaths = trimmedPaths
        
        return trimmedPaths


    " get the trimmed version of child and parent paths that are overlapping in some fashion "
    def getOverlapDeparture(self, parentPathID, childPathID, paths):

        "Assumption:  one section of the medial axis is closely aligned with the path "        
            
        plotIter = False
        print "getOverlapDeparture():"
        
        path1 = paths[parentPathID]
        path2 = paths[childPathID]
        
        
        isExist1 = False
        isInterior1 = False
        departurePoint1 = 0
        isExist2 = False
        isInterior2 = False
        departurePoint2 = 0


        branchNodeID = self.pathClasses[childPathID]["branchNodeID"]
        poseOrigin = Pose(self.nodeHash[branchNodeID].getGlobalGPACPose())
        #localJunctionPoint = self.pathClasses[childPathID]["localJunctionPose"]
        
        
        #poseOrigin2 = Pose(self.nodeHash[branchNodeID].getEstPose())
        #globalJunctionPoint = poseOrigin2.convertLocalToGlobal(localJunctionPoint)
        
        globalJunctionPoint = self.getGlobalJunctionPose(childPathID)
        
        " orienting the medial axis of the branch node correctly "
        #branchDir = self.pathParents[childPathID][3]
        
        """
        branchHull, branchMedial = computeHullAxis(branchNodeID, self.nodeHash[branchNodeID], tailCutOff = False)

        branchMedialSpline = SplineFit(branchMedial, smooth=0.1)
        branchPoints = branchMedialSpline.getUniformSamples()
        
        globalBranchPoints = []
        for p in branchPoints:
            result = poseOrigin.convertLocalToGlobal(p)
            globalBranchPoints.append(result)
        
        #if not branchDir:
        #    " backward "
        #    globalBranchPoints.reverse()
            
        globalBranchSpline = SplineFit(globalBranchPoints, smooth = 0.1)
        """
            

        " return exception if we receive an invalid path "        
        if len(path1) == 0:
            print "path1 has zero length"
            raise
            #return departurePoint1, isInterior1, isExist1, departurePoint2, isInterior2, isExist2


        " make sure the overlap of both paths are oriented the same way "
        path1Spline = SplineFit(path1, smooth=0.1)
        path2Spline = SplineFit(path2, smooth=0.1)

        path2Reverse = deepcopy(path2)
        path2Reverse.reverse()
        path2SplineReverse = SplineFit(path2Reverse, smooth=0.1)
        
        overlapMatch = []
        angleSum1 = 0.0
        angleSum2 = 0.0
        for i in range(0,len(path1)):
            p_1 = path1[i]
            p_2, j, minDist = gen_icp.findClosestPointInA(path2, p_1)
            if minDist < 0.1:
                overlapMatch.append((i,j,minDist))

                pathU1 = path1Spline.findU(p_1)    
                pathU2 = path2Spline.findU(p_2)    
                pathU2_R = path2SplineReverse.findU(p_2)    

                pathVec1 = path1Spline.getUVector(pathU1)
                pathVec2 = path2Spline.getUVector(pathU2)
                pathVec2_R = path2SplineReverse.getUVector(pathU2_R)

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
        print "angle sums for orienting child path with parent path:", angleSum1, angleSum2

        if angleSum1 > angleSum2:
            orientedPath2 = path2Reverse
        else:
            orientedPath2 = path2
            
        path2Spline = SplineFit(orientedPath2, smooth=0.1)



        " for each point on the child path, find its closest pair on the parent path "
        
        pathPoints1 = path1Spline.getUniformSamples()
        pathPoints2 = path2Spline.getUniformSamples()

        distances = []
        indices = []
        for i in range(0,len(pathPoints2)):
            p_2 = pathPoints2[i]
            p_1, i_1, minDist = gen_icp.findClosestPointInA(pathPoints1, p_2)
            
            " keep the distance information "
            distances.append(minDist)
            " and the associated index of the point on the parent path "
            indices.append(i_1)

        " match distances of the tip points of child path "        
        maxFront = distances[0]
        maxBack = distances[-1]
        print "match distances of tip points:", maxFront, maxBack

        " walk back from tip point until we have a non-monotic increase in match distance "
        " this becomes our departure point "
        
        "TODO:  why is there a 3-offset in this comparison? "
        currI = 1
        try:
            while distances[currI+3] < maxFront:
                maxFront = distances[currI]
                currI += 1
        except:
            pass
        
        " departure index on child path "
        frontDepI = currI
        
        " departure index of child path and match distance "
        frontPoint = [frontDepI, distances[frontDepI]]

        " FIXME:  index out of bounds case "
        currI = 2
        try:
            while distances[-currI-3] < maxBack:
                maxBack = distances[-currI]
                currI += 1
        except:
            pass

        " departure index on child path "
        backDepI = len(distances) - currI

        " departure index of child path and match distance "
        backPoint = [backDepI, distances[backDepI]]

        print "lengths of parent and child paths:", len(pathPoints1), len(pathPoints2)
        print "non monotonic departures:", maxFront, maxBack, frontPoint, backPoint


        "reset to the tip match distance "
        maxFront = distances[0]
        maxBack = distances[-1]
        
        
        " count the number of times a matched point on the parent path is used "
        foo = indices[0:frontDepI+1]
        d1 = {}
        for i in set(foo):
            d1[i] = foo.count(i)

        foo = indices[backDepI:]
        d2 = {}
        for i in set(foo):
            d2[i] = foo.count(i)

        " select the match point that is used the most on the parent path "
        max1 = max(d1, key=d1.get)
        max2 = max(d2, key=d2.get)


        " NOTE:  This guarantees detection of a departure.  We've assumed that there exists a"
        " departure between these two paths "
        
        " compute two candidate departure points "
        
        " NOTE:  departure point on parent path is not used here "
        if True:

            #departurePoint1 = pathPoints1[max1]
            isExist1 = True

            " if the highest match point is either of the tips of the parent path, then this is an external departure "
            if max1 == 0 or max1 == len(pathPoints1)-1:
                isInterior1 = False
            else:
                isInterior1 = True

        if True:
            
            #departurePoint2 = pathPoints1[max2]
            isExist2 = True

            " if the highest match point is either of the tips of the parent path, then this is an external departure "
            if max2 == 0 or max2 == len(pathPoints1)-1:
                isInterior2 = False
            else:
                isInterior2 = True

        print "isExist1 =", isExist1, "isInterior1 =", isInterior1
        print "isExist2 =", isExist2, "isInterior2 =", isInterior2

        " sum of closest points on front and back "
        " select the one with minimal cost "        
        
        " now we compare our candidates to the known direction and location of the branching point "
        angleSum1 = 0.0
        overlapSum1 = 0.0
        matchCount1 = 0
        
        " path section for our front departure hypothesis "
        pathSec1 = pathPoints2[:frontDepI+1]
        pathSec1.reverse()
        
        """
        if len(pathSec1) > 5:
            pathSec1Spline = SplineFit(pathSec1, smooth = 0.1)        
            
            
            for i in range(0,len(pathSec1)):
                p_1 = pathSec1[i]
                p_2, j, minDist = gen_icp.findClosestPointInA(globalBranchPoints, p_1)
                if minDist < 0.1:
                    overlapMatch.append((i,j,minDist))
                    matchCount1 += 1
                    
                    overlapSum1 += minDist
    
                    pathU1 = pathSec1Spline.findU(p_1)    
                    pathUB = globalBranchSpline.findU(p_2)    
    
                    pathVec1 = pathSec1Spline.getUVector(pathU1)
                    pathVecB = globalBranchSpline.getUVector(pathUB)
    
                    val1 = pathVec1[0]*pathVecB[0] + pathVec1[1]*pathVecB[1]
                    if val1 > 1.0:
                        val1 = 1.0
                    elif val1 < -1.0:
                        val1 = -1.0
                    ang1 = acos(val1)
                    
                    angleSum1 += ang1

            if matchCount1 > 0:
                angleSum1 /= matchCount1
                overlapSum1 /= matchCount1
        """

        angleSum2 = 0.0
        overlapSum2 = 0.0
        matchCount2 = 0
        
        " path section for our back departure hypothesis "
        pathSec2 = pathPoints2[backDepI:]

        """
        if len(pathSec2) > 5:
            pathSec2Spline = SplineFit(pathSec2, smooth = 0.1)        
    
    
            for i in range(0,len(pathSec2)):
                p_1 = pathSec2[i]
                p_2, j, minDist = gen_icp.findClosestPointInA(globalBranchPoints, p_1)
                if minDist < 0.1:
                    overlapMatch.append((i,j,minDist))
                    matchCount2 += 1
                    
                    overlapSum2 += minDist
    
                    pathU2 = pathSec2Spline.findU(p_1)    
                    pathUB = globalBranchSpline.findU(p_2)    
    
                    pathVec2 = pathSec2Spline.getUVector(pathU2)
                    pathVecB = globalBranchSpline.getUVector(pathUB)
    
                    val2 = pathVec2[0]*pathVecB[0] + pathVec2[1]*pathVecB[1]
                    if val2 > 1.0:
                        val2 = 1.0
                    elif val2 < -1.0:
                        val2 = -1.0
                    ang2 = acos(val2)
                    
                    angleSum2 += ang2

            if matchCount2 > 0:
                angleSum2 /= matchCount2
                overlapSum2 /= matchCount2
        """
        print "pathSec1 hypothesis angle and overlap sum and match count:", angleSum1, overlapSum1, matchCount1
        print "pathSec2 hypothesis angle and overlap sum and match count:", angleSum2, overlapSum2, matchCount2

        " distance of departure point from known junction point "
        #juncDist1 = -1
        #juncDist2 = -1
        #if len(pathSec1) > 0:
        p0 = pathSec1[0]
        juncDist1 = sqrt((globalJunctionPoint[0]-p0[0])**2 + (globalJunctionPoint[1]-p0[1])**2)

        #if len(pathSec2) > 0:
        p0 = pathSec2[0]
        juncDist2 = sqrt((globalJunctionPoint[0]-p0[0])**2 + (globalJunctionPoint[1]-p0[1])**2)

        print "pathSec1 hypothesis discrepancy distance:", juncDist1
        print "pathSec2 hypothesis discrepancy distance:", juncDist2

        if plotIter:
            pylab.clf()
            xP = []
            yP = []
            for p in path2:
                xP.append(p[0])
                yP.append(p[1])
            pylab.plot(xP,yP, color=(0.5,0.5,1.0))
    
            if True:    
                P1 = pathPoints2[0]
                P2 = pathPoints2[frontDepI]
                
                P3 = pathPoints2[backDepI]
                P4 = pathPoints2[-1]

                " 0, frontDepI "
                xP = [P1[0], P2[0]]
                yP = [P1[1], P2[1]]
                pylab.scatter(xP,yP, color='b')        

                " backDepI, -1 "
                xP = [P3[0], P4[0]]
                yP = [P3[1], P4[1]]
                pylab.scatter(xP,yP, color='g')        
    
                
            pylab.scatter([globalJunctionPoint[0]],[globalJunctionPoint[1]], color='r')        
    
            xP = []
            yP = []
            for p in path1:
                xP.append(p[0])
                yP.append(p[1])
            pylab.plot(xP,yP, color=(1.0,0.5,0.5))
    
            #xP = []
            #yP = []
            #for p in globalBranchPoints:
            #    xP.append(p[0])
            #    yP.append(p[1])
            #pylab.plot(xP,yP, color='g')
    
            pylab.xlim(-5,10)
            pylab.ylim(-8,8)
            pylab.title("%d %d %d %d %d %3.2f %3.2f %3.2f %d %3.2f %3.2f %3.2f" % ( isExist1, isExist2, isInterior1, isInterior2, matchCount1, overlapSum1, angleSum1, juncDist1, matchCount2, overlapSum2, angleSum2, juncDist2))
            pylab.savefig("trimDeparture_%04u.png" % self.pathPlotCount)
            
            self.pathPlotCount += 1

        secP1 = []
        secP2 = []

        """
        cases
        1) near-zero length of path section, never departs, not the path, high juncDist
        2) near-zero length, no match count, low juncDist, correct path
        3) high path length, high match count, low juncDist, incorrect path 
        4) no match count, high junc dist, incorrect
        5) no match count, low junc dist, correct
        """
        
        if juncDist1 < juncDist2:
            
            " FIXME: if [0,frontDepI] has high overlap compared to total length, then we reject it "
            #if matchCount1 > 0 and frontDepI :
            #    pass
            
            
            secP1 = pathPoints2[0]
            secP2 = pathPoints2[frontDepI]
        else:
            " FIXME: if [backDepI,len(distances)-1] has high overlap compared to total length, then we reject it "
            secP1 = pathPoints2[backDepI]
            secP2 = pathPoints2[len(distances)-1]

        if len(secP1) == 0:
            print "no departures found"
            raise
        

            
        
        return secP1, secP2



    " get the trimmed version of child and parent paths that are overlapping in some fashion "
    def getBranchPoint(self, parentPathID, childPathID, paths):

        "Assumption:  one section of the medial axis is closely aligned with the path "        
            
        plotIter = False
        print "getOverlapDeparture():"
        
        path1 = paths[parentPathID]
        path2 = paths[childPathID]
        
        
        isExist1 = False
        isInterior1 = False
        departurePoint1 = 0
        isExist2 = False
        isInterior2 = False
        departurePoint2 = 0


        branchNodeID = self.pathClasses[childPathID]["branchNodeID"]
        poseOrigin = Pose(self.nodeHash[branchNodeID].getGlobalGPACPose())
        localJunctionPoint = self.pathClasses[childPathID]["localJunctionPose"]
        
        
        poseOrigin2 = Pose(self.nodeHash[branchNodeID].getEstPose())
        globalJunctionPoint = poseOrigin2.convertLocalToGlobal(localJunctionPoint)
        
        
        " orienting the medial axis of the branch node correctly "
        #branchDir = self.pathParents[childPathID][3]
        
        branchHull, branchMedial = computeHullAxis(branchNodeID, self.nodeHash[branchNodeID], tailCutOff = False)

        branchMedialSpline = SplineFit(branchMedial, smooth=0.1)
        branchPoints = branchMedialSpline.getUniformSamples()
        
        globalBranchPoints = []
        for p in branchPoints:
            result = poseOrigin.convertLocalToGlobal(p)
            globalBranchPoints.append(result)
        
        #if not branchDir:
        #    " backward "
        #    globalBranchPoints.reverse()
            
        globalBranchSpline = SplineFit(globalBranchPoints, smooth = 0.1)

            

        " return exception if we receive an invalid path "        
        if len(path1) == 0:
            print "path1 has zero length"
            raise
            #return departurePoint1, isInterior1, isExist1, departurePoint2, isInterior2, isExist2


        " make sure the overlap of both paths are oriented the same way "
        path1Spline = SplineFit(path1, smooth=0.1)
        path2Spline = SplineFit(path2, smooth=0.1)

        path2Reverse = deepcopy(path2)
        path2Reverse.reverse()
        path2SplineReverse = SplineFit(path2Reverse, smooth=0.1)
        
        overlapMatch = []
        angleSum1 = 0.0
        angleSum2 = 0.0
        for i in range(0,len(path1)):
            p_1 = path1[i]
            p_2, j, minDist = gen_icp.findClosestPointInA(path2, p_1)
            if minDist < 0.1:
                overlapMatch.append((i,j,minDist))

                pathU1 = path1Spline.findU(p_1)    
                pathU2 = path2Spline.findU(p_2)    
                pathU2_R = path2SplineReverse.findU(p_2)    

                pathVec1 = path1Spline.getUVector(pathU1)
                pathVec2 = path2Spline.getUVector(pathU2)
                pathVec2_R = path2SplineReverse.getUVector(pathU2_R)

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
        print "angle sums for orienting child path with parent path:", angleSum1, angleSum2

        if angleSum1 > angleSum2:
            orientedPath2 = path2Reverse
        else:
            orientedPath2 = path2
            
        path2Spline = SplineFit(orientedPath2, smooth=0.1)



        " for each point on the child path, find its closest pair on the parent path "
        
        pathPoints1 = path1Spline.getUniformSamples()
        pathPoints2 = path2Spline.getUniformSamples()

        distances = []
        indices = []
        for i in range(0,len(pathPoints2)):
            p_2 = pathPoints2[i]
            p_1, i_1, minDist = gen_icp.findClosestPointInA(pathPoints1, p_2)
            
            " keep the distance information "
            distances.append(minDist)
            " and the associated index of the point on the parent path "
            indices.append(i_1)

        " match distances of the tip points of child path "        
        maxFront = distances[0]
        maxBack = distances[-1]
        print "match distances of tip points:", maxFront, maxBack

        " walk back from tip point until we have a non-monotic increase in match distance "
        " this becomes our departure point "
        
        "TODO:  why is there a 3-offset in this comparison? "
        currI = 1
        try:
            while distances[currI+3] < maxFront:
                maxFront = distances[currI]
                currI += 1
        except:
            pass
        
        " departure index on child path "
        frontDepI = currI
        
        " departure index of child path and match distance "
        frontPoint = [frontDepI, distances[frontDepI]]

        " FIXME:  index out of bounds case "
        currI = 2
        try:
            while distances[-currI-3] < maxBack:
                maxBack = distances[-currI]
                currI += 1
        except:
            pass

        " departure index on child path "
        backDepI = len(distances) - currI

        " departure index of child path and match distance "
        backPoint = [backDepI, distances[backDepI]]

        print "lengths of parent and child paths:", len(pathPoints1), len(pathPoints2)
        print "non monotonic departures:", maxFront, maxBack, frontPoint, backPoint


        "reset to the tip match distance "
        maxFront = distances[0]
        maxBack = distances[-1]
        
        
        " count the number of times a matched point on the parent path is used "
        foo = indices[0:frontDepI+1]
        d1 = {}
        for i in set(foo):
            d1[i] = foo.count(i)

        foo = indices[backDepI:]
        d2 = {}
        for i in set(foo):
            d2[i] = foo.count(i)

        " select the match point that is used the most on the parent path "
        max1 = max(d1, key=d1.get)
        max2 = max(d2, key=d2.get)


        " NOTE:  This guarantees detection of a departure.  We've assumed that there exists a"
        " departure between these two paths "
        
        " compute two candidate departure points "
        
        " NOTE:  departure point on parent path is not used here "
        if True:

            #departurePoint1 = pathPoints1[max1]
            isExist1 = True

            " if the highest match point is either of the tips of the parent path, then this is an external departure "
            if max1 == 0 or max1 == len(pathPoints1)-1:
                isInterior1 = False
            else:
                isInterior1 = True

        if True:
            
            #departurePoint2 = pathPoints1[max2]
            isExist2 = True

            " if the highest match point is either of the tips of the parent path, then this is an external departure "
            if max2 == 0 or max2 == len(pathPoints1)-1:
                isInterior2 = False
            else:
                isInterior2 = True

        print "isExist1 =", isExist1, "isInterior1 =", isInterior1
        print "isExist2 =", isExist2, "isInterior2 =", isInterior2

        " sum of closest points on front and back "
        " select the one with minimal cost "        
        
        " now we compare our candidates to the known direction and location of the branching point "
        angleSum1 = 0.0
        overlapSum1 = 0.0
        matchCount1 = 0
        
        " path section for our front departure hypothesis "
        pathSec1 = pathPoints2[:frontDepI+1]
        pathSec1.reverse()
        
        if len(pathSec1) > 5:
            pathSec1Spline = SplineFit(pathSec1, smooth = 0.1)        
            
            
            for i in range(0,len(pathSec1)):
                p_1 = pathSec1[i]
                p_2, j, minDist = gen_icp.findClosestPointInA(globalBranchPoints, p_1)
                if minDist < 0.1:
                    overlapMatch.append((i,j,minDist))
                    matchCount1 += 1
                    
                    overlapSum1 += minDist
    
                    pathU1 = pathSec1Spline.findU(p_1)    
                    pathUB = globalBranchSpline.findU(p_2)    
    
                    pathVec1 = pathSec1Spline.getUVector(pathU1)
                    pathVecB = globalBranchSpline.getUVector(pathUB)
    
                    val1 = pathVec1[0]*pathVecB[0] + pathVec1[1]*pathVecB[1]
                    if val1 > 1.0:
                        val1 = 1.0
                    elif val1 < -1.0:
                        val1 = -1.0
                    ang1 = acos(val1)
                    
                    angleSum1 += ang1

            if matchCount1 > 0:
                angleSum1 /= matchCount1
                overlapSum1 /= matchCount1
        

        angleSum2 = 0.0
        overlapSum2 = 0.0
        matchCount2 = 0
        
        " path section for our back departure hypothesis "
        pathSec2 = pathPoints2[backDepI:]

        if len(pathSec2) > 5:
            pathSec2Spline = SplineFit(pathSec2, smooth = 0.1)        
    
    
            for i in range(0,len(pathSec2)):
                p_1 = pathSec2[i]
                p_2, j, minDist = gen_icp.findClosestPointInA(globalBranchPoints, p_1)
                if minDist < 0.1:
                    overlapMatch.append((i,j,minDist))
                    matchCount2 += 1
                    
                    overlapSum2 += minDist
    
                    pathU2 = pathSec2Spline.findU(p_1)    
                    pathUB = globalBranchSpline.findU(p_2)    
    
                    pathVec2 = pathSec2Spline.getUVector(pathU2)
                    pathVecB = globalBranchSpline.getUVector(pathUB)
    
                    val2 = pathVec2[0]*pathVecB[0] + pathVec2[1]*pathVecB[1]
                    if val2 > 1.0:
                        val2 = 1.0
                    elif val2 < -1.0:
                        val2 = -1.0
                    ang2 = acos(val2)
                    
                    angleSum2 += ang2

            if matchCount2 > 0:
                angleSum2 /= matchCount2
                overlapSum2 /= matchCount2

        print "pathSec1 hypothesis angle and overlap sum and match count:", angleSum1, overlapSum1, matchCount1
        print "pathSec2 hypothesis angle and overlap sum and match count:", angleSum2, overlapSum2, matchCount2

        " distance of departure point from known junction point "
        #juncDist1 = -1
        #juncDist2 = -1
        #if len(pathSec1) > 0:
        p0 = pathSec1[0]
        juncDist1 = sqrt((globalJunctionPoint[0]-p0[0])**2 + (globalJunctionPoint[1]-p0[1])**2)

        #if len(pathSec2) > 0:
        p0 = pathSec2[0]
        juncDist2 = sqrt((globalJunctionPoint[0]-p0[0])**2 + (globalJunctionPoint[1]-p0[1])**2)

        print "pathSec1 hypothesis discrepancy distance:", juncDist1
        print "pathSec2 hypothesis discrepancy distance:", juncDist2

        if plotIter:
            pylab.clf()
            xP = []
            yP = []
            for p in path2:
                xP.append(p[0])
                yP.append(p[1])
            pylab.plot(xP,yP, color=(0.5,0.5,1.0))
    
            if True:    
                P1 = pathPoints2[0]
                P2 = pathPoints2[frontDepI]
                
                P3 = pathPoints2[backDepI]
                P4 = pathPoints2[-1]

                " 0, frontDepI "
                xP = [P1[0], P2[0]]
                yP = [P1[1], P2[1]]
                pylab.scatter(xP,yP, color='b')        

                " backDepI, -1 "
                xP = [P3[0], P4[0]]
                yP = [P3[1], P4[1]]
                pylab.scatter(xP,yP, color='g')        
    
                
            pylab.scatter([globalJunctionPoint[0]],[globalJunctionPoint[1]], color='r')        
    
            xP = []
            yP = []
            for p in path1:
                xP.append(p[0])
                yP.append(p[1])
            pylab.plot(xP,yP, color=(1.0,0.5,0.5))
    
            #xP = []
            #yP = []
            #for p in globalBranchPoints:
            #    xP.append(p[0])
            #    yP.append(p[1])
            #pylab.plot(xP,yP, color='g')
    
            pylab.xlim(-5,10)
            pylab.ylim(-8,8)
            pylab.title("%d %d %d %d %d %3.2f %3.2f %3.2f %d %3.2f %3.2f %3.2f" % ( isExist1, isExist2, isInterior1, isInterior2, matchCount1, overlapSum1, angleSum1, juncDist1, matchCount2, overlapSum2, angleSum2, juncDist2))
            pylab.savefig("trimDeparture_%04u.png" % self.pathPlotCount)
            
            self.pathPlotCount += 1

        secP1 = []
        secP2 = []

        """
        cases
        1) near-zero length of path section, never departs, not the path, high juncDist
        2) near-zero length, no match count, low juncDist, correct path
        3) high path length, high match count, low juncDist, incorrect path 
        4) no match count, high junc dist, incorrect
        5) no match count, low junc dist, correct
        """
        
        if juncDist1 < juncDist2:
            
            " FIXME: if [0,frontDepI] has high overlap compared to total length, then we reject it "
            #if matchCount1 > 0 and frontDepI :
            #    pass
            
            
            secP1 = pathPoints2[0]
            secP2 = pathPoints2[frontDepI]
        else:
            " FIXME: if [backDepI,len(distances)-1] has high overlap compared to total length, then we reject it "
            secP1 = pathPoints2[backDepI]
            secP2 = pathPoints2[len(distances)-1]

        if len(secP1) == 0:
            print "no departures found"
            raise
        


        minI_1 = 0        
        minI_2 = 0
        minDist_1 = 1e100
        minDist_2 = 1e100        
        for i in range(len(path2)):
            pnt = path2[i]
            dist1 = sqrt((pnt[0]-secP1[0])**2 + (pnt[1]-secP1[1])**2)
            dist2 = sqrt((pnt[0]-secP2[0])**2 + (pnt[1]-secP2[1])**2)
        
            if dist1 < minDist_1:
                minDist_1 = dist1
                minI_1 = i

            if dist2 < minDist_2:
                minDist_2 = dist2
                minI_2 = i

        term0 = path2[0]
        termN = path2[-1]

        " smallest distance is the terminal point "
        dist0_1 = sqrt((term0[0]-secP1[0])**2 + (term0[1]-secP1[1])**2)
        distN_1 = sqrt((termN[0]-secP1[0])**2 + (termN[1]-secP1[1])**2)
        dist0_2 = sqrt((term0[0]-secP2[0])**2 + (term0[1]-secP2[1])**2)
        distN_2 = sqrt((termN[0]-secP2[0])**2 + (termN[1]-secP2[1])**2)
        print "terminal distances:", dist0_1, distN_1, dist0_2, distN_2
        
        distList = [dist0_1, distN_1, dist0_2, distN_2]
        
        minI = -1
        minDist = 1e100
        for i in range(len(distList)):
            if distList[i] < minDist:
                minDist = distList[i]
                minI = i
        
        #junctionNodeID = self.pathClasses[childPathID]["branchNodeID"]
        #localJunctionPoint = self.pathClasses[childPathID]["localJunctionPose"]
        
        #poseOrigin = Pose(self.nodeHash[junctionNodeID].getEstPose())
        #junctionPoint = poseOrigin.convertLocalToGlobal(localJunctionPoint)
        
        junctionPoint = globalJunctionPoint
        
        
        minDist2 = 1e100
        juncI = 0        
        for i in range(len(path2)):
            pnt = path2[i]
            dist = sqrt((pnt[0]-junctionPoint[0])**2 + (pnt[1]-junctionPoint[1])**2)
        
            if dist < minDist2:
                minDist2 = dist
                juncI = i
         
         
        
        print "len(path1):", len(path1)
        print "len(path2):", len(path2)
        print "juncI:", juncI
        print "minDist:", minDist_1, minDist_2
        
        " now we have closest point to departure point. "
        " Which side is the departing side? "    
        #sum1 = self.getOverlapSum(path1, path2[0:minI2+1] + [junctionPoint])
        #sum2 = self.getOverlapSum(path1, [junctionPoint] + path2[minI2:])
        

                
        if minI == 0:
            "secP1 is terminal 0"
            index = juncI-10
            #index = juncI
            if index < 1:
                index = 1
            #newPath2 = path2[:index+1] + [junctionPoint]
            
            #self.pathClasses[childPathID]["globalJunctionPose"] = path2[index]

            
            newPath2 = path2[:index+1]
            newPath2.reverse()
        
        elif minI == 1:
            "secP1 is terminal N"
            index = juncI+10
            #index = juncI
            if index >= len(path2)-1:
                
                " ensure at least 2 elements in path "
                index = len(path2)-2

            #self.pathClasses[childPathID]["globalJunctionPose"] = path2[index]

            #newPath2 = [junctionPoint] + path2[index:]
            newPath2 = path2[index:]

            
        elif minI == 2:
            "secP2 is terminal 0"
            index = juncI-10
            #index = juncI
            if index < 1:
                index = 1

            #self.pathClasses[childPathID]["globalJunctionPose"] = path2[index]

            #newPath2 = path2[:index+1] + [junctionPoint]
            newPath2 = path2[:index+1]
            newPath2.reverse()
            
        elif minI == 3:
            "secP2 is terminal N"
            index = juncI+10
            #index = juncI
            if index >= len(path2)-1:
                " ensure at least 2 elements in path "
                index = len(path2)-2
            
            #self.pathClasses[childPathID]["globalJunctionPose"] = path2[index]
            
            #newPath2 = [junctionPoint] + path2[index:]
            newPath2 = path2[index:]
        
        else:
            print "no terminal found"
            raise
        
        
        " extrapolate newPath2 to intersect the parent path "
        
        
        
        
        
        
        
        #if minI_1 > minI_2:
        #    path2[minI_2:minI_1+1]
        #else:
        #    path2[minI_1:minI_2+1]
        
        
        #if startI == 0:
        #    newPath2 = path2[startI:endI+1] + [junctionPoint]
        #    newPath2.reverse()
        #else:
        #    newPath2 = [junctionPoint] + path2[startI:endI+1]
            

        #print "Sums:", sum1, sum2
        
        #if sum1 > sum2:
        #    index = minI2+1-10
        #    if index < 1:
        #        index = 1
        #    newPath2 = path2[0:index] + [junctionPoint]
        #    newPath2.reverse()
        #else:
        #    index = minI2+10
        #    if index >= len(path2):
        #        index = len(path2)-1
        #    newPath2 = [junctionPoint] + path2[index:]

        " convert path so that the points are uniformly distributed "
        
        
        max_spacing = 0.08
        #max_spacing = 0.04
        #max_spacing = 0.1
        newPath3 = []

        " make sure path is greater than 5 points "
        while len(newPath3) <= 5:

            #print "len(newPath3) =", len(newPath3), "<= 5"

            max_spacing /= 2
            print "max_spacing =", max_spacing
            
            newPath3 = [copy(newPath2[0])]
                                
            for i in range(len(newPath2)-1):
                p0 = newPath2[i]
                p1 = newPath2[(i+1)]
                dist = sqrt((p0[0]-p1[0])**2 + (p0[1]-p1[1])**2)
    
                vec = [p1[0]-p0[0], p1[1]-p0[1]]
                vec[0] /= dist
                vec[1] /= dist
                
                
                if dist > max_spacing:
                    " cut into pieces max_spacing length or less "
                    numCount = int(floor(dist / max_spacing))
                    
                    for j in range(1, numCount+1):
                        newP = [j*max_spacing*vec[0] + p0[0], j*max_spacing*vec[1] + p0[1]]
                        newPath3.append(newP)

                newPath3.append(copy(p1))            
        
        #trimmedPaths.append(deepcopy(newPath3))
        #trimmedPaths[pathID] = deepcopy(newPath3)






        leafPath = deepcopy(newPath3)
        
        frontVec = [0.,0.]
        backVec = [0.,0.]
        indic = range(3)
        indic.reverse()
        
        for i in indic:
            if i+2 < len(leafPath):
                p1 = leafPath[i+2]
                p2 = leafPath[i]
                vec = [p2[0]-p1[0], p2[1]-p1[1]]
                frontVec[0] += vec[0]
                frontVec[1] += vec[1]
    
    
        frontMag = math.sqrt(frontVec[0]*frontVec[0] + frontVec[1]*frontVec[1])
    
        frontVec[0] /= frontMag
        frontVec[1] /= frontMag

    
        newP1 = (leafPath[0][0] + frontVec[0]*10, leafPath[0][1] + frontVec[1]*10)
    
        leafPath.insert(0,newP1)
        
        medial2 = deepcopy(leafPath)

        " take the long length segments at tips of medial axis"
        edge1 = medial2[0:2]
        
        frontVec = [edge1[0][0]-edge1[1][0], edge1[0][1]-edge1[1][1]]
        frontMag = math.sqrt(frontVec[0]*frontVec[0] + frontVec[1]*frontVec[1])
        
        frontVec[0] /= frontMag
        frontVec[1] /= frontMag
        
        " make a smaller version of these edges "
        newP1 = (edge1[1][0] + frontVec[0]*2, edge1[1][1] + frontVec[1]*2)

        edge1 = [newP1, edge1[1]]

        " find the intersection points with the hull "
        interPoints = []
        for k in range(len(path1)-1):
            hullEdge = [path1[k],path1[k+1]]
            isIntersect1, point1 = Intersect(edge1, hullEdge)
            if isIntersect1:
                interPoints.append(point1)
                break

        
        " replace the extended edges with a termination point at the hull edge "            
        medial2 = medial2[1:]
        
        if isIntersect1:
            medial2.insert(0, point1)


        #juncAng = acos(-frontVec[0])
        

        #frontVec[0], frontVec[1]

        juncAng = acos(-frontVec[0])
        if -frontVec[1] < 0.0:
            juncAng = -juncAng
        
        
        globJuncPose = [medial2[0][0], medial2[0][1], juncAng]
        
        " last junction point is the intersection point "
        
        return globJuncPose


    def pathTermVisited(self, pathID):
        
        self.pathTermsVisited[pathID] = True

    def getPathTermsVisited(self):
        return self.pathTermsVisited

    def resetTerms(self):
        
        for k, v in self.pathTermsVisited.iteritems():
            self.pathTermsVisited[k] = False

    def getPathTerms(self):

        terms = {}

        pathIDs = self.getPathIDs()

        for pathID in pathIDs:

            pID, pathI = self.topDict["t%u" % (pathID+1)]
            path = self.trimmedPaths[pID]
            terms[pathID] = [path[pathI][0],path[pathI][1], 0.0]
            
        print "returning path terms:", terms

        return terms

    def getClosestPath(self, locPose):
        pass
    
        pathIDs = self.getPathIDs()
        
        minDist = 1e100
        minPathID = 0
        for pathID in pathIDs:
            
            pathSeg = self.trimmedPaths[pathID]
            
            for p in pathSeg:
                dist = sqrt((p[0]-locPose[0])**2+(p[1]-locPose[1])**2)
                
                if dist < minDist:
                    minDist = dist
                    minPathID = pathID
        
        return minPathID

    def getOrderedOverlappingPaths(self, nodeID):
         
        ordered1 = self.findSplice2(nodeID)
        return ordered1
        
    def getOrderedOverlappingPaths_old(self, nodeID):

        overlappedPaths = []
        
        "FIXME:  add departure detection and junction continuity checking "
        "        to reject paths that we are only incident to "
    
        print "getOrderedOverlappingPaths(", nodeID, ")"
    
        overlapSums = {}

        pathIDs = self.getPathIDs()    
        for pathID in pathIDs:
            sum1 = self.getOverlapCondition(self.trimmedPaths[pathID], nodeID)
            
            print "overlap sum", pathID, "=", sum1
            
            if sum1 <= 1e10:
                overlappedPaths.append(pathID)
                overlapSums[pathID] = sum1

        print "overlapSums:", overlapSums

        orderedPathIDs = self.getPathOrdering(nodeID, overlappedPaths)

        print "orderedPathIDs:", orderedPathIDs

        departures = []
        interiors = []
        depPoints = []
        contig = []

        for pathID in orderedPathIDs:
            departurePoint1, depAngle2, isInterior1, isExist1, discDist1, departurePoint2, depAngle2, isInterior2, isExist2, discDist2, contigFrac, overlapSum = self.getDeparturePoint(self.trimmedPaths[pathID], nodeID)
            departures.append([isExist1,isExist2])
            interiors.append([isInterior1, isInterior2])
            depPoints.append([departurePoint1, departurePoint2])
            contig.append((contigFrac, overlapSum))
       
        print "departures:", departures
        print "interiors:", interiors
            
        " determine which paths are leaves "
        pathIDs = self.getPathIDs()
        isAParent = {}
        for k in pathIDs:
            isAParent[k] = False
        for k in orderedPathIDs:
            currPath = self.getPath(k)
            currParent = currPath["parentID"]
            if currParent != None:
                isAParent[currParent] = True


        print "isAParent:", isAParent
        
        toBeRemoved = []
        seqSize = len(orderedPathIDs)
        
        " TODO: is this necessary still? "
        " check for proper parent-child relationship and "
        " that there is interior departure between paths "

        
        
        for k in range(seqSize-1):
            pathID1 = orderedPathIDs[k]
            pathID2 = orderedPathIDs[k+1]

            " there should be a interior departure between these paths on the parent path "
            if self.getPath(pathID1)["parentID"] == pathID2:
                pass
                #parentPath = pathID2
                
            elif self.getPath(pathID2)["parentID"] == pathID1:
                pass
                #parentPath = pathID1

                    
            else:
                print "ordered path IDs could not determine relationship", pathID1, pathID2, orderedPathIDs
                raise
        
        """        
        for k in range(seqSize-1):
            pathID1 = orderedPathIDs[k]
            pathID2 = orderedPathIDs[k+1]

            " there should be a interior departure between these paths on the parent path "
            if self.pathParents[pathID1][0] == pathID2:
                parentPath = pathID2

                " if the departure doesn't show, remove the least overlapped path section "
                if not interiors[k+1][0]:
                    if overlapSums[pathID1] > overlapSums[pathID2]:
                        print "remove pathID1", pathID1
                        #toBeRemoved.append(k)
                    else:
                        print "remove pathID2", pathID2
                        #toBeRemoved.append(k+1)
                    
                    if k != 0 and k != seqSize-1:
                        print "paths in the middle violate a departure constraint!"
                        print pathID1, pathID2, orderedPathIDs
                        print departures, interiors
                        print overlapSums
                        print toBeRemoved
                        raise
                
                
            elif self.pathParents[pathID2][0] == pathID1:
                parentPath = pathID1
                if not interiors[k][1]:
                    if overlapSums[pathID1] > overlapSums[pathID2]:
                        print "remove pathID1", pathID1
                        #toBeRemoved.append(k)
                    else:
                        print "remove pathID2", pathID2
                        #toBeRemoved.append(k+1)

                    if k != 0 and k != seqSize-1:
                        print "paths in the middle violate a departure constraint!"
                        print pathID1, pathID2, orderedPathIDs
                        print departures, interiors
                        print overlapSums
                        print toBeRemoved
                        raise
                    
            else:
                print "ordered path IDs could not determine relationship", pathID1, pathID2, orderedPathIDs
                raise
        """
        
        print "toBeRemoved:", toBeRemoved

        finalPathIDs = []
        for k in range(seqSize):
            
            if toBeRemoved.count(k) == 0:
                finalPathIDs.append(orderedPathIDs[k])
            
        print "finalPathIDs:", finalPathIDs    
        ordered1 = self.findSplice(nodeID)
        return ordered1
        
        return finalPathIDs
    
    
    def getPathPath(self, startPathID, endPathID):
        " given start pathID and end pathID, what is the series of paths I need to follow to get from start to end "

        if startPathID == endPathID:
            return [startPathID]
        
        pathGraph = graph.graph()

        pathIDs = self.getPathIDs()

        for pathID in pathIDs:
            pathGraph.add_node(pathID, [])
        
        for pathID in pathIDs:    
            parentPathID = self.getPath(pathID)["parentID"]
            if parentPathID != None:
                pathGraph.add_edge(pathID, parentPathID)
        
        shortestPathSpanTree, shortestDist = pathGraph.shortest_path(endPathID)
        nextPathID = shortestPathSpanTree[startPathID]

        splicedPaths = [startPathID, nextPathID]
        while nextPathID != endPathID:
            nextPathID = shortestPathSpanTree[nextPathID]
            splicedPaths.append(nextPathID)

        return splicedPaths

    def getPathOrdering(self, nodeID, pathIDs):

        plotIter = False

        """
        for pathID in pathIDs:
            
            path = self.paths[pathID]
            hull = self.hulls[pathID]
            overlapCost1 = self.getOverlapCondition(path, nodeID1)
            overlapCost2 = self.getOverlapCondition(path, nodeID2)


            departurePoint1, isInterior1, isExist1 = self.getDeparturePoint(path, nodeID1, hull)
            departurePoint2, isInterior2, isExist2 = self.getDeparturePoint(path, nodeID2, hull)
        """

        node2 = self.nodeHash[nodeID]

        hull2, medial2 = computeHullAxis(nodeID, node2, tailCutOff = False)

        estPose2 = node2.getGlobalGPACPose()        
            
        #minMatchDist2 = 0.2
        #minMatchDist2 = 0.5
        #minMatchDist2 = 2.0
            
        " set the initial guess "
        poseOrigin = Pose(estPose2)
                        
        localizedPaths = []
        for pathID in pathIDs:
            path = self.trimmedPaths[pathID]
            
            localPath = []
            for pnt in path:
                localPath.append(poseOrigin.convertGlobalToLocal(pnt))
            
            localizedPaths.append(localPath)
            
        " find minimum distance and its associated path "                    
        #supportSpline = SplineFit(localSupport, smooth=0.1)        
        #vecPoints1 = supportSpline.getUniformSamples()
        #supportPoints = gen_icp.addGPACVectorCovariance(vecPoints1,high_var=0.05, low_var = 0.001)
            
        medialSpline2 = SplineFit(medial2, smooth=0.1)
        vecPoints2 = medialSpline2.getUniformSamples()
        points2 = gen_icp.addGPACVectorCovariance(vecPoints2,high_var=0.05, low_var = 0.001)

        " transformed points without associated covariance "
        poly2 = []
        for p in points2:
            poly2.append([p[0],p[1]])            
        
        pathSelect = []
        #minMatchDist2 = 0.5
        minMatchDist2 = 1.0

        for i in range(len(vecPoints2)):

            minDist = 1e100
            minPathID = -1
            
            p_2 = vecPoints2[i]

            for j in range(len(localizedPaths)):


                path = localizedPaths[j]
                medialSpline1 = SplineFit(path, smooth=0.1)
                vecPoints1 = medialSpline1.getUniformSamples()
                supportPoints = gen_icp.addGPACVectorCovariance(vecPoints1,high_var=0.05, low_var = 0.001)
    
                " for every transformed point of A, find it's closest neighbor in B "
                try:
                    p_1, minI, minDist_I = gen_icp.findClosestPointWithAngle(vecPoints1, p_2, math.pi/8.0)
        
                    if minDist_I <= minMatchDist2:
                        C2 = points2[i][2]
                        C1 = supportPoints[minI][2]

                        ax = supportPoints[minI][0]
                        ay = supportPoints[minI][1]        
                        bx = points2[i][0]
                        by = points2[i][1]
                
                        c11 = C1[0][0]
                        c12 = C1[0][1]
                        c21 = C1[1][0]
                        c22 = C1[1][1]
                                
                        b11 = C2[0][0]
                        b12 = C2[0][1]
                        b21 = C2[1][0]
                        b22 = C2[1][1]    
                        
                        dist = gen_icp.computeMatchErrorP([0.0,0.0,0.0], [ax,ay], [bx,by], [c11,c12,c21,c22], [b11,b12,b21,b22])

                        if dist < minDist:
                            minDist = dist
                            minPathID = pathIDs[j]
                        #" we store the untransformed point, but the transformed covariance of the A point "
                        #support_pairs.append([points2[i],supportPoints[minI],C2,C1])

                
                except:
                    #raise
                    pass

                #print "minPathID:", minPathID

            if minPathID != -1:
                pathSelect.append(minPathID)
        
        print "pathSelect:", pathSelect
        
        " find the average position of each path ID"
        avgIDPosition = {}
        avgIDCount = {}
        for pathID in pathIDs:
            avgIDPosition[pathID] = 0
            avgIDCount[pathID] = 0
            
        for i in range(len(pathSelect)):
            avgIDPosition[pathSelect[i]] += i
            avgIDCount[pathSelect[i]] += 1


        for pathID in pathIDs:
            if avgIDCount[pathID] != 0:
                avgIDPosition[pathID] /= avgIDCount[pathID]
            else:
                del avgIDPosition[pathID]
            
        " now sort out the order "
        orderedPaths = sorted(avgIDPosition, key=avgIDPosition.__getitem__)        
        
        #print "orderedPaths:", orderedPaths
        
        
        if plotIter:
            pylab.clf()
            
            xP = []
            yP = []
            
            for p in poly2:
                xP.append(p[0])
                yP.append(p[1])
            pylab.plot(xP,yP)
            
            
            for i in range(len(localizedPaths)):
                xP = []
                yP = []
                path = localizedPaths[i]
                for p in path:
                    xP.append(p[0])
                    yP.append(p[1])
                
                pylab.plot(xP,yP)
    
            pylab.title("%d: %s %s" % (nodeID, repr(pathIDs), repr(orderedPaths)))
            pylab.savefig("orderedPath_%04u.png" % self.orderCount)
            
            self.orderCount += 1
        
        return orderedPaths

    def getPathDeparture(self, path1, path2, pathID1, pathID2, termPathIDs1 = [], termPathIDs2 = [], plotIter = False):
        
        isExist1 = False
        isInterior1 = False
        departurePoint1 = 0
        angle1 = 0.0
        isExist2 = False
        isInterior2 = False
        departurePoint2 = 0
        angle2 = 0.0
        
        if len(path1) == 0:
            return departurePoint1, angle1, isInterior1, isExist1, 0.0, departurePoint2, angle2, isInterior2, isExist2, 0.0
        
        #node2 = self.nodeHash[nodeID]

        #hull2, medial2 = computeHullAxis(nodeID, node2, tailCutOff = False)
        
        #estPose2 = node2.getGlobalGPACPose()        
        
        "Assumption:  one section of the medial axis is closely aligned with the path "
        #poseOrigin = Pose(estPose2)
        
        
        pathSpline2 = SplineFit(path2, smooth=0.1)
        pathPoints2 = pathSpline2.getUniformSamples()
        print "path2:", path2[0], pathPoints2[0], path2[-1], pathPoints2[-1]

        pathSpline1 = SplineFit(path1, smooth=0.1)
        pathPoints1 = pathSpline1.getUniformSamples()
        print "path1:", path1[0], pathPoints1[0], path1[-1], pathPoints1[-1]


        " tip angles "
        angSum1 = 0.0
        angSum2 = 0.0
        angs1 = []
        angs2 = []
        phi1 = normalizeAngle(pathPoints2[0][2])
        phi2 = normalizeAngle(pathPoints2[-1][2])
        for i in range(10):
            ang1 = normalizeAngle(pathPoints2[i][2]-phi1)
            ang2 = normalizeAngle(pathPoints2[-i-1][2]-phi2)
            angSum1 += ang1
            angSum2 += ang2
            
            angs1.append(ang1+phi1)
            angs2.append(ang2+phi2)

        angle1 = angSum1 / 10.0 + phi1
        angle2 = angSum2 / 10.0 + phi2

        " invert one angle so opposite tips have opposite angles "
        angle1 = normalizeAngle(angle1 + pi)

        #print "ang1:", angle1, angs1
        #print "ang2:", angle2, angs2
        #print "diff:", diffAngle(angle1, angle2)

        distances = []
        indices = []
        for i in range(0,len(pathPoints2)):
            p_2 = pathPoints2[i]
            #p_1, minDist = gen_icp.findClosestPointInB(pathPoints, p_2, [0.0,0.0,0.0])
            p_1, i_1, minDist = gen_icp.findClosestPointInA(pathPoints1, p_2)
            distances.append(minDist)
            indices.append(i_1)
        
        
        " Compute the front and back departure points by finding the inflection point on the distance curve "
        " these indices become frontDepI and backDepI respectively "
        maxFront = distances[0]
        maxBack = distances[-1]

        currI = 1
        try:
            while distances[currI+3] < maxFront:
                maxFront = distances[currI]
                currI += 1
        except:
            pass
        
        frontDepI = currI
        frontPoint = [frontDepI, distances[frontDepI]]

        " FIXME:  index out of bounds case "
        currI = 2
        try:
            while distances[-currI-3] < maxBack:
                maxBack = distances[-currI]
                currI += 1
        except:
            pass

        backDepI = len(distances) - currI
        backPoint = [backDepI, distances[backDepI]]

        "reset to the maximum distances "
        maxFront = distances[0]
        maxBack = distances[-1]
        
        " perform a line fit on the distance curves from the departure point to the tip"    
        #(ar1,br1)= scipy.polyfit(range(0,frontDepI+1),distances[0:frontDepI+1],1)
        #(ar2,br2)= scipy.polyfit(range(backDepI,len(distances)),distances[backDepI:],1)

        #t1 = range(0,frontDepI+1)
        #xr1 = scipy.polyval([ar1,br1],t1)            
        
        #t2 = range(backDepI,len(distances))
        #xr2 = scipy.polyval([ar2,br2],t2)            
        
        " compute the average of distances between departure and termination "
        #frontAvg = 0.0
        #for i in range(0,frontDepI+1):
        #    frontAvg += distances[i]
        #frontAvg /= (frontDepI+1)

        #backAvg = 0.0
        #for i in range(backDepI,len(distances)):
        #    backAvg += distances[i]
        #backAvg /= len(distances)-backDepI

        " for all the points matched from the local curve to the path curve "
        " count the number of times they are matched "
        foo = indices[0:frontDepI+1]
        d1 = {}
        for i in set(foo):
            d1[i] = foo.count(i)

        foo = indices[backDepI:]
        d2 = {}
        for i in set(foo):
            d2[i] = foo.count(i)

        " find the point that has the most matches "
        max1 = max(d1, key=d1.get)
        max2 = max(d2, key=d2.get)



        arr1 = numpy.array(deepcopy(indices[0:frontDepI+1]))
        arr2 = numpy.array(deepcopy(indices[backDepI:]))
        matchVar1 = arr1.var()
        matchVar2 = arr2.var()

        " compute the selected point index by average instead of by maximum"
        #frontI = 0.0
        #for i in range(0,frontDepI+1):
        #    frontI += indices[i]
        #frontI /= (frontDepI+1)

        #backI = 0.0
        #for i in range(backDepI,len(distances)):
        #    backI += indices[i]
        #backI /= len(distances)-backDepI

        tipMatch1 = pathPoints1[indices[0]]

        " discrepancy distance between tip closest point and average departure point "
        dist1 = sqrt((tipMatch1[0]-pathPoints1[max1][0])**2 + (tipMatch1[1] - pathPoints1[max1][1])**2)

        tipMatch2 = pathPoints1[indices[-1]]

        " discrepancy distance between tip closest point and average departure point "
        dist2 = sqrt((tipMatch2[0]-pathPoints1[max2][0])**2 + (tipMatch2[1] - pathPoints1[max2][1])**2)


        DEP_THRESH = 0.1
        #DEP_THRESH = 0.3

        if maxFront > DEP_THRESH:
            departurePoint1 = pathPoints1[max1]
            isExist1 = True



            if max1 == 0 or max1 == len(pathPoints1)-1:
                isInterior1 = False
            else:
                isInterior1 = True

            # Test purposes only
            #if nodeID == 4 or nodeID == 5:
            #    isInterior1 = True


        if maxBack > DEP_THRESH:
            departurePoint2 = pathPoints1[max2]
            isExist2 = True


            if max2 == 0 or max2 == len(pathPoints1)-1:
                isInterior2 = False
            else:
                isInterior2 = True

        " sum of closest points on front and back "
        " select the one with minimal cost "
        
        if False:
            pylab.clf()
            xP = range(len(pathPoints2))
            yP = distances
            pylab.plot(xP,yP, color ='b')
            
            if maxFront > 0.5:
                xP = [frontPoint[0]]
                yP = [frontPoint[1]]
                pylab.scatter(xP,yP, color='b')
    
            if maxBack > 0.5:
                xP = [backPoint[0]]
                yP = [backPoint[1]]
                pylab.scatter(xP,yP, color='b')
            
            pylab.xlim(0,200)
            pylab.ylim(0,2)
            #pylab.title("%1.2f,%1.2f,%1.2f,%1.2f,%1.2f,%1.2f,%d,%d,%d,%d" % (frontSum,backSum,frontAvg,backAvg,frontI,backI,max1,max2,d1[max1],d2[max2]))
            pylab.title("%1.2f %1.2f %d %d %d" % ( maxFront, maxBack, len(pathPoints1), max1, max2))
            pylab.savefig("distances_%04u.png" % self.pathPlotCount)

        if plotIter:
            pylab.clf()
            xP = []
            yP = []
            for p in pathPoints2:
                xP.append(p[0])
                yP.append(p[1])
            pylab.plot(xP,yP, color='b')
    
            if True:    
                xP = [pathPoints1[max1][0]]
                yP = [pathPoints1[max1][1]]
                pylab.scatter(xP,yP, color='b')        
    
                xP = [tipMatch1[0]]
                yP = [tipMatch1[1]]
                pylab.scatter(xP,yP, color='r')        
    
    
            if True:
                xP = [pathPoints1[max2][0]]
                yP = [pathPoints1[max2][1]]
                pylab.scatter(xP,yP, color='b')        
    
                xP = [tipMatch2[0]]
                yP = [tipMatch2[1]]
                pylab.scatter(xP,yP, color='r')        
    
    
            xP = []
            yP = []
            for p in pathPoints1:
                xP.append(p[0])
                yP.append(p[1])
            pylab.plot(xP,yP, color='r')
            
            pylab.xlim(-5,10)
            pylab.ylim(-8,8)
            #pylab.title("nodeID %d: %d %d" % (nodeID, isInterior, isExist))
            pylab.title("%1.2f, %1.2f, %1.2f, %1.2f, %1.2f, %1.2f, %1.2f, %1.2f, [%d,%d] [%d,%d], (%d,%d)" % (maxFront, maxBack, dist1, dist2, matchVar1, matchVar2, angle1, angle2, isExist1, isExist2, isInterior1, isInterior2, pathID1, pathID2))
            #pylab.title("%s,%s,%s,%s" % (repr(pathPoints1[0]),repr(pathPoints1[-1]),repr(pathPoints2[0]),repr(pathPoints2[-1])))
            #pylab.title("(%1.2f,%1.2f),(%1.2f,%1.2f),(%1.2f,%1.2f),(%1.2f,%1.2f)" % (pathPoints1[0][0],pathPoints1[0][1],pathPoints1[-1][0],pathPoints1[-1][1],pathPoints2[0][0],pathPoints2[0][1],pathPoints2[-1][0],pathPoints2[-1][1]))
            pylab.savefig("pathDeparture_%04u.png" % self.pathPlotCount)
            
            self.pathPlotCount += 1
        
        
        #return departurePoint, isInterior, isExist
        return departurePoint1, angle1, isInterior1, isExist1, dist1, departurePoint2, angle2, isInterior2, isExist2, dist2

    def getPathDeparture2(self, path1, path2, pathID1, pathID2, termPathIDs1 = [], termPathIDs2 = []):
        
        isExist1 = False
        isInterior1 = False
        departurePoint1 = 0
        angle1 = 0.0
        isExist2 = False
        isInterior2 = False
        departurePoint2 = 0
        angle2 = 0.0
        
        if len(path1) == 0:
            return departurePoint1, angle1, isInterior1, isExist1, 0.0, departurePoint2, angle2, isInterior2, isExist2, 0.0
        
        #node2 = self.nodeHash[nodeID]

        #hull2, medial2 = computeHullAxis(nodeID, node2, tailCutOff = False)
        
        #estPose2 = node2.getGlobalGPACPose()        
        
        "Assumption:  one section of the medial axis is closely aligned with the path "
        #poseOrigin = Pose(estPose2)
        
        
        pathSpline2 = SplineFit(path2, smooth=0.1)
        pathPoints2 = pathSpline2.getUniformSamples()
        print "path2:", path2[0], pathPoints2[0], path2[-1], pathPoints2[-1]

        pathSpline1 = SplineFit(path1, smooth=0.1)
        pathPoints1 = pathSpline1.getUniformSamples()
        print "path1:", path1[0], pathPoints1[0], path1[-1], pathPoints1[-1]


        " tip angles "
        angSum1 = 0.0
        angSum2 = 0.0
        angs1 = []
        angs2 = []
        phi1 = normalizeAngle(pathPoints2[0][2])
        phi2 = normalizeAngle(pathPoints2[-1][2])
        for i in range(10):
            ang1 = normalizeAngle(pathPoints2[i][2]-phi1)
            ang2 = normalizeAngle(pathPoints2[-i-1][2]-phi2)
            angSum1 += ang1
            angSum2 += ang2
            
            angs1.append(ang1+phi1)
            angs2.append(ang2+phi2)

        angle1 = angSum1 / 10.0 + phi1
        angle2 = angSum2 / 10.0 + phi2

        " invert one angle so opposite tips have opposite angles "
        angle1 = normalizeAngle(angle1 + pi)

        #print "ang1:", angle1, angs1
        #print "ang2:", angle2, angs2
        #print "diff:", diffAngle(angle1, angle2)

        distances = []
        indices = []
        for i in range(0,len(pathPoints2)):
            p_2 = pathPoints2[i]
            #p_1, minDist = gen_icp.findClosestPointInB(pathPoints, p_2, [0.0,0.0,0.0])
            p_1, i_1, minDist = gen_icp.findClosestPointInA(pathPoints1, p_2)
            distances.append(minDist)
            indices.append(i_1)
        
        
        " Compute the front and back departure points by finding the inflection point on the distance curve "
        " these indices become frontDepI and backDepI respectively "
        maxFront = distances[0]
        maxBack = distances[-1]

        currI = 1
        try:
            while distances[currI+3] < maxFront:
                maxFront = distances[currI]
                currI += 1
        except:
            pass
        
        frontDepI = currI
        frontPoint = [frontDepI, distances[frontDepI]]

        " FIXME:  index out of bounds case "
        currI = 2
        try:
            while distances[-currI-3] < maxBack:
                maxBack = distances[-currI]
                currI += 1
        except:
            pass

        backDepI = len(distances) - currI
        backPoint = [backDepI, distances[backDepI]]

        "reset to the maximum distances "
        maxFront = distances[0]
        maxBack = distances[-1]
        
        " perform a line fit on the distance curves from the departure point to the tip"    
        #(ar1,br1)= scipy.polyfit(range(0,frontDepI+1),distances[0:frontDepI+1],1)
        #(ar2,br2)= scipy.polyfit(range(backDepI,len(distances)),distances[backDepI:],1)

        #t1 = range(0,frontDepI+1)
        #xr1 = scipy.polyval([ar1,br1],t1)            
        
        #t2 = range(backDepI,len(distances))
        #xr2 = scipy.polyval([ar2,br2],t2)            
        
        " compute the average of distances between departure and termination "
        #frontAvg = 0.0
        #for i in range(0,frontDepI+1):
        #    frontAvg += distances[i]
        #frontAvg /= (frontDepI+1)

        #backAvg = 0.0
        #for i in range(backDepI,len(distances)):
        #    backAvg += distances[i]
        #backAvg /= len(distances)-backDepI

        " for all the points matched from the local curve to the path curve "
        " count the number of times they are matched "
        foo = indices[0:frontDepI+1]
        d1 = {}
        for i in set(foo):
            d1[i] = foo.count(i)

        foo = indices[backDepI:]
        d2 = {}
        for i in set(foo):
            d2[i] = foo.count(i)

        " find the point that has the most matches "
        max1 = max(d1, key=d1.get)
        max2 = max(d2, key=d2.get)



        arr1 = numpy.array(deepcopy(indices[0:frontDepI+1]))
        arr2 = numpy.array(deepcopy(indices[backDepI:]))
        matchVar1 = arr1.var()
        matchVar2 = arr2.var()

        " compute the selected point index by average instead of by maximum"
        #frontI = 0.0
        #for i in range(0,frontDepI+1):
        #    frontI += indices[i]
        #frontI /= (frontDepI+1)

        #backI = 0.0
        #for i in range(backDepI,len(distances)):
        #    backI += indices[i]
        #backI /= len(distances)-backDepI

        tipMatch1 = pathPoints1[indices[0]]

        " discrepancy distance between tip closest point and average departure point "
        dist1 = sqrt((tipMatch1[0]-pathPoints1[max1][0])**2 + (tipMatch1[1] - pathPoints1[max1][1])**2)

        tipMatch2 = pathPoints1[indices[-1]]

        " discrepancy distance between tip closest point and average departure point "
        dist2 = sqrt((tipMatch2[0]-pathPoints1[max2][0])**2 + (tipMatch2[1] - pathPoints1[max2][1])**2)


        DEP_THRESH = 0.1
        #DEP_THRESH = 0.3

        if maxFront > DEP_THRESH:
            departurePoint1 = pathPoints1[max1]
            isExist1 = True



            if max1 == 0 or max1 == len(pathPoints1)-1:
                isInterior1 = False
            else:
                isInterior1 = True

            # Test purposes only
            #if nodeID == 4 or nodeID == 5:
            #    isInterior1 = True


        if maxBack > DEP_THRESH:
            departurePoint2 = pathPoints1[max2]
            isExist2 = True


            if max2 == 0 or max2 == len(pathPoints1)-1:
                isInterior2 = False
            else:
                isInterior2 = True

        " sum of closest points on front and back "
        " select the one with minimal cost "
        
        if False:
            pylab.clf()
            xP = range(len(pathPoints2))
            yP = distances
            pylab.plot(xP,yP, color ='b')
            
            if maxFront > 0.5:
                xP = [frontPoint[0]]
                yP = [frontPoint[1]]
                pylab.scatter(xP,yP, color='b')
    
            if maxBack > 0.5:
                xP = [backPoint[0]]
                yP = [backPoint[1]]
                pylab.scatter(xP,yP, color='b')
            
            pylab.xlim(0,200)
            pylab.ylim(0,2)
            #pylab.title("%1.2f,%1.2f,%1.2f,%1.2f,%1.2f,%1.2f,%d,%d,%d,%d" % (frontSum,backSum,frontAvg,backAvg,frontI,backI,max1,max2,d1[max1],d2[max2]))
            pylab.title("%1.2f %1.2f %d %d %d" % ( maxFront, maxBack, len(pathPoints1), max1, max2))
            pylab.savefig("distances_%04u.png" % self.pathPlotCount)

        if False:
            pylab.clf()
            xP = []
            yP = []
            for p in pathPoints2:
                xP.append(p[0])
                yP.append(p[1])
            pylab.plot(xP,yP, color='b')
    
            if True:    
                xP = [pathPoints1[max1][0]]
                yP = [pathPoints1[max1][1]]
                pylab.scatter(xP,yP, color='b')        
    
                xP = [tipMatch1[0]]
                yP = [tipMatch1[1]]
                pylab.scatter(xP,yP, color='r')        
    
    
            if True:
                xP = [pathPoints1[max2][0]]
                yP = [pathPoints1[max2][1]]
                pylab.scatter(xP,yP, color='b')        
    
                xP = [tipMatch2[0]]
                yP = [tipMatch2[1]]
                pylab.scatter(xP,yP, color='r')        
    
    
            xP = []
            yP = []
            for p in pathPoints1:
                xP.append(p[0])
                yP.append(p[1])
            pylab.plot(xP,yP, color='r')
            
            pylab.xlim(-5,10)
            pylab.ylim(-8,8)
            #pylab.title("nodeID %d: %d %d" % (nodeID, isInterior, isExist))
            pylab.title("%1.2f, %1.2f, %1.2f, %1.2f, %1.2f, %1.2f, %1.2f, %1.2f, [%d,%d] [%d,%d], (%d,%d)" % (maxFront, maxBack, dist1, dist2, matchVar1, matchVar2, angle1, angle2, isExist1, isExist2, isInterior1, isInterior2, pathID1, pathID2))
            #pylab.title("%s,%s,%s,%s" % (repr(pathPoints1[0]),repr(pathPoints1[-1]),repr(pathPoints2[0]),repr(pathPoints2[-1])))
            #pylab.title("(%1.2f,%1.2f),(%1.2f,%1.2f),(%1.2f,%1.2f),(%1.2f,%1.2f)" % (pathPoints1[0][0],pathPoints1[0][1],pathPoints1[-1][0],pathPoints1[-1][1],pathPoints2[0][0],pathPoints2[0][1],pathPoints2[-1][0],pathPoints2[-1][1]))
            pylab.savefig("pathDeparture2_%04u.png" % self.pathPlotCount2)
            
            self.pathPlotCount2 += 1
        
        
        #return departurePoint, isInterior, isExist
        return departurePoint1, angle1, isInterior1, isExist1, dist1, departurePoint2, angle2, isInterior2, isExist2, dist2


    def getDeparturePoint(self, currPath, nodeID, plotIter = False):
        
        isExist1 = False
        isInterior1 = False
        departurePoint1 = 0
        angle1 = 0.0
        isExist2 = False
        isInterior2 = False
        departurePoint2 = 0
        angle2 = 0.0
        
        if len(currPath) == 0:
            return departurePoint1, angle1, isInterior1, isExist1, 0.0, departurePoint2, angle2, isInterior2, isExist2, 0.0
        
        node2 = self.nodeHash[nodeID]
        
        hull2, medial2 = computeHullAxis(nodeID, node2, tailCutOff = False)
        
        estPose2 = node2.getGlobalGPACPose()        
        
        "Assumption:  one section of the medial axis is closely aligned with the path "
        poseOrigin = Pose(estPose2)
        
        medialSpline2 = SplineFit(medial2, smooth=0.1)
        points2 = medialSpline2.getUniformSamples(interpAngle=True)

        points2_offset = []
        for p in points2:
            result = poseOrigin.convertLocalOffsetToGlobal(p)
            points2_offset.append(result)


        globalMedialSpline = SplineFit(points2_offset, smooth=0.1)
        globalMedialPoints = globalMedialSpline.getUniformSamples()



        pathSpline = SplineFit(currPath, smooth=0.1)
        pathPoints = pathSpline.getUniformSamples()

        currPathReverse = deepcopy(currPath)
        currPathReverse.reverse()
        pathSplineReverse = SplineFit(currPathReverse, smooth=0.1)
        pathPointsReverse = pathSplineReverse.getUniformSamples()

        overlapMatch = []
        angleSum1 = 0.0
        angleSum2 = 0.0
        for i in range(0,len(points2_offset)):
            p_1 = points2_offset[i]
            p_2, j, minDist = gen_icp.findClosestPointInA(pathPoints, p_1)
            if minDist < 0.5:
                overlapMatch.append((i,j,minDist))

                pathU1 = globalMedialSpline.findU(p_1)    
                pathU2 = pathSpline.findU(p_2)    
                pathU2_R = pathSplineReverse.findU(p_2)    

                pathVec1 = globalMedialSpline.getUVector(pathU1)
                pathVec2 = pathSpline.getUVector(pathU2)
                pathVec2_R = pathSplineReverse.getUVector(pathU2_R)

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
        print "getDeparturePoint:", i, "angleSum1 =", angleSum1, "angleSum2 =", angleSum2
        if angleSum1 > angleSum2:
            pathPoints = pathPointsReverse
            pathSpline = pathSplineReverse
                    


        
        
        
        
        
        
        
        
  

        " tip angles "
        angSum1 = 0.0
        angSum2 = 0.0
        angs1 = []
        angs2 = []
        phi1 = normalizeAngle(points2_offset[0][2])
        phi2 = normalizeAngle(points2_offset[-1][2])
        for i in range(10):
            ang1 = normalizeAngle(points2_offset[i][2]-phi1)
            ang2 = normalizeAngle(points2_offset[-i-1][2]-phi2)
            angSum1 += ang1
            angSum2 += ang2
            
            angs1.append(ang1+phi1)
            angs2.append(ang2+phi2)

        angle1 = angSum1 / 10.0 + phi1
        angle2 = angSum2 / 10.0 + phi2

        " invert one angle so opposite tips have opposite angles "
        angle1 = normalizeAngle(angle1 + pi)

        print "ang1:", angle1, angs1
        print "ang2:", angle2, angs2
        print "diff:", diffAngle(angle1, angle2)
        
        distSum = 0.0
        contigCount = 0
        maxContig = 0
        distances = []
        indices = []
        for i in range(0,len(points2_offset)):
            p_2 = points2_offset[i]
            #p_1, minDist = gen_icp.findClosestPointInB(pathPoints, p_2, [0.0,0.0,0.0])
            p_1, i_1, minDist = gen_icp.findClosestPointInA(pathPoints, p_2)
            distances.append(minDist)
            indices.append(i_1)
            distSum += minDist
            
            if minDist < 0.2:
                contigCount += 1
                if contigCount > maxContig:
                    maxContig = contigCount
            else:
                contigCount = 0
        
        overlapSum = distSum / float(len(points2_offset))
        
        print "maxContig,overlapSum:", maxContig, overlapSum
        
        " Compute the front and back departure points by finding the inflection point on the distance curve "
        " these indices become frontDepI and backDepI respectively "
        maxFront = distances[0]
        maxBack = distances[-1]

        currI = 1
        try:
            while distances[currI+3] < maxFront:
                maxFront = distances[currI]
                currI += 1
        except:
            pass
        
        frontDepI = currI
        frontPoint = [frontDepI, distances[frontDepI]]


        frontAngleRefI = frontDepI
        while distances[frontAngleRefI] < 0.1:
            frontAngleRefI -= 1
            
            if frontAngleRefI < 0:
                frontAngleRefI = 0
                break
                        
        forePathIndex = indices[frontAngleRefI]
        forePathAngle = pathPoints[forePathIndex][2]
        forePathAngle = normalizeAngle(forePathAngle + pi)

        foreDiffAngle = diffAngle(angle1, forePathAngle)
        
        newFrontDepI = 4
        while fabs(diffAngle(normalizeAngle(points2_offset[newFrontDepI][2]+pi), forePathAngle)) > pi/6.0:
            
            if newFrontDepI < frontDepI:
                newFrontDepI += 1
            else:
                break        
        
        print "front:", nodeID, frontDepI, distances[frontDepI], frontAngleRefI, forePathIndex, forePathAngle, angle1, foreDiffAngle, newFrontDepI, diffAngle(normalizeAngle(points2_offset[newFrontDepI][2]+pi), forePathAngle) 
        

        " FIXME:  index out of bounds case "
        currI = 2
        try:
            while distances[-currI-3] < maxBack:
                maxBack = distances[-currI]
                currI += 1
        except:
            pass

        backDepI = len(distances) - currI
        backPoint = [backDepI, distances[backDepI]]

        backAngleRefI = backDepI
        while distances[backAngleRefI] < 0.1:
            backAngleRefI += 1
            
            if backAngleRefI >= len(distances):
                backAngleRefI = len(distances)-1
                break
        
        backPathIndex = indices[backAngleRefI]
        backPathAngle = pathPoints[backPathIndex][2]
        
        backDiffAngle = diffAngle(angle2, backPathAngle)
        
        newBackDepI = len(distances)-5
        while fabs(diffAngle(points2_offset[newBackDepI][2], backPathAngle)) > pi/6.0:
            
            if newBackDepI > backDepI:
                newBackDepI -= 1
            else:
                break        

        print "back:", nodeID, backDepI, distances[backDepI], backAngleRefI, backPathIndex, backPathAngle, angle2, backDiffAngle, newBackDepI, diffAngle(points2_offset[newBackDepI][2], backPathAngle) 


        foreAngDiffs = []
        #for i in range(0,frontDepI+1):
        for i in range(0,len(points2_offset)):
            foreAngDiffs.append(fabs(diffAngle(normalizeAngle(points2_offset[i][2] + pi), forePathAngle)))

        backAngDiffs = []
        #for i in range(backDepI,len(points2_offset)):
        for i in range(0,len(points2_offset)):
            backAngDiffs.append(fabs(diffAngle(points2_offset[i][2], backPathAngle)))



        "reset to the maximum distances "
        maxFront = distances[0]
        maxBack = distances[-1]
        
        " perform a line fit on the distance curves from the departure point to the tip"    
        #(ar1,br1)= scipy.polyfit(range(0,frontDepI+1),distances[0:frontDepI+1],1)
        #(ar2,br2)= scipy.polyfit(range(backDepI,len(distances)),distances[backDepI:],1)

        #t1 = range(0,frontDepI+1)
        #xr1 = scipy.polyval([ar1,br1],t1)            
        
        #t2 = range(backDepI,len(distances))
        #xr2 = scipy.polyval([ar2,br2],t2)            
        
        " compute the average of distances between departure and termination "
        #frontAvg = 0.0
        #for i in range(0,frontDepI+1):
        #    frontAvg += distances[i]
        #frontAvg /= (frontDepI+1)

        #backAvg = 0.0
        #for i in range(backDepI,len(distances)):
        #    backAvg += distances[i]
        #backAvg /= len(distances)-backDepI

        " for all the points matched from the local curve to the path curve "
        " count the number of times they are matched "
        foo = indices[0:frontDepI+1]
        d1 = {}
        for i in set(foo):
            d1[i] = foo.count(i)

        foo = indices[backDepI:]
        d2 = {}
        for i in set(foo):
            d2[i] = foo.count(i)

        " find the point that has the most matches "
        max1 = max(d1, key=d1.get)
        max2 = max(d2, key=d2.get)



        arr1 = numpy.array(deepcopy(indices[0:frontDepI+1]))
        arr2 = numpy.array(deepcopy(indices[backDepI:]))
        matchVar1 = arr1.var()
        matchVar2 = arr2.var()

        " compute the selected point index by average instead of by maximum"
        #frontI = 0.0
        #for i in range(0,frontDepI+1):
        #    frontI += indices[i]
        #frontI /= (frontDepI+1)

        #backI = 0.0
        #for i in range(backDepI,len(distances)):
        #    backI += indices[i]
        #backI /= len(distances)-backDepI

        tipMatch1 = pathPoints[indices[0]]

        " discrepancy distance between tip closest point and average departure point "
        dist1 = sqrt((tipMatch1[0]-pathPoints[max1][0])**2 + (tipMatch1[1] - pathPoints[max1][1])**2)

        tipMatch2 = pathPoints[indices[-1]]

        " discrepancy distance between tip closest point and average departure point "
        dist2 = sqrt((tipMatch2[0]-pathPoints[max2][0])**2 + (tipMatch2[1] - pathPoints[max2][1])**2)


        DEP_THRESH = 0.3

        if maxFront > DEP_THRESH:
            departurePoint1 = pathPoints[max1]
            
            departurePoint1 = points2_offset[newFrontDepI]
            
            isExist1 = True

            if max1 == 0 or max1 == len(pathPoints)-1:
                isInterior1 = False
            else:
                isInterior1 = True

        if maxBack > DEP_THRESH:
            departurePoint2 = pathPoints[max2]

            departurePoint2 = points2_offset[newBackDepI]

            isExist2 = True
            
            if max2 == 0 or max2 == len(pathPoints)-1:
                isInterior2 = False
            else:
                isInterior2 = True

        """
        #p_1, depI1, pDist1 = gen_icp.findClosestPointInA(points2_offset, pathPoints[max1])
        #p_2, depI2, pDist2 = gen_icp.findClosestPointInA(points2_offset, pathPoints[max2])

        #p_1, depI1, pDist1 = gen_icp.findClosestPointInA(medial2, points2[depI1])
        #p_2, depI2, pDist2 = gen_icp.findClosestPointInA(medial2, points2[depI2])
        """
            
        """
        if maxFront > 0.5 and maxFront > maxBack:
            departurePoint = pathPoints[max1]
            isExist = True

            if max1 == 0 or max1 == len(pathPoints)-1:
                isInterior = False
            else:
                isInterior = True

        if maxBack > 0.5 and maxBack > maxFront:
            departurePoint = pathPoints[max2]
            isExist = True

            if max2 == 0 or max2 == len(pathPoints)-1:
                isInterior = False
            else:
                isInterior = True
        """    
            
        " sum of closest points on front and back "
        " select the one with minimal cost "
        
        if False:
            pylab.clf()
            xP = range(len(points2_offset))
            yP = distances
            pylab.plot(xP,yP, color ='b')
            
            yP = foreAngDiffs
            pylab.plot(xP,yP, color ='r')

            yP = backAngDiffs
            pylab.plot(xP,yP, color ='g')
            
            if maxFront > 0.5:
                xP = [frontPoint[0]]
                yP = [frontPoint[1]]
                pylab.scatter(xP,yP, color='b')
    
            if maxBack > 0.5:
                xP = [backPoint[0]]
                yP = [backPoint[1]]
                pylab.scatter(xP,yP, color='b')
            
            pylab.xlim(0,200)
            pylab.ylim(0,2)
            #pylab.title("%1.2f,%1.2f,%1.2f,%1.2f,%1.2f,%1.2f,%d,%d,%d,%d" % (frontSum,backSum,frontAvg,backAvg,frontI,backI,max1,max2,d1[max1],d2[max2]))
            pylab.title("nodeID %d: %1.2f %1.2f %d %d %d %d %1.2f" % (nodeID, maxFront, maxBack, len(pathPoints), max1, max2, maxContig, overlapSum))
            pylab.savefig("distances_%04u.png" % self.pathPlotCount)

        if plotIter:
            pylab.clf()
            xP = []
            yP = []
            for p in points2_offset:
                xP.append(p[0])
                yP.append(p[1])
            pylab.plot(xP,yP, color='b')
    
            if True:    
                xP = [pathPoints[max1][0]]
                yP = [pathPoints[max1][1]]
                pylab.scatter(xP,yP, color='b')        
    
                xP = [tipMatch1[0]]
                yP = [tipMatch1[1]]
                pylab.scatter(xP,yP, color='r')        

                xP = [points2_offset[frontDepI][0]]
                yP = [points2_offset[frontDepI][1]]
                pylab.scatter(xP,yP, color='g')        

                xP = [pathPoints[forePathIndex][0]]
                yP = [pathPoints[forePathIndex][1]]
                pylab.scatter(xP,yP, color='y')        

                xP = [points2_offset[newFrontDepI][0]]
                yP = [points2_offset[newFrontDepI][1]]
                pylab.scatter(xP,yP, color='m')        
    
    
            if True:
                xP = [pathPoints[max2][0]]
                yP = [pathPoints[max2][1]]
                pylab.scatter(xP,yP, color='b')        
    
                xP = [tipMatch2[0]]
                yP = [tipMatch2[1]]
                pylab.scatter(xP,yP, color='r')        

                xP = [points2_offset[backDepI][0]]
                yP = [points2_offset[backDepI][1]]
                pylab.scatter(xP,yP, color='g')        

                xP = [pathPoints[backPathIndex][0]]
                yP = [pathPoints[backPathIndex][1]]
                pylab.scatter(xP,yP, color='y')        

                xP = [points2_offset[newBackDepI][0]]
                yP = [points2_offset[newBackDepI][1]]
                pylab.scatter(xP,yP, color='m')        

    
            xP = []
            yP = []
            for p in currPath:
                xP.append(p[0])
                yP.append(p[1])
            pylab.plot(xP,yP, color='r')
            
            pylab.xlim(-5,10)
            pylab.ylim(-8,8)
            #pylab.title("nodeID %d: %d %d" % (nodeID, isInterior, isExist))
            pylab.title("%d: %1.2f, %1.2f, %1.2f, %1.2f, %1.2f, %1.2f, %1.2f, %1.2f, [%d,%d] [%d,%d]" % (nodeID, maxFront, maxBack, dist1, dist2, matchVar1, matchVar2, angle1, angle2, isExist1, isExist2, isInterior1, isInterior2))
            pylab.savefig("departure_%04u.png" % self.pathPlotCount)
            
            self.pathPlotCount += 1
        
        print "departure %d: %1.2f, %1.2f, %1.2f, %1.2f, %1.2f, %1.2f, %1.2f, %1.2f, [%d,%d] [%d,%d], [%d,%d]" % (nodeID, maxFront, maxBack, dist1, dist2, matchVar1, matchVar2, angle1, angle2, isExist1, isExist2, isInterior1, isInterior2, frontDepI, backDepI)
        
        
        " if the medial axis does not overlap the path contiguously enough, mark a high discrepancy "
        #if float(maxContig)/float(len(points2_offset)) < 0.4:
        #    dist1 = 1e100
        #    dist2 = 1e100
        
        maxContig, overlapSum
        contigFrac = float(maxContig)/float(len(points2_offset))
        print "returning:", contigFrac, overlapSum
        #return departurePoint, isInterior, isExist
        #return departurePoint1, angle1, isInterior1, isExist1, dist1, departurePoint2, angle2, isInterior2, isExist2, dist2
        return departurePoint1, angle1, isInterior1, isExist1, dist1, departurePoint2, angle2, isInterior2, isExist2, dist2, contigFrac, overlapSum
        

    " returns the endpoints of a node medial axis overlap path "
    def getOverlapBoundaries(self, currPath, nodeID):
        
        isExist1 = False
        isInterior1 = False
        departurePoint1 = 0
        angle1 = 0.0
        isExist2 = False
        isInterior2 = False
        departurePoint2 = 0
        angle2 = 0.0
        
        if len(currPath) == 0:
            return departurePoint1, angle1, isInterior1, isExist1, 0.0, departurePoint2, angle2, isInterior2, isExist2, 0.0
        
        node2 = self.nodeHash[nodeID]
        
        hull2, medial2 = computeHullAxis(nodeID, node2, tailCutOff = False)
        
        estPose2 = node2.getGlobalGPACPose()        
        
        "Assumption:  one section of the medial axis is closely aligned with the path "
        poseOrigin = Pose(estPose2)
        
        medialSpline2 = SplineFit(medial2, smooth=0.1)
        points2 = medialSpline2.getUniformSamples()

        pathSpline = SplineFit(currPath, smooth=0.1)
        pathPoints = pathSpline.getUniformSamples()

        
        points2_offset = []
        for p in points2:
            result = poseOrigin.convertLocalOffsetToGlobal(p)
            points2_offset.append(result)


        " tip angles "
        angSum1 = 0.0
        angSum2 = 0.0
        angs1 = []
        angs2 = []
        phi1 = normalizeAngle(points2_offset[0][2])
        phi2 = normalizeAngle(points2_offset[-1][2])
        for i in range(10):
            ang1 = normalizeAngle(points2_offset[i][2]-phi1)
            ang2 = normalizeAngle(points2_offset[-i-1][2]-phi2)
            angSum1 += ang1
            angSum2 += ang2
            
            angs1.append(ang1+phi1)
            angs2.append(ang2+phi2)

        angle1 = angSum1 / 10.0 + phi1
        angle2 = angSum2 / 10.0 + phi2

        " invert one angle so opposite tips have opposite angles "
        angle1 = normalizeAngle(angle1 + pi)

        print "ang1:", angle1, angs1
        print "ang2:", angle2, angs2
        print "diff:", diffAngle(angle1, angle2)

        distances = []
        indices = []
        for i in range(0,len(points2_offset)):
            p_2 = points2_offset[i]
            #p_1, minDist = gen_icp.findClosestPointInB(pathPoints, p_2, [0.0,0.0,0.0])
            p_1, i_1, minDist = gen_icp.findClosestPointInA(pathPoints, p_2)
            distances.append(minDist)
            indices.append(i_1)
        
        
        " Compute the front and back departure points by finding the inflection point on the distance curve "
        " these indices become frontDepI and backDepI respectively "
        maxFront = distances[0]
        maxBack = distances[-1]

        currI = 1
        try:
            while distances[currI+3] < maxFront:
                maxFront = distances[currI]
                currI += 1
        except:
            pass
        
        frontDepI = currI
        frontPoint = [frontDepI, distances[frontDepI]]

        " FIXME:  index out of bounds case "
        currI = 2
        try:
            while distances[-currI-3] < maxBack:
                maxBack = distances[-currI]
                currI += 1
        except:
            pass

        backDepI = len(distances) - currI
        backPoint = [backDepI, distances[backDepI]]

        "reset to the maximum distances "
        maxFront = distances[0]
        maxBack = distances[-1]
        
        " for all the points matched from the local curve to the path curve "
        " count the number of times they are matched "
        foo = indices[0:frontDepI+1]
        d1 = {}
        for i in set(foo):
            d1[i] = foo.count(i)

        foo = indices[backDepI:]
        d2 = {}
        for i in set(foo):
            d2[i] = foo.count(i)

        " find the point that has the most matches "
        max1 = max(d1, key=d1.get)
        max2 = max(d2, key=d2.get)

        arr1 = numpy.array(deepcopy(indices[0:frontDepI+1]))
        arr2 = numpy.array(deepcopy(indices[backDepI:]))
        matchVar1 = arr1.var()
        matchVar2 = arr2.var()


        tipMatch1 = pathPoints[indices[0]]

        " discrepancy distance between tip closest point and average departure point "
        dist1 = sqrt((tipMatch1[0]-pathPoints[max1][0])**2 + (tipMatch1[1] - pathPoints[max1][1])**2)

        tipMatch2 = pathPoints[indices[-1]]

        " discrepancy distance between tip closest point and average departure point "
        dist2 = sqrt((tipMatch2[0]-pathPoints[max2][0])**2 + (tipMatch2[1] - pathPoints[max2][1])**2)


        DEP_THRESH = 0.3

        if maxFront > DEP_THRESH:
            departurePoint1 = pathPoints[max1]
            isExist1 = True

            if max1 == 0 or max1 == len(pathPoints)-1:
                isInterior1 = False
            else:
                isInterior1 = True

        if maxBack > DEP_THRESH:
            departurePoint2 = pathPoints[max2]
            isExist2 = True


            if max2 == 0 or max2 == len(pathPoints)-1:
                isInterior2 = False
            else:
                isInterior2 = True
        
        " sum of closest points on front and back "
        " select the one with minimal cost "
        
        if False:
            pylab.clf()
            xP = range(len(points2_offset))
            yP = distances
            pylab.plot(xP,yP, color ='b')
            
            if maxFront > 0.5:
                xP = [frontPoint[0]]
                yP = [frontPoint[1]]
                pylab.scatter(xP,yP, color='b')
    
            if maxBack > 0.5:
                xP = [backPoint[0]]
                yP = [backPoint[1]]
                pylab.scatter(xP,yP, color='b')
            
            pylab.xlim(0,200)
            pylab.ylim(0,2)
            #pylab.title("%1.2f,%1.2f,%1.2f,%1.2f,%1.2f,%1.2f,%d,%d,%d,%d" % (frontSum,backSum,frontAvg,backAvg,frontI,backI,max1,max2,d1[max1],d2[max2]))
            pylab.title("nodeID %d: %1.2f %1.2f %d %d %d" % (nodeID, maxFront, maxBack, len(pathPoints), max1, max2))
            pylab.savefig("distances_%04u.png" % self.pathPlotCount)

        if False:
            pylab.clf()
            xP = []
            yP = []
            for p in points2_offset:
                xP.append(p[0])
                yP.append(p[1])
            pylab.plot(xP,yP, color='b')
    
            if True:    
                xP = [pathPoints[max1][0]]
                yP = [pathPoints[max1][1]]
                pylab.scatter(xP,yP, color='b')        
    
                xP = [tipMatch1[0]]
                yP = [tipMatch1[1]]
                pylab.scatter(xP,yP, color='r')        
    
    
            if True:
                xP = [pathPoints[max2][0]]
                yP = [pathPoints[max2][1]]
                pylab.scatter(xP,yP, color='b')        
    
                xP = [tipMatch2[0]]
                yP = [tipMatch2[1]]
                pylab.scatter(xP,yP, color='r')        
    
    
            xP = []
            yP = []
            for p in currPath:
                xP.append(p[0])
                yP.append(p[1])
            pylab.plot(xP,yP, color='r')
            
            pylab.xlim(-5,10)
            pylab.ylim(-8,8)
            #pylab.title("nodeID %d: %d %d" % (nodeID, isInterior, isExist))
            pylab.title("%d: %1.2f, %1.2f, %1.2f, %1.2f, %1.2f, %1.2f, %1.2f, %1.2f, [%d,%d] [%d,%d]" % (nodeID, maxFront, maxBack, dist1, dist2, matchVar1, matchVar2, angle1, angle2, isExist1, isExist2, isInterior1, isInterior2))
            pylab.savefig("departure_%04u.png" % self.pathPlotCount)
            
            self.pathPlotCount += 1
        
        print "departure %d: %1.2f, %1.2f, %1.2f, %1.2f, %1.2f, %1.2f, %1.2f, %1.2f, [%d,%d] [%d,%d]" % (nodeID, maxFront, maxBack, dist1, dist2, matchVar1, matchVar2, angle1, angle2, isExist1, isExist2, isInterior1, isInterior2)
        
        #return departurePoint, isInterior, isExist
        return departurePoint1, angle1, isInterior1, isExist1, dist1, frontDepI, departurePoint2, angle2, isInterior2, isExist2, dist2, backDepI
        

    def getPathOverlapCondition2(self, path1, path2, pathID1, pathID2):

        plotIter = False

        if len(path1) <= 5 or len(path2) <= 5:
            return 0
    
        medial2 = path2
               
        minMatchDist2 = 0.2
        #minMatchDist2 = 0.5
        #minMatchDist2 = 1.0
        #minMatchDist2 = 2.0
                    
    
        supportSpline = SplineFit(path1, smooth=0.1)        
        vecPoints1 = supportSpline.getUniformSamples()
        supportPoints = gen_icp.addGPACVectorCovariance(vecPoints1,high_var=0.05, low_var = 0.001)
        
        medialSpline2 = SplineFit(medial2, smooth=0.1)
        vecPoints2 = medialSpline2.getUniformSamples()
        points2 = gen_icp.addGPACVectorCovariance(vecPoints2,high_var=0.05, low_var = 0.001)

        " transformed points without associated covariance "
        poly2 = []
        for p in points2:
            poly2.append([p[0],p[1]])            
        
        " get the circles and radii "
        #radius2, center2 = gen_icp.computeEnclosingCircle(points2)
                
        support_pairs = []
        #for i in range(len(poly2)):
        for i in range(len(vecPoints2)):
            
            p_2 = vecPoints2[i]
    
            " for every transformed point of A, find it's closest neighbor in B "
            #_1, minDist = gen_icp.findClosestPointInB(supportPoints, p_2, [0.0,0.0,0.0])

            try:
                p_1, minI, minDist = gen_icp.findClosestPointWithAngle(vecPoints1, p_2, math.pi/8.0)
    
                #if gen_icp.isInCircle(p_1, radius2, center2):
    
                if minDist <= minMatchDist2:
                    C2 = points2[i][2]
                    C1 = supportPoints[minI][2]
                    #C1 = p_1[2]
    
                    " we store the untransformed point, but the transformed covariance of the A point "
                    #support_pairs.append([points2[i],p_1,C2,C1])
                    support_pairs.append([points2[i],supportPoints[minI],C2,C1])
            except:
                pass

        cost = 0.0
        if len(support_pairs) == 0:
            cost = 1e100
        else:
            vals = []
            sum1 = 0.0
            for pair in support_pairs:
        
                a = pair[0]
                b = pair[1]
                Ca = pair[2]
                Cb = pair[3]
        
                ax = a[0]
                ay = a[1]        
                bx = b[0]
                by = b[1]
        
                c11 = Ca[0][0]
                c12 = Ca[0][1]
                c21 = Ca[1][0]
                c22 = Ca[1][1]
                        
                b11 = Cb[0][0]
                b12 = Cb[0][1]
                b21 = Cb[1][0]
                b22 = Cb[1][1]    
            
                val = gen_icp.computeMatchErrorP([0.0,0.0,0.0], [ax,ay], [bx,by], [c11,c12,c21,c22], [b11,b12,b21,b22])
                
                vals.append(val)
                sum1 += val
                
            cost = sum1 / len(support_pairs)
            

        if plotIter:
            pylab.clf()
            xP = []
            yP = []
            for p in supportPoints:
                xP.append(p[0])
                yP.append(p[1])
            pylab.plot(xP,yP, color='b')
    
            xP = []
            yP = []
            for p in poly2:
                xP.append(p[0])
                yP.append(p[1])
            pylab.plot(xP,yP, color='r')
            
            for pair in support_pairs:
                p1 = pair[0]
                p2 = pair[1]
                xP = [p1[0],p2[0]]
                yP = [p1[1],p2[1]]
                pylab.plot(xP,yP)
            
            pylab.xlim(-5,10)
            pylab.ylim(-8,8)
            pylab.title("%d %d cost = %f, count = %d" % (pathID1, pathID2, cost, len(support_pairs)))
            pylab.savefig("pathOverlapCost2_%04u.png" % self.overlapPlotCount2)
            self.overlapPlotCount2 += 1


        return len(support_pairs)


    def getPathOverlapCondition(self, path1, path2, pathID1, pathID2, plotIter = False):

        #plotIter = True

        if len(path1) == 0:
            return 1e100
    
        medial2 = path2
               
        minMatchDist2 = 0.2
        #minMatchDist2 = 0.5
        #minMatchDist2 = 1.0
        #minMatchDist2 = 2.0
                    
    
        supportSpline = SplineFit(path1, smooth=0.1)        
        vecPoints1 = supportSpline.getUniformSamples()
        supportPoints = gen_icp.addGPACVectorCovariance(vecPoints1,high_var=0.05, low_var = 0.001)
            
        medialSpline2 = SplineFit(medial2, smooth=0.1)
        vecPoints2 = medialSpline2.getUniformSamples()
        points2 = gen_icp.addGPACVectorCovariance(vecPoints2,high_var=0.05, low_var = 0.001)

        " transformed points without associated covariance "
        poly2 = []
        for p in points2:
            poly2.append([p[0],p[1]])            
        
        " get the circles and radii "
        #radius2, center2 = gen_icp.computeEnclosingCircle(points2)
                
        support_pairs = []
        #for i in range(len(poly2)):
        for i in range(len(vecPoints2)):
            
            p_2 = vecPoints2[i]
    
            " for every transformed point of A, find it's closest neighbor in B "
            #_1, minDist = gen_icp.findClosestPointInB(supportPoints, p_2, [0.0,0.0,0.0])

            try:
                p_1, minI, minDist = gen_icp.findClosestPointWithAngle(vecPoints1, p_2, math.pi/8.0)
    
                #if gen_icp.isInCircle(p_1, radius2, center2):
    
                if minDist <= minMatchDist2:
                    C2 = points2[i][2]
                    C1 = supportPoints[minI][2]
                    #C1 = p_1[2]
    
                    " we store the untransformed point, but the transformed covariance of the A point "
                    #support_pairs.append([points2[i],p_1,C2,C1])
                    support_pairs.append([points2[i],supportPoints[minI],C2,C1])
            except:
                pass

        cost = 0.0
        if len(support_pairs) == 0:
            cost = 1e100
        else:
            vals = []
            sum1 = 0.0
            for pair in support_pairs:
        
                a = pair[0]
                b = pair[1]
                Ca = pair[2]
                Cb = pair[3]
        
                ax = a[0]
                ay = a[1]        
                bx = b[0]
                by = b[1]
        
                c11 = Ca[0][0]
                c12 = Ca[0][1]
                c21 = Ca[1][0]
                c22 = Ca[1][1]
                        
                b11 = Cb[0][0]
                b12 = Cb[0][1]
                b21 = Cb[1][0]
                b22 = Cb[1][1]    
            
                val = gen_icp.computeMatchErrorP([0.0,0.0,0.0], [ax,ay], [bx,by], [c11,c12,c21,c22], [b11,b12,b21,b22])
                
                vals.append(val)
                sum1 += val
                
            cost = sum1 / len(support_pairs)
            

        if plotIter:
            pylab.clf()
            xP = []
            yP = []
            for p in supportPoints:
                xP.append(p[0])
                yP.append(p[1])
            pylab.plot(xP,yP, color='b')
    
            xP = []
            yP = []
            for p in poly2:
                xP.append(p[0])
                yP.append(p[1])
            pylab.plot(xP,yP, color='r')
            
            for pair in support_pairs:
                p1 = pair[0]
                p2 = pair[1]
                xP = [p1[0],p2[0]]
                yP = [p1[1],p2[1]]
                pylab.plot(xP,yP)
            
            pylab.xlim(-5,10)
            pylab.ylim(-8,8)
            pylab.title("%d %d cost = %f, count = %d" % (pathID1, pathID2, cost, len(support_pairs)))
            pylab.savefig("pathOverlapCost_%04u.png" % self.overlapPlotCount)
            self.overlapPlotCount += 1
        
        if len(support_pairs) == 0:
            return 1e100

        return cost


    def getOverlapCondition(self, supportLine, nodeID, plotIter = False):

        #plotIter = False

        if len(supportLine) == 0:
            return 1e100

        node2 = self.nodeHash[nodeID]
  
        hull2, medial2 = computeHullAxis(nodeID, node2, tailCutOff = False)
        #hull2, medial2 = computeHullAxis(nodeID, node2, tailCutOff = True)
    
        estPose2 = node2.getGlobalGPACPose()        
        
        minMatchDist2 = 0.2
        #minMatchDist2 = 0.5
        #minMatchDist2 = 1.0
        #minMatchDist2 = 2.0
            
        " set the initial guess "
        poseOrigin = Pose(estPose2)
        
        localSupport = []
        for pnt in supportLine:
            localSupport.append(poseOrigin.convertGlobalToLocal(pnt))
    
        supportSpline = SplineFit(localSupport, smooth=0.1)        
        vecPoints1 = supportSpline.getUniformSamples()
        supportPoints = gen_icp.addGPACVectorCovariance(vecPoints1,high_var=0.05, low_var = 0.001)
            
        medialSpline2 = SplineFit(medial2, smooth=0.1)
        vecPoints2 = medialSpline2.getUniformSamples()
        points2 = gen_icp.addGPACVectorCovariance(vecPoints2,high_var=0.05, low_var = 0.001)

        " transformed points without associated covariance "
        poly2 = []
        for p in points2:
            poly2.append([p[0],p[1]])            
        
        " get the circles and radii "
        radius2, center2 = gen_icp.computeEnclosingCircle(points2)
                
        support_pairs = []
        #for i in range(len(poly2)):
        for i in range(len(vecPoints2)):
            
            p_2 = vecPoints2[i]
    
            " for every transformed point of A, find it's closest neighbor in B "
            #_1, minDist = gen_icp.findClosestPointInB(supportPoints, p_2, [0.0,0.0,0.0])

            try:
                p_1, minI, minDist = gen_icp.findClosestPointWithAngle(vecPoints1, p_2, math.pi/8.0)
    
                #if gen_icp.isInCircle(p_1, radius2, center2):
    
                if minDist <= minMatchDist2:
                    C2 = points2[i][2]
                    C1 = supportPoints[minI][2]
                    #C1 = p_1[2]
    
                    " we store the untransformed point, but the transformed covariance of the A point "
                    #support_pairs.append([points2[i],p_1,C2,C1])
                    support_pairs.append([points2[i],supportPoints[minI],C2,C1])
            except:
                pass

        cost = 0.0
        if len(support_pairs) == 0:
            cost = 1e100
        else:
            vals = []
            sum1 = 0.0
            for pair in support_pairs:
        
                a = pair[0]
                b = pair[1]
                Ca = pair[2]
                Cb = pair[3]
        
                ax = a[0]
                ay = a[1]        
                bx = b[0]
                by = b[1]
        
                c11 = Ca[0][0]
                c12 = Ca[0][1]
                c21 = Ca[1][0]
                c22 = Ca[1][1]
                        
                b11 = Cb[0][0]
                b12 = Cb[0][1]
                b21 = Cb[1][0]
                b22 = Cb[1][1]    
            
                val = gen_icp.computeMatchErrorP([0.0,0.0,0.0], [ax,ay], [bx,by], [c11,c12,c21,c22], [b11,b12,b21,b22])
                
                vals.append(val)
                sum1 += val
                
            cost = sum1 / len(support_pairs)
            

        if plotIter:
            pylab.clf()
            xP = []
            yP = []
            for p in supportPoints:
                xP.append(p[0])
                yP.append(p[1])
            pylab.plot(xP,yP, color='b')
    
            xP = []
            yP = []
            for p in poly2:
                xP.append(p[0])
                yP.append(p[1])
            pylab.plot(xP,yP, color='r')
            
            for pair in support_pairs:
                p1 = pair[0]
                p2 = pair[1]
                xP = [p1[0],p2[0]]
                yP = [p1[1],p2[1]]
                pylab.plot(xP,yP)
                    
            pylab.xlim(-5,10)
            pylab.ylim(-8,8)
            pylab.title("nodeID %d, cost = %f, count = %d" % (nodeID, cost, len(support_pairs)))
            pylab.savefig("overlapCost_%04u.png" % self.overlapPlotCount)
            self.overlapPlotCount += 1
        
        if len(support_pairs) == 0:
            return 1e100

        return cost


    def getOverlapCondition2(self, supportLine, nodeID, plotIter = False):

        #plotIter = True

        if len(supportLine) <= 5:
            return 0

        node2 = self.nodeHash[nodeID]
  
        hull2, medial2 = computeHullAxis(nodeID, node2, tailCutOff = False)
        #hull2, medial2 = computeHullAxis(nodeID, node2, tailCutOff = True)
    
        estPose2 = node2.getGlobalGPACPose()        
        
        #minMatchDist2 = 0.2
        minMatchDist2 = 0.5
        #minMatchDist2 = 1.0
        #minMatchDist2 = 2.0
            
        " set the initial guess "
        poseOrigin = Pose(estPose2)
        
        localSupport = []
        for pnt in supportLine:
            localSupport.append(poseOrigin.convertGlobalToLocal(pnt))
    
        supportSpline = SplineFit(localSupport, smooth=0.1)        
        vecPoints1 = supportSpline.getUniformSamples()
        supportPoints = gen_icp.addGPACVectorCovariance(vecPoints1,high_var=0.05, low_var = 0.001)
            
        medialSpline2 = SplineFit(medial2, smooth=0.1)
        vecPoints2 = medialSpline2.getUniformSamples()
        points2 = gen_icp.addGPACVectorCovariance(vecPoints2,high_var=0.05, low_var = 0.001)

        " transformed points without associated covariance "
        poly2 = []
        for p in points2:
            poly2.append([p[0],p[1]])            
        
        " get the circles and radii "
        radius2, center2 = gen_icp.computeEnclosingCircle(points2)
                
        support_pairs = []
        #for i in range(len(poly2)):
        for i in range(len(vecPoints2)):
            
            p_2 = vecPoints2[i]
    
            " for every transformed point of A, find it's closest neighbor in B "
            #_1, minDist = gen_icp.findClosestPointInB(supportPoints, p_2, [0.0,0.0,0.0])

            try:
                p_1, minI, minDist = gen_icp.findClosestPointWithAngle(vecPoints1, p_2, math.pi/8.0)
    
                #if gen_icp.isInCircle(p_1, radius2, center2):
    
                if minDist <= minMatchDist2:
                    C2 = points2[i][2]
                    C1 = supportPoints[minI][2]
                    #C1 = p_1[2]
    
                    " we store the untransformed point, but the transformed covariance of the A point "
                    #support_pairs.append([points2[i],p_1,C2,C1])
                    support_pairs.append([points2[i],supportPoints[minI],C2,C1])
            except:
                pass

        cost = 0.0
        if len(support_pairs) == 0:
            cost = 1e100
        else:
            vals = []
            sum1 = 0.0
            for pair in support_pairs:
        
                a = pair[0]
                b = pair[1]
                Ca = pair[2]
                Cb = pair[3]
        
                ax = a[0]
                ay = a[1]        
                bx = b[0]
                by = b[1]
        
                c11 = Ca[0][0]
                c12 = Ca[0][1]
                c21 = Ca[1][0]
                c22 = Ca[1][1]
                        
                b11 = Cb[0][0]
                b12 = Cb[0][1]
                b21 = Cb[1][0]
                b22 = Cb[1][1]    
            
                val = gen_icp.computeMatchErrorP([0.0,0.0,0.0], [ax,ay], [bx,by], [c11,c12,c21,c22], [b11,b12,b21,b22])
                
                vals.append(val)
                sum1 += val
                
            cost = sum1 / len(support_pairs)
            

        if plotIter:
            pylab.clf()
            xP = []
            yP = []
            for p in supportPoints:
                xP.append(p[0])
                yP.append(p[1])
            pylab.plot(xP,yP, color='b')
    
            xP = []
            yP = []
            for p in poly2:
                xP.append(p[0])
                yP.append(p[1])
            pylab.plot(xP,yP, color='r')
            
            for pair in support_pairs:
                p1 = pair[0]
                p2 = pair[1]
                xP = [p1[0],p2[0]]
                yP = [p1[1],p2[1]]
                pylab.plot(xP,yP)
                    
            pylab.xlim(-5,10)
            pylab.ylim(-8,8)
            pylab.title("nodeID %d, cost = %f, count = %d" % (nodeID, cost, len(support_pairs)))
            pylab.savefig("overlapCost_%04u.png" % self.overlapPlotCount)
            self.overlapPlotCount += 1
        
        return len(support_pairs)

        
    def checkUniqueBranch(self, parentPathID, nodeID1, depAngle, depPoint):


        #return True
        print "checkUniqueBranch(", parentPathID, ",", nodeID1, ",", depAngle, ",", depPoint, ")"

        " check cartesian distance to similar junction points from parent path "

        " BEING MORE PERMISSIVE IN CREATING NEW BRANCHES BECAUSE WE CAN MERGE LATER "

        " cartesian distance "
        #DISC_THRESH = 1.0
        DISC_THRESH = 0.5

        " 60 degree threshold "
        #ANG_THRESH = 1.047
        #ANG_THRESH = 0.7
        ANG_THRESH = 0.523 # pi/6
        
        pathIDs = self.getPathIDs()
        
        for pathID in pathIDs:
            
            path = self.getPath(pathID)
            
            if path["parentID"] == parentPathID:
                
                junctionNodeID = path["branchNodeID"]
                #localJunctionPoint = path["localJunctionPose"]
                " check dist "

                #poseOrigin = Pose(self.nodeHash[junctionNodeID].getEstPose())
                #junctionPoint = poseOrigin.convertLocalOffsetToGlobal(localJunctionPoint)
                
                junctionPoint = self.getGlobalJunctionPose(pathID)
                
                dist = sqrt((depPoint[0]-junctionPoint[0])**2 + (depPoint[1]-junctionPoint[1])**2 )

                " check difference of tangential angle "
                angDiff = diffAngle(depAngle, junctionPoint[2])

                print "result = ", junctionNodeID, dist, angDiff
        
                isNodeFeatureless = self.nodeHash[nodeID1].getIsFeatureless()
                print "node", nodeID1, "is featureless =", isNodeFeatureless
        
                if dist < DISC_THRESH:
                    if fabs(angDiff) < ANG_THRESH:
                        
                        print "DUPLICATE junction of", junctionPoint, "rejecting with differences of", dist, angDiff
                        return False, pathID
                
                " Explicit limitation that we wont be able to detect a T-Junction coming from left or right"
                if isNodeFeatureless:
                    print "REJECT new branch because the branching node is featureless"
                    return False, -1
        
        return True, -1



    def determineBranch(self, nodeID1, nodeID2, frontExist1, frontExist2, frontInterior1, frontInterior2, depAngle1, depAngle2, depPoint1, depPoint2, parentPathID1, parentPathID2, dirFlag):
        
        """
        1)  both nodes departing
            - are same departure  ( equal departing angle and are on top of each other)  1 new path
            - are different departures ( different departing angle and are not necessarily on top of each other) 2 new paths
            
        2) neither node is departing
            - are on same path
            
        3) one node is departing only
            - is only one departure ( other node has no departure, or has external departure but not same angle )
            - is a false departure  ( falseness happens if bad map data, ignore )
            - are both departing ( should have at least external departure on both, compare departure angle )
        
        Now, we have the departing angles, but we need to compare them with each other. 
        
        dirFlag specifies which node we want to check branches on.  We do not want to branch with the anchored end of a probe sweep   
        """

        print "determineBranch(", nodeID1, ",", nodeID2, ",", frontExist1, ",", frontExist2, ",", frontInterior1, ",", frontInterior2, ",", depAngle1, ",", depAngle2, ",", depPoint1, ",", depPoint2, ",", parentPathID1, ",", parentPathID2, ",", dirFlag, ")"
        
        foreTerm1 = frontInterior1 and frontExist1
        foreTerm2 = frontInterior2 and frontExist2

        frontAngDiff = diffAngle(depAngle1, depAngle2)

        pathBranchIDs = [-1, -1]
        isBranch = [False, False]
        isNew = [False, False]
        
        " 60 degree threshold "
        ANG_THRESH = 1.047
        #ANG_THRESH = 0.3

        
        if foreTerm1 or foreTerm2:

            if foreTerm1 and foreTerm2:
                " both departing case: check if same or different "
                if fabs(frontAngDiff) < ANG_THRESH:
                    " are on top of each other and are branching to the same path "

                    isUnique1, duplicatePathID1 = self.checkUniqueBranch(parentPathID1, nodeID1, depAngle1, depPoint1)
                    isUnique2, duplicatePathID2 = self.checkUniqueBranch(parentPathID2, nodeID2, depAngle2, depPoint2)

                    if isUnique1 and isUnique2:
                        
                        pathID = parentPathID1
                        branchNodeID = nodeID1
                        globalJunctionPoint = depPoint1
                        depAng = depAngle1
                        junctionAug = [globalJunctionPoint[0], globalJunctionPoint[1], depAng]

                        poseOrigin = Pose(self.nodeHash[branchNodeID].getEstPose())
                        
                        #newPathID = self.addNewPath(pathID, branchNodeID, poseOrigin.convertGlobalPoseToLocal(junctionAug))
                        newPathID = self.addPath(pathID, branchNodeID, poseOrigin.convertGlobalPoseToLocal(junctionAug))

                        pathBranchIDs[0] = newPathID
                        pathBranchIDs[1] = newPathID

                        isBranch[0] = True
                        isBranch[1] = True
                        isNew[0] = True
                        isNew[1] = True
                        
                        print "foreA"

                    if isUnique1 and not isUnique2:
                        if dirFlag == 0:
                            " foreTerm1 has unique departure "    
                            pathID = parentPathID1
                            branchNodeID = nodeID1
                            globalJunctionPoint = depPoint1
                            depAng = depAngle1
                            junctionAug = [globalJunctionPoint[0], globalJunctionPoint[1], depAng]
    
                            poseOrigin = Pose(self.nodeHash[branchNodeID].getEstPose())
                            #newPathID = self.addNewPath(pathID, branchNodeID, poseOrigin.convertGlobalPoseToLocal(junctionAug))
                            newPathID = self.addPath(pathID, branchNodeID, poseOrigin.convertGlobalPoseToLocal(junctionAug))
                            pathBranchIDs[0] = newPathID
                            isBranch[0] = True
                            isNew[0] = True

                    if not isUnique1 and isUnique2:
                        if dirFlag == 1:
                            " foreTerm2 has unique departure "    
                            pathID = parentPathID2
                            branchNodeID = nodeID2
                            globalJunctionPoint = depPoint2
                            depAng = depAngle2
                            junctionAug = [globalJunctionPoint[0], globalJunctionPoint[1], depAng]
    
                            poseOrigin = Pose(self.nodeHash[branchNodeID].getEstPose())
                            #newPathID = self.addNewPath(pathID, branchNodeID, poseOrigin.convertGlobalPoseToLocal(junctionAug))
                            newPathID = self.addPath(pathID, branchNodeID, poseOrigin.convertGlobalPoseToLocal(junctionAug))
                            pathBranchIDs[1] = newPathID
                            isBranch[1] = True
                            isNew[1] = True
                        
                    
                    if duplicatePathID1 != -1:
                        pathBranchIDs[0] = duplicatePathID1
                        isBranch[0] = True
                        

                    if duplicatePathID2 != -1:
                        pathBranchIDs[1] = duplicatePathID2
                        isBranch[1] = True
                        
                    
                else:
                    " these are branching to separate new paths "

                    isUnique1, duplicatePathID1 = self.checkUniqueBranch(parentPathID1, nodeID1, depAngle1, depPoint1)
                    isUnique2, duplicatePathID2 = self.checkUniqueBranch(parentPathID2, nodeID2, depAngle2, depPoint2)

                    if isUnique1 and isUnique2:

                        if dirFlag == 0:
                            pathID = parentPathID1
                            branchNodeID = nodeID1
                            globalJunctionPoint = depPoint1
                            depAng = depAngle1
                            junctionAug = [globalJunctionPoint[0], globalJunctionPoint[1], depAng]
    
                            poseOrigin = Pose(self.nodeHash[branchNodeID].getEstPose())
                            #newPathID1 = self.addNewPath(pathID, branchNodeID, poseOrigin.convertGlobalPoseToLocal(junctionAug))
                            newPathID1 = self.addPath(pathID, branchNodeID, poseOrigin.convertGlobalPoseToLocal(junctionAug))
                            pathBranchIDs[0] = newPathID1
                            isBranch[0] = True
                            isNew[0] = True


                        if dirFlag == 1:
                            pathID = parentPathID2
                            branchNodeID = nodeID2
                            globalJunctionPoint = depPoint2
                            depAng = depAngle2
                            junctionAug = [globalJunctionPoint[0], globalJunctionPoint[1], depAng]
    
                            poseOrigin = Pose(self.nodeHash[branchNodeID].getEstPose())
                            #newPathID2 = self.addNewPath(pathID, branchNodeID, poseOrigin.convertGlobalPoseToLocal(junctionAug))
                            newPathID2 = self.addPath(pathID, branchNodeID, poseOrigin.convertGlobalPoseToLocal(junctionAug))
                            pathBranchIDs[1] = newPathID2
                            isBranch[1] = True
                            isNew[1] = True

                        print "foreB"

                    if isUnique1 and not isUnique2:
                        if dirFlag == 0:
                            " foreTerm1 has unique departure "    
                            pathID = parentPathID1
                            branchNodeID = nodeID1
                            globalJunctionPoint = depPoint1
                            depAng = depAngle1
                            junctionAug = [globalJunctionPoint[0], globalJunctionPoint[1], depAng]
    
                            poseOrigin = Pose(self.nodeHash[branchNodeID].getEstPose())
                            #newPathID = self.addNewPath(pathID, branchNodeID, poseOrigin.convertGlobalPoseToLocal(junctionAug))
                            newPathID = self.addPath(pathID, branchNodeID, poseOrigin.convertGlobalPoseToLocal(junctionAug))
                            pathBranchIDs[0] = newPathID
                            isBranch[0] = True
                            isNew[0] = True

                    if not isUnique1 and isUnique2:
                        if dirFlag == 1:
                            " foreTerm2 has unique departure "    
                            pathID = parentPathID2
                            branchNodeID = nodeID2
                            globalJunctionPoint = depPoint2
                            depAng = depAngle2
                            junctionAug = [globalJunctionPoint[0], globalJunctionPoint[1], depAng]
    
                            poseOrigin = Pose(self.nodeHash[branchNodeID].getEstPose())
                            #newPathID = self.addNewPath(pathID, branchNodeID, poseOrigin.convertGlobalPoseToLocal(junctionAug))
                            newPathID = self.addPath(pathID, branchNodeID, poseOrigin.convertGlobalPoseToLocal(junctionAug))
                            pathBranchIDs[1] = newPathID
                            isBranch[1] = True
                            isNew[1] = True

                    if duplicatePathID1 != -1:
                        pathBranchIDs[0] = duplicatePathID1
                        isBranch[0] = True
                        
                    if duplicatePathID2 != -1:
                        pathBranchIDs[1] = duplicatePathID2
                        isBranch[1] = True
                        
            
            elif foreTerm1 and not foreTerm2:

                isUnique1, duplicatePathID1 = self.checkUniqueBranch(parentPathID1, nodeID1, depAngle1, depPoint1)
                    
                if isUnique1:
                       
                    if frontExist2 and fabs(frontAngDiff) < ANG_THRESH:
                        " both have same departure "
                        pathID = parentPathID1
                        branchNodeID = nodeID1
                        globalJunctionPoint = depPoint1
                        depAng = depAngle1
                        junctionAug = [globalJunctionPoint[0], globalJunctionPoint[1], depAng]


                        poseOrigin = Pose(self.nodeHash[branchNodeID].getEstPose())
                        #newPathID = self.addNewPath(pathID, branchNodeID, poseOrigin.convertGlobalPoseToLocal(junctionAug))
                        newPathID = self.addPath(pathID, branchNodeID, poseOrigin.convertGlobalPoseToLocal(junctionAug))
                        pathBranchIDs[0] = newPathID
                        pathBranchIDs[1] = newPathID
                        isBranch[0] = True
                        isBranch[1] = True
                        isNew[0] = True
                        isNew[1] = True

                        print "foreC"
                        
                    else:
                        if dirFlag == 0:
                            " foreTerm1 has unique departure "    
                            pathID = parentPathID1
                            branchNodeID = nodeID1
                            globalJunctionPoint = depPoint1
                            depAng = depAngle1
                            junctionAug = [globalJunctionPoint[0], globalJunctionPoint[1], depAng]
    
                            poseOrigin = Pose(self.nodeHash[branchNodeID].getEstPose())
                            #newPathID = self.addNewPath(pathID, branchNodeID, poseOrigin.convertGlobalPoseToLocal(junctionAug))
                            newPathID = self.addPath(pathID, branchNodeID, poseOrigin.convertGlobalPoseToLocal(junctionAug))
                            pathBranchIDs[0] = newPathID
                            isBranch[0] = True
                            isNew[0] = True

                        print "foreD"
                        
                if duplicatePathID1 != -1:
                    pathBranchIDs[0] = duplicatePathID1
                    isBranch[0] = True
 
            elif foreTerm2 and not foreTerm1:

                isUnique2, duplicatePathID2 = self.checkUniqueBranch(parentPathID2, nodeID2, depAngle2, depPoint2)

                if isUnique2:
                    if frontExist1 and fabs(frontAngDiff) < ANG_THRESH:

                        " both have same departure "
                        pathID = parentPathID2
                        branchNodeID = nodeID2
                        globalJunctionPoint = depPoint2
                        depAng = depAngle2
                        junctionAug = [globalJunctionPoint[0], globalJunctionPoint[1], depAng]

                        poseOrigin = Pose(self.nodeHash[branchNodeID].getEstPose())
                        #newPathID = self.addNewPath(pathID, branchNodeID, poseOrigin.convertGlobalPoseToLocal(junctionAug))
                        newPathID = self.addPath(pathID, branchNodeID, poseOrigin.convertGlobalPoseToLocal(junctionAug))
                        pathBranchIDs[0] = newPathID
                        pathBranchIDs[1] = newPathID
                        isBranch[0] = True
                        isBranch[1] = True
                        isNew[0] = True
                        isNew[1] = True

                        print "foreE"

                    else:
                        
                        if dirFlag == 1:
                            " foreTerm2 has unique departure "    
                            pathID = parentPathID2
                            branchNodeID = nodeID2
                            globalJunctionPoint = depPoint2
                            depAng = depAngle2
                            junctionAug = [globalJunctionPoint[0], globalJunctionPoint[1], depAng]
    
                            poseOrigin = Pose(self.nodeHash[branchNodeID].getEstPose())
                            #newPathID = self.addNewPath(pathID, branchNodeID, poseOrigin.convertGlobalPoseToLocal(junctionAug))
                            newPathID = self.addPath(pathID, branchNodeID, poseOrigin.convertGlobalPoseToLocal(junctionAug))
                            pathBranchIDs[1] = newPathID
                            isBranch[1] = True
                            isNew[1] = True

                        print "foreF"

                if duplicatePathID2 != -1:
                    pathBranchIDs[1] = duplicatePathID2
                    isBranch[1] = True

            else:
                print "foreG"
                " no new departure, just add nodes to the leaf paths "
                pathID = parentPathID1
                branchNodeID = nodeID1
                globalJunctionPoint = depPoint1
                
                if parentPathID1 != parentPathID2:
                    print "departing path IDs of two colocated nodes are not the same"
                    raise
        
        
        return isBranch, pathBranchIDs, isNew
        
        
    def splicePathIDs(self, pathIDs):
        
        if len(pathIDs) == 0:
            return []
        
        if len(pathIDs) == 1:
            return [self.trimmedPaths[pathIDs[0]]]

        " Assumption:  the pathIDs are connected paths with parent-child relations "


        " find the root "
        " pick any node, and go up the tree until we hit the root "
        currPathID = pathIDs[0]
        
        while True:
            #print "currPathID =", currPathID
            parent = self.getPath(currPathID)
            thisParentID = parent["parentID"]
            
            if thisParentID == None:
                break
            
            if pathIDs.count(thisParentID) == 0:
                break
            
            currPathID = thisParentID

        " currPathID is the root path "
        rootPathID = currPathID
        
            
        " if both terminal paths are child paths, 1 resultant spliced path "
        " if one terminal path is the root, then 2 resultant spliced paths "
        rightPathID = pathIDs[0]
        leftPathID = pathIDs[-1]

        newPaths = []        
        if rightPathID == rootPathID or leftPathID == rootPathID:
            
            print "one of paths is root make 2 resultant spliced paths"
            
            if rightPathID != rootPathID:
                childPathID = rightPathID
            else:
                childPathID = leftPathID
            
            " find path 1 from rootPathID to childPathID "
            " find path 2 from rootPathID to childPathID "
            rootPath = self.trimmedPaths[rootPathID]
            childPath = self.trimmedPaths[childPathID]


            for join in self.joins:
                if join[0][0] == childPathID:
                    childJunctionPoint = join[2]
                    childI = join[0][1]

            if childI < abs(childI-len(childPath)-1):
                childTermI = len(childPath)-1
            else:
                childTermI = 0

            "find path between: (childPathID, childTermI) to (rootPathID, 0) and (rootPathID, len(rootPath)-1)"
            shortestPathSpanTree, shortestDist = self.pathGraph.shortest_path((childPathID, childTermI))
            startNode1 = shortestPathSpanTree[(rootPathID, 0)]
            startNode2 = shortestPathSpanTree[(rootPathID, len(rootPath)-1)]
            
            #print "len(shortestPathSpanTree) =", len(shortestPathSpanTree)
            #print "shorestPathSpanTree:", shortestPathSpanTree

            currNode = startNode1
            splicedPath1 = []
            while currNode != (childPathID, childTermI):
                #print currNode, "!=", (childPathID, childTermI)
                splicedPath1.append(self.pathGraph.get_node_attributes(currNode))
                currNode = shortestPathSpanTree[currNode]
            splicedPath1.append(self.pathGraph.get_node_attributes(currNode))

            currNode = startNode2
            splicedPath2 = []
            while currNode != (childPathID, childTermI):
                #print currNode, "!=", (childPathID, childTermI)
                splicedPath2.append(self.pathGraph.get_node_attributes(currNode))
                currNode = shortestPathSpanTree[currNode]
            splicedPath2.append(self.pathGraph.get_node_attributes(currNode))

            #print "splicedPath1:", splicedPath1
            #print "splicedPath2:", splicedPath2
            
            newPaths.append(splicedPath1)
            newPaths.append(splicedPath2)
            
        else:
            "find entire path from rightPathID to leftPathID "
            " start from point farthest from child's junction point "

            #print "both paths are childs,  make 1 resultant spliced paths"

            rightPath = self.trimmedPaths[rightPathID]
            leftPath = self.trimmedPaths[leftPathID]
        
            " find the junction of this child's to its parent "    
            #joins.append([(pathID, minI1),(parentPathID, minI2), junctionPoint])
            
            for join in self.joins:
                if join[0][0] == rightPathID:
                    rightJunctionPoint = join[2]
                    rightI = join[0][1]

                if join[0][0] == leftPathID:
                    leftJunctionPoint = join[2]
                    leftI = join[0][1]
                    
            
            if leftI < abs(leftI-len(leftPath)-1):
                leftTermI = len(leftPath)-1
            else:
                leftTermI = 0

            if rightI < abs(rightI-len(rightPath)-1):
                rightTermI = len(rightPath)-1
            else:
                rightTermI = 0
                
            "find path between: (rightPathID, rightTermI) to (leftPathID, rightTermI) "
            shortestPathSpanTree, shortestDist = self.pathGraph.shortest_path((rightPathID, rightTermI))
            currNode = shortestPathSpanTree[(leftPathID, leftTermI)]

            splicedPath = []
            while currNode != (rightPathID, rightTermI):
                #print currNode, "!=", (rightPathID, rightTermI)
                splicedPath.append(self.pathGraph.get_node_attributes(currNode))
                currNode = shortestPathSpanTree[currNode]
            splicedPath.append(self.pathGraph.get_node_attributes(currNode))

            #print "len(shortestPathSpanTree) =", len(shortestPathSpanTree)
            #print "splicedPath:", splicedPath
            #print "shorestPathSpanTree:", shortestPathSpanTree

            newPaths.append(splicedPath)





        return newPaths

    def getNearestPathPoint(self, originPoint):

        " find closest point on each path "
        minI1 = 0
        minJ1 = 0
        minDist1 = 1e100
        
        pathIDs = self.getPathIDs()    
        for pathID in pathIDs:
            path = self.trimmedPaths[pathID]
            for j in range(len(path)):
                p0 = path[j]
                dist = sqrt((originPoint[0]-p0[0])**2 + (originPoint[1]-p0[1])**2)
                
                if dist < minDist1:
                    minDist1 = dist
                    minI1 = pathID
                    minJ1 = j
        
        newPoint = self.trimmedPaths[minI1][minJ1]
        
        return newPoint
        
    def computeNavigationPath(self, startPose, endPose):
        
        " get the trimmed paths "
        
        " find closest point on each path "
        minI1 = 0
        minJ1 = 0
        minDist1 = 1e100
        minI2 = 0
        minJ2 = 0
        minDist2 = 1e100
        
        pathIDs = self.getPathIDs()    
        for pathID in pathIDs:
            #for i in range(len(self.trimmedPaths)):
            path = self.trimmedPaths[pathID]
            for j in range(len(path)):
                p0 = path[j]
                dist = sqrt((startPose[0]-p0[0])**2 + (startPose[1]-p0[1])**2)
                
                if dist < minDist1:
                    minDist1 = dist
                    minI1 = pathID
                    minJ1 = j

                dist = sqrt((endPose[0]-p0[0])**2 + (endPose[1]-p0[1])**2)

                if dist < minDist2:
                    minDist2 = dist
                    minI2 = pathID
                    minJ2 = j
        
        " select shortest one "
        
        " get origin and target path "
        
        " get splice of origin and target path "
        " splice minJ1 and minJ2 "
        orderPathIDs = self.getPathPath(minI1, minI2)        
        
        splicedPaths = self.splicePathIDs(orderPathIDs)
        
        print len(splicedPaths), "spliced paths for path finding "
        
        " again, find closest point of start and end pose "
        
        minTotalDist = 3e100
        minSplicePathID = 0
        spliceIndex = (0,1)
        for i in range(len(splicedPaths)):

            minDist1 = 1e100
            minDist2 = 1e100

            minJ1 = 0
            minJ2 = 0

            
            path = splicedPaths[i]

            for j in range(len(path)):
                p0 = path[j]
                dist1 = sqrt((startPose[0]-p0[0])**2 + (startPose[1]-p0[1])**2)
                
                if dist1 < minDist1:
                    minDist1 = dist1
                    minJ1 = j

                dist2 = sqrt((endPose[0]-p0[0])**2 + (endPose[1]-p0[1])**2)

                if dist2 < minDist2:
                    minDist2 = dist2
                    minJ2 = j

            
            totalDist = minDist1 + minDist2
            if totalDist < minTotalDist:
                minTotalDist = totalDist
                minSplicePathID = i
                spliceIndex = (minJ1,minJ2)

            print i, path[0], path[-1], minJ1, minJ2, minDist1, minDist2, totalDist, minTotalDist, minSplicePathID, spliceIndex
        
        
        " return splicePath[startIndex:endIndex+1]"
        if spliceIndex[0] < spliceIndex[1]:
            #return splicedPaths[minSplicePathID][spliceIndex[0]:spliceIndex[1]+1]
            return splicedPaths[minSplicePathID]
        else:
            #path = splicedPaths[minSplicePathID][spliceIndex[1]:spliceIndex[0]+1]
            path = splicedPaths[minSplicePathID]
            path.reverse()
            return path

            
