#!/usr/bin/python

"""
2D Generalized ICP

WHAT:

This is a python implementation of the generalized ICP algorithm described in:
"Generalized-ICP", Aleksandr V. Segal, Dirk Haehnel, and Sebastian Thrun, In Robotics: Science and Systems 2009

Paper available on Aleksander Segal's personal website:
http://www.stanford.edu/~avsegal/

TO INSTALL:

This code requires 3 modules to be installed on your system.  Numpy, Scipy, and the Python NIPALS PCA module by Henning Risvik.

http://www.scipy.org/kk

http://numpy.scipy.org/

http://folk.uio.no/henninri/pca_module/


TO RUN:

To run the example program, simply execute this file from the command-line with the following command:
python gen_icp.py

Copy the main section to your own code to have a working implementation and include this file as a module.

This code is provided under the GNU General Public License, Version 3 

Copyright 2009, by Jacob Everist
jacob.everist@gmail.com

http://jacobeverist.com/gen_icp

"""

import sys
import numpy
import scipy 
import scipy.linalg
import scipy.optimize
import math
import pylab
import pca_module
import functions
from subprocess import *
import traceback
from math import cos, sin, pi, ceil
import cProfile

from matplotlib.patches import Circle

from copy import copy
from copy import deepcopy

from Pose import Pose
from StableCurve import StableCurve
from SplineFit import SplineFit

import Image
#from medialaxis import computeMedialAxis
#import graph

from icp import shapeCostC
from icp import computeMatchErrorP
from icp import matchPairs

import time

import multiprocessing as processing
import ctypes, os

numIterations = 0
globalPlotCount = 0

fig = pylab.figure()

import nelminICP


" pyCUDA imports "
try:
    #import pycuda.driver as drv
    #import pycuda.tools
    #import pycuda.autoinit
    #import pycuda.gpuarray as gpuarray
    #from pycuda.compiler import SourceModule
    #import cProfile
    
    threadsPerBlock = 256
    matchStr = """
    #define THREAD_COUNT %d
    
    __global__ void computeMatchErrorG(const float* A, const float* offset, float *sumVal, int N)
    {
    
        __shared__ float cacheVal[THREAD_COUNT];
    
        int tid = threadIdx.x + blockIdx.x * blockDim.x;
        int cacheIndex = threadIdx.x;
    
        float xd, yd, theta = 0.0;
        float ax = 0.0;
        float ay = 0.0;
        float bx = 0.0;
        float by1 = 0.0;
        float tx, ty, dx, dy = 0.0;
        float r11, r12, r21, r22 = 0.0;
        float c11, c12, c21, c22 = 0.0;
        float b11, b12, b21, b22 = 0.0;
        float res11, res12, res21, res22, resDet = 0.0;
        float q11, q12, q21, q22 = 0.0;
        float errVal = 0.0;
    
    
        float tempVal = 0.0;
    
        while (tid < N) {
    
            xd = offset[0];
            yd = offset[1];
            theta = offset[2];
    
            ax = A[tid*12];
            ay = A[tid*12 + 1];
            bx = A[tid*12 + 2];
            by1 = A[tid*12 + 3];
            c11 = A[tid*12 + 4];
            c12 = A[tid*12 + 5];
            c21 = A[tid*12 + 6];
            c22 = A[tid*12 + 7];
            b11 = A[tid*12 + 8];
            b12 = A[tid*12 + 9];
            b21 = A[tid*12 + 10];
            b22 = A[tid*12 + 11];
        
            tx = ax*cos(theta) - ay*sin(theta) + xd;
            ty = ax*sin(theta) + ay*cos(theta) + yd;
            dx = bx - tx;
            dy = by1 - ty;
        
            r11 = cos(theta);
            r12 = -sin(theta);
            r21 = sin(theta);
            r22 = r11;
        
            res11 = b11 + r11*(c11*r11 + c12*r12) + r12*(c21*r11 + c22*r12);
            res12 = b12 + r11*(c11*r21 + c12*r22) + r12*(c21*r21 + c22*r22);
            res21 = b21 + r21*(c11*r11 + c12*r12) + r22*(c21*r11 + c22*r12);
            res22 = b22 + r21*(c11*r21 + c12*r22) + r22*(c21*r21 + c22*r22);
            
            resDet = res22*res11 - res12*res21;
            
            q11 = res22/resDet;
            q12 = -res12/resDet;
            q21 = -res21/resDet;
            q22 = res11/resDet;
        
            errVal = dx*(dx*q11 + dy*q12) + dy*(dx*q21 + dy*q22);
        
            
            //return errVal;
            tempVal += errVal;
    
            tid += blockDim.x * gridDim.x;
        }
    
        cacheVal[cacheIndex] = tempVal;
    
        __syncthreads();
    
    
    
        int i = blockDim.x/2;
        while (i != 0) {
            if (cacheIndex < i) {
                cacheVal[cacheIndex] += cacheVal[cacheIndex + i];
            }
    
            __syncthreads();
            i /= 2;
        }
    
        // thread 0 returns the value 
        if (cacheIndex == 0) {
            sumVal[blockIdx.x] = cacheVal[0];
        }
        //tid = threadIdx.x + blockIdx.x * blockDim.x;
        //sumVal[tid] = cacheVal[threadIdx.x];
    }
    
    """ # % threadsPerBlock

    #mod3 = SourceModule(matchStr)
    #match_them = mod3.get_function("computeMatchErrorG")
except:
    isGPU = False
else:
    isGPU = True

def __num_processors():
    if os.name == 'nt': # Windows
        return int(os.getenv('NUMBER_OF_PROCESSORS'))
    else: # glibc (Linux, *BSD, Apple)
        get_nprocs = ctypes.cdll.libc.get_nprocs
        get_nprocs.restype = ctypes.c_int
        get_nprocs.argtypes = []
        return get_nprocs()

def __do_nothing(data, someArg1, someArg2, someArg3):
    
    knn = [someArg1]
    
    for p in data:
        for p2 in data:
            pass
    
    return copy(data)

def __remote_process(rank, qin, qout, someArg1, someArg2, someArg3):

    while 1:
        # read input queue (block until data arrives)
        nc, data = qin.get()
        # process data
        knn = __do_nothing(data, nc, someArg2, someArg3)
        # write to output queue
        qout.put((nc,knn))

def __remote_ICP(rank, qin, qout, globalPath, medial):

    while 1:
        # read input queue (block until data arrives)
        nc, args = qin.get()
        #print rank, "received", nc, args
        # process data
        #knn = __do_nothing(data, nc, someArg2, someArg3)
        results = []
        for arg in args:
            resultPose, lastCost = globalOverlapICP_GPU2(arg, globalPath, medial)
            results.append([resultPose, lastCost])
        # write to output queue
        qout.put((nc,results))


def batchGlobalICP(globalPath, medial, args):

 

    ndata = len(args)
    nproc = __num_processors()
    
    print "nproc =", nproc
    
    # compute chunk size
    #chunk_size = ceil(float(ndata) / float(nproc))
    #chunk_size = int(chunk_size)
    chunk_size = ndata / nproc
    chunk_size = 2 if chunk_size < 2 else chunk_size
    
    print "chunk_size =", chunk_size
    print "max_size =", ndata/chunk_size
    
    # set up a pool of processes
    qin = processing.Queue(maxsize=ndata/chunk_size)
    qout = processing.Queue(maxsize=ndata/chunk_size)
    pool = [processing.Process(target=__remote_ICP,
                args=(rank, qin, qout, globalPath, medial))
                    for rank in range(nproc)]
    for p in pool: p.start()
    
    # put data chunks in input queue
    cur, nc = 0, 0
    while 1:
        #_data = data[:,cur:cur+chunk_size]
        _data = args[cur:cur+chunk_size]
        if len(_data) == 0: break
        print "put(", (nc,_data), ")"
        qin.put((nc,_data))
        print "DONE"
        cur += chunk_size
        nc += 1
    
    # read output queue
    knn = []
    while len(knn) < nc:
        knn += [qout.get()]
    # avoid race condition
    _knn = [n for i,n in sorted(knn)]
    knn = []
    for tmp in _knn:
        knn += tmp
    # terminate workers
    for p in pool: p.terminate()
    return knn


def knn_search():

    """ find the K nearest neighbours for data points in data,
        using an O(n log n) kd-tree, exploiting all logical
        processors on the computer """

    #param = data.shape[0]
    # build kdtree
    #tree = kdtree(data.copy(), leafsize=leafsize)

    
    dataFoo = range(0,101)
    
    data = dataFoo
    ndata = len(data)
    nproc = __num_processors()
    
    print "nproc =", nproc
    
    # compute chunk size
    chunk_size = len(data) / nproc
    chunk_size = 100 if chunk_size < 100 else chunk_size
    
    someArg1 = 5
    someArg2 = range(20)
    someArg3 = "foobar maybe"
    
    
    # set up a pool of processes
    qin = processing.Queue(maxsize=ndata/chunk_size)
    qout = processing.Queue(maxsize=ndata/chunk_size)
    pool = [processing.Process(target=__remote_process,
                args=(rank, qin, qout, someArg1, someArg2, someArg3))
                    for rank in range(nproc)]
    for p in pool: p.start()
    
    # put data chunks in input queue
    cur, nc = 0, 0
    while 1:
        #_data = data[:,cur:cur+chunk_size]
        _data = data[cur:cur+chunk_size]
        if len(_data) == 0: break
        qin.put((nc,_data))
        cur += chunk_size
        nc += 1
    
    # read output queue
    knn = []
    while len(knn) < nc:
        knn += [qout.get()]
    # avoid race condition
    _knn = [n for i,n in sorted(knn)]
    knn = []
    for tmp in _knn:
        knn += tmp
    # terminate workers
    for p in pool: p.terminate()
    return knn


def matchPairs2(points, points_offset, globalPoints, minMatchDist):
    
    def findClosestPointInA(a_trans, b):
    
        #a_p, a_i, minDist = findClosestPointInA(a_trans, b_p, [0.0,0.0,0.0])
        
        minDist = 1e100
        minPoint = None
        min_i = 0
    
        for i in range(len(a_trans)):
            p = a_trans[i]
    
            dist = math.sqrt((p[0]-b[0])**2 + (p[1]-b[1])**2)
            if dist < minDist:
                minPoint = copy(p)
                minDist = dist
                min_i = i
    
        if minPoint != None:
            return minPoint, min_i, minDist
        else:
            raise

    def computeEnclosingCircle(a_data):
            
        maxA = 0.0
        a_p1 = []
        a_p2 = []
        for i in range(len(a_data)):
            p1 = a_data[i]
    
            for j in range(i+1, len(a_data)):
                p2 = a_data[j]
    
                dist = math.sqrt((p1[0]-p2[0])**2 + (p1[1]-p2[1])**2)
    
                if dist > maxA:
                    maxA = dist
                    a_p1 = copy(p1)
                    a_p2 = copy(p2)
    
        radiusA = maxA/2.0
    
        centerA = [(a_p1[0] + a_p2[0])/2.0, (a_p1[1] + a_p2[1])/2.0]
        
        return radiusA, centerA

    def isInCircle(p, radius, center):
        
        dist = math.sqrt((p[0] - center[0]) ** 2 + (p[1] - center[1]) ** 2)
        if dist < radius:
            return True
    
        return False
    
    match_pairs = []
    
    " get the circles and radii "
    radius2, center2 = computeEnclosingCircle(points_offset)
    
    for i in range(len(globalPoints)):
        p_1 = [globalPoints[i][0],globalPoints[i][1]]
        #p_1 = localPoly[i]
    
        if isInCircle(p_1, radius2, center2):
    
            " for every point of B, find it's closest neighbor in transformed A "
            p_2, i_2, minDist = findClosestPointInA(points_offset, p_1)
    
            if minDist <= minMatchDist:
    
                " add to the list of match pairs less than 1.0 distance apart "
                " keep A points and covariances untransformed "
                C2 = points[i_2][2]
                
                C1 = globalPoints[i][2]
                
                " we store the untransformed point, but the transformed covariance of the A point "
                match_pairs.append([points[i_2],globalPoints[i],C2,C1])
                


    return match_pairs

def getLongestPath(node, currSum, currPath, tree, isVisited, nodePath, nodeSum):

    if isVisited[node]:
        return

    isVisited[node] = 1
    nodePath[node] = copy(currPath) + [node]

    if nodeSum[node] < currSum:
        nodeSum[node] = currSum

    for childNode in tree[node]:
        getLongestPath(childNode, currSum + 1, nodePath[node], tree, isVisited, nodePath, nodeSum)

    isVisited[node] = 0

def getMedialAxis2(hull):

    minX = 1e100
    maxX = -1e100
    minY = 1e100
    maxY = -1e100
    for p in hull:
        if p[0] > maxX:
            maxX = p[0]
        if p[0] < minX:
            minX = p[0]
        if p[1] > maxY:
            maxY = p[1]
        if p[1] < minY:
            minY = p[1]


    PIXELSIZE = 0.05
    #mapSize = 0.15*40 + 2.0 + 2.0
    #mapSize = 5.0
    #mapSize = max(maxX-minX,maxY-minY) + 10*PIXELSIZE
    mapSize = 2*max(max(maxX,math.fabs(minX)),max(maxY,math.fabs(minY))) + 1
    pixelSize = PIXELSIZE
    #print mapSize, mapSize/pixelSize
    numPixel = int(2.0*mapSize / pixelSize + 1.0)
    divPix = math.floor((2.0*mapSize/pixelSize)/mapSize)

    #print "numPixel = ", numPixel

    def realToGrid(point):
        indexX = int(math.floor(point[0]*divPix)) + numPixel/2 + 1
        indexY = int(math.floor(point[1]*divPix)) + numPixel/2 + 1
        return indexX, indexY

    def gridToReal(indices):
        i = indices[0]
        j = indices[1]
        point = [(i - numPixel/2 - 1)*(pixelSize/2.0) + pixelSize/2.0, (j - numPixel/2 - 1)*(pixelSize/2.0) + pixelSize/2.0]
        return point


    #print "hull:", hull

    gridHull = []
    #for p in hull:
    for i in range(len(hull)):
        p = hull[i]
        #if abs(p[0]) >= mapSize or abs(p[1]) >= mapSize:
        #    print "point:", p
        gridHull.append(realToGrid(p))

    #print "gridHull:", gridHull
    
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



    #print "minMax:", minX, maxX, minY, maxY
            
    xRange = range(minX, maxX + 1)
    yRange = range(minY, maxY + 1)

    #print "xRange:", xRange
    #print "yRange:", yRange
    interior = []

    for i in xRange:
        for j in yRange:
            if functions.point_inside_polygon(i, j, gridHull):
                interior.append((i,j))


    inputImg = Image.new('L', (numPixel,numPixel), 255)
    imga = inputImg.load()
    
    #for p in gridHull:
    for i in range(len(gridHull)):
        p = gridHull[i]
        #print p, hull[i]
        imga[p[0],p[1]] = 0
    
    for p in interior:
        imga[p[0],p[1]] = 0
        #i, j = realToGrid(p)
        #imga[i,j] = 0

    resultImg = Image.new('L', inputImg.size)
    resultImg = computeMedialAxis(inputImg, resultImg)
    imgA = resultImg.load()

    points = []
    for i in range(1, inputImg.size[0]-1):
        for j in range(1, inputImg.size[1]-1):
            if imgA[i,j] == 0:
                points.append((i,j))

    medialGraph = graph.graph()
    for p in points:
        medialGraph.add_node(p, [])




    builtGraph = {}
    for i in range(2, inputImg.size[0]-2):
        for j in range(2, inputImg.size[1]-2):
            if imgA[i,j] == 0:
                builtGraph[(i,j)] = []
                for k in range(i-1, i+2):
                    for l in range(j-1, j+2):
                        if imgA[k,l] == 0:
                            builtGraph[(i,j)].append((k,l))
                            medialGraph.add_edge((i,j), (k,l))
                            


    mst = medialGraph.minimal_spanning_tree()

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

    leaves = []
    for k, v in uni_mst.items():
        if len(v) == 1:
            leaves.append(k)


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
            #print k, v
            if v > maxDist:
                maxNode = k
                maxDist = v

        #print leaf, "-", maxNode, maxDist

        if maxDist > maxPairDist:
            maxPairDist = maxDist
            maxPair = [leaf, maxNode]
            maxPath = nodePath[maxNode]


    frontVec = [0.,0.]
    backVec = [0.,0.]
    indic = range(6)
    indic.reverse()

    for i in indic:
        p1 = maxPath[i+2]
        p2 = maxPath[i]
        vec = [p2[0]-p1[0], p2[1]-p1[1]]
        frontVec[0] += vec[0]
        frontVec[1] += vec[1]

        p1 = maxPath[-i-3]
        p2 = maxPath[-i-1]
        vec = [p2[0]-p1[0], p2[1]-p1[1]]
        backVec[0] += vec[0]
        backVec[1] += vec[1]

    frontMag = math.sqrt(frontVec[0]*frontVec[0] + frontVec[1]*frontVec[1])
    backMag = math.sqrt(backVec[0]*backVec[0] + backVec[1]*backVec[1])

    frontVec[0] /= frontMag
    frontVec[1] /= frontMag
    backVec[0] /= backMag
    backVec[1] /= backMag

    newP1 = (maxPath[0][0] + frontVec[0]*500, maxPath[0][1] + frontVec[1]*500)
    newP2 = (maxPath[-1][0] + backVec[0]*500, maxPath[-1][1] + backVec[1]*500)

    maxPath.insert(0,newP1)
    maxPath.append(newP2)


    " convert path points to real "
    realPath = []
    for p in maxPath:
        realPath.append(gridToReal(p))


    return realPath



def computeSeparation(bbox2, uniform2_trans, hull):

    interPoints = []
    
    " start from 0 and find first intersection with bbox "
    firstFound = False
    subSec1 = []
    bbEdge1 = []
    for j in range(len(uniform2_trans)-1):
        p1 = uniform2_trans[j]
        p2 = uniform2_trans[j+1]
        
        edge2 = [p1,p2]

        for i in range(len(bbox2)):
            p1 = bbox2[i]
            p2 = bbox2[(i+1) % 4]
            
            edge1 = [p1,p2]
        
            result, point = functions.Intersect(edge1, edge2)
            if result:
                firstFound = True
                interPoints.append(point)
                subSec1 = uniform2_trans[j+1:]
                subSec1.insert(0,point)
                bbEdge1 = [i,i+1]
                break
            
        if firstFound:
            break
        
    
    " start from -1 and find first intersection with bbox"
    secondFound = False
    subSec2 = []
    bbEdge2 = []
    for j in range(len(subSec1)-1):
        p1 = subSec1[-1-j]
        p2 = subSec1[-2-j]
        
        edge2 = [p1,p2]

        for i in range(len(bbox2)):
            p1 = bbox2[i]
            p2 = bbox2[(i+1) % 4]
            
            edge1 = [p1,p2]
        
            result, point = functions.Intersect(edge1, edge2)
            if result:
                secondFound = True
                interPoints.append(point)
                
                subSec2 = subSec1[0:-1-j]
                subSec2.append(point)
                
                bbEdge2 = [i,i+1]
                break
            
        if secondFound:
            break

    " split the bounding box in 2 "
    halfBox1 = []
    for i in range(len(bbox2)):
        halfBox1.append(bbox2[(bbEdge1[1]+i) % 4])
        
        if (bbEdge1[1]+i) % 4 == bbEdge2[0]:
            break
    
    subSec3 = deepcopy(subSec2)
    subSec3.reverse()
    halfBox1 += subSec3


    " split the bounding box in 2 "
    halfBox2 = []
    for i in range(len(bbox2)):
        halfBox2.append(bbox2[(bbEdge2[1]+i) % 4])
        
        if (bbEdge2[1]+i) % 4 == bbEdge1[0]:
            break
    
    halfBox2 += subSec2


    " separated hull points now have a break in the sequence "
    halfHullA = []
    halfHullB = []
    lastA = True
    lastB = False
    for p in hull:
        
        inA = False
        inB = False
        
        if functions.point_inside_polygon(p[0], p[1], halfBox1):
            #halfHullA.append(p)
            inA = True

        if functions.point_inside_polygon(p[0], p[1], halfBox2):
            #halfHullB.append(p)
            inB = True
        
        #print "hulls:", inA, inB
        if inA and inB or not inA and not inB:
            #print "BOTH hulls"
            
            #inA = lastA
            #inB = lastB
            inA = False
            inB = False

        if inA:
            halfHullA.append(p)
            
        if inB:
            halfHullB.append(p)
            
        lastA = inA
        lastB = inB

        #if not inA and not inB:
            #print "NEITHER hulls"
            #print "box:", bbox2

            #raise ValueError

    #print "points found:", len(interPoints)    
    if len(interPoints) != 2:
        print "interPoints:", interPoints
        print "box:", bbox2
        raise

    return halfHullA, halfHullB
    #return halfBox1, halfBox2

    #return a_data_A, a_data_B

def computeSeparation2(bbox, medialAxis, hull):

    interPoints = []

    foreIntersection = []
    backIntersection = []
    bbEdge1 = []
    bbEdge2 = []
    
    for i in range(len(bbox)):
        p1 = bbox[i]
        p2 = bbox[(i+1) % 4]
        
        edge1 = [p1,p2]
        edge2 = [medialAxis[0],medialAxis[1]]
    
        result, point = functions.Intersect(edge1, edge2)
        if result:
            foreIntersection = point
            bbEdge1 = [i,(i+1) % 4]

        edge1 = [p1,p2]
        edge2 = [medialAxis[-1],medialAxis[-2]]
    
        result, point = functions.Intersect(edge1, edge2)
        if result:
            backIntersection = point
            bbEdge2 = [i,(i+1) % 4]
            
    if len(foreIntersection) == 0 or len(backIntersection) == 0:
        raise

    medialSubSec = medialAxis[1:-2]

    halfBox1 = [foreIntersection] + medialSubSec + [backIntersection]
    halfBox2 = [foreIntersection] + medialSubSec + [backIntersection]
    
    indStart = bbEdge2[1]
    indTerm = bbEdge1[0]
    for i in range(4):
        halfBox1.append(bbox[(indStart + i) % 4])
    
        if (indStart + i) % 4 == indTerm:
            break
    
    indStart = bbEdge2[0]
    indTerm = bbEdge1[1]
    for i in range(4):
        halfBox2.append(bbox[(indStart - i) % 4])
    
        if (indStart - i) % 4 == indTerm:
            break

    
    #halfBox2.reverse()

    " separated hull points now have a break in the sequence "
    halfHullA = []
    halfHullB = []
    for p in hull:
        
        inA = False
        inB = False
        
        if functions.point_inside_polygon(p[0], p[1], halfBox1):
            inA = True

        if functions.point_inside_polygon(p[0], p[1], halfBox2):
            inB = True
        
        if inA and inB or not inA and not inB:
            print "BOTH or NEITHER hulls"
            
            inA = False
            inB = False

        if inA:
            halfHullA.append(p)
            
        if inB:
            halfHullB.append(p)

    return halfHullA, halfHullB


def computeBBox(points):
    
    xMax = -1e100
    xMin = 1e100
    yMax = -1e100
    yMin = 1e100
    
    
    for p in points:
        
        if p[0] > xMax:
            xMax = p[0]
        
        if p[0] < xMin:
            xMin = p[0]
            
        if p[1] > yMax:
            yMax = p[1]
        
        if p[1] < yMin:
            yMin = p[1]


    " pad the edges so not points on are on the boundary "
    xMax = xMax + 0.4
    xMin = xMin - 0.4
    yMax = yMax + 0.4
    yMin = yMin - 0.4
    
    bbox = [[xMin,yMin], [xMin,yMax], [xMax,yMax], [xMax,yMin]]
            
    return bbox



def dispP(p, offset):
    
    xd = offset[0]
    yd = offset[1]
    theta = offset[2]

    px = p[0]
    py = p[1]
    
    tx = px*math.cos(theta) - py*math.sin(theta) + xd
    ty = px*math.sin(theta) + py*math.cos(theta) + yd
    
    p_off = [tx, ty]

    return p_off

" displace the point by the offset plus modify it's covariance "
def dispPointWithAngle(p, offset):
    
    xd = offset[0]
    yd = offset[1]
    theta = offset[2]
    
    px = p[0]
    py = p[1]
    
    tx = px*math.cos(theta) - py*math.sin(theta) + xd
    ty = px*math.sin(theta) + py*math.cos(theta) + yd
    
    p_off = [tx, ty]

    Cv = p[2]
    
    r11 = math.cos(theta)
    r12 = -math.sin(theta)
    r21 = math.sin(theta)
    r22 = r11

    c11 = Cv[0][0]
    c12 = Cv[0][1]
    c21 = Cv[1][0]
    c22 = Cv[1][1]
    
    res11 = r11*(c11*r11 + c12*r12) + r12*(c21*r11 + c22*r12)
    res12 = r11*(c11*r21 + c12*r22) + r12*(c21*r21 + c22*r22)
    res21 = r21*(c11*r11 + c12*r12) + r22*(c21*r11 + c22*r12)
    res22 = r21*(c11*r21 + c12*r22) + r22*(c21*r21 + c22*r22)

    Ca = [[res11, res12], [res21, res22]]
    p_off.append(Ca)

    pa = p[3]
    p_off.append(functions.normalizeAngle(pa + theta))
    
    return p_off

" displace the point by the offset plus modify it's covariance "
def dispPoint(p, offset):
    
    xd = offset[0]
    yd = offset[1]
    theta = offset[2]

    #T = numpy.matrix([    [math.cos(theta), -math.sin(theta), xd],
    #        [math.sin(theta), math.cos(theta), yd],
    #        [0.0, 0.0, 1.0]
    #        ])
    
    px = p[0]
    py = p[1]
    
    tx = px*math.cos(theta) - py*math.sin(theta) + xd
    ty = px*math.sin(theta) + py*math.cos(theta) + yd
    
    #p_hom = numpy.matrix([[p[0]],[p[1]],[1.0]])
    #temp = T*p_hom
    #p_off = [temp[0,0],temp[1,0]]
    p_off = [tx, ty]

    Cv = p[2]
    
    r11 = math.cos(theta)
    r12 = -math.sin(theta)
    r21 = math.sin(theta)
    r22 = r11

    c11 = Cv[0][0]
    c12 = Cv[0][1]
    c21 = Cv[1][0]
    c22 = Cv[1][1]
    
    res11 = r11*(c11*r11 + c12*r12) + r12*(c21*r11 + c22*r12)
    res12 = r11*(c11*r21 + c12*r22) + r12*(c21*r21 + c22*r22)
    res21 = r21*(c11*r11 + c12*r12) + r22*(c21*r11 + c22*r12)
    res22 = r21*(c11*r21 + c12*r22) + r22*(c21*r21 + c22*r22)

    Ca = [[res11, res12], [res21, res22]]
    p_off.append(Ca)
    
    #print 
    #print p_off

    return p_off
    
    
    xd = offset[0]
    yd = offset[1]
    theta = offset[2]

    T = numpy.matrix([    [math.cos(theta), -math.sin(theta), xd],
            [math.sin(theta), math.cos(theta), yd],
            [0.0, 0.0, 1.0]
            ])

    p_hom = numpy.matrix([[p[0]],[p[1]],[1.0]])
    temp = T*p_hom
    p_off = [temp[0,0],temp[1,0]]

    R = numpy.matrix([    [math.cos(theta), -math.sin(theta)],
        [math.sin(theta), math.cos(theta)] ])

    Cv = p[2]
    Cv = numpy.matrix([[Cv[0][0],Cv[0][1]], [Cv[1][0], Cv[1][1]]])

    Ca = R * Cv * numpy.transpose(R)
    p_off.append(Ca)

    #print p_off

    return p_off

" displace the point by the offset only.  No covariance "
def dispOffset(p, offset):
    xd = offset[0]
    yd = offset[1]
    theta = offset[2]

    px = p[0]
    py = p[1]
    
    #p_off = [0.0,0.0]
    #p_off[0] = px*math.cos(theta) - py*math.sin(theta) + xd
    #p_off[1] = px*math.sin(theta) + py*math.cos(theta) + yd

    #p_off = [0.0,0.0]
    p_off = [px*math.cos(theta) - py*math.sin(theta) + xd, px*math.sin(theta) + py*math.cos(theta) + yd]
    
    return p_off

    #print
    #print p_off

    #T = numpy.matrix([    [math.cos(theta), -math.sin(theta), xd],
    #        [math.sin(theta), math.cos(theta), yd],
    #        [0.0, 0.0, 1.0]
    #        ])

    #p_hom = numpy.matrix([[p[0]],[p[1]],[1.0]])
    #temp = T*p_hom
    #p_off = [temp[0,0],temp[1,0]]

    #print p_off

    #return p_off

def dispOffsetMany(points, offset):
    xd = offset[0]
    yd = offset[1]
    theta = offset[2]

    T = numpy.matrix([    [math.cos(theta), -math.sin(theta), xd],
            [math.sin(theta), math.cos(theta), yd],
            [0.0, 0.0, 1.0]
            ])

    results = []

    for p in points:
        p_hom = numpy.matrix([[p[0]],[p[1]],[1.0]])
        temp = T*p_hom
        p_off = [temp[0,0],temp[1,0]]
        results.append(p_off)

    return results

def disp(ai, bi, T):

    temp = T*ai
    result = bi-temp    

    result[2] = 1.0

    return result

"""
def computeMatchErrorC(offset, a, b, Ca, Cb):
    #sum1 += computeMatchErrorC([xd, yd, theta], [ax,ay], [bx,by], [c11,c12,c21,c22], [b11,b12,b21,b22])

    xd = offset[0]
    yd = offset[1]
    theta = offset[2]
    
    ax = a[0]
    ay = a[1]
    
    bx = b[0]
    by = b[1]

    tx = ax*math.cos(theta) - ay*math.sin(theta) + xd
    ty = ax*math.sin(theta) + ay*math.cos(theta) + yd
    dx = bx - tx
    dy = by - ty

    r11 = math.cos(theta)
    r12 = -math.sin(theta)
    r21 = math.sin(theta)
    r22 = r11

    c11 = Ca[0]
    c12 = Ca[1]
    c21 = Ca[2]
    c22 = Ca[3]
    
    b11 = Cb[0]
    b12 = Cb[1]
    b21 = Cb[2]
    b22 = Cb[3]
    
    res11 = b11 + r11*(c11*r11 + c12*r12) + r12*(c21*r11 + c22*r12)
    res12 = b12 + r11*(c11*r21 + c12*r22) + r12*(c21*r21 + c22*r22)
    res21 = b21 + r21*(c11*r11 + c12*r12) + r22*(c21*r11 + c22*r12)
    res22 = b22 + r21*(c11*r21 + c12*r22) + r22*(c21*r21 + c22*r22)
    
    resDet = res22*res11 - res12*res21
    
    q11 = res22/resDet
    q12 = -res12/resDet
    q21 = -res21/resDet
    q22 = res11/resDet

    errVal = dx*(dx*q11 + dy*q12) + dy*(dx*q21 + dy*q22)

    return errVal
"""

def computeMatchError(offset, a, b, Ca, Cb):


    #print "received:", offset
    #xd = offset[0]
    #yd = offset[1]
    #theta = offset[2]
    
    #T = numpy.matrix([    [math.cos(theta), -math.sin(theta), xd],
    #        [math.sin(theta), math.cos(theta), yd],
    #        [0.0, 0.0, 1.0]
    #        ])

    #R = numpy.matrix([    [math.cos(theta), -math.sin(theta)],
    #    [math.sin(theta), math.cos(theta)] ])

    #Cv = Ca

    #d_vec = disp(numpy.concatenate((a,numpy.matrix([1.0]))), numpy.concatenate((b,numpy.matrix([1.0]))), T)
    #ai = numpy.concatenate((a,numpy.matrix([1.0])))
    #bi = numpy.concatenate((b,numpy.matrix([1.0])))
    #temp = T*ai
    #d_vec = bi-temp
    #d_vec[2] = 1.0    

    #res = Cb + R * Cv * numpy.transpose(R)

    #invMat = scipy.linalg.inv(res)

    # add homogeneous dimension back
    #invMat = numpy.concatenate((invMat,numpy.matrix([[0.0],[0.0]])), 1)
    #invMat = numpy.concatenate((invMat,numpy.matrix([0.0,0.0,0.0])))

    #error = numpy.transpose(d_vec)*invMat*d_vec
    #errVal = error[0,0]

    #return errVal

    xd = offset[0]
    yd = offset[1]
    theta = offset[2]
    
    #ax = a[0,0]
    #ay = a[1,0]
    ax = a[0]
    ay = a[1]
    
    #bx = b[0,0]
    #by = b[1,0]
    bx = b[0]
    by = b[1]

    #print ax, ay, xd, yd, theta
    tx = ax*math.cos(theta) - ay*math.sin(theta) + xd
    ty = ax*math.sin(theta) + ay*math.cos(theta) + yd
    dx = bx - tx
    dy = by - ty

    r11 = math.cos(theta)
    r12 = -math.sin(theta)
    r21 = math.sin(theta)
    r22 = r11

    c11 = Ca[0][0]
    c12 = Ca[0][1]
    c21 = Ca[1][0]
    c22 = Ca[1][1]
    
    b11 = Cb[0][0]
    b12 = Cb[0][1]
    b21 = Cb[1][0]
    b22 = Cb[1][1]
    
    res11 = b11 + r11*(c11*r11 + c12*r12) + r12*(c21*r11 + c22*r12)
    res12 = b12 + r11*(c11*r21 + c12*r22) + r12*(c21*r21 + c22*r22)
    res21 = b21 + r21*(c11*r11 + c12*r12) + r22*(c21*r11 + c22*r12)
    res22 = b22 + r21*(c11*r21 + c12*r22) + r22*(c21*r21 + c22*r22)
    
    resDet = res22*res11 - res12*res21
    
    q11 = res22/resDet
    q12 = -res12/resDet
    q21 = -res21/resDet
    q22 = res11/resDet

    errVal = dx*(dx*q11 + dy*q12) + dy*(dx*q21 + dy*q22)

    return errVal


def computeMatchErrorSimple(offset, a, b, Ca, Cb):

    xd = offset[0]
    yd = offset[1]
    theta = offset[2]

    T = numpy.matrix([    [math.cos(theta), -math.sin(theta), xd],
            [math.sin(theta), math.cos(theta), yd],
            [0.0, 0.0, 1.0]
            ])

    #R = [[math.cos(theta), -math.sin(theta)], [math.sin(theta), math.cos(theta)]]
    
    R = numpy.matrix([    [math.cos(theta), -math.sin(theta)],
        [math.sin(theta), math.cos(theta)] ])

    r00 = math.cos(theta)
    r01 = -math.sin(theta)
    r10 = math.sin(theta)
    r11 = math.cos(theta)

    #[high_var, 0.0],
    #[0.0, high_var]
    
    #Cv = Ca
    #v00 = c[0,0]
    #v01 = c[0,1]
    #v10 = c[1,0]
    #v11 = c[1,1]

    #d_vec = disp(numpy.concatenate((a,numpy.matrix([1.0]))), numpy.concatenate((b,numpy.matrix([1.0]))), T)

    ai = numpy.concatenate((a,numpy.matrix([1.0])))
    bi = numpy.concatenate((b,numpy.matrix([1.0])))
    temp = T*ai
    d_vec = bi-temp
    d_vec[2] = 1.0    

    res = Cb + R * Cv * numpy.transpose(R)

    # remove the 3rd homogeneous dimension for inverse
    invMat = scipy.linalg.inv(res)

    # add homogeneous dimension back
    invMat = numpy.concatenate((invMat,numpy.matrix([[0.0],[0.0]])), 1)
    invMat = numpy.concatenate((invMat,numpy.matrix([0.0,0.0,0.0])))

    error = numpy.transpose(d_vec)*invMat*d_vec
    errVal = error[0,0]

    return errVal

def findLocalNormal(pnt,points):

    x_list = []
    y_list = []

    pnt_count = 0
    pnts_copy = deepcopy(points)

    while pnt_count < 10:
        # select 3 closest points
        minDist = 1e100
        minPoint = []
        for p in pnts_copy:
            dist = math.sqrt((p[0]-pnt[0])**2 + (p[1]-pnt[1])**2)
        
            if dist < minDist:
                minDist = dist
                minPoint = p

        x_list.append(minPoint[0])
        y_list.append(minPoint[1])

        pnts_copy.remove(minPoint)
        pnt_count += 1

    x_list.append(pnt[0])
    y_list.append(pnt[1])

    cov_a = scipy.cov(x_list,y_list)

    loadings = []

    # NOTE:  seems to create opposing colinear vectors if data is colinear, not orthogonal vectors

    try:
        scores, loadings, E = pca_module.nipals_mat(cov_a, 2, 0.000001, False)

    except:
        raise

    if len(loadings) < 2:
        raise

    # return the second vector returned from PCA because this has the least variance (orthogonal to plane)
    return loadings[1]

def computeVectorCovariance(vec,x_var,y_var):

    c11 = x_var
    c12 = c21 = 0.0
    c22 = y_var

    mag = math.sqrt(vec[0]**2 + vec[1]**2)
    normVec = [vec[0]/mag, vec[1]/mag]

    if normVec[1] == 0:
        r22 = r11 = 1.0
        r12 = r21 = 0.0

    else:
        B = -1 / (normVec[1] + normVec[0]**2/normVec[1])
        A = -normVec[0]*B/normVec[1]
        #R = [[A, -B], [B, A]]

        r22 = r11 = A
        r12 = -B
        r21 = B

    res11 = r11*(c11*r11 + c12*r21) + r21*(c21*r11 + c22*r21)
    res12 = r11*(c11*r12 + c12*r22) + r21*(c21*r12 + c22*r22)
    res21 = r12*(c11*r11 + c12*r21) + r22*(c21*r11 + c22*r21)
    res22 = r12*(c11*r12 + c12*r22) + r22*(c21*r12 + c22*r22)
    
    Ca = [[res11, res12], [res21, res22]]
    
    return Ca


def findClosestPointInHull2(hull2_trans, p_1):

    #a_p, a_i, minDist = findClosestPointInA(a_trans, b_p, [0.0,0.0,0.0])

    minDist = 1e100
    minPoint = None
    min_i = 0

    for i in range(len(hull2_trans)):
        p = hull2_trans[i]

        dist = math.sqrt((p[0]-p_1[0])**2 + (p[1]-p_1[1])**2)
        if dist < minDist:
            minPoint = copy(p)
            minDist = dist
            min_i = i

    if minPoint != None:
        return minPoint, min_i, minDist
    else:
        raise


def findClosestPointInA(a_trans, b):

    #a_p, a_i, minDist = findClosestPointInA(a_trans, b_p, [0.0,0.0,0.0])

    minDist = 1e100
    minPoint = None
    min_i = 0

    for i in range(len(a_trans)):
        p = a_trans[i]

        dist = math.sqrt((p[0]-b[0])**2 + (p[1]-b[1])**2)
        if dist < minDist:
            minPoint = copy(p)
            minDist = dist
            min_i = i

    if minPoint != None:
        return minPoint, min_i, minDist
    else:
        raise

def findClosestPointWithAngle(points1, pnt, angleThresh):

    ax = pnt[0]
    ay = pnt[1]    
    a_theta = pnt[2]
        
    minDist = 1e100
    minPoint = None
    minI = 0

    #for p in points1:
    for i in range(len(points1)):

        p = points1[i]
        
        angDiff = abs(functions.normalizeAngle(a_theta-p[2]))
        #print "angDiff =", angDiff
        if angleThresh > angDiff or (math.pi - angleThresh) < angDiff:
    
            dist = math.sqrt((p[0]-ax)**2 + (p[1]-ay)**2)
            if dist < minDist:
                minPoint = copy(p)
                minDist = dist
                minI = i

    if minPoint != None:
        return minPoint, minI, minDist
    else:
        raise


# for point T*a, find the closest point b in B
def findClosestPointInB(b_data, a, offset):

    #xd = offset[0]
    #yd = offset[1]
    #theta = offset[2]

    #T = numpy.matrix([    [math.cos(theta), -math.sin(theta), xd],
    #        [math.sin(theta), math.cos(theta), yd],
    #        [0.0, 0.0, 1.0]
    #        ])


    #a_hom = numpy.matrix([[a[0]],[a[1]],[1.0]])
    #temp = T*a_hom
    #a_off = [temp[0,0],temp[1,0]]

    xd = offset[0]
    yd = offset[1]
    theta = offset[2]

    ax = a[0]
    ay = a[1]    
    a_off = [ax*math.cos(theta) - ay*math.sin(theta) + xd, ax*math.sin(theta) + ay*math.cos(theta) + yd]
    
    #print a_off
    
    minDist = 1e100
    minPoint = None

    for p in b_data:

        dist = math.sqrt((p[0]-a_off[0])**2 + (p[1]-a_off[1])**2)
        if dist < minDist:
            minPoint = copy(p)
            minDist = dist


    if minPoint != None:
        return minPoint, minDist
    else:
        raise


# for point T*a, find the closest point b in B
def findClosestPointInHull1(hull1, p_2, offset):

    xd = offset[0]
    yd = offset[1]
    theta = offset[2]

    ax = p_2[0]
    ay = p_2[1]    
    a_off = [ax*math.cos(theta) - ay*math.sin(theta) + xd, ax*math.sin(theta) + ay*math.cos(theta) + yd]
    
    minDist = 1e100
    minPoint = None

    for p in hull1:

        dist = math.sqrt((p[0]-a_off[0])**2 + (p[1]-a_off[1])**2)
        if dist < minDist:
            minPoint = copy(p)
            minDist = dist


    if minPoint != None:
        return minPoint, minDist
    else:
        raise

def cost_func2(offset, match_pairs, a_data_raw = [], polyB = [], circles = []):
    global numIterations
    global fig
    
    errors = []
    #    print "offset =", offset
    for pair in match_pairs:

        #print "pair:", pair
        a = pair[0]
        b = pair[1]
        
        errors.append(computeMatchError(offset, a, b, pair[2], pair[3]))
        
    #return 0,0,0
    #print "returning", errors
    return errors


def cornerMatchQuality(currAng, match_pairs, point1, point2, ang1, ang2, thresh):
    global numIterations


    def computeOffset(point1, point2, ang1, ang2):
    
        " corner points and orientations "
        corner1Pose = [point1[0], point1[1], ang1]
        corner2Pose = [point2[0], point2[1], ang2]
        
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

    offset = computeOffset(point1, point2, ang1, ang2 + currAng)

    vals = []
    sum1 = 0.0
    for pair in match_pairs:

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
    
        val = computeMatchErrorP(offset, [ax,ay], [bx,by], [c11,c12,c21,c22], [b11,b12,b21,b22])
        vals.append(val)
        sum1 += val

    """
    for pair in match_pairs:

        #print "pair:", pair
        a = pair[0]
        b = pair[1]
        
        val = computeMatchError(offset, a, b, pair[2], pair[3])
        vals.append(val)
        sum += val
    """

    upCount = 0
    downCount = 0
    for val in vals:
        if val >= thresh:
            upCount += 1
        else:
            downCount += 1
    
    #print "match errors:", vals        
    #print "match errors:", upCount, "> threshold,", downCount, "< threshold"
        
    return upCount, downCount, vals

def medialOverlapInPlaceCostFunc(params, match_pairs, medialSpline1, medialSpline2, u2, u1):
    global numIterations

    def computeOffset(point1, point2, ang1, ang2):
    
        " corner points and orientations "
        overlapPoint1Pose = [point1[0], point1[1], ang1]
        overlapPoint2Pose = [point2[0], point2[1], ang2]
        
        " convert the desired intersection point on curve 1 into global coordinates "
        poseProfile1 = Pose([0.0,0.0,0.0])
        
        " now convert this point into a pose, and perform the inverse transform using corner2Pose "
        desGlobalPose2 = Pose(overlapPoint1Pose)
        
        " perform inverse offset from the destination pose "
        negCurve2Pose = desGlobalPose2.doInverse(overlapPoint2Pose)
        
        " relative pose between pose 1 and pose 2 to make corners coincide and same angle "
        resPose2 = desGlobalPose2.convertLocalOffsetToGlobal(negCurve2Pose)
        localOffset = poseProfile1.convertGlobalPoseToLocal(resPose2)
        
        return [localOffset[0], localOffset[1], localOffset[2]]

    currU = u2
    currAng = params[0]
    #print "trying:", currU, currAng
    
    if currU < 0.1:
        return 1e5 * (0.15-currU)
    if currU > 0.9:
        return 1e5 * (currU - 0.85)

    " compute the point and angles from the paramters "
    #pose1 = medialSpline1.getUVecSet([u1, u1+0.02])[0]
    #pose2 = medialSpline2.getUVecSet([currU, currU + 0.02])[0]
    if u1+0.02 > 1.0:
        pose1 = medialSpline1.getUVecSet([0.98, 1.0])[0]
    elif u1 < 0.0:
        pose1 = medialSpline1.getUVecSet([0.0, 0.02])[0]
    else:
        pose1 = medialSpline1.getUVecSet([u1, u1+0.02])[0]

    if currU+0.02 > 1.0:
        pose2 = medialSpline2.getUVecSet([0.98, 1.0])[0]
    elif currU < 0.0:
        pose2 = medialSpline2.getUVecSet([0.0, 0.02])[0]
    else:
        pose2 = medialSpline2.getUVecSet([currU, currU + 0.02])[0]

    point1 = [pose1[0],pose1[1]]
    point2 = [pose2[0],pose2[1]]
    ang1 = pose1[2]
    ang2 = pose2[2]

    offset = computeOffset(point1, point2, ang1, ang2 + currAng)

    vals = []
    sum1 = 0.0
    for pair in match_pairs:

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
    
        val = computeMatchErrorP(offset, [ax,ay], [bx,by], [c11,c12,c21,c22], [b11,b12,b21,b22])
        
        #val += a[3] - (b[3] + offset[2])
        
        vals.append(val)
        sum1 += val
        
    return sum1

def medialOverlapTransformCostFunc(params, match_pairs, arcDist1):
    
    currArcDist = params[0]
    currAngOffset = params[1]

    """
    sum = 0.0
    for pair in match_pairs:
        p1 = pair[0]
        p2 = pair[1]
        
        dist = math.sqrt((p1[0] + arcDist1 - currArcDist - p2[0])**2+(p1[1] + currAngOffset - p2[1])**2)
        sum += dist

    print "currArcDist,currAngOffset:", currArcDist, currAngOffset
    print "match count:", len(match_pairs), sum

    return sum
    """
    
    vals = []
    sum1 = 0.0
    for pair in match_pairs:

        a = pair[0]
        b = pair[1]
        Ca = a[2]
        Cb = b[2]

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
    
        val = computeMatchErrorP([arcDist1 - currArcDist, currAngOffset, 0.0], [ax,ay], [bx,by], [c11,c12,c21,c22], [b11,b12,b21,b22])        
        vals.append(val)
        sum1 += val

    print "currArcDist,currAngOffset:", currArcDist, currAngOffset
    print "match count:", len(match_pairs), sum1
        
    return sum1



    newPoints = []
    for p in tPoints2_trans:
        newPoints.append((p[0]+arcDist1-currArcDist, p[1] + currAngOffset))

    sum1 = 0.0
    count = 0
    for p in newPoints:


        thisDist = p[0]
        index = math.floor(thisDist*10)
        try:
            pointSet = gPointHash[index]

        except:
            val = 0.0
            
        else:
            if len(pointSet) == 0:
                raise
            
            minDist = 1e100
            minI = -1
            for i in range(len(pointSet)):
                
                currDist = abs(p[0]-pointSet[i][0])
                if currDist < minDist:
                    minDist = currDist
                    minI = i
            
            val = abs(p[1] - pointSet[minI][1])

            sum1 += val
            count += 1
    
    print "match count:", count, sum1

    if count != 0:
        sum1 /= count
        return sum1

    return 1e100

def medialOverlapCostFunc_GPU(params, input_GPU, N, medialSpline1, medialSpline2, uHigh, uLow, u1):
    global numIterations

    def computeOffset(point1, point2, ang1, ang2):
    
        " corner points and orientations "
        overlapPoint1Pose = [point1[0], point1[1], ang1]
        overlapPoint2Pose = [point2[0], point2[1], ang2]
        
        " convert the desired intersection point on curve 1 into global coordinates "
        poseProfile1 = Pose([0.0,0.0,0.0])
        
        " now convert this point into a pose, and perform the inverse transform using corner2Pose "
        desGlobalPose2 = Pose(overlapPoint1Pose)
        
        " perform inverse offset from the destination pose "
        negCurve2Pose = desGlobalPose2.doInverse(overlapPoint2Pose)
        
        " relative pose between pose 1 and pose 2 to make corners coincide and same angle "
        resPose2 = desGlobalPose2.convertLocalOffsetToGlobal(negCurve2Pose)
        localOffset = poseProfile1.convertGlobalPoseToLocal(resPose2)
        
        return [localOffset[0], localOffset[1], localOffset[2]]

    currU = params[0]
    currAng = params[1]

    if currU < 0.0:
        return 1e5 * (0.15-currU)
    if currU > 1.0:
        return 1e5 * (currU - 0.85)

    if currU < uLow:
        return 1e5 * (uLow+0.05-currU)
    if currU > uHigh:
        return 1e5 * (currU - uHigh+0.05)

    " compute the point and angles from the parameters "
    if u1+0.02 > 1.0:
        pose1 = medialSpline1.getUVecSet([0.98, 1.0])[0]
    elif u1 < 0.0:
        pose1 = medialSpline1.getUVecSet([0.0, 0.02])[0]
    else:
        pose1 = medialSpline1.getUVecSet([u1, u1+0.02])[0]

    if currU+0.02 > 1.0:
        pose2 = medialSpline2.getUVecSet([0.98, 1.0])[0]
    elif currU < 0.0:
        pose2 = medialSpline2.getUVecSet([0.0, 0.02])[0]
    else:
        pose2 = medialSpline2.getUVecSet([currU, currU + 0.02])[0]
        

    point1 = [pose1[0],pose1[1]]
    point2 = [pose2[0],pose2[1]]
    ang1 = pose1[2]
    ang2 = pose2[2]

    offset = computeOffset(point1, point2, ang1, ang2 + currAng)

    cost2 = getCost_GPU(offset, input_GPU, N)

    return cost2


def convertToGPU(match_pairs):
    numPairs = len(match_pairs)
    bufSize = numPairs * 12
    inputArray = numpy.arange(bufSize).astype(numpy.float32)
    saveStr = ""

    for k in range(numPairs):
        pair = match_pairs[k]
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

        inputArray[k*12] = ax
        inputArray[k*12 + 1] = ay
        inputArray[k*12 + 2] = bx
        inputArray[k*12 + 3] = by
        inputArray[k*12 + 4] = c11
        inputArray[k*12 + 5] = c12
        inputArray[k*12 + 6] = c21
        inputArray[k*12 + 7] = c22
        inputArray[k*12 + 8] = b11
        inputArray[k*12 + 9] = b12
        inputArray[k*12 + 10] = b21
        inputArray[k*12 + 11] = b22

    input_gpu = gpuarray.to_gpu(inputArray)
    return input_gpu

def getCost_GPU(offset, inputGPU, numPairs):

    #numPairs = len(match_pairs)
    blocksPerGrid = (numPairs + threadsPerBlock - 1) / threadsPerBlock

    bufSize = numPairs * 12
    #inputArray = numpy.arange(bufSize).astype(numpy.float32)
    offsetArray = numpy.array(offset).astype(numpy.float32)
    sumVals = numpy.zeros(blocksPerGrid).astype(numpy.float32)
    #sumVals = numpy.zeros(numPairs).astype(numpy.float32)

    #computeMatchErrorG(inputArray, offsetArray, float *sumVal, len(match_pairs))
    
    #input_gpu = gpuarray.to_gpu(inputArray)

    match_them(
        inputGPU, drv.In(offsetArray), drv.Out(sumVals), numpy.int32(numPairs),
        block=(threadsPerBlock, 1, 1),
        grid=(blocksPerGrid,1)
        )    

    #print "sumVals:", sumVals
    totalSum = 0.0
    for val in sumVals:
        totalSum += val

    return totalSum


def medialOverlapCostFunc(params, match_pairs, poses_1, poses_2, uHigh, uLow, u1):
    global numIterations

    def computeOffset(point1, point2, ang1, ang2):
    
        " corner points and orientations "
        overlapPoint1Pose = [point1[0], point1[1], ang1]
        overlapPoint2Pose = [point2[0], point2[1], ang2]
        
        #print "overlapPoint1Pose =", overlapPoint1Pose
        #print "overlapPoint2Pose =", overlapPoint2Pose
        
        " convert the desired intersection point on curve 1 into global coordinates "
        poseProfile1 = Pose([0.0,0.0,0.0])
        
        " now convert this point into a pose, and perform the inverse transform using corner2Pose "
        desGlobalPose2 = Pose(overlapPoint1Pose)
        
        " perform inverse offset from the destination pose "
        negCurve2Pose = desGlobalPose2.doInverse(overlapPoint2Pose)
         
        #print "negCurve2Pose =", negCurve2Pose  
        " relative pose between pose 1 and pose 2 to make corners coincide and same angle "
        resPose2 = desGlobalPose2.convertLocalOffsetToGlobal(negCurve2Pose)
        #print "resPose2 =", resPose2  
        localOffset = poseProfile1.convertGlobalPoseToLocal(resPose2)
        
        return [localOffset[0], localOffset[1], localOffset[2]]

    currU = params[0]
    currAng = params[1]
    #print "currU:", currU
    #print "currAng:", currAng
    #print "u1:", u1
    #print "trying:", currU, currAng
    
    #if currU < 0.1:
    #    return 1e5 * (0.15-currU)
    #if currU > 0.9:
    #    return 1e5 * (currU - 0.85)

    if currU < 0.0:
        return 1e5 * (0.15-currU)
    if currU > 1.0:
        return 1e5 * (currU - 0.85)

    if currU < uLow:
        return 1e5 * (uLow+0.05-currU)
    if currU > uHigh:
        return 1e5 * (currU - uHigh+0.05)

    " compute the point and angles from the parameters "
    #u1_iter = u1 + 0.02
    #currU_iter = currU + 0.02
    #if u1_iter <= 1.0:
    #    pose1 = medialSpline1.getUVecSet([u1, u1_iter])[0]
    #else:
    #    pose1 = medialSpline1.getUVecSet([u1-0.02, u1])[0]
        
    #if currU_iter <= 1.0:
    #    pose2 = medialSpline2.getUVecSet([currU, currU_iter])[0]
    #else:
    #    pose2 = medialSpline2.getUVecSet([currU - 0.02, currU])[0]

    if u1 >= 1.0:
        pose1 = poses_1[-1]
    elif u1 < 0.0:
        pose1 = poses_1[0]
    else:
        pose1 = poses_1[int(u1*100)]

    if currU >= 1.0:
        pose2 = poses_2[-1]
    elif currU < 0.0:
        pose2 = poses_2[0]
    else:
        pose2 = poses_2[int(currU*100)]
        
    #if u1+0.02 > 1.0:
    #    pose1 = medialSpline1.getUVecSet([0.98, 1.0])[0]
    #elif u1 < 0.0:
    #    pose1 = medialSpline1.getUVecSet([0.0, 0.02])[0]
    #else:
    #    pose1 = medialSpline1.getUVecSet([u1, u1+0.02])[0]

    #if currU+0.02 > 1.0:
    #    pose2 = medialSpline2.getUVecSet([0.98, 1.0])[0]
    #elif currU < 0.0:
    #    pose2 = medialSpline2.getUVecSet([0.0, 0.02])[0]
    #else:
    #    pose2 = medialSpline2.getUVecSet([currU, currU + 0.02])[0]
        

    point1 = [pose1[0],pose1[1]]
    point2 = [pose2[0],pose2[1]]
    ang1 = pose1[2]
    ang2 = pose2[2]

    offset = computeOffset(point1, point2, ang1, ang2 + currAng)

    #print "offset =", offset

    #print "offset:", repr(offset)
    #print "match_pairs:", repr(match_pairs)
    #exit()
    vals = []
    sum1 = 0.0
    for pair in match_pairs:

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
    
        val = computeMatchErrorP(offset, [ax,ay], [bx,by], [c11,c12,c21,c22], [b11,b12,b21,b22])
        
        #val += a[3] - (b[3] + offset[2])
        
        vals.append(val)
        sum1 += val
        
    return sum1


def overlapHistogram(params, match_pairs, medialSpline1, medialSpline2, u1):
    global numIterations

    def computeOffset(point1, point2, ang1, ang2):
    
        " corner points and orientations "
        overlapPoint1Pose = [point1[0], point1[1], ang1]
        overlapPoint2Pose = [point2[0], point2[1], ang2]
        
        " convert the desired intersection point on curve 1 into global coordinates "
        poseProfile1 = Pose([0.0,0.0,0.0])
        
        " now convert this point into a pose, and perform the inverse transform using corner2Pose "
        desGlobalPose2 = Pose(overlapPoint1Pose)
        
        " perform inverse offset from the destination pose "
        negCurve2Pose = desGlobalPose2.doInverse(overlapPoint2Pose)
        
        " relative pose between pose 1 and pose 2 to make corners coincide and same angle "
        resPose2 = desGlobalPose2.convertLocalOffsetToGlobal(negCurve2Pose)
        localOffset = poseProfile1.convertGlobalPoseToLocal(resPose2)
        
        return [localOffset[0], localOffset[1], localOffset[2]]

    currU = params[0]
    currAng = params[1]
    #print "trying:", currU, currAng
    
    " compute the point and angles from the paramters "
    #pose1 = medialSpline1.getUVecSet([u1, u1+0.02])[0]
    #pose2 = medialSpline2.getUVecSet([currU, currU + 0.02])[0]
    if u1+0.02 > 1.0:
        pose1 = medialSpline1.getUVecSet([0.98, 1.0])[0]
    elif u1 < 0.0:
        pose1 = medialSpline1.getUVecSet([0.0, 0.02])[0]
    else:
        pose1 = medialSpline1.getUVecSet([u1, u1+0.02])[0]

    if currU+0.02 > 1.0:
        pose2 = medialSpline2.getUVecSet([0.98, 1.0])[0]
    elif currU < 0.0:
        pose2 = medialSpline2.getUVecSet([0.0, 0.02])[0]
    else:
        pose2 = medialSpline2.getUVecSet([currU, currU + 0.02])[0]

    point1 = [pose1[0],pose1[1]]
    point2 = [pose2[0],pose2[1]]
    ang1 = pose1[2]
    ang2 = pose2[2]

    offset = computeOffset(point1, point2, ang1, ang2 + currAng)

    vals = []
    sum1 = 0.0
    for pair in match_pairs:

        a = pair[0]
        b = pair[1]
        Ca = pair[2]
        Cb = pair[3]

        ax = a[0]
        ay = a[1]        
        bx = b[0]
        by = b[1]

        xd = offset[0]
        yd = offset[1]
        theta = offset[2]
    
        tx = ax*math.cos(theta) - ay*math.sin(theta) + xd
        ty = ax*math.sin(theta) + ay*math.cos(theta) + yd
        dx = bx - tx
        dy = by - ty
        
        val = math.sqrt(dx*dx + dy*dy)
    
        vals.append(val)
    
    lowCount = 0
    midCount = 0
    highCount = 0
    for val in vals:
        if val < 0.2:
            lowCount += 1
        elif val >= 0.2 and val <= 0.7:
            midCount += 1
        else:
            highCount += 1
    
    #print "histogram:", lowCount, midCount, highCount    
    return lowCount, midCount, highCount
    #return sum


def cornerCostFunc(currAng, match_pairs, point1, point2, ang1, ang2, a_data_raw = [], polyB = [], circles = [], isPrint = False):
    global numIterations


    def computeOffset(point1, point2, ang1, ang2):
    
        " corner points and orientations "
        corner1Pose = [point1[0], point1[1], ang1]
        corner2Pose = [point2[0], point2[1], ang2]
        
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

    offset = computeOffset(point1, point2, ang1, ang2 + currAng)

    vals = []
    sum1 = 0.0
    for pair in match_pairs:

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
    
        val = computeMatchErrorP(offset, [ax,ay], [bx,by], [c11,c12,c21,c22], [b11,b12,b21,b22])
        vals.append(val)
        sum1 += val

    if isPrint:
        thresh = 0.01
        upCount = 0
        downCount = 0
        for val in vals:
            if val >= thresh:
                upCount += 1
            else:
                downCount += 1
                
        print "match errors:", upCount, "> threshold,", downCount, "< threshold"
        
    return sum1


def shapeCost(offset, match_pairs):
    global numIterations
    global fig
    
    sum1 = 0.0
    for pair in match_pairs:

        #print "pair:", repr(pair)
        a = pair[0]
        b = pair[1]
        
        sum1 += computeMatchError(offset, a, b, pair[2], pair[3])
        
    return sum1

"""
def shapeCostC(offset, match_pairs):
    global numIterations
    global fig


    sum1 = 0.0
    for pair in match_pairs:

        #print "pair:", pair
        a = pair[0]
        b = pair[1]
        Ca = pair[2]
        Cb = pair[3]

        xd = offset[0]
        yd = offset[1]
        theta = offset[2]
        
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

        sum1 += computeMatchErrorC([xd, yd, theta], [ax,ay], [bx,by], [c11,c12,c21,c22], [b11,b12,b21,b22])
        
    return sum1
"""

def cost_func(offset, match_pairs, a_data_raw = [], polyB = [], circles = []):
    global numIterations
    global fig
    
    sum1 = 0.0
    
    for pair in match_pairs:

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
    
        sum1 += computeMatchErrorP(offset, [ax,ay], [bx,by], [c11,c12,c21,c22], [b11,b12,b21,b22])
    
    
    """
    for pair in match_pairs:

        #print "pair:", pair
        a = pair[0]
        b = pair[1]
        
        sum += computeMatchError(offset, a, b, pair[2], pair[3])
    """
    
    return sum1


def addDistanceFromOriginCovariance(points, tan_var=0.1, perp_var=0.01):

    for p in points:

        # the covariance matrix that represents increasing uncertainty as distance from origin increases
        pVec = p
        dist = math.sqrt(pVec[0]**2 + pVec[1]**2)
        
        # 1. get the orthogonal vector, rotation 90 degrees
        rotVec = [0.,0.]
        rotVec[0] = -pVec[1]
        rotVec[1] = pVec[0]

        # 2. normalize vector    
        rotVec[0] /= dist
        rotVec[1] /= dist

        # 3. computeVectorCovariance() with tan_var, perp_var proportional to dist, and tan_var larger
        C = computeVectorCovariance(rotVec,dist*tan_var,dist*perp_var)

        #print "point:", p
        #print "adding:", C
        
        if len(p) <= 2:
            p.append(C)
        else:
            p[2][0][0] += C[0][0]
            p[2][0][1] += C[0][1]
            p[2][1][0] += C[1][0]
            p[2][1][1] += C[1][1]
            
            #p[2] += C

def addGPACVectorCovariance(points, high_var=1.0, low_var=0.001):

    newPoints = []

    for p in points:
        #C = numpy.matrix([    [high_var, 0.0],
        #        [0.0, high_var]
        #        ])
        C = [ [high_var, 0.], [0., high_var]]

        # the covariance matrix that enforces the point-to-plane constraint
        normVec = [math.cos(p[2]+math.pi/2.0), math.sin(p[2]+math.pi/2.0)]        
        C = computeVectorCovariance(normVec,low_var,high_var)


        newPoints.append([p[0], p[1], C])

    return newPoints

def addGPACVectorCovarianceWithAngle(points, high_var=1.0, low_var=0.001):

    newPoints = []

    for p in points:
        #C = numpy.matrix([    [high_var, 0.0],
        #        [0.0, high_var]
        #        ])
        C = [ [high_var, 0.], [0., high_var]]

        # the covariance matrix that enforces the point-to-plane constraint
        normVec = [math.cos(p[2]+math.pi/2.0), math.sin(p[2]+math.pi/2.0)]        
        C = computeVectorCovariance(normVec,low_var,high_var)


        newPoints.append([p[0], p[1], C, p[2]])

    return newPoints

def addPointToLineCovariance(points, high_var=1.0, low_var=0.001):

    for p in points:
        C = [ [high_var, 0.], [0., high_var]]
        #C = numpy.matrix([    [high_var, 0.0],
        #        [0.0, high_var]
        #        ])

        try:
            # the covariance matrix that enforces the point-to-plane constraint
            normVec = findLocalNormal(p,points)
            C = computeVectorCovariance(normVec,low_var,high_var)

        except:
            pass

        if len(p) <= 2:
            p.append(C)
        else:
            p[2][0][0] += C[0][0]
            p[2][0][1] += C[0][1]
            p[2][1][0] += C[1][0]
            p[2][1][1] += C[1][1]
            
            #p[2] += C
    #return newPoints

" check if a point is contained in a circle "
def isInCircle(p, radius, center):
    
    dist = math.sqrt((p[0] - center[0]) ** 2 + (p[1] - center[1]) ** 2)
    if dist < radius:
        return True

    return False

" check if a point is contained in the polygon, use radius and center to quicken check "
def isValid(p, radius, center, poly):
    
    dist = math.sqrt((p[0] - center[0]) ** 2 + (p[1] - center[1]) ** 2)
    if dist < radius:
        if not functions.point_inside_polygon(p[0],p[1],poly):
            return True

    return False
    
    
def isValidPast(p, pastCircles, poly):
    
    for circle in pastCircles:
        
        radius, center = circle
        dist = math.sqrt((p[0] - center[0]) ** 2 + (p[1] - center[1]) ** 2)
    
        if dist < radius:
            if not functions.point_inside_polygon(p[0],p[1],poly):
                return True

    return False

def filterVertices(offset, radius, a_trans, b_data):

    circleA = [[offset[0],offset[1]], radius]
    circleB = [[0.0,0.0], radius]

    filteredAPoints = []
    filteredBPoints = []
    refilteredAPoints = []
    refilteredBPoints = []
    
    for p in b_data:
        dist = math.sqrt((p[0] - offset[0]) ** 2 + (p[1] - offset[1]) ** 2)
        
        if dist < radius:
            filteredBPoints.append(copy(p))

    for p in a_trans:
        dist = math.sqrt((p[0] - 0.0) ** 2 + (p[1] - 0.0) ** 2)
        
        if dist < radius:
            filteredAPoints.append(copy(p))

    for p in filteredBPoints:
        if not functions.point_inside_polygon(p[0],p[1],a_trans):
            refilteredBPoints.append(copy(p))

    for p in filteredAPoints:
        if not functions.point_inside_polygon(p[0],p[1],b_data):
            refilteredAPoints.append(copy(p))

    return refilteredAPoints, refilteredBPoints

def computeEnclosingCircle(a_data):
        
    maxA = 0.0
    a_p1 = []
    a_p2 = []
    for i in range(len(a_data)):
        p1 = a_data[i]

        for j in range(i+1, len(a_data)):
            p2 = a_data[j]

            dist = math.sqrt((p1[0]-p2[0])**2 + (p1[1]-p2[1])**2)

            if dist > maxA:
                maxA = dist
                a_p1 = copy(p1)
                a_p2 = copy(p2)

    radiusA = maxA/2.0

    centerA = [(a_p1[0] + a_p2[0])/2.0, (a_p1[1] + a_p2[1])/2.0]
    
    return radiusA, centerA

def gen_ICP(pastPose, targetPose, pastHull, targetHull, pastCircles, costThresh = 0.004, minMatchDist = 2.0, plotIter = False):

    global numIterations
    
    lastCost = 1e100
    
    startIteration = numIterations

    estPoseOrigin = pastPose
    
    #print "estPoseOrigin =", estPoseOrigin
    
    " set the initial guess "
    estPose2 = targetPose
    #print "estPose2 =", estPose2
    
    poseOrigin = Pose(estPoseOrigin)
    
    offset = poseOrigin.convertGlobalPoseToLocal(estPose2)
    #print "offset =", offset

    #exit()
    " transform the past poses "
    a_data_raw = targetHull
    a_data = []
    for p in a_data_raw:
        result = dispPoint(p, offset)
        #print offset, p, result
        #print
        
        a_data.append(result)
    
    #print "targetHull =", targetHull
    #print "a_trans =", a_trans
    

    
    polyB = []        
    for p in pastHull:
        polyB.append([p[0],p[1]])    

    while True:
        " find the matching pairs "
        match_pairs = []

        a_data_raw = targetHull
        b_data = pastHull            

        " transform the target Hull with the latest offset "
        a_data = []
        for p in a_data_raw:
            result = dispPoint(p, offset)
            a_data.append(result)

        " transformed points without associated covariance "
        polyA = []
        for p in a_data:
            polyA.append([p[0],p[1]])    
        
        " get the circles and radii "
        radiusA, centerA = computeEnclosingCircle(a_data)
        #radiusB, centerB = computeEnclosingCircle(pastHull)
        
        if True:
            for i in range(len(a_data)):
                a_p = polyA[i]
    
                #if isValid(a_p, radiusB, centerB, polyB):
                if isValidPast(a_p, pastCircles, polyB):
    
                    " for every transformed point of A, find it's closest neighbor in B "
                    b_p, minDist = findClosestPointInB(b_data, a_p, [0.0,0.0,0.0])
        
                    if minDist <= minMatchDist:
            
                        " add to the list of match pairs less than 1.0 distance apart "
                        " keep A points and covariances untransformed "
                        Ca = a_data_raw[i][2]
                        Cb = b_p[2]
        
                        " we store the untransformed point, but the transformed covariance of the A point "
                        match_pairs.append([a_data_raw[i],b_p,Ca,Cb])

        if True:
            for i in range(len(b_data)):
                b_p = polyB[i]
        
                if isValid(b_p, radiusA, centerA, polyA):
            
                    #print "selected", b_p, "in circle", radiusA, centerA, "with distance"
                    " for every point of B, find it's closest neighbor in transformed A "
                    a_p, a_i, minDist = findClosestPointInA(a_data, b_p)
        
                    if minDist <= minMatchDist:
            
                        " add to the list of match pairs less than 1.0 distance apart "
                        " keep A points and covariances untransformed "
                        Ca = a_data_raw[a_i][2]
                        
                        Cb = b_data[i][2]
        
                        " we store the untransformed point, but the transformed covariance of the A point "
                        match_pairs.append([a_data_raw[a_i],b_p,Ca,Cb])


        " plot the initial configuration "
        if startIteration == numIterations:

            numIterations += 1

            pylab.clf()
            pylab.axes()
            match_global = []
            
            for pair in match_pairs:
                p1 = pair[0]
                p2 = pair[1]
                
                p1_o = dispOffset(p1, offset)
                #p2_o = dispOffset(p2, offset)

                
                p1_g = poseOrigin.convertLocalToGlobal(p1_o)
                p2_g = poseOrigin.convertLocalToGlobal(p2)
                match_global.append([p1_g,p2_g])
            
            draw_matches(match_global, [0.0,0.0,0.0])
            
            xP = []
            yP = []
            for b in polyB:
                p1 = poseOrigin.convertLocalToGlobal(b)

                xP.append(p1[0])    
                yP.append(p1[1])
            
            p1 = poseOrigin.convertLocalToGlobal(polyB[0])
            xP.append(p1[0])    
            yP.append(p1[1])
            
            pylab.plot(xP,yP,linewidth=1, color=(0.0,0.0,1.0))

            xP = []
            yP = []
            for b in a_data_raw:
                p = [b[0],b[1]]
                p = dispOffset(p,offset)
                
                p1 = poseOrigin.convertLocalToGlobal(p)
                xP.append(p1[0])    
                yP.append(p1[1])
            
            p = [a_data_raw[0][0],a_data_raw[0][1]]
            p = dispOffset(p,offset)

            p1 = poseOrigin.convertLocalToGlobal(p)
            xP.append(p1[0])    
            yP.append(p1[1])

            pylab.plot(xP,yP,linewidth=1, color=(0.0,0.0,1.0))

            plotEnv()        
            
            pylab.xlim(-4,9)
            pylab.ylim(-3,9)
            #pylab.ylim(-7,5)
            pylab.savefig("ICP_plot_%04u.png" % numIterations)
            pylab.clf()                    
            
        if False:
            for pair in match_pairs:
                a = pair[0]
                b = pair[1]
            
                pylab.plot([a[0],b[0]], [a[1],b[1]],color=(1.0,0.0,0.0))
            
            pylab.show()
            
        #allCircles = deepcopy(pastCircles)
        #allCircles.append([radiusA,centerA])

        allCircles = [[radiusA,centerA]]

        # optimize the match error for the current list of match pairs
        newOffset = scipy.optimize.fmin(cost_func, offset, [match_pairs, a_data_raw, polyB, allCircles])
    
        # get the current cost
        newCost = cost_func(newOffset, match_pairs)
    
        # check for convergence condition, different between last and current cost is below threshold
        if abs(lastCost - newCost) < costThresh or (numIterations - startIteration) > 10:
            offset = newOffset
            lastCost = newCost
            break
    
        numIterations += 1

        offset = newOffset

        " reduce the minMatch distance for each step down to a floor value "
        minMatchDist /= 2
        if minMatchDist < 0.25:
            minMatchDist = 0.25
    
        # optionally draw the position of the points in current transform
        if plotIter:
            pylab.clf()
            pylab.axes()
            match_global = []
            
            #match_pairs.append([a_data_raw[a_i],b_p,Ca,Cb])
            
            for pair in match_pairs:
                p1 = pair[0]
                p2 = pair[1]
                
                p1_o = dispOffset(p1, offset)
                #p2_o = dispOffset(p2, offset)

                
                p1_g = poseOrigin.convertLocalToGlobal(p1_o)
                p2_g = poseOrigin.convertLocalToGlobal(p2)
                match_global.append([p1_g,p2_g])
            
            #draw_matches(match_pairs, offset)
            draw_matches(match_global, [0.0,0.0,0.0])
            
            xP = []
            yP = []
            for b in polyB:
                p1 = poseOrigin.convertLocalToGlobal(b)

                xP.append(p1[0])    
                yP.append(p1[1])
                #xP.append(b[0])    
                #yP.append(b[1])
            
            p1 = poseOrigin.convertLocalToGlobal(polyB[0])
            xP.append(p1[0])    
            yP.append(p1[1])
            #xP.append(polyB[0][0])    
            #yP.append(polyB[0][1])
            
            pylab.plot(xP,yP,linewidth=1, color=(0.0,0.0,1.0))

            xP = []
            yP = []
            for b in a_data_raw:
                p = [b[0],b[1]]
                p = dispOffset(p,offset)
                
                p1 = poseOrigin.convertLocalToGlobal(p)
                xP.append(p1[0])    
                yP.append(p1[1])
                
                #xP.append(p[0])    
                #yP.append(p[1])
            
            p = [a_data_raw[0][0],a_data_raw[0][1]]
            p = dispOffset(p,offset)

            p1 = poseOrigin.convertLocalToGlobal(p)
            xP.append(p1[0])    
            yP.append(p1[1])

            #xP.append(p[0])    
            #yP.append(p[1])
            
            pylab.plot(xP,yP,linewidth=1, color=(0.0,0.0,1.0))

            #cir = Circle( (0,0), radius=0.5)
            #a.add_patch(cir)

            plotEnv()        
            
            #for circle in pastCircles:
            #    radius, center = circle
            #    cir = pylab.Circle((center[0],center[1]), radius=radius,  fc='r')
            #    pylab.gca().add_patch(cir)
                
            #cir = pylab.Circle((centerA[0],centerA[1]), radius=radiusA,  fc='r')
            #pylab.gca().add_patch(cir)
            #cir = pylab.Circle((.5,.5), radius=0.25, alpha =.2, fc='b')
            #pylab.gca().add_patch(cir)

            #pylab.xlim(-4.5,4.5)
            #pylab.ylim(-4,4)
            #pylab.xlim(-3,6)
            #pylab.ylim(-4,4)
            #pylab.xlim(-4,7)
            #pylab.ylim(-7,3)
            #pylab.xlim(-4,9)
            #pylab.ylim(-7,5)
            pylab.xlim(-4,9)
            pylab.ylim(-3,9)
            #pylab.axis('equal')
            pylab.savefig("ICP_plot_%04u.png" % numIterations)
            pylab.clf()            
                        
        # save the current offset and cost
        offset = newOffset
        lastCost = newCost

    return offset

def gen_ICP2(estPose1, offset, pastHull, targetHull, pastCircles, costThresh = 0.004, minMatchDist = 2.0, plotIter = False, n1 = 0, n2 = 0):

    global numIterations
    
    lastCost = 1e100
    
    startIteration = numIterations

    " set the initial guess "
    poseOrigin = Pose(estPose1)
    
    #offset = poseOrigin.convertGlobalPoseToLocal(estPose2)
    #print "offset =", offset

    #exit()
    " transform the past poses "
    a_data_raw = targetHull
    a_data = []
    for p in a_data_raw:
        result = dispPoint(p, offset)
        #print offset, p, result
        #print
        
        a_data.append(result)
    
    #print "targetHull =", targetHull
    #print "a_trans =", a_trans
    

    
    polyB = []        
    for p in pastHull:
        polyB.append([p[0],p[1]])    

    while True:
        " find the matching pairs "
        match_pairs = []

        a_data_raw = targetHull
        b_data = pastHull            

        " transform the target Hull with the latest offset "
        a_data = []
        for p in a_data_raw:
            result = dispPoint(p, offset)
            a_data.append(result)

        " transformed points without associated covariance "
        polyA = []
        for p in a_data:
            polyA.append([p[0],p[1]])    
        
        " get the circles and radii "
        radiusA, centerA = computeEnclosingCircle(a_data)
        #radiusB, centerB = computeEnclosingCircle(pastHull)
        
        if True:
            for i in range(len(a_data)):
                a_p = polyA[i]
    
                #if isValid(a_p, radiusB, centerB, polyB):
                if isValidPast(a_p, pastCircles, polyB):
    
                    " for every transformed point of A, find it's closest neighbor in B "
                    b_p, minDist = findClosestPointInB(b_data, a_p, [0.0,0.0,0.0])
        
                    if minDist <= minMatchDist:
            
                        " add to the list of match pairs less than 1.0 distance apart "
                        " keep A points and covariances untransformed "
                        Ca = a_data_raw[i][2]
                        Cb = b_p[2]
        
                        " we store the untransformed point, but the transformed covariance of the A point "
                        match_pairs.append([a_data_raw[i],b_p,Ca,Cb])

        if True:
            for i in range(len(b_data)):
                b_p = polyB[i]
        
                if isValid(b_p, radiusA, centerA, polyA):
            
                    #print "selected", b_p, "in circle", radiusA, centerA, "with distance"
                    " for every point of B, find it's closest neighbor in transformed A "
                    a_p, a_i, minDist = findClosestPointInA(a_data, b_p)
        
                    if minDist <= minMatchDist:
            
                        " add to the list of match pairs less than 1.0 distance apart "
                        " keep A points and covariances untransformed "
                        Ca = a_data_raw[a_i][2]
                        
                        Cb = b_data[i][2]
        
                        " we store the untransformed point, but the transformed covariance of the A point "
                        match_pairs.append([a_data_raw[a_i],b_p,Ca,Cb])


        " plot the initial configuration "
        if plotIter and startIteration == numIterations:

            numIterations += 1

            pylab.clf()
            pylab.axes()
            match_global = []
            
            for pair in match_pairs:
                p1 = pair[0]
                p2 = pair[1]
                
                p1_o = dispOffset(p1, offset)
                #p2_o = dispOffset(p2, offset)

                
                p1_g = poseOrigin.convertLocalToGlobal(p1_o)
                p2_g = poseOrigin.convertLocalToGlobal(p2)
                match_global.append([p1_g,p2_g])
            
            draw_matches(match_global, [0.0,0.0,0.0])
            
            xP = []
            yP = []
            for b in polyB:
                p1 = poseOrigin.convertLocalToGlobal(b)

                xP.append(p1[0])    
                yP.append(p1[1])
            
            p1 = poseOrigin.convertLocalToGlobal(polyB[0])
            xP.append(p1[0])    
            yP.append(p1[1])
            
            pylab.plot(xP,yP,linewidth=1, color=(0.0,0.0,1.0))

            xP = []
            yP = []
            for b in a_data_raw:
                p = [b[0],b[1]]
                p = dispOffset(p,offset)
                
                p1 = poseOrigin.convertLocalToGlobal(p)
                xP.append(p1[0])    
                yP.append(p1[1])
            
            p = [a_data_raw[0][0],a_data_raw[0][1]]
            p = dispOffset(p,offset)

            p1 = poseOrigin.convertLocalToGlobal(p)
            xP.append(p1[0])    
            yP.append(p1[1])

            pylab.plot(xP,yP,linewidth=1, color=(0.0,0.0,1.0))

            plotEnv()        
            
            pylab.xlim(-4,9)
            pylab.ylim(-3,9)
            pylab.title("%u %u" % (n1, n2))
            #pylab.ylim(-7,5)
            pylab.savefig("ICP_plot_%04u.png" % numIterations)
            pylab.clf()                    
            
        if False:
            for pair in match_pairs:
                a = pair[0]
                b = pair[1]
            
                pylab.plot([a[0],b[0]], [a[1],b[1]],color=(1.0,0.0,0.0))
            
            pylab.show()
            
        #allCircles = deepcopy(pastCircles)
        #allCircles.append([radiusA,centerA])

        allCircles = [[radiusA,centerA]]

        # optimize the match error for the current list of match pairs
        newOffset = scipy.optimize.fmin(cost_func, offset, [match_pairs, a_data_raw, polyB, allCircles], disp = 0)
    
    
    
        # get the current cost
        newCost = cost_func(newOffset, match_pairs)
    
        # check for convergence condition, different between last and current cost is below threshold
        if abs(lastCost - newCost) < costThresh or (numIterations - startIteration) > 10:
            offset = newOffset
            lastCost = newCost
            break
    
        numIterations += 1

        offset = newOffset

        " reduce the minMatch distance for each step down to a floor value "
        minMatchDist /= 2
        if minMatchDist < 0.25:
            minMatchDist = 0.25
    
        # optionally draw the position of the points in current transform
        if plotIter:
            pylab.clf()
            pylab.axes()
            match_global = []
            
            #match_pairs.append([a_data_raw[a_i],b_p,Ca,Cb])
            
            for pair in match_pairs:
                p1 = pair[0]
                p2 = pair[1]
                
                p1_o = dispOffset(p1, offset)
                #p2_o = dispOffset(p2, offset)

                
                p1_g = poseOrigin.convertLocalToGlobal(p1_o)
                p2_g = poseOrigin.convertLocalToGlobal(p2)
                match_global.append([p1_g,p2_g])
            
            #draw_matches(match_pairs, offset)
            draw_matches(match_global, [0.0,0.0,0.0])
            
            xP = []
            yP = []
            for b in polyB:
                p1 = poseOrigin.convertLocalToGlobal(b)

                xP.append(p1[0])    
                yP.append(p1[1])
                #xP.append(b[0])    
                #yP.append(b[1])
            
            p1 = poseOrigin.convertLocalToGlobal(polyB[0])
            xP.append(p1[0])    
            yP.append(p1[1])
            #xP.append(polyB[0][0])    
            #yP.append(polyB[0][1])
            
            pylab.plot(xP,yP,linewidth=1, color=(0.0,0.0,1.0))

            xP = []
            yP = []
            for b in a_data_raw:
                p = [b[0],b[1]]
                p = dispOffset(p,offset)
                
                p1 = poseOrigin.convertLocalToGlobal(p)
                xP.append(p1[0])    
                yP.append(p1[1])
                
                #xP.append(p[0])    
                #yP.append(p[1])
            
            p = [a_data_raw[0][0],a_data_raw[0][1]]
            p = dispOffset(p,offset)

            p1 = poseOrigin.convertLocalToGlobal(p)
            xP.append(p1[0])    
            yP.append(p1[1])

            #xP.append(p[0])    
            #yP.append(p[1])
            
            pylab.plot(xP,yP,linewidth=1, color=(0.0,0.0,1.0))

            #cir = Circle( (0,0), radius=0.5)
            #a.add_patch(cir)

            plotEnv()        
            
            #for circle in pastCircles:
            #    radius, center = circle
            #    cir = pylab.Circle((center[0],center[1]), radius=radius,  fc='r')
            #    pylab.gca().add_patch(cir)
                
            #cir = pylab.Circle((centerA[0],centerA[1]), radius=radiusA,  fc='r')
            #pylab.gca().add_patch(cir)
            #cir = pylab.Circle((.5,.5), radius=0.25, alpha =.2, fc='b')
            #pylab.gca().add_patch(cir)

            #pylab.xlim(-4.5,4.5)
            #pylab.ylim(-4,4)
            #pylab.xlim(-3,6)
            #pylab.ylim(-4,4)
            #pylab.xlim(-4,7)
            #pylab.ylim(-7,3)
            #pylab.xlim(-4,9)
            #pylab.ylim(-7,5)
            pylab.title("%u %u" % (n1, n2))
            pylab.xlim(-4,9)
            pylab.ylim(-3,9)
            #pylab.axis('equal')
            pylab.savefig("ICP_plot_%04u.png" % numIterations)
            pylab.clf()            
                        
        # save the current offset and cost
        offset = newOffset
        lastCost = newCost


    if True:
        pylab.clf()
        pylab.axes()
        match_global = []
        
        for pair in match_pairs:
            p1 = pair[0]
            p2 = pair[1]
            
            p1_o = dispOffset(p1, offset)

            p1_g = poseOrigin.convertLocalToGlobal(p1_o)
            p2_g = poseOrigin.convertLocalToGlobal(p2)
            match_global.append([p1_g,p2_g])
        
        draw_matches(match_global, [0.0,0.0,0.0])
        
        xP = []
        yP = []
        for b in polyB:
            p1 = poseOrigin.convertLocalToGlobal(b)

            xP.append(p1[0])    
            yP.append(p1[1])
        
        p1 = poseOrigin.convertLocalToGlobal(polyB[0])
        xP.append(p1[0])    
        yP.append(p1[1])
        
        pylab.plot(xP,yP,linewidth=1, color=(0.0,0.0,1.0))

        xP = []
        yP = []
        for b in a_data_raw:
            p = [b[0],b[1]]
            p = dispOffset(p,offset)
            
            p1 = poseOrigin.convertLocalToGlobal(p)
            xP.append(p1[0])    
            yP.append(p1[1])
            
        p = [a_data_raw[0][0],a_data_raw[0][1]]
        p = dispOffset(p,offset)

        p1 = poseOrigin.convertLocalToGlobal(p)
        xP.append(p1[0])    
        yP.append(p1[1])

        pylab.plot(xP,yP,linewidth=1, color=(1.0,0.0,0.0))

        pylab.title("%u %u" % (n1, n2))
        pylab.xlim(-4,9)
        pylab.ylim(-3,9)
        pylab.savefig("ICP_plot_%04u.png" % numIterations)
        pylab.clf()            



    offset[2] =  functions.normalizeAngle(offset[2])
    return offset, lastCost


            
def shapeICP2(estPose1, estPose2, hull1, hull2, medial1, medial2, posture1_unstable, posture1_stable, posture2_unstable, posture2_stable, costThresh = 0.004, minMatchDist = 2.0, plotIter = False, n1 = 0, n2 = 0):
    
    ACCEPTED = 0
    REJECTED_CONTAIN = 1
    
    status = ACCEPTED
    
    poly1 = []
    for p in hull1:
        poly1.append([p[0],p[1]])
    poly2 = []
    for p in hull2:
        poly2.append([p[0],p[1]])

    count1_unstable = 0
    count1_stable = 0
    for p in posture1_unstable:
        if functions.point_inside_polygon(p[0],p[1],poly1):
            count1_unstable += 1
    for p in posture1_stable:
        if functions.point_inside_polygon(p[0],p[1],poly1):
            count1_stable += 1

    posture1 = posture1_stable
    if count1_stable < 38:
        if count1_unstable > count1_stable:
            posture1 = posture2_unstable

    count2_unstable = 0
    count2_stable = 0
    for p in posture2_unstable:
        if functions.point_inside_polygon(p[0],p[1],poly2):
            count2_unstable += 1
    for p in posture2_stable:
        if functions.point_inside_polygon(p[0],p[1],poly2):
            count2_stable += 1

    posture2 = posture2_stable
    if count2_stable < 38:
        if count2_unstable > count2_stable:
            posture2 = posture2_unstable

    curve1 = StableCurve(posture1)
    curve2 = StableCurve(posture2)
    uniform1 = curve1.getPlot()
    uniform2 = curve2.getPlot()

    headPoint1, tailPoint1 = curve1.getTips()
    headPoint2, tailPoint2 = curve2.getTips()
    
    stablePoints1 = curve1.getPoints()
    stablePoints2 = curve2.getPoints()
    
    radius, center = computeEnclosingCircle(hull1)
    circle1 = [radius, center]
    

    def computeConstraintOrientation(offset, p1_stable, p2_stable):

        def dispOffset(p, offset):
            xd = offset[0]
            yd = offset[1]
            theta = offset[2]
        
            px = p[0]
            py = p[1]
            pa = p[2]
            
            newAngle = functions.normalizeAngle(pa + theta)
            
            p_off = [px*math.cos(theta) - py*math.sin(theta) + xd, px*math.sin(theta) + py*math.cos(theta) + yd, newAngle]
            
            return p_off

        
        #node1 = self.nodeHash[n1]
        #node2 = self.nodeHash[n2]
                                
        #curve1 = node1.getStableGPACurve()
        #curve2 = node2.getStableGPACurve()

        curve1 = SplineFit(p1_stable, smooth = 0.5, kp = 2)
        curve2 = SplineFit(p2_stable, smooth = 0.5, kp = 2)

        points1 = curve1.getUniformSamples()
        points2 = curve2.getUniformSamples()

        " transform the new pose "
        points2_offset = []
        for p in points2:
            result = dispOffset(p, offset)        
            points2_offset.append(result)
        
        poly1 = []
        for p in points1:
            poly1.append([p[0],p[1]])    

        poly2 = []
        for p in points2_offset:
            poly2.append([p[0],p[1]])    
        
        " find the matching pairs "
        match_pairs = []
            
        " get the circles and radii "
        radius2, center2 = computeEnclosingCircle(points2_offset)
        radius1, center1 = computeEnclosingCircle(points1)

        minMatchDist = 0.05
            
        for i in range(len(points2_offset)):
            p_2 = poly2[i]

            if isValid(p_2, radius1, center1, poly1):

                " for every transformed point of A, find it's closest neighbor in B "
                p_1, minDist = findClosestPointInHull1(points1, p_2, [0.0,0.0,0.0])
    
                if minDist <= minMatchDist:

                    " we store the untransformed point, but the transformed covariance of the A point "
                    match_pairs.append([points2_offset[i],p_1])

        for i in range(len(points1)):
            p_1 = poly1[i]
    
            if isValid(p_1, radius2, center2, poly2):
        
                #print "selected", b_p, "in circle", radiusA, centerA, "with distance"
                " for every point of B, find it's closest neighbor in transformed A "
                p_2, i_2, minDist = findClosestPointInHull2(points2_offset, p_1)
    
                if minDist <= minMatchDist:
    
                    " we store the untransformed point, but the transformed covariance of the A point "
                    match_pairs.append([points2_offset[i_2],points1[i]])
                    
        " if we made no matches, default to regular orientation "
        if len(match_pairs) == 0:
            return True, 0.0

        distances = []
        angSum = 0
        for i in range(len(match_pairs)):
            pair = match_pairs[i]
            p1 = pair[0]
            p2 = pair[1]
            
            angDist = functions.normalizeAngle(p1[2]-p2[2])
            distances.append(angDist)
            angSum += angDist

        
        angDiff = angSum / len(distances)


        if math.fabs(angDiff) >= 0.7:
            return False, angDiff
        else:
            return True, angDiff

    firstGuess = [0.0,0.0,0.0]

    isForward, angDiff = computeConstraintOrientation(firstGuess, posture1_stable, posture2_stable)
    
    termPoints = (headPoint1,tailPoint1,headPoint2,tailPoint2)
    pastCircles = [circle1]
    offset = firstGuess


    global numIterations
    
    lastCost = 1e100
    
    startIteration = numIterations

    " set the initial guess "
    poseOrigin = Pose(estPose1)
    
    " transform the past poses "
    hull2_trans = []
    for p in hull2:
        result = dispPoint(p, offset)
        hull2_trans.append(result)
        
    poly1 = []        
    for p in hull1:
        poly1.append([p[0],p[1]])    
    bbox1 = computeBBox(poly1)

    poly2 = []        
    for p in hull2:
        poly2.append([p[0],p[1]])        
    bbox2 = computeBBox(poly2)

    

    breakOut = False

    #medial1 = getMedialAxis(hull1)
    #medial2 = getMedialAxis(hull2)

    

    " check if the medial axis is ordered correctly "
    #distFore1 = (medial1[1][0]-uniform1[1][0])**2 + (medial1[1][1]-uniform1[1][1])**2
    #distBack1 = (medial1[-2][0]-uniform1[-2][0])**2 + (medial1[-2][1]-uniform1[-2][1])**2

    #if distFore1 > 1 and distBack1 > 1:
    #    medial1.reverse()

    #distFore2 = (medial2[1][0]-uniform2[1][0])**2 + (medial2[1][1]-uniform2[1][1])**2
    #distBack2 = (medial2[-2][0]-uniform2[-2][0])**2 + (medial2[-2][1]-uniform2[-2][1])**2

    #if distFore2 > 1 and distBack2 > 1:
    #    medial2.reverse()



    halfHull1A, halfHull1B = computeSeparation2(bbox1, medial1, hull1)
    poly1A = []        
    for p in halfHull1A:
        poly1A.append([p[0],p[1]])    
    poly1B = []        
    for p in halfHull1B:
        poly1B.append([p[0],p[1]])    

    halfHull2A, halfHull2B = computeSeparation2(bbox2, medial2, hull2)
    poly2A = []        
    for p in halfHull2A:
        poly2A.append([p[0],p[1]])    
    poly2B = []        
    for p in halfHull2B:
        poly2B.append([p[0],p[1]])    


    while True:
        " find the matching pairs "
        match_pairs = []

        #b_data = hull1            

        " transform the target Hull with the latest offset "
        hull2_trans = []
        for p in hull2:
            result = dispPoint(p, offset)
            hull2_trans.append(result)

        " transformed points without associated covariance "
        poly2_trans = []
        for p in hull2_trans:
            poly2_trans.append([p[0],p[1]])    

        bbox2_trans = computeBBox(poly2_trans)

        uniform2_trans = []
        for b in uniform2:
            p = dispOffset(b, offset)
            uniform2_trans.append(p)

        medial2_trans = []
        for b in medial2:
            p = dispOffset(b, offset)
            medial2_trans.append(p)

        halfHull2A_trans = []
        for b in halfHull2A:
            p = dispOffset(b, offset)
            halfHull2A_trans.append(p)

        halfHull2B_trans = []
        for b in halfHull2B:
            p = dispOffset(b, offset)
            halfHull2B_trans.append(p)

        poly2A_trans = []        
        for p in halfHull2A_trans:
            poly2A_trans.append([p[0],p[1]])    
        poly2B_trans = []        
        for p in halfHull2B_trans:
            poly2B_trans.append([p[0],p[1]])    


        """
        try:    

            #halfHull2A_trans, halfHull2B_trans = computeSeparation(bbox2_trans, uniform2_trans, hull2_trans)
            halfHull2A_trans, halfHull2B_trans = computeSeparation(bbox2_trans, medial2_trans, medial2_trans)

            poly2A_trans = []        
            for p in halfHull2A_trans:
                poly2A_trans.append([p[0],p[1]])    
            poly2B_trans = []        
            for p in halfHull2B_trans:
                poly2B_trans.append([p[0],p[1]])    

        except ValueError:
            print "Exception caught"
            breakOut = True
        """

        #print "processing continued"
        
        " get the circles and radii "
        radius2, center2 = computeEnclosingCircle(hull2_trans)
        #radiusB, centerB = computeEnclosingCircle(pastHull)


        if isForward:
            for i in range(len(halfHull2A_trans)):
                a_p = poly2A_trans[i]
    
                if isValidPast(a_p, pastCircles, poly1A):
    
                    " for every transformed point of A, find it's closest neighbor in B "
                    #b_p, minDist = findClosestPointInB(hull1, a_p, [0.0,0.0,0.0])
                    b_p, minDist = findClosestPointInHull1(halfHull1A, a_p, [0.0,0.0,0.0])
        
                    if minDist <= minMatchDist:
            
                        " add to the list of match pairs less than 1.0 distance apart "
                        " keep A points and covariances untransformed "
                        Ca = halfHull2A[i][2]
                        Cb = b_p[2]
        
                        " we store the untransformed point, but the transformed covariance of the A point "
                        match_pairs.append([halfHull2A[i],b_p,Ca,Cb])
    
            for i in range(len(halfHull1A)):
                b_p = poly1A[i]
        
                if isValid(b_p, radius2, center2, poly2A_trans):
            
                    #print "selected", b_p, "in circle", radiusA, centerA, "with distance"
                    " for every point of B, find it's closest neighbor in transformed A "
                    #a_p, a_i, minDist = findClosestPointInA(hull2_trans, b_p)
                    a_p, a_i, minDist = findClosestPointInHull2(halfHull2A_trans, b_p)
        
                    if minDist <= minMatchDist:
            
                        " add to the list of match pairs less than 1.0 distance apart "
                        " keep A points and covariances untransformed "
                        Ca = halfHull2A[a_i][2]
                        Cb = halfHull1A[i][2]
        
                        " we store the untransformed point, but the transformed covariance of the A point "
                        match_pairs.append([halfHull2A[a_i],b_p,Ca,Cb])
    
            for i in range(len(halfHull2B_trans)):
                a_p = poly2B_trans[i]
    
                if isValidPast(a_p, pastCircles, poly1B):
    
                    " for every transformed point of A, find it's closest neighbor in B "
                    #b_p, minDist = findClosestPointInB(hull1, a_p, [0.0,0.0,0.0])
                    b_p, minDist = findClosestPointInHull1(halfHull1B, a_p, [0.0,0.0,0.0])
        
                    if minDist <= minMatchDist:
            
                        " add to the list of match pairs less than 1.0 distance apart "
                        " keep A points and covariances untransformed "
                        Ca = halfHull2B[i][2]
                        Cb = b_p[2]
        
                        " we store the untransformed point, but the transformed covariance of the A point "
                        match_pairs.append([halfHull2B[i],b_p,Ca,Cb])
    
            for i in range(len(halfHull1B)):
                b_p = poly1B[i]
        
                if isValid(b_p, radius2, center2, poly2B_trans):
            
                    #print "selected", b_p, "in circle", radiusA, centerA, "with distance"
                    " for every point of B, find it's closest neighbor in transformed A "
                    #a_p, a_i, minDist = findClosestPointInA(hull2_trans, b_p)
                    a_p, a_i, minDist = findClosestPointInHull2(halfHull2B_trans, b_p)
        
                    if minDist <= minMatchDist:
            
                        " add to the list of match pairs less than 1.0 distance apart "
                        " keep A points and covariances untransformed "
                        Ca = halfHull2B[a_i][2]
                        Cb = halfHull1B[i][2]
        
                        " we store the untransformed point, but the transformed covariance of the A point "
                        match_pairs.append([halfHull2B[a_i],b_p,Ca,Cb])
        else:

            for i in range(len(halfHull2A_trans)):
                a_p = poly2A_trans[i]
    
                if isValidPast(a_p, pastCircles, poly1B):
    
                    " for every transformed point of A, find it's closest neighbor in B "
                    #b_p, minDist = findClosestPointInB(hull1, a_p, [0.0,0.0,0.0])
                    b_p, minDist = findClosestPointInHull1(halfHull1B, a_p, [0.0,0.0,0.0])
        
                    if minDist <= minMatchDist:
            
                        " add to the list of match pairs less than 1.0 distance apart "
                        " keep A points and covariances untransformed "
                        Ca = halfHull2A[i][2]
                        Cb = b_p[2]
        
                        " we store the untransformed point, but the transformed covariance of the A point "
                        match_pairs.append([halfHull2A[i],b_p,Ca,Cb])
    
            for i in range(len(halfHull1B)):
                b_p = poly1B[i]
        
                if isValid(b_p, radius2, center2, poly2A_trans):
            
                    #print "selected", b_p, "in circle", radiusA, centerA, "with distance"
                    " for every point of B, find it's closest neighbor in transformed A "
                    #a_p, a_i, minDist = findClosestPointInA(hull2_trans, b_p)
                    a_p, a_i, minDist = findClosestPointInHull2(halfHull2A_trans, b_p)
        
                    if minDist <= minMatchDist:
            
                        " add to the list of match pairs less than 1.0 distance apart "
                        " keep A points and covariances untransformed "
                        Ca = halfHull2A[a_i][2]
                        Cb = halfHull1B[i][2]
        
                        " we store the untransformed point, but the transformed covariance of the A point "
                        match_pairs.append([halfHull2A[a_i],b_p,Ca,Cb])
    
            for i in range(len(halfHull2B_trans)):
                a_p = poly2B_trans[i]
    
                if isValidPast(a_p, pastCircles, poly1A):
    
                    " for every transformed point of A, find it's closest neighbor in B "
                    #b_p, minDist = findClosestPointInB(hull1, a_p, [0.0,0.0,0.0])
                    b_p, minDist = findClosestPointInHull1(halfHull1A, a_p, [0.0,0.0,0.0])
        
                    if minDist <= minMatchDist:
            
                        " add to the list of match pairs less than 1.0 distance apart "
                        " keep A points and covariances untransformed "
                        Ca = halfHull2B[i][2]
                        Cb = b_p[2]
        
                        " we store the untransformed point, but the transformed covariance of the A point "
                        match_pairs.append([halfHull2B[i],b_p,Ca,Cb])
    
            for i in range(len(halfHull1A)):
                b_p = poly1A[i]
        
                if isValid(b_p, radius2, center2, poly2B_trans):
            
                    #print "selected", b_p, "in circle", radiusA, centerA, "with distance"
                    " for every point of B, find it's closest neighbor in transformed A "
                    #a_p, a_i, minDist = findClosestPointInA(hull2_trans, b_p)
                    a_p, a_i, minDist = findClosestPointInHull2(halfHull2B_trans, b_p)
        
                    if minDist <= minMatchDist:
            
                        " add to the list of match pairs less than 1.0 distance apart "
                        " keep A points and covariances untransformed "
                        Ca = halfHull2B[a_i][2]
                        Cb = halfHull1A[i][2]
        
                        " we store the untransformed point, but the transformed covariance of the A point "
                        match_pairs.append([halfHull2B[a_i],b_p,Ca,Cb])
            

        """        
        for i in range(len(hull2_trans)):
            a_p = poly2_trans[i]

            #if isValid(a_p, radiusB, centerB, poly1):
            if isValidPast(a_p, pastCircles, poly1):

                " for every transformed point of A, find it's closest neighbor in B "
                #b_p, minDist = findClosestPointInB(hull1, a_p, [0.0,0.0,0.0])
                b_p, minDist = findClosestPointInHull1(hull1, a_p, [0.0,0.0,0.0])
    
                if minDist <= minMatchDist:
        
                    " add to the list of match pairs less than 1.0 distance apart "
                    " keep A points and covariances untransformed "
                    Ca = hull2[i][2]
                    Cb = b_p[2]
    
                    " we store the untransformed point, but the transformed covariance of the A point "
                    match_pairs.append([hull2[i],b_p,Ca,Cb])

        for i in range(len(hull1)):
            b_p = poly1[i]
    
            if isValid(b_p, radius2, center2, poly2_trans):
        
                #print "selected", b_p, "in circle", radiusA, centerA, "with distance"
                " for every point of B, find it's closest neighbor in transformed A "
                #a_p, a_i, minDist = findClosestPointInA(hull2_trans, b_p)
                a_p, a_i, minDist = findClosestPointInHull2(hull2_trans, b_p)
    
                if minDist <= minMatchDist:
        
                    " add to the list of match pairs less than 1.0 distance apart "
                    " keep A points and covariances untransformed "
                    Ca = hull2[a_i][2]
                    Cb = hull1[i][2]
    
                    " we store the untransformed point, but the transformed covariance of the A point "
                    match_pairs.append([hull2[a_i],b_p,Ca,Cb])
        """
        
        allCircles = [[radius2,center2]]

        # optimize the match error for the current list of match pairs
        #newOffset = scipy.optimize.fmin(cost_func, offset, [match_pairs, hull2, poly1, allCircles], disp = 0)
        newOffset = scipy.optimize.fmin(shapeCostC, offset, [match_pairs], disp = 0)
    
        
        # get the current cost
        #newCost = cost_func(newOffset, match_pairs)
        newCost = shapeCostC(newOffset, match_pairs)
    
        # check for convergence condition, different between last and current cost is below threshold
        if abs(lastCost - newCost) < costThresh or (numIterations - startIteration) > 10:
            offset = newOffset
            lastCost = newCost
            break
    
        numIterations += 1

        offset = newOffset

        " reduce the minMatch distance for each step down to a floor value "
        minMatchDist /= 2
        if minMatchDist < 0.25:
            minMatchDist = 0.25
                            
        # save the current offset and cost
        offset = newOffset
        lastCost = newCost

    headPoint1 = termPoints[0]
    tailPoint1 = termPoints[1]
    headPoint2 = termPoints[2]
    tailPoint2 = termPoints[3]
        
    headPoint2_trans = dispOffset(headPoint2, offset)
    tailPoint2_trans = dispOffset(tailPoint2, offset)

    poly2_trans = []
    for b in hull2:
        p = [b[0],b[1]]
        p = dispOffset(p,offset)
        poly2_trans.append([p[0],p[1]])

    uniform2_trans = []
    for b in uniform2:
        p = dispOffset(b, offset)
        uniform2_trans.append(p)

    medial2_trans = []
    for b in medial2:
        p = dispOffset(b, offset)
        medial2_trans.append(p)

    halfHull2A_trans = []
    for b in halfHull2A:
        p = dispOffset(b, offset)
        halfHull2A_trans.append(p)

    halfHull2B_trans = []
    for b in halfHull2B:
        p = dispOffset(b, offset)
        halfHull2B_trans.append(p)

    poly2A_trans = []        
    for p in halfHull2A_trans:
        poly2A_trans.append([p[0],p[1]])    
    poly2B_trans = []        
    for p in halfHull2B_trans:
        poly2B_trans.append([p[0],p[1]])    
        
    stablePoints2_trans = []
    for b in stablePoints2:
        p = dispOffset(b, offset)
        stablePoints2_trans.append(p)
    
    inForeCount1 = 0
    inBackCount1 = 0
    inForeCount2 = 0
    inBackCount2 = 0

    #fore1 = stablePoints1[0:20]
    #back1 = stablePoints1[20:]
    #fore2 = stablePoints2_trans[0:20]
    #back2 = stablePoints2_trans[20:]

    print "len(medial1) =", len(medial1)

    fore1 = medial1[15:35]
    back1 = medial1[-36:-16]
    fore2 = medial2_trans[15:35]
    back2 = medial2_trans[-36:-16]



    for p in fore1:
        if functions.point_inside_polygon(p[0],p[1],poly2_trans):
            inForeCount1 += 1

    for p in back1:
        if functions.point_inside_polygon(p[0],p[1],poly2_trans):
            inBackCount1 += 1

    for p in fore2:
        if functions.point_inside_polygon(p[0],p[1],poly1):
            inForeCount2 += 1

    for p in back2:
        if functions.point_inside_polygon(p[0],p[1],poly1):
            inBackCount2 += 1
            

    res1 = functions.point_inside_polygon(headPoint2_trans[0],headPoint2_trans[1],poly1)
    res2 = functions.point_inside_polygon(tailPoint2_trans[0],tailPoint2_trans[1],poly1)
    res3 = functions.point_inside_polygon(headPoint1[0],headPoint1[1],poly2_trans)
    res4 = functions.point_inside_polygon(tailPoint1[0],tailPoint1[1],poly2_trans)

    " we want to capture the situation where one half of pose 1 is contained in pose 2 and vice versa simultaneously "
    if ((inForeCount1 >= 18) or (inBackCount1 >= 18)) and ((inForeCount2 >= 18) or (inBackCount2 >= 18)):
        pass
    else:
        lastCost = 1e100
        status = REJECTED_CONTAIN
    
    #if (res1 or res2) and (res3 or res4):
    #    pass
    #else:
    #    lastCost = 1e100

    if True:
        pylab.clf()
        pylab.axes()
        match_global = []
        
        for pair in match_pairs:
            p1 = pair[0]
            p2 = pair[1]
            
            p1_o = dispOffset(p1, offset)

            p1_g = poseOrigin.convertLocalToGlobal(p1_o)
            p2_g = poseOrigin.convertLocalToGlobal(p2)
            match_global.append([p1_g,p2_g])
        
        draw_matches(match_global, [0.0,0.0,0.0])
        
        """
        xP = []
        yP = []
        for b in poly1:
            p1 = poseOrigin.convertLocalToGlobal(b)

            xP.append(p1[0])    
            yP.append(p1[1])
        
        p1 = poseOrigin.convertLocalToGlobal(poly1[0])
        xP.append(p1[0])    
        yP.append(p1[1])
        
        pylab.plot(xP,yP,linewidth=1, color=(0.0,0.0,1.0))
        """


        """
        xP = []
        yP = []
        for p in uniform1:
            p1 = poseOrigin.convertLocalToGlobal(p)
            xP.append(p1[0])    
            yP.append(p1[1])            
        pylab.plot(xP,yP,linewidth=1, color=(0.0,0.0,1.0))
        """


        xP = []
        yP = []
        for p in medial1:
            p1 = poseOrigin.convertLocalToGlobal(p)
            xP.append(p1[0])    
            yP.append(p1[1])            
        pylab.plot(xP,yP,linewidth=1, color=(0.0,0.0,0.5))

        
        xP = []
        yP = []
        for b in bbox1:
            p1 = poseOrigin.convertLocalToGlobal(b)

            xP.append(p1[0])    
            yP.append(p1[1])
        
        p1 = poseOrigin.convertLocalToGlobal(bbox1[0])
        xP.append(p1[0])    
        yP.append(p1[1])
        pylab.plot(xP,yP,linewidth=1, color=(0.0,0.0,1.0))
        
        """
        xP = []
        yP = []
        for p in fore1:
            p1 = poseOrigin.convertLocalToGlobal(p)
            xP.append(p1[0])    
            yP.append(p1[1])            
        pylab.plot(xP,yP,linewidth=2, color=(0.0,0.0,1.0))
        xP = []
        yP = []
        for p in back1:
            p1 = poseOrigin.convertLocalToGlobal(p)
            xP.append(p1[0])    
            yP.append(p1[1])            
        pylab.plot(xP,yP,linewidth=2, color=(0.0,0.0,1.0))
        """

        p1 = poseOrigin.convertLocalToGlobal(headPoint1)
        p2 = poseOrigin.convertLocalToGlobal(tailPoint1)
        xP = [p1[0],p2[0]]
        yP = [p1[1],p2[1]]

        #pylab.scatter(xP,yP,linewidth=3, color=(0.0,0.0,1.0))        
        
        try:    
            #bPoints1, bPoints2 = computeSeparation(bbox1, uniform1, hull1)
            #bPoints1, bPoints2 = computeSeparation(bbox1, medial1, hull1)


            bPoints1 = halfHull1A
            bPoints2 = halfHull1B

            xP = []
            yP = []
            for p in bPoints1:
                p1 = poseOrigin.convertLocalToGlobal(p)
                xP.append(p1[0])    
                yP.append(p1[1])    

            pylab.scatter(xP,yP,linewidth=1, color=(0.0,0.0,1.0))

            xP = []
            yP = []
            for p in bPoints2:
                p1 = poseOrigin.convertLocalToGlobal(p)
                xP.append(p1[0])    
                yP.append(p1[1])    

            pylab.scatter(xP,yP,linewidth=1, color=(0.0,0.5,1.0))
        except Exception:
            traceback.print_exc()
            pass


        #for b in hull2:
        #    p = [b[0],b[1]]
        #    p = dispOffset(p,offset)
        #    poly2_trans.append([p[0],p[1]])

        """
        xP = []
        yP = []
        for p in poly2_trans:
            
            p1 = poseOrigin.convertLocalToGlobal(p)
            xP.append(p1[0])    
            yP.append(p1[1])
            
        p = poly2_trans[0]        
        p1 = poseOrigin.convertLocalToGlobal(p)
        xP.append(p1[0])    
        yP.append(p1[1])

        pylab.plot(xP,yP,linewidth=1, color=(1.0,0.0,0.0))
        """

        """
        xP = []
        yP = []
        for p in uniform2_trans:
            p1 = poseOrigin.convertLocalToGlobal(p)
            xP.append(p1[0])    
            yP.append(p1[1])    
        pylab.plot(xP,yP,linewidth=1, color=(1.0,0.0,0.0))
        """

        xP = []
        yP = []
        for p in medial2_trans:
            p1 = poseOrigin.convertLocalToGlobal(p)
            xP.append(p1[0])    
            yP.append(p1[1])    
        pylab.plot(xP,yP,linewidth=1, color=(0.5,0.0,0.0))

        xP = []
        yP = []
        for b in bbox2_trans:
            p1 = poseOrigin.convertLocalToGlobal(b)

            xP.append(p1[0])    
            yP.append(p1[1])
        
        p1 = poseOrigin.convertLocalToGlobal(bbox2_trans[0])
        xP.append(p1[0])    
        yP.append(p1[1])
        pylab.plot(xP,yP,linewidth=1, color=(1.0,0.0,0.0))

        """
        xP = []
        yP = []
        for p in fore2:
            p1 = poseOrigin.convertLocalToGlobal(p)
            xP.append(p1[0])    
            yP.append(p1[1])    
        pylab.plot(xP,yP,linewidth=2, color=(1.0,0.0,0.0))
        xP = []
        yP = []
        for p in back2:
            p1 = poseOrigin.convertLocalToGlobal(p)
            xP.append(p1[0])    
            yP.append(p1[1])    
        pylab.plot(xP,yP,linewidth=2, color=(1.0,0.0,0.0))
        """

        p1 = poseOrigin.convertLocalToGlobal(headPoint2_trans)
        p2 = poseOrigin.convertLocalToGlobal(tailPoint2_trans)
        
        xP = [p1[0],p2[0]]
        yP = [p1[1],p2[1]]
        #pylab.scatter(xP,yP,linewidth=3, color=(1.0,0.0,0.0))

        try:    
            #bPoints1, bPoints2 = computeSeparation(bbox2_trans, uniform2_trans, hull2_trans)
            #bPoints1, bPoints2 = computeSeparation(bbox2_trans, medial2_trans, hull2_trans)
            bPoints1 = halfHull2A_trans
            bPoints2 = halfHull2B_trans 

            xP = []
            yP = []
            for p in bPoints1:
                p1 = poseOrigin.convertLocalToGlobal(p)
                xP.append(p1[0])    
                yP.append(p1[1])    

            pylab.scatter(xP,yP,linewidth=1, color=(1.0,0.5,0.0))

            xP = []
            yP = []
            for p in bPoints2:
                p1 = poseOrigin.convertLocalToGlobal(p)
                xP.append(p1[0])    
                yP.append(p1[1])    

            pylab.scatter(xP,yP,linewidth=1, color=(1.0,0.0,0.0))
        except Exception:
            traceback.print_exc()
            pass


        pylab.title("%u %u, cost = %f, fore1,back1,fore2,back2 = %d,%d,%d,%d" % (n1, n2, lastCost, inForeCount1, inBackCount1, inForeCount2, inBackCount2))
        #pylab.xlim(-4,9)
        #pylab.ylim(-3,9)
        pylab.xlim(estPose1[0]-4, estPose1[0]+4)                    
        pylab.ylim(estPose1[1]-3, estPose1[1]+3)



        if lastCost > 1.6:
            pylab.savefig("badICP_plot_%04u_%04u.png" % (n1,n2))
        else:
            pylab.savefig("goodICP_plot_%04u_%04u.png" % (n1,n2))
        
        #if lastCost > 1.6:
        #    pylab.savefig("badICP_plot_%04u.png" % numIterations)
        #else:
        #    pylab.savefig("goodICP_plot_%04u.png" % numIterations)
            
        pylab.clf()

    offset[2] =  functions.normalizeAngle(offset[2])

    #res1 = functions.point_inside_polygon(headPoint2_trans[0],headPoint2_trans[1],poly1)
    #res2 = functions.point_inside_polygon(tailPoint2_trans[0],tailPoint2_trans[1],poly1)
    #res3 = functions.point_inside_polygon(headPoint1[0],headPoint1[1],poly2_trans)
    #res4 = functions.point_inside_polygon(tailPoint1[0],tailPoint1[1],poly2_trans)

    #if (res1 or res2) and (res3 or res4):

    #if breakOut:
    #    raise

    return offset, lastCost, status
    
        
    
    
def overlapICP(estPose1, gndOffset, initGuess, hull1, hull2, medialPoints1, medialPoints2, rootPose1, rootPose2, inPlace = False, plotIter = False, n1 = 0, n2 = 0, uRange = 1.5):

    global numIterations

    def computeOffset(point1, point2, ang1, ang2):
    
        " corner points and orientations "
        corner1Pose = [point1[0], point1[1], ang1]
        corner2Pose = [point2[0], point2[1], ang2]
        
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

    u1 = initGuess[0]
    u2 = initGuess[1]
    
    #uHigh = u2 + 1.0
    #uLow = u2 - 1.0
    #uHigh = u2 + 0.08
    #uLow = u2 - 0.08
    #uHigh = u2 + uRange
    #uLow = u2 - uRange


    
    currU = u2
    currAng = initGuess[2]
    medialSpline1 = SplineFit(medialPoints1, smooth=0.1)
    medialSpline2 = SplineFit(medialPoints2, smooth=0.1)

    uSet = [i*0.01 for i in range(100)]
    poses_1 = medialSpline1.getUVecSet(uSet)
    poses_2 = medialSpline2.getUVecSet(uSet)



    if u1 >= 1.0:
        pose1 = poses_1[-1]
    elif u1 < 0.0:
        pose1 = poses_1[0]
    else:
        pose1 = poses_1[int(u1*100)]

    if currU >= 1.0:
        pose2 = poses_2[-1]
    elif currU < 0.0:
        pose2 = poses_2[0]
    else:
        pose2 = poses_2[int(currU*100)]

    
    uHigh = medialSpline2.getUOfDist(u2, uRange, distIter = 0.001)
    uLow = medialSpline2.getUOfDist(u2, -uRange, distIter = 0.001)

    if inPlace:
        uHigh = u2 + 0.08
        uLow = u2 - 0.08


    #if u1+0.02 > 1.0:
    #    pose1 = medialSpline1.getUVecSet([0.98, 1.0])[0]
    #elif u1 < 0.0:
    #    pose1 = medialSpline1.getUVecSet([0.0, 0.02])[0]
    #else:
    #    pose1 = medialSpline1.getUVecSet([u1, u1+0.02])[0]

    #if currU+0.02 > 1.0:
    #    pose2 = medialSpline2.getUVecSet([0.98, 1.0])[0]
    #elif currU < 0.0:
    #    pose2 = medialSpline2.getUVecSet([0.0, 0.02])[0]
    #else:
    #    pose2 = medialSpline2.getUVecSet([currU, currU + 0.02])[0]

    point1 = [pose1[0],pose1[1]]
    point2 = [pose2[0],pose2[1]]
    ang1 = pose1[2]
    ang2 = pose2[2]

    offset = computeOffset(point1, point2, ang1, ang2 + currAng)


    " transform the past poses "
    a_data_raw = hull2
    a_data = []
    for p in a_data_raw:
        result = dispOffset(p, offset)        
        a_data.append(result)    

    polyB = []        
    for p in hull1:
        polyB.append([p[0],p[1]])    
    
    costThresh = 0.004
    minMatchDist = 2.0
    lastCost = 1e100
    
    startIteration = numIterations

    " set the initial guess "
    poseOrigin = Pose(estPose1)
    
    " sample a series of points from the medial curves "
    " augment points with point-to-line covariances "
    " treat the points with the point-to-line constraint "
    #points1 = addGPACVectorCovariance(medialSpline1.getUniformSamples(),high_var=0.1, low_var = 0.05)
    #points2 = addGPACVectorCovariance(medialSpline2.getUniformSamples(),high_var=0.1, low_var = 0.05)
    points1 = addGPACVectorCovariance(medialSpline1.getUniformSamples(),high_var=0.05, low_var = 0.001)
    points2 = addGPACVectorCovariance(medialSpline2.getUniformSamples(),high_var=0.05, low_var = 0.001)
    #points1 = addGPACVectorCovariance(medialSpline1.getUniformSamples(),high_var=0.1, low_var = 0.0005)
    #points2 = addGPACVectorCovariance(medialSpline2.getUniformSamples(),high_var=0.1, low_var = 0.0005)
    #points1 = addGPACVectorCovarianceWithAngle(medialSpline1.getUniformSamples(),high_var=0.1, low_var = 0.001)
    #points2 = addGPACVectorCovarianceWithAngle(medialSpline2.getUniformSamples(),high_var=0.1, low_var = 0.001)

    " transform pose 2 by initial offset guess "    
    " transform the new pose "
    points2_offset = []
    for p in points2:
        result = dispPoint(p, offset)        
        #result = dispPointWithAngle(p, offset)
        points2_offset.append(result)


    
    poly1 = []        
    for p in points1:
        poly1.append([p[0],p[1]])    

    #if inPlace:
    if False:
        offset[2] =  functions.normalizeAngle(offset[2])    

        #if plotIter:
        if True:
            
            " set the origin of pose 1 "
            poseOrigin = Pose(estPose1)
    
            #pylab.clf()
            #pylab.axes()
            
            #pylab.subplot(121)
            fig1 = pylab.figure(1)
            fig1.clf()

            axes2 = fig1.add_subplot(111)

            xP = []
            yP = []
            for b in poly1:
                p1 = poseOrigin.convertLocalToGlobal(b)
    
                xP.append(p1[0])    
                yP.append(p1[1])
    
            axes2.plot(xP,yP,linewidth=1, color=(0.0,0.0,1.0))
    
            xP = []
            yP = []
            for b in points2:
                p = [b[0],b[1]]
                p = dispOffset(p,offset)
                
                p1 = poseOrigin.convertLocalToGlobal(p)
                xP.append(p1[0])    
                yP.append(p1[1])

            axes2.plot(xP,yP,linewidth=1, color=(1.0,0.0,0.0))

            xP = []
            yP = []
            for b in points2:
                p = [b[0],b[1]]
                p = dispOffset(p,gndOffset)
                
                p1 = poseOrigin.convertLocalToGlobal(p)
                xP.append(p1[0])    
                yP.append(p1[1])

            axes2.plot(xP,yP,linewidth=1, color=(1.0,0.5,0.5))

            point2_trans = dispOffset(point2, offset)
            p1 = poseOrigin.convertLocalToGlobal(point1)
            p2 = poseOrigin.convertLocalToGlobal(point2_trans)

            xP = [p1[0],p2[0]]
            yP = [p1[1],p2[1]]
            axes2.scatter(xP,yP,color='b')


            gpac2_trans = dispOffset(rootPose1, offset)
            gpac1 = poseOrigin.convertLocalToGlobal(rootPose2)
            gpac2 = poseOrigin.convertLocalToGlobal(gpac2_trans)

            xP = [gpac1[0],gpac2[0]]
            yP = [gpac1[1],gpac2[1]]
            axes2.scatter(xP,yP,color='r')


            xP = []
            yP = []
            for b in polyB:
                p1 = poseOrigin.convertLocalToGlobal(b)

                xP.append(p1[0])    
                yP.append(p1[1])
            
            axes2.plot(xP,yP,linewidth=1, color=(0.0,0.0,1.0))


            xP = []
            yP = []
            for b in a_data_raw:
                p = [b[0],b[1]]
                p = dispOffset(p,offset)
                
                p1 = poseOrigin.convertLocalToGlobal(p)
                xP.append(p1[0])    
                yP.append(p1[1])

            axes2.plot(xP,yP,linewidth=1, color=(1.0,0.0,0.0))

            xP = []
            yP = []
            for b in a_data_raw:
                p = [b[0],b[1]]
                p = dispOffset(p,gndOffset)
                
                p1 = poseOrigin.convertLocalToGlobal(p)
                xP.append(p1[0])    
                yP.append(p1[1])

            axes2.plot(xP,yP,linewidth=1, color=(1.0,0.5,0.5))

            plotEnv(axes2)        
            
            
            axes2.set_title("%s %s, ang = %1.3f" % (repr(n1), repr(n2), currAng))
            axes2.set_xlim(estPose1[0]-4, estPose1[0]+4)                    
            axes2.set_ylim(estPose1[1]-3, estPose1[1]+3)
            

            """
            axes1 = fig1.add_subplot(121)
            axes1.grid(True)
            
            
            

            
            #pylab.subplot(122)
            #pylab.xlim(-4,4)
            #pylab.ylim(-4,4)
            
            tPoints1 = medialSpline1.getTransformCurve()
            #tPoints2 = medialSpline2.getTransformCurve()
            medialPoints2_gnd = []
            for b in medialPoints2:
                p = [b[0],b[1]]
                p = dispOffset(p,gndOffset)
                
                #p1 = poseOrigin.convertLocalToGlobal(p)
                medialPoints2_gnd.append(p)

            medialSpline2_gnd = SplineFit(medialPoints2_gnd, smooth=0.1)
            tPoints2_gnd = medialSpline2_gnd.getTransformCurve()

            medialPoints2_trans = []
            for b in medialPoints2:
                p = [b[0],b[1]]
                p = dispOffset(p,offset)
                
                #p1 = poseOrigin.convertLocalToGlobal(p)
                medialPoints2_trans.append(p)

            medialSpline2_trans = SplineFit(medialPoints2_trans, smooth=0.1)
            tPoints2_trans = medialSpline2_trans.getTransformCurve()



            dist1 = medialSpline1.dist(0.0, u1, iter = 0.005)    

            xP = []
            yP = []
            for p in tPoints1:
                xP.append(p[0]+2.0)
                yP.append(p[1])
            axes1.plot(xP,yP, color='b')

            dist2 = medialSpline1.dist(0.0, currU, iter = 0.005)                

            xP = []
            yP = []
            for p in tPoints2_trans:
                xP.append(p[0]+dist1-dist2+2.0)
                yP.append(p[1])
            axes1.plot(xP,yP, color='r')

            #xP = []
            #yP = []
            #for p in tPoints2_gnd:
            #    xP.append(p[0])
            #    yP.append(p[1])
            #axes1.plot(xP,yP, color=(1.0,0.5,0.5))
            
            axes1.set_xlabel("distance")
            axes1.set_ylabel("angle (radians)")
            axes1.set_xlim(0,10)
            axes1.set_ylim(-4,4)        
            
            #fig1.set_size_inches(12,6)

            """
            #fig1 = pylab.gcf()
            #fig1.set_size_inches(12,6)
            fig1.savefig("ICP_plot_%04u_0.png" % numIterations)
            #fig1.set_size_inches(8,6)
            pylab.clf()            
            
            pylab.figure(1)
    
            
            numIterations += 1
    
        
        #return offset, (0,0,0)


    while True:
        
        " find the matching pairs "
        match_pairs = []

        " transform the target Hull with the latest offset "
        points2_offset = []
        for p in points2:
            result = dispPoint(p, offset)        
            #result = dispPointWithAngle(p, offset)
            points2_offset.append(result)
        
        " transformed points without associated covariance "
        poly2 = []
        for p in points2_offset:
            poly2.append([p[0],p[1]])    
        
        " get the circles and radii "
        radius2, center2 = computeEnclosingCircle(points2_offset)
        radius1, center1 = computeEnclosingCircle(points1)
    
        for i in range(len(points2_offset)):
            p_2 = poly2[i]

            if isInCircle(p_2, radius1, center1):

                " for every transformed point of A, find it's closest neighbor in B "
                p_1, minDist = findClosestPointInB(points1, p_2, [0.0,0.0,0.0])
    
                if minDist <= minMatchDist:
        
                    " add to the list of match pairs less than 1.0 distance apart "
                    " keep A points and covariances untransformed "
                    C2 = points2[i][2]
                    C1 = p_1[2]
    
                    " we store the untransformed point, but the transformed covariance of the A point "
                    match_pairs.append([points2[i],p_1,C2,C1])

        for i in range(len(points1)):
            p_1 = poly1[i]
    
            if isInCircle(p_1, radius2, center2):
        
                #print "selected", b_p, "in circle", radiusA, centerA, "with distance"
                " for every point of B, find it's closest neighbor in transformed A "
                p_2, i_2, minDist = findClosestPointInA(points2_offset, p_1)
    
                if minDist <= minMatchDist:
        
                    " add to the list of match pairs less than 1.0 distance apart "
                    " keep A points and covariances untransformed "
                    C2 = points2[i_2][2]
                    
                    C1 = points1[i][2]
    
                    " we store the untransformed point, but the transformed covariance of the A point "
                    #match_pairs.append([points2[i_2],p_1,C2,C1])
                    match_pairs.append([points2[i_2],points1[i],C2,C1])

        #if plotIter and startIteration == numIterations:
        if plotIter:
            
            " set the origin of pose 1 "
            poseOrigin = Pose(estPose1)
    
            pylab.clf()
            pylab.axes()
            match_global = []
            
            for pair in match_pairs:
                p1 = pair[0]
                p2 = pair[1]
                
                p1_o = dispOffset(p1, offset)
                
                p1_g = poseOrigin.convertLocalToGlobal(p1_o)
                p2_g = poseOrigin.convertLocalToGlobal(p2)
                match_global.append([p1_g,p2_g])

            draw_matches(match_global, [0.0,0.0,0.0])

            
            xP = []
            yP = []
            for b in poly1:
                p1 = poseOrigin.convertLocalToGlobal(b)
    
                xP.append(p1[0])    
                yP.append(p1[1])
    
            pylab.plot(xP,yP,linewidth=1, color=(0.0,0.0,1.0))

            #center1Global = poseOrigin.convertLocalToGlobal(center1)
            #cir = pylab.Circle((center1Global[0],center1Global[1]), radius=radius1,  ec='b', fc='none')
            #pylab.gca().add_patch(cir)

    
            xP = []
            yP = []
            for b in points2:
                p = [b[0],b[1]]
                p = dispOffset(p,offset)
                
                p1 = poseOrigin.convertLocalToGlobal(p)
                xP.append(p1[0])    
                yP.append(p1[1])

            pylab.plot(xP,yP,linewidth=1, color=(1.0,0.0,0.0))

            xP = []
            yP = []
            for b in points2:
                p = [b[0],b[1]]
                p = dispOffset(p,gndOffset)
                
                p1 = poseOrigin.convertLocalToGlobal(p)
                xP.append(p1[0])    
                yP.append(p1[1])

            pylab.plot(xP,yP,linewidth=1, color=(1.0,0.5,0.5))

            #center2Global = poseOrigin.convertLocalToGlobal(center2)
            #cir = pylab.Circle((center2Global[0],center2Global[1]), radius=radius2,  ec='r', fc='none')
            #pylab.gca().add_patch(cir)


            point2_trans = dispOffset(point2, offset)
            p1 = poseOrigin.convertLocalToGlobal(point1)
            p2 = poseOrigin.convertLocalToGlobal(point2_trans)

            xP = [p1[0],p2[0]]
            yP = [p1[1],p2[1]]
            pylab.scatter(xP,yP,color='b')


            gpac2_trans = dispOffset(rootPose1, offset)
            gpac1 = poseOrigin.convertLocalToGlobal(rootPose2)
            gpac2 = poseOrigin.convertLocalToGlobal(gpac2_trans)
            #gpac2_trans = dispOffset([0.0,0.0], offset)
            #gpac1 = poseOrigin.convertLocalToGlobal([0.0,0.0])
            #gpac2 = poseOrigin.convertLocalToGlobal(gpac2_trans)

            xP = [gpac1[0],gpac2[0]]
            yP = [gpac1[1],gpac2[1]]
            pylab.scatter(xP,yP,color='r')


            xP = []
            yP = []
            for b in polyB:
                p1 = poseOrigin.convertLocalToGlobal(b)

                xP.append(p1[0])    
                yP.append(p1[1])
            
            pylab.plot(xP,yP,linewidth=1, color=(0.0,0.0,1.0))


            xP = []
            yP = []
            for b in a_data_raw:
                p = [b[0],b[1]]
                p = dispOffset(p,offset)
                
                p1 = poseOrigin.convertLocalToGlobal(p)
                xP.append(p1[0])    
                yP.append(p1[1])

            pylab.plot(xP,yP,linewidth=1, color=(1.0,0.0,0.0))

            xP = []
            yP = []
            for b in a_data_raw:
                p = [b[0],b[1]]
                p = dispOffset(p,gndOffset)
                
                p1 = poseOrigin.convertLocalToGlobal(p)
                xP.append(p1[0])    
                yP.append(p1[1])

            pylab.plot(xP,yP,linewidth=1, color=(1.0,0.5,0.5))

            """
            xP = []
            yP = []
            for b in medialPoints1:
                p1 = poseOrigin.convertLocalToGlobal(b)

                xP.append(p1[0])    
                yP.append(p1[1])

            pylab.scatter(xP,yP,color=(0.0,0.0,1.0))

            xP = []
            yP = []
            for b in medialPoints2:
                p = [b[0],b[1]]
                p = dispOffset(p,offset)
                
                p1 = poseOrigin.convertLocalToGlobal(p)
                xP.append(p1[0])    
                yP.append(p1[1])
                
            pylab.scatter(xP,yP,color=(1.0,0.0,0.0))
            """
            plotEnv()        
            
            #pylab.title("%u %u, u = %f" % (n1, n2, currU))
            lowCount, midCount, highCount = overlapHistogram([currU,currAng], match_pairs, medialSpline1, medialSpline2, u1)            
            pylab.title("%s %s, ang = %1.3f, hist = %d %d %d" % (repr(n1), repr(n2), currAng, lowCount, midCount, highCount))
            pylab.xlim(estPose1[0]-4, estPose1[0]+4)                    
            pylab.ylim(estPose1[1]-3, estPose1[1]+3)
            #pylab.xlim(-3,3)
            #pylab.ylim(-3,3)
            #pylab.axis('equal')
            pylab.savefig("ICP_plot_%04u_0.png" % numIterations)
            pylab.clf()            
        

        #if inPlace:
        #    newParam = scipy.optimize.fmin(medialOverlapInPlaceCostFunc, [currAng], [match_pairs, medialSpline1, medialSpline2, currU, u1], disp = 0)
        #    newU = currU
        #    newAng = newParam[0]
        #else:
        #    pass
        
        
        newParam = scipy.optimize.fmin(medialOverlapCostFunc, [currU,currAng], [match_pairs, poses_1, poses_2, uHigh, uLow, u1], disp = 0)
        #newParam = scipy.optimize.fmin(medialOverlapCostFunc, [currU,currAng], [match_pairs, medialSpline1, medialSpline2, uHigh, uLow, u1], disp = 0)
        newU = newParam[0]
        newAng = newParam[1]
        
        # get the current cost
        #if inPlace:
        #    newCost = medialOverlapInPlaceCostFunc([newAng], match_pairs, medialSpline1, medialSpline2, newU, u1)
        #else:
        newCost = medialOverlapCostFunc([newU, newAng], match_pairs, poses_1, poses_2, uHigh, uLow, u1)
        
        " set the current parameters "
        currAng = functions.normalizeAngle(newAng)
        currU = newU
        
        #print "iteration", numIterations
        if inPlace:
            print "currU =", currU, "currAng =", currAng, "newCost =", newCost

        " compute offset from newU and newAng"

        if u1+0.02 > 1.0:
            pose1 = medialSpline1.getUVecSet([0.98, 1.0])[0]
        elif u1 < 0.0:
            pose1 = medialSpline1.getUVecSet([0.0, 0.02])[0]
        else:
            pose1 = medialSpline1.getUVecSet([u1, u1+0.02])[0]
    
        if currU+0.02 > 1.0:
            pose2 = medialSpline2.getUVecSet([0.98, 1.0])[0]
        elif currU < 0.0:
            pose2 = medialSpline2.getUVecSet([0.0, 0.02])[0]
        else:
            pose2 = medialSpline2.getUVecSet([currU, currU + 0.02])[0]

        #pose1 = medialSpline1.getUVecSet([u1, u1+0.02])[0]
        #pose2 = medialSpline2.getUVecSet([currU, currU+0.02])[0]
    
        point1 = [pose1[0],pose1[1]]
        point2 = [pose2[0],pose2[1]]
        ang1 = pose1[2]
        ang2 = pose2[2]
    
        newOffset = computeOffset(point1, point2, ang1, ang2 + currAng)

        # save the current offset and cost
        offset = newOffset

    
        # check for convergence condition, different between last and current cost is below threshold
        if abs(lastCost - newCost) < costThresh or (numIterations - startIteration) > 10:
            #print "breaking:", abs(lastCost - newCost) < costThresh, (numIterations - startIteration) > 10
            #print "lastCost =", lastCost, "newCost =", newCost, "costThresh =", costThresh
            break



        " update after check for termination condition "
        lastCost = newCost

        # optionally draw the position of the points in current transform
        if plotIter:
            
            " set the origin of pose 1 "
            poseOrigin = Pose(estPose1)
    
            pylab.clf()
            pylab.axes()
            match_global = []
            
            for pair in match_pairs:
                p1 = pair[0]
                p2 = pair[1]
                
                p1_o = dispOffset(p1, offset)
                
                p1_g = poseOrigin.convertLocalToGlobal(p1_o)
                p2_g = poseOrigin.convertLocalToGlobal(p2)
                match_global.append([p1_g,p2_g])
            
            draw_matches(match_global, [0.0,0.0,0.0])
            
            xP = []
            yP = []
            for b in poly1:
                p1 = poseOrigin.convertLocalToGlobal(b)

                xP.append(p1[0])    
                yP.append(p1[1])

            pylab.plot(xP,yP,linewidth=1, color=(0.0,0.0,1.0))

            #center1Global = poseOrigin.convertLocalToGlobal(center1)
            #cir = pylab.Circle((center1Global[0],center1Global[1]), radius=radius1,  ec='b', fc='none')
            #pylab.gca().add_patch(cir)

            xP = []
            yP = []
            for b in points2:
                p = [b[0],b[1]]
                p = dispOffset(p,offset)
                
                p1 = poseOrigin.convertLocalToGlobal(p)
                xP.append(p1[0])    
                yP.append(p1[1])
            
            pylab.plot(xP,yP,linewidth=1, color=(1.0,0.0,0.0))

            xP = []
            yP = []
            for b in points2:
                p = [b[0],b[1]]
                p = dispOffset(p,gndOffset)
                
                p1 = poseOrigin.convertLocalToGlobal(p)
                xP.append(p1[0])    
                yP.append(p1[1])

            pylab.plot(xP,yP,linewidth=1, color=(1.0,0.5,0.5))

            #center2Global = poseOrigin.convertLocalToGlobal(center2)
            #cir = pylab.Circle((center2Global[0],center2Global[1]), radius=radius2,  ec='r', fc='none')
            #pylab.gca().add_patch(cir)

            point2_trans = dispOffset(point2, offset)
            p1 = poseOrigin.convertLocalToGlobal(point1)
            p2 = poseOrigin.convertLocalToGlobal(point2_trans)

            #xP = [point1[0],point2_trans[0]]
            #yP = [point1[1],point2_trans[1]]
            xP = [p1[0],p2[0]]
            yP = [p1[1],p2[1]]

            pylab.scatter(xP,yP,color='b')

            gpac2_trans = dispOffset(rootPose1, offset)
            gpac1 = poseOrigin.convertLocalToGlobal(rootPose2)
            gpac2 = poseOrigin.convertLocalToGlobal(gpac2_trans)
            #gpac2_trans = dispOffset([0.0,0.0], offset)
            #gpac1 = poseOrigin.convertLocalToGlobal([0.0,0.0])
            #gpac2 = poseOrigin.convertLocalToGlobal(gpac2_trans)

            xP = [gpac1[0],gpac2[0]]
            yP = [gpac1[1],gpac2[1]]
            pylab.scatter(xP,yP,color='r')

            xP = []
            yP = []
            for b in polyB:
                p1 = poseOrigin.convertLocalToGlobal(b)

                xP.append(p1[0])    
                yP.append(p1[1])
            
            pylab.plot(xP,yP,linewidth=1, color=(0.0,0.0,1.0))

            xP = []
            yP = []
            for b in a_data_raw:
                p = [b[0],b[1]]
                p = dispOffset(p,offset)
                
                p1 = poseOrigin.convertLocalToGlobal(p)
                xP.append(p1[0])    
                yP.append(p1[1])
            
            pylab.plot(xP,yP,linewidth=1, color=(1.0,0.0,0.0))

            xP = []
            yP = []
            for b in a_data_raw:
                p = [b[0],b[1]]
                p = dispOffset(p,gndOffset)
                
                p1 = poseOrigin.convertLocalToGlobal(p)
                xP.append(p1[0])    
                yP.append(p1[1])

            pylab.plot(xP,yP,linewidth=1, color=(1.0,0.5,0.5))

            
            """
            xP = []
            yP = []
            for b in medialPoints1:
                p1 = poseOrigin.convertLocalToGlobal(b)

                xP.append(p1[0])    
                yP.append(p1[1])

            pylab.scatter(xP,yP,color=(0.0,0.0,1.0))

            xP = []
            yP = []
            for b in medialPoints2:
                p = [b[0],b[1]]
                p = dispOffset(p,offset)
                
                p1 = poseOrigin.convertLocalToGlobal(p)
                xP.append(p1[0])    
                yP.append(p1[1])
                
            pylab.scatter(xP,yP,color=(1.0,0.0,0.0))
            """
            plotEnv()        
            
            lowCount, midCount, highCount = overlapHistogram([currU,currAng], match_pairs, medialSpline1, medialSpline2, u1)            
            #pylab.title("%u %u, u1 = %1.3f, u2 = %1.3f, hist = %d %d %d" % (n1, n2, u1, currU, lowCount, midCount, highCount))
            #pylab.title("%u %u, u1 = %1.3f, u2 = %1.3f, ang = %1.3f, hist = %d %d %d" % (n1, n2, u1, currU, currAng, lowCount, midCount, highCount))
            pylab.title("%s %s, ang = %1.3f, hist = %d %d %d" % (repr(n1), repr(n2), currAng, lowCount, midCount, highCount))
            #pylab.title("%u %u, %u, u1 = %1.3f, u2 = %1.3f" % (n1, n2, len(match_pairs), u1, currU))
            pylab.xlim(estPose1[0]-4, estPose1[0]+4)                    
            pylab.ylim(estPose1[1]-3, estPose1[1]+3)
            #pylab.axis('equal')
            pylab.savefig("ICP_plot_%04u_1.png" % numIterations)
            pylab.clf()            

        
        numIterations += 1

        " reduce the minMatch distance for each step down to a floor value "
        #minMatchDist /= 2
        #if minMatchDist < 0.25:
        #    minMatchDist = 0.25
    
    if False:

        " set the origin of pose 1 "
        poseOrigin = Pose(estPose1)

        fig1 = pylab.figure(1)
        fig1.clf()

        axes2 = fig1.add_subplot(111)

        xP = []
        yP = []
        for b in poly1:
            p1 = poseOrigin.convertLocalToGlobal(b)

            xP.append(p1[0])    
            yP.append(p1[1])

        axes2.plot(xP,yP,linewidth=1, color=(0.0,0.0,1.0))

        xP = []
        yP = []
        for b in points2:
            p = [b[0],b[1]]
            p = dispOffset(p,offset)
            
            p1 = poseOrigin.convertLocalToGlobal(p)
            xP.append(p1[0])    
            yP.append(p1[1])

        axes2.plot(xP,yP,linewidth=1, color=(1.0,0.0,0.0))

        xP = []
        yP = []
        for b in points2:
            p = [b[0],b[1]]
            p = dispOffset(p,gndOffset)
            
            p1 = poseOrigin.convertLocalToGlobal(p)
            xP.append(p1[0])    
            yP.append(p1[1])

        axes2.plot(xP,yP,linewidth=1, color=(1.0,0.5,0.5))

        point2_trans = dispOffset(point2, offset)
        p1 = poseOrigin.convertLocalToGlobal(point1)
        p2 = poseOrigin.convertLocalToGlobal(point2_trans)

        xP = [p1[0],p2[0]]
        yP = [p1[1],p2[1]]
        axes2.scatter(xP,yP,color='b')


        gpac2_trans = dispOffset(rootPose1, offset)
        gpac1 = poseOrigin.convertLocalToGlobal(rootPose2)
        gpac2 = poseOrigin.convertLocalToGlobal(gpac2_trans)

        xP = [gpac1[0],gpac2[0]]
        yP = [gpac1[1],gpac2[1]]
        axes2.scatter(xP,yP,color='r')


        xP = []
        yP = []
        for b in polyB:
            p1 = poseOrigin.convertLocalToGlobal(b)

            xP.append(p1[0])    
            yP.append(p1[1])
        
        axes2.plot(xP,yP,linewidth=1, color=(0.0,0.0,1.0))


        xP = []
        yP = []
        for b in a_data_raw:
            p = [b[0],b[1]]
            p = dispOffset(p,offset)
            
            p1 = poseOrigin.convertLocalToGlobal(p)
            xP.append(p1[0])    
            yP.append(p1[1])

        axes2.plot(xP,yP,linewidth=1, color=(1.0,0.0,0.0))

        xP = []
        yP = []
        for b in a_data_raw:
            p = [b[0],b[1]]
            p = dispOffset(p,gndOffset)
            
            p1 = poseOrigin.convertLocalToGlobal(p)
            xP.append(p1[0])    
            yP.append(p1[1])

        axes2.plot(xP,yP,linewidth=1, color=(1.0,0.5,0.5))

        plotEnv(axes2)        
        
        axes2.set_title("%s %s, ang = %1.3f" % (repr(n1), repr(n2), currAng))
        axes2.set_xlim(estPose1[0]-4, estPose1[0]+4)                    
        axes2.set_ylim(estPose1[1]-3, estPose1[1]+3)
        
        """
        axes1 = fig1.add_subplot(121)
        axes1.grid(True)

        tPoints1 = medialSpline1.getTransformCurve()
        #tPoints2 = medialSpline2.getTransformCurve()
        medialPoints2_gnd = []
        for b in medialPoints2:
            p = [b[0],b[1]]
            p = dispOffset(p,gndOffset)
            
            #p1 = poseOrigin.convertLocalToGlobal(p)
            medialPoints2_gnd.append(p)

        medialSpline2_gnd = SplineFit(medialPoints2_gnd, smooth=0.1)
        tPoints2_gnd = medialSpline2_gnd.getTransformCurve()

        medialPoints2_trans = []
        for b in medialPoints2:
            p = [b[0],b[1]]
            p = dispOffset(p,offset)
            
            #p1 = poseOrigin.convertLocalToGlobal(p)
            medialPoints2_trans.append(p)

        medialSpline2_trans = SplineFit(medialPoints2_trans, smooth=0.1)
        tPoints2_trans = medialSpline2_trans.getTransformCurve()

        dist1 = medialSpline1.dist(0.0, u1, iter = 0.005)    

        xP = []
        yP = []
        for p in tPoints1:
            xP.append(p[0]+2.0)
            yP.append(p[1])
        axes1.plot(xP,yP, color='b')

        dist2 = medialSpline1.dist(0.0, currU, iter = 0.005)                

        xP = []
        yP = []
        for p in tPoints2_trans:
            xP.append(p[0]+dist1-dist2+2.0)
            yP.append(p[1])
        axes1.plot(xP,yP, color='r')



        #xP = []
        #yP = []
        #for p in tPoints2_gnd:
        #    xP.append(p[0])
        #    yP.append(p[1])
        #axes1.plot(xP,yP, color=(1.0,0.5,0.5))
        
        axes1.set_xlabel("distance")
        axes1.set_ylabel("angle (radians)")
        axes1.set_xlim(0,10)
        axes1.set_ylim(-4,4)        
        """
        fig1.savefig("ICP_plot_%04u_0.png" % numIterations)
        pylab.clf()            
        pylab.figure(1)
        """
        " set the origin of pose 1 "
        poseOrigin = Pose(estPose1)

        pylab.clf()
        pylab.axes()
        match_global = []
        
        for pair in match_pairs:
            p1 = pair[0]
            p2 = pair[1]
            
            p1_o = dispOffset(p1, offset)
            
            p1_g = poseOrigin.convertLocalToGlobal(p1_o)
            p2_g = poseOrigin.convertLocalToGlobal(p2)
            match_global.append([p1_g,p2_g])
        
        draw_matches(match_global, [0.0,0.0,0.0])
        
        xP = []
        yP = []
        for b in poly1:
            p1 = poseOrigin.convertLocalToGlobal(b)

            xP.append(p1[0])    
            yP.append(p1[1])

        pylab.plot(xP,yP,linewidth=1, color=(0.0,0.0,1.0))

        xP = []
        yP = []
        for b in points2:
            p = [b[0],b[1]]
            p = dispOffset(p,offset)
            
            p1 = poseOrigin.convertLocalToGlobal(p)
            xP.append(p1[0])    
            yP.append(p1[1])
        
        pylab.plot(xP,yP,linewidth=1, color=(1.0,0.0,0.0))

        xP = []
        yP = []
        for b in points2:
            p = [b[0],b[1]]
            p = dispOffset(p,gndOffset)
            
            p1 = poseOrigin.convertLocalToGlobal(p)
            xP.append(p1[0])    
            yP.append(p1[1])

        pylab.plot(xP,yP,linewidth=1, color=(1.0,0.5,0.5))


        point2_trans = dispOffset(point2, offset)
        p1 = poseOrigin.convertLocalToGlobal(point1)
        p2 = poseOrigin.convertLocalToGlobal(point2_trans)


        xP = [p1[0],p2[0]]
        yP = [p1[1],p2[1]]

        pylab.scatter(xP,yP,color='b')

        gpac2_trans = dispOffset(rootPose1, offset)
        gpac1 = poseOrigin.convertLocalToGlobal(rootPose2)
        gpac2 = poseOrigin.convertLocalToGlobal(gpac2_trans)

        xP = [gpac1[0],gpac2[0]]
        yP = [gpac1[1],gpac2[1]]
        pylab.scatter(xP,yP,color='r')

        xP = []
        yP = []
        for b in polyB:
            p1 = poseOrigin.convertLocalToGlobal(b)

            xP.append(p1[0])    
            yP.append(p1[1])
        
        pylab.plot(xP,yP,linewidth=1, color=(0.0,0.0,1.0))

        xP = []
        yP = []
        for b in a_data_raw:
            p = [b[0],b[1]]
            p = dispOffset(p,offset)
            
            p1 = poseOrigin.convertLocalToGlobal(p)
            xP.append(p1[0])    
            yP.append(p1[1])
        
        pylab.plot(xP,yP,linewidth=1, color=(1.0,0.0,0.0))

        xP = []
        yP = []
        for b in a_data_raw:
            p = [b[0],b[1]]
            p = dispOffset(p,gndOffset)
            
            p1 = poseOrigin.convertLocalToGlobal(p)
            xP.append(p1[0])    
            yP.append(p1[1])

        pylab.plot(xP,yP,linewidth=1, color=(1.0,0.5,0.5))

        
        plotEnv()        
        
        lowCount, midCount, highCount = overlapHistogram([currU,currAng], match_pairs, medialSpline1, medialSpline2, u1)            
        pylab.title("%s %s, ang = %1.3f, hist = %d %d %d" % (repr(n1), repr(n2), currAng, lowCount, midCount, highCount))
        pylab.xlim(estPose1[0]-4, estPose1[0]+4)                    
        pylab.ylim(estPose1[1]-3, estPose1[1]+3)
        pylab.savefig("ICP_plot_%04u_1.png" % numIterations)
        pylab.clf()        
        """    
        numIterations += 1


        
    histogram = overlapHistogram([currU,currAng], match_pairs, medialSpline1, medialSpline2, u1)            


    offset[2] =  functions.normalizeAngle(offset[2])    
    return offset, histogram


def overlapICP_GPU2(estPose1, gndOffset, initGuess, hull1, hull2, medialPoints1, medialPoints2, rootPose1, rootPose2, inPlace = False, plotIter = False, n1 = 0, n2 = 0, uRange = 1.5):

    global numIterations

    def computeOffset(point1, point2, ang1, ang2):
    
        " corner points and orientations "
        corner1Pose = [point1[0], point1[1], ang1]
        corner2Pose = [point2[0], point2[1], ang2]
        
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

    u1 = initGuess[0]
    u2 = initGuess[1]
    
    currU = u2
    currAng = initGuess[2]
    medialSpline1 = SplineFit(medialPoints1, smooth=0.1)
    medialSpline2 = SplineFit(medialPoints2, smooth=0.1)


    uSet = [i*0.01 for i in range(100)]
    poses_1 = medialSpline1.getUVecSet(uSet)
    poses_2 = medialSpline2.getUVecSet(uSet)

    uHigh = medialSpline2.getUOfDist(u2, uRange, distIter = 0.001)
    uLow = medialSpline2.getUOfDist(u2, -uRange, distIter = 0.001)

    if inPlace:
        uHigh = u2 + 0.08
        uLow = u2 - 0.08

    if u1 >= 1.0:
        pose1 = poses_1[-1]
    elif u1 < 0.0:
        pose1 = poses_1[0]
    else:
        pose1 = poses_1[int(u1*100)]

    if currU >= 1.0:
        pose2 = poses_2[-1]
    elif currU < 0.0:
        pose2 = poses_2[0]
    else:
        pose2 = poses_2[int(currU*100)]

    point1 = [pose1[0],pose1[1]]
    point2 = [pose2[0],pose2[1]]
    ang1 = pose1[2]
    ang2 = pose2[2]

    offset = computeOffset(point1, point2, ang1, ang2 + currAng)


    " transform the past poses "
    a_data_raw = hull2
    a_data = []
    for p in a_data_raw:
        result = dispOffset(p, offset)        
        a_data.append(result)    

    polyB = []        
    for p in hull1:
        polyB.append([p[0],p[1]])    
    
    costThresh = 0.004
    minMatchDist = 2.0
    lastCost = 1e100
    
    startIteration = numIterations

    " set the initial guess "
    poseOrigin = Pose(estPose1)
    
    " sample a series of points from the medial curves "
    " augment points with point-to-line covariances "
    " treat the points with the point-to-line constraint "

    points1 = addGPACVectorCovariance(medialSpline1.getUniformSamples(),high_var=0.05, low_var = 0.001)
    points2 = addGPACVectorCovariance(medialSpline2.getUniformSamples(),high_var=0.05, low_var = 0.001)

    " transform pose 2 by initial offset guess "    
    " transform the new pose "
    points2_offset = []
    for p in points2:
        result = dispPoint(p, offset)        
        #result = dispPointWithAngle(p, offset)
        points2_offset.append(result)


    
    poly1 = []        
    for p in points1:
        poly1.append([p[0],p[1]])    


    while True:
        
        " find the matching pairs "
        match_pairs = []

        " transform the target Hull with the latest offset "
        points2_offset = []
        for p in points2:
            result = dispPoint(p, offset)        
            #result = dispPointWithAngle(p, offset)
            points2_offset.append(result)
        
        " transformed points without associated covariance "
        poly2 = []
        for p in points2_offset:
            poly2.append([p[0],p[1]])    
        
        " get the circles and radii "
        radius2, center2 = computeEnclosingCircle(points2_offset)
        radius1, center1 = computeEnclosingCircle(points1)
    
        for i in range(len(points2_offset)):
            p_2 = poly2[i]

            if isInCircle(p_2, radius1, center1):

                " for every transformed point of A, find it's closest neighbor in B "
                p_1, minDist = findClosestPointInB(points1, p_2, [0.0,0.0,0.0])
    
                if minDist <= minMatchDist:
        
                    " add to the list of match pairs less than 1.0 distance apart "
                    " keep A points and covariances untransformed "
                    C2 = points2[i][2]
                    C1 = p_1[2]
    
                    " we store the untransformed point, but the transformed covariance of the A point "
                    match_pairs.append([points2[i],p_1,C2,C1])

        for i in range(len(points1)):
            p_1 = poly1[i]
    
            if isInCircle(p_1, radius2, center2):
        
                #print "selected", b_p, "in circle", radiusA, centerA, "with distance"
                " for every point of B, find it's closest neighbor in transformed A "
                p_2, i_2, minDist = findClosestPointInA(points2_offset, p_1)
    
                if minDist <= minMatchDist:
        
                    " add to the list of match pairs less than 1.0 distance apart "
                    " keep A points and covariances untransformed "
                    C2 = points2[i_2][2]
                    
                    C1 = points1[i][2]
    
                    " we store the untransformed point, but the transformed covariance of the A point "
                    #match_pairs.append([points2[i_2],p_1,C2,C1])
                    match_pairs.append([points2[i_2],points1[i],C2,C1])




        flatMatchPairs = []
        for pair in match_pairs:
            p1 = pair[0]
            p2 = pair[1]
            C1 = pair[2]
            C2 = pair[3]

            flatMatchPairs.append(p1[0])
            flatMatchPairs.append(p1[1])
            flatMatchPairs.append(C1[0][0])
            flatMatchPairs.append(C1[0][1])
            flatMatchPairs.append(C1[1][0])
            flatMatchPairs.append(C1[1][1])
            flatMatchPairs.append(p2[0])
            flatMatchPairs.append(p2[1])
            flatMatchPairs.append(C2[0][0])
            flatMatchPairs.append(C2[0][1])
            flatMatchPairs.append(C2[1][0])
            flatMatchPairs.append(C2[1][1])
        c_poses_1 = [item for sublist in poses_1 for item in sublist]
        c_poses_2 = [item for sublist in poses_2 for item in sublist]
        newParam, newCost = nelminICP.ICPmin(flatMatchPairs, len(match_pairs), initGuess, c_poses_1, c_poses_2, len(poses_1))


        #newParam = scipy.optimize.fmin(medialOverlapCostFunc, [currU,currAng], [match_pairs, medialSpline1, medialSpline2, uHigh, uLow, u1], disp = 0)
        newU = newParam[0]
        newAng = newParam[1]
        
        
        
        
        """
        # get the current cost
        #if inPlace:
        #    newCost = medialOverlapInPlaceCostFunc([newAng], match_pairs, medialSpline1, medialSpline2, newU, u1)
        #else:
        newCost = medialOverlapCostFunc([newU, newAng], match_pairs, medialSpline1, medialSpline2, uHigh, uLow, u1)
        """
        
        " set the current parameters "
        currAng = functions.normalizeAngle(newAng)
        currU = newU
        
        #print "iteration", numIterations
        if inPlace:
            print "currU =", currU, "currAng =", currAng, "newCost =", newCost

        " compute offset from newU and newAng"

        if u1+0.02 > 1.0:
            pose1 = medialSpline1.getUVecSet([0.98, 1.0])[0]
        elif u1 < 0.0:
            pose1 = medialSpline1.getUVecSet([0.0, 0.02])[0]
        else:
            pose1 = medialSpline1.getUVecSet([u1, u1+0.02])[0]
    
        if currU+0.02 > 1.0:
            pose2 = medialSpline2.getUVecSet([0.98, 1.0])[0]
        elif currU < 0.0:
            pose2 = medialSpline2.getUVecSet([0.0, 0.02])[0]
        else:
            pose2 = medialSpline2.getUVecSet([currU, currU + 0.02])[0]

        #pose1 = medialSpline1.getUVecSet([u1, u1+0.02])[0]
        #pose2 = medialSpline2.getUVecSet([currU, currU+0.02])[0]
    
        point1 = [pose1[0],pose1[1]]
        point2 = [pose2[0],pose2[1]]
        ang1 = pose1[2]
        ang2 = pose2[2]
    
        newOffset = computeOffset(point1, point2, ang1, ang2 + currAng)

        # save the current offset and cost
        offset = newOffset

    
        # check for convergence condition, different between last and current cost is below threshold
        if abs(lastCost - newCost) < costThresh or (numIterations - startIteration) > 10:
            #print "breaking:", abs(lastCost - newCost) < costThresh, (numIterations - startIteration) > 10
            #print "lastCost =", lastCost, "newCost =", newCost, "costThresh =", costThresh
            break



        " update after check for termination condition "
        lastCost = newCost
        
        numIterations += 1

        " reduce the minMatch distance for each step down to a floor value "
        #minMatchDist /= 2
        #if minMatchDist < 0.25:
        #    minMatchDist = 0.25
    
    if False:

        " set the origin of pose 1 "
        poseOrigin = Pose(estPose1)

        fig1 = pylab.figure(1)
        fig1.clf()

        axes2 = fig1.add_subplot(111)

        xP = []
        yP = []
        for b in poly1:
            p1 = poseOrigin.convertLocalToGlobal(b)

            xP.append(p1[0])    
            yP.append(p1[1])

        axes2.plot(xP,yP,linewidth=1, color=(0.0,0.0,1.0))

        xP = []
        yP = []
        for b in points2:
            p = [b[0],b[1]]
            p = dispOffset(p,offset)
            
            p1 = poseOrigin.convertLocalToGlobal(p)
            xP.append(p1[0])    
            yP.append(p1[1])

        axes2.plot(xP,yP,linewidth=1, color=(1.0,0.0,0.0))

        xP = []
        yP = []
        for b in points2:
            p = [b[0],b[1]]
            p = dispOffset(p,gndOffset)
            
            p1 = poseOrigin.convertLocalToGlobal(p)
            xP.append(p1[0])    
            yP.append(p1[1])

        axes2.plot(xP,yP,linewidth=1, color=(1.0,0.5,0.5))

        point2_trans = dispOffset(point2, offset)
        p1 = poseOrigin.convertLocalToGlobal(point1)
        p2 = poseOrigin.convertLocalToGlobal(point2_trans)

        xP = [p1[0],p2[0]]
        yP = [p1[1],p2[1]]
        axes2.scatter(xP,yP,color='b')


        gpac2_trans = dispOffset(rootPose1, offset)
        gpac1 = poseOrigin.convertLocalToGlobal(rootPose2)
        gpac2 = poseOrigin.convertLocalToGlobal(gpac2_trans)

        xP = [gpac1[0],gpac2[0]]
        yP = [gpac1[1],gpac2[1]]
        axes2.scatter(xP,yP,color='r')


        xP = []
        yP = []
        for b in polyB:
            p1 = poseOrigin.convertLocalToGlobal(b)

            xP.append(p1[0])    
            yP.append(p1[1])
        
        axes2.plot(xP,yP,linewidth=1, color=(0.0,0.0,1.0))


        xP = []
        yP = []
        for b in a_data_raw:
            p = [b[0],b[1]]
            p = dispOffset(p,offset)
            
            p1 = poseOrigin.convertLocalToGlobal(p)
            xP.append(p1[0])    
            yP.append(p1[1])

        axes2.plot(xP,yP,linewidth=1, color=(1.0,0.0,0.0))

        xP = []
        yP = []
        for b in a_data_raw:
            p = [b[0],b[1]]
            p = dispOffset(p,gndOffset)
            
            p1 = poseOrigin.convertLocalToGlobal(p)
            xP.append(p1[0])    
            yP.append(p1[1])

        axes2.plot(xP,yP,linewidth=1, color=(1.0,0.5,0.5))

        plotEnv(axes2)        
        
        axes2.set_title("%s %s, ang = %1.3f" % (repr(n1), repr(n2), currAng))
        axes2.set_xlim(estPose1[0]-4, estPose1[0]+4)                    
        axes2.set_ylim(estPose1[1]-3, estPose1[1]+3)
        fig1.savefig("ICP_plot_%04u_0.png" % numIterations)
        pylab.clf()            
        pylab.figure(1)

        numIterations += 1


        
    histogram = overlapHistogram([currU,currAng], match_pairs, medialSpline1, medialSpline2, u1)            


    offset[2] =  functions.normalizeAngle(offset[2])    
    return offset, histogram



def globalOverlapICP(initGuess, globalPath, medialPoints,  plotIter = False, n1 = 0, n2 = 0):

    global numIterations
    global globalPlotCount
    
    def computeOffset(point1, point2, ang1, ang2):
    
        " corner points and orientations "
        corner1Pose = [point1[0], point1[1], ang1]
        corner2Pose = [point2[0], point2[1], ang2]
        
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

    u1 = initGuess[0]
    u2 = initGuess[1]

    uHigh = u2 + 0.2
    uLow = u2 - 0.2

    currU = u2
    currAng = initGuess[2]
    
    globalSpline = SplineFit(globalPath, smooth=0.1)
    medialSpline = SplineFit(medialPoints, smooth=0.1)


    uSet = [i*0.01 for i in range(100)]
    poses_1 = globalSpline.getUVecSet(uSet)
    poses_2 = medialSpline.getUVecSet(uSet)

    if u1 >= 1.0:
        pose1 = poses_1[-1]
    elif u1 < 0.0:
        pose1 = poses_1[0]
    else:
        pose1 = poses_1[int(u1*100)]

    if currU >= 1.0:
        pose2 = poses_2[-1]
    elif currU < 0.0:
        pose2 = poses_2[0]
    else:
        pose2 = poses_2[int(currU*100)]

    #if u1+0.02 > 1.0:
    #    pose1 = globalSpline.getUVecSet([0.98, 1.0])[0]
    #elif u1 < 0.0:
    #    pose1 = globalSpline.getUVecSet([0.0, 0.02])[0]
    #else:
    #    pose1 = globalSpline.getUVecSet([u1, u1+0.02])[0]

    #if currU+0.02 > 1.0:
    #    pose2 = medialSpline.getUVecSet([0.98, 1.0])[0]
    #elif currU < 0.0:
    #    pose2 = medialSpline.getUVecSet([0.0, 0.02])[0]
    #else:
    #    pose2 = medialSpline.getUVecSet([currU, currU + 0.02])[0]
    

    #pose2 = medialSpline.getUVecSet([currU, currU+0.02])[0]
    #pose1 = globalSpline.getUVecSet([u1, u1 + 0.02])[0]

    #pathOrigin = Pose(pose1)
    #localPath = []
    #for p in globalPath:
    #    localPath.append(pathOrigin.convertGlobalToLocal(p))        

    #localSpline = SplineFit(localPath, smooth=0.1)
    #pose1 = localSpline.getUVecSet([u1, u1 + 0.02])[0]
        
    point1 = [pose1[0],pose1[1]]
    point2 = [pose2[0],pose2[1]]
    ang1 = pose1[2]
    ang2 = pose2[2]

    currPose = computeOffset(point1, point2, ang1, ang2 + currAng)

    " transform the new pose "
    poseOrigin = Pose(currPose)
    #a_data_raw = hull
    #a_data = []
    #for p in a_data_raw:
    #    #result = poseOrigin.convertLocalToGlobal(p)
    #    result = dispOffset(p, currPose)        
    #    a_data.append(result)    
    
    costThresh = 0.004
    minMatchDist = 2.0
    lastCost = 1e100

    startIteration = numIterations

    " set the initial guess "
    poseOrigin = Pose(currPose)

    " sample a series of points from the medial curves "
    " augment points with point-to-line covariances "
    " treat the points with the point-to-line constraint "
    globalPoints = addGPACVectorCovariance(globalSpline.getUniformSamples(),high_var=0.05, low_var = 0.05)
    #localPoints = addGPACVectorCovariance(localSpline.getUniformSamples(),high_var=0.05, low_var = 0.001)
    points = addGPACVectorCovariance(medialSpline.getUniformSamples(),high_var=0.05, low_var = 0.05)

    " transform pose 2 by initial offset guess "    
    " transform the new pose "
    points_offset = []
    #for p in points:
    #    result = dispPoint(p, offset)        
    #    points_offset.append(result)

    globalPoly = []        
    for p in globalPoints:
        globalPoly.append([p[0],p[1]])    

    #localPoly = []        
    #for p in localPoints:
    #    localPoly.append([p[0],p[1]])    

    while True:
        
        " find the matching pairs "
        match_pairs = []

        " transform the target Hull with the latest offset "
        points_offset = []
        for p in points:
            result = dispPoint(p, currPose)        
            points_offset.append(result)
        
        " transformed points without associated covariance "
        posePoly = []
        for p in points_offset:
            posePoly.append([p[0],p[1]])    
        
        " get the circles and radii "
        radius2, center2 = computeEnclosingCircle(points_offset)
    
        """
        for i in range(len(points_offset)):
            p_2 = posePoly[i]

            " for every transformed point of A, find it's closest neighbor in B "
            p_1, minDist = findClosestPointInB(globalPoints, p_2, [0.0,0.0,0.0])
            #p_1, minDist = findClosestPointInB(localPoints, p_2, [0.0,0.0,0.0])

            if minDist <= minMatchDist:
    
                " add to the list of match pairs less than 1.0 distance apart "
                " keep A points and covariances untransformed "
                C2 = points[i][2]
                C1 = p_1[2]

                " we store the untransformed point, but the transformed covariance of the A point "
                match_pairs.append([points[i],p_1,C2,C1])
        """

        for i in range(len(globalPoints)):
            p_1 = globalPoly[i]
            #p_1 = localPoly[i]
    
            if isInCircle(p_1, radius2, center2):
        
                " for every point of B, find it's closest neighbor in transformed A "
                p_2, i_2, minDist = findClosestPointInA(points_offset, p_1)
    
                if minDist <= minMatchDist:
        
                    " add to the list of match pairs less than 1.0 distance apart "
                    " keep A points and covariances untransformed "
                    C2 = points[i_2][2]
                    
                    C1 = globalPoints[i][2]
                    #C1 = localPoints[i][2]
    
                    " we store the untransformed point, but the transformed covariance of the A point "
                    match_pairs.append([points[i_2],globalPoints[i],C2,C1])
                    #match_pairs.append([points[i_2],localPoints[i],C2,C1])

        #if plotIter and startIteration == numIterations:
        if plotIter:
            

            #center2Global = poseOrigin.convertLocalToGlobal(center2)
            #cir = pylab.Circle((center2Global[0],center2Global[1]), radius=radius2,  ec='r', fc='none')
            #pylab.gca().add_patch(cir)


            #point_trans = dispOffset(point2, offset)
            #p1 = point1
            #p2 = poseOrigin.convertLocalToGlobal(point_trans)

            #xP = [p1[0],p2[0]]
            #yP = [p1[1],p2[1]]
            #pylab.scatter(xP,yP,color='b')


            #gpac2_trans = dispOffset(rootPose1, offset)
            #gpac1 = poseOrigin.convertLocalToGlobal(rootPose2)
            #gpac2 = poseOrigin.convertLocalToGlobal(gpac2_trans)
            #xP = [gpac1[0],gpac2[0]]
            #yP = [gpac1[1],gpac2[1]]
            #pylab.scatter(xP,yP,color='r')

            
            """
            xP = []
            yP = []
            for p1 in globalPath:
                xP.append(p1[0])    
                yP.append(p1[1])
            
            pylab.plot(xP,yP,linewidth=1, color=(0.0,0.0,1.0))
            """

            fig1 = pylab.figure(1)
            fig1.clf()

            axes2 = fig1.add_subplot(111)

            match_global = []
            
            for pair in match_pairs:
                p1 = pair[0]
                p2 = pair[1]
                
                p1_o = dispOffset(p1, currPose)
                
                p1_g = p1_o
                p2_g = p2
                match_global.append([p1_g,p2_g])

            draw_matches(match_global, [0.0,0.0,0.0], axes2)

            
            xP = []
            yP = []
            for b in globalPoly:
                p1 = b
                xP.append(p1[0])    
                yP.append(p1[1])
    
            axes2.plot(xP,yP,linewidth=1, color=(0.0,0.0,1.0))

            #axes2.scatter([globalPoly[0][0]],[globalPoly[0][1]],color=(0.0,0.0,1.0))
            axes2.scatter([pose1[0]],[pose1[1]],color=(0.0,0.0,1.0))

            xP = []
            yP = []
            for b in points:
                p = [b[0],b[1]]
                p1 = dispOffset(p,currPose)
                
                xP.append(p1[0])    
                yP.append(p1[1])

            axes2.plot(xP,yP,linewidth=1, color=(1.0,0.0,0.0))

            u0 = points[0]
            u0L = dispOffset(u0,currPose)

            #axes2.scatter([u0L[0]],[u0L[1]],color=(1.0,0.0,0.0))

            #xP = []
            #yP = []
            #for b in a_data_raw:
            #    p = [b[0],b[1]]
            #    p1 = dispOffset(p,currPose)
            #    
            #    #p1 = poseOrigin.convertLocalToGlobal(p)
            #    xP.append(p1[0])    
            #    yP.append(p1[1])

            #axes2.plot(xP,yP,linewidth=1, color=(1.0,0.0,0.0))



            
            cost = medialOverlapCostFunc([currU, currAng], match_pairs, poses_1, poses_2, uHigh, uLow, u1)            

            plotEnv(axes2)                    
            axes2.set_title("%u, u1 = %1.3f, u2 = %1.3f, ang = %1.3f, cost = %f" % (n1, u1, currU, currAng, cost))

            #axes2.set_title("%s %s, ang = %1.3f" % (repr(n1), repr(n2), currAng))
            axes2.set_xlim(currPose[0]-4, currPose[0]+4)                    
            axes2.set_ylim(currPose[1]-3, currPose[1]+3)
            
            """
            axes1 = fig1.add_subplot(121)
            axes1.grid(True)
            
            gPoints = globalSpline.getTransformCurve()
            angOffset = -gPoints[0][1]
            gPoints = globalSpline.getTransformCurve(angOffset)

            dist1 = globalSpline.dist(0.0, u1, iter = 0.005)                

            xP = []
            yP = []
            for p in gPoints:
                xP.append(p[0])
                yP.append(p[1])
            axes1.plot(xP,yP, color='b')

            medialPoints_trans = []
            for b in medialPoints:
                p = [b[0],b[1]]
                p = dispOffset(p,currPose)
                
                medialPoints_trans.append(p)

            medialSpline_trans = SplineFit(medialPoints_trans, smooth=0.1)
            tPoints2_trans = medialSpline_trans.getTransformCurve(offset = angOffset, compareAngle = gPoints[0][1])




            dist2 = medialSpline_trans.dist(0.0, currU, iter = 0.005)    

            xP = []
            yP = []
            for p in tPoints2_trans:
                xP.append(p[0]+dist1-dist2)
                yP.append(p[1])
            axes1.plot(xP,yP, color='r')

            
            axes1.set_xlabel("distance")
            axes1.set_ylabel("angle (radians)")
            axes1.set_xlim(-2,10)
            axes1.set_ylim(-6,6)        
            """
            fig1.savefig("ICP_plot_%04u_0.png" % numIterations)
            pylab.clf()            
            
            pylab.figure(1)            
            


        #newParam = scipy.optimize.fmin(medialOverlapCostFunc, [currU,currAng], [match_pairs, localSpline, medialSpline, uHigh, uLow, u1], disp = 0)
        newParam = scipy.optimize.fmin(medialOverlapCostFunc, [currU,currAng], [match_pairs, poses_1, poses_2, uHigh, uLow, u1], disp = 0)
        newU = newParam[0]
        newAng = newParam[1]
        
        #newCost = medialOverlapCostFunc([newU, newAng], match_pairs, localSpline, medialSpline, uHigh, uLow, u1)
        newCost = medialOverlapCostFunc([newU, newAng], match_pairs, poses_1, poses_2, uHigh, uLow, u1)

        " set the current parameters "
        currAng = functions.normalizeAngle(newAng)
        currU = newU
        
        " compute offset from newU and newAng"
        #if u1 >= 0.98:
        #    pose1 = globalSpline.getUVecSet([1.0-0.01, 1.0])[0]
        #else:                
        #    pose1 = globalSpline.getUVecSet([u1, u1+0.02])[0]
        #pose1 = localSpline.getUVecSet([u1, u1+0.02])[0]
        #if currU >= 0.98:
        #    pose2 = medialSpline.getUVecSet([1.0 - 0.01,1.0])[0]
        #else:
        #    pose2 = medialSpline.getUVecSet([currU, currU+0.02])[0]


        if u1+0.02 > 1.0:
            pose1 = globalSpline.getUVecSet([0.98, 1.0])[0]
        elif u1 < 0.0:
            pose1 = globalSpline.getUVecSet([0.0, 0.02])[0]
        else:
            pose1 = globalSpline.getUVecSet([u1, u1+0.02])[0]
    
        if currU+0.02 > 1.0:
            pose2 = medialSpline.getUVecSet([0.98, 1.0])[0]
        elif currU < 0.0:
            pose2 = medialSpline.getUVecSet([0.0, 0.02])[0]
        else:
            pose2 = medialSpline.getUVecSet([currU, currU + 0.02])[0]
            


        point1 = [pose1[0],pose1[1]]
        point2 = [pose2[0],pose2[1]]
        ang1 = pose1[2]
        ang2 = pose2[2]
    
        currPose = computeOffset(point1, point2, ang1, ang2 + currAng)
    
        #point1 = [pose1[0],pose1[1]]
        #point2 = [pose2[0],pose2[1]]
        #ang1 = pose1[2]
        #ang2 = pose2[2]
    
        #newOffset = computeOffset(point1, point2, ang1, ang2 + currAng)

        # save the current offset and cost
        #offset = newOffset


        # check for convergence condition, different between last and current cost is below threshold
        
        isTerminate = False
        if abs(lastCost - newCost) < costThresh or (numIterations - startIteration) > 10:
            isTerminate = True
            print "terminating globalOverlap:", lastCost, newCost, costThresh, numIterations-startIteration
        
        #if abs(lastCost - newCost) < costThresh or (numIterations - startIteration) > 10:
        #    break



        " update after check for termination condition "
        lastCost = newCost

        # optionally draw the position of the points in current transform
        if plotIter:
            

            fig1 = pylab.figure(1)
            fig1.clf()

            axes2 = fig1.add_subplot(111)

            match_global = []
            
            for pair in match_pairs:
                p1 = pair[0]
                p2 = pair[1]
                
                p1_o = dispOffset(p1, currPose)
                
                p1_g = p1_o
                p2_g = p2
                match_global.append([p1_g,p2_g])

            draw_matches(match_global, [0.0,0.0,0.0], axes2)

            
            xP = []
            yP = []
            for b in globalPoly:
                p1 = b
                xP.append(p1[0])    
                yP.append(p1[1])
    
            axes2.plot(xP,yP,linewidth=1, color=(0.0,0.0,1.0))

            #axes2.scatter([globalPoly[0][0]],[globalPoly[0][1]],color=(0.0,0.0,1.0))
            axes2.scatter([pose1[0]],[pose1[1]],color=(0.0,0.0,1.0))

            xP = []
            yP = []
            for b in points:
                p = [b[0],b[1]]
                p1 = dispOffset(p,currPose)
                
                xP.append(p1[0])    
                yP.append(p1[1])

            axes2.plot(xP,yP,linewidth=1, color=(1.0,0.0,0.0))

            u0 = points[0]
            u0L = dispOffset(u0,currPose)

            #axes2.scatter([u0L[0]],[u0L[1]],color=(1.0,0.0,0.0))

            #xP = []
            #yP = []
            #for b in a_data_raw:
            #    p = [b[0],b[1]]
            #    p1 = dispOffset(p,currPose)
            #    
            #    #p1 = poseOrigin.convertLocalToGlobal(p)
            #    xP.append(p1[0])    
            #    yP.append(p1[1])

            #axes2.plot(xP,yP,linewidth=1, color=(1.0,0.0,0.0))



            
            cost = medialOverlapCostFunc([currU, currAng], match_pairs, globalSpline, medialSpline, uHigh, uLow, u1)            

            plotEnv(axes2)                    
            axes2.set_title("%u, u1 = %1.3f, u2 = %1.3f, ang = %1.3f, cost = %f" % (n1, u1, currU, currAng, cost))

            #axes2.set_title("%s %s, ang = %1.3f" % (repr(n1), repr(n2), currAng))
            axes2.set_xlim(currPose[0]-4, currPose[0]+4)                    
            axes2.set_ylim(currPose[1]-3, currPose[1]+3)
            
            """

            axes1 = fig1.add_subplot(121)
            axes1.grid(True)
            
            gPoints = globalSpline.getTransformCurve()
            angOffset = -gPoints[0][1]
            gPoints = globalSpline.getTransformCurve(angOffset)

            dist1 = globalSpline.dist(0.0, u1, iter = 0.005)                

            xP = []
            yP = []
            for p in gPoints:
                xP.append(p[0])
                yP.append(p[1])
            axes1.plot(xP,yP, color='b')

            medialPoints_trans = []
            for b in medialPoints:
                p = [b[0],b[1]]
                p = dispOffset(p,currPose)
                
                medialPoints_trans.append(p)

            medialSpline_trans = SplineFit(medialPoints_trans, smooth=0.1)
            #tPoints2_trans = medialSpline_trans.getTransformCurve(angOffset)
            tPoints2_trans = medialSpline_trans.getTransformCurve(offset = angOffset, compareAngle = gPoints[0][1])

            dist2 = medialSpline_trans.dist(0.0, currU, iter = 0.005)    

            xP = []
            yP = []
            for p in tPoints2_trans:
                xP.append(p[0]+dist1-dist2)
                yP.append(p[1])
            axes1.plot(xP,yP, color='r')

            
            axes1.set_xlabel("distance")
            axes1.set_ylabel("angle (radians)")
            axes1.set_xlim(-2,10)
            axes1.set_ylim(-6,6)        
            """
            fig1.savefig("ICP_plot_%04u_0.png" % numIterations)
            pylab.clf()            
            
            pylab.figure(1)            

        
        numIterations += 1
        
        if isTerminate:
            break

    " draw final position "
    if True:
        
        " set the origin of pose 1 "
        poseOrigin = Pose(currPose)
        
        pylab.clf()
        pylab.axes()
        match_global = []
        
        for pair in match_pairs:
            p1 = pair[0]
            p2 = pair[1]
            
            p1_o = dispOffset(p1, currPose)
            
            p1_g = p1_o
            p2_g = p2
            match_global.append([p1_g,p2_g])

        draw_matches(match_global, [0.0,0.0,0.0])

        
        xP = []
        yP = []
        for b in globalPoly:
            p1 = b
            xP.append(p1[0])    
            yP.append(p1[1])

        pylab.plot(xP,yP,linewidth=1, color=(0.0,0.0,1.0))

        pylab.scatter([pose1[0]],[pose1[1]],color=(0.0,0.0,1.0))

        xP = []
        yP = []
        for b in points:
            p = [b[0],b[1]]
            p1 = dispOffset(p,currPose)
            
            #p1 = poseOrigin.convertLocalToGlobal(p)
            xP.append(p1[0])    
            yP.append(p1[1])

        pylab.plot(xP,yP,linewidth=1, color=(1.0,0.0,0.0))

        u0 = points[0]
        u0L = dispOffset(u0,currPose)

       # pylab.scatter([u0L[0]],[u0L[1]],color=(1.0,0.0,0.0))

        #xP = []
        #yP = []
        #for b in a_data_raw:
        #    p = [b[0],b[1]]
        #    p1 = dispOffset(p,currPose)
        #    
        #    #p1 = poseOrigin.convertLocalToGlobal(p)
        #    xP.append(p1[0])    
        #    yP.append(p1[1])

        pylab.plot(xP,yP,linewidth=1, color=(1.0,0.0,0.0))

        plotEnv()        
        pylab.title("Py: %u, u1 = %1.3f, u2 = %1.3f, ang = %1.3f, cost = %f" % (n1, u1, currU, currAng, newCost))

        pylab.xlim(currPose[0]-4, currPose[0]+4)                    
        pylab.ylim(currPose[1]-3, currPose[1]+3)
        pylab.savefig("ICP_plot_%04u_1.png" % globalPlotCount)
        pylab.clf()
        " save inputs "
        saveFile = ""
        saveFile += "initGuess = " + repr(initGuess) + "\n"
        saveFile += "globalPath = " + repr(globalPath) + "\n"
        saveFile += "medialPoints = " + repr(medialPoints) + "\n"

        f = open("icpInputSave_%04u.txt" % globalPlotCount, 'w')
        f.write(saveFile)
        f.close()        
            
        globalPlotCount += 1
        
    #histogram = overlapHistogram([currU,currAng], match_pairs, medialSpline1, medialSpline2, u1)            


    #offset[2] =  functions.normalizeAngle(offset[2])    
    
    return currPose, newCost
    #return offset
    #return offset, histogram


def globalOverlapICP_GPU(initGuess, globalPath, medialPoints, poses_1, poses_2, plotIter = False, n1 = 0, n2 = 0):

    global numIterations
    global globalPlotCount
    
    def computeOffset(point1, point2, ang1, ang2):
    
        " corner points and orientations "
        corner1Pose = [point1[0], point1[1], ang1]
        corner2Pose = [point2[0], point2[1], ang2]
        
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

    
    #saveStr = ""
    #saveStr += "initGuess = " + repr(initGuess) + '\n'
    #saveStr += "poses_1 = " + repr([item for sublist in poses_1 for item in sublist]) + '\n'
    #saveStr += "poses_2 = " + repr([item for sublist in poses_2 for item in sublist]) + '\n'
    #saveStr += "numPoses = " + repr(len(poses_1)) + '\n'

    u1 = initGuess[0]
    u2 = initGuess[1]

    uHigh = u2 + 0.2
    uLow = u2 - 0.2

    currU = u2
    currAng = initGuess[2]
    
    globalSpline = SplineFit(globalPath, smooth=0.1)
    medialSpline = SplineFit(medialPoints, smooth=0.1)

    if u1 >= 1.0:
        pose1 = poses_1[-1]
    elif u1 < 0.0:
        pose1 = poses_1[0]
    else:
        pose1 = poses_1[int(u1*100)]

    if currU >= 1.0:
        pose2 = poses_2[-1]
    elif currU < 0.0:
        pose2 = poses_2[0]
    else:
        pose2 = poses_2[int(currU*100)]
                    
    point1 = [pose1[0],pose1[1]]
    point2 = [pose2[0],pose2[1]]
    ang1 = pose1[2]
    ang2 = pose2[2]

    currPose = computeOffset(point1, point2, ang1, ang2 + currAng)

    " transform the new pose "
    poseOrigin = Pose(currPose)
    
    costThresh = 0.004
    minMatchDist = 2.0
    lastCost = 1e100
    
    startIteration = numIterations

    " set the initial guess "
    poseOrigin = Pose(currPose)
    
    " sample a series of points from the medial curves "
    " augment points with point-to-line covariances "
    " treat the points with the point-to-line constraint "
    globalPoints = addGPACVectorCovariance(globalSpline.getUniformSamples(),high_var=0.05, low_var = 0.05)
    #localPoints = addGPACVectorCovariance(localSpline.getUniformSamples(),high_var=0.05, low_var = 0.001)
    points = addGPACVectorCovariance(medialSpline.getUniformSamples(),high_var=0.05, low_var = 0.05)

    " transform pose 2 by initial offset guess "    
    " transform the new pose "
    points_offset = []

    globalPoly = []        
    for p in globalPoints:
        globalPoly.append([p[0],p[1]])    


    while True:
        
        " find the matching pairs "
        match_pairs = []

        " transform the target Hull with the latest offset "
        points_offset = []
        for p in points:
            result = dispPoint(p, currPose)        
            points_offset.append(result)
        
        " transformed points without associated covariance "
        posePoly = []
        for p in points_offset:
            posePoly.append([p[0],p[1]])    
        
        " get the circles and radii "
        radius2, center2 = computeEnclosingCircle(points_offset)

        for i in range(len(globalPoints)):
            p_1 = globalPoly[i]
            #p_1 = localPoly[i]
    
            if isInCircle(p_1, radius2, center2):
        
                " for every point of B, find it's closest neighbor in transformed A "
                p_2, i_2, minDist = findClosestPointInA(points_offset, p_1)
    
                if minDist <= minMatchDist:
        
                    " add to the list of match pairs less than 1.0 distance apart "
                    " keep A points and covariances untransformed "
                    C2 = points[i_2][2]
                    
                    C1 = globalPoints[i][2]
    
                    " we store the untransformed point, but the transformed covariance of the A point "
                    match_pairs.append([points[i_2],globalPoints[i],C1,C2])

        input_GPU = convertToGPU(match_pairs)

        #oldCost = medialOverlapCostFunc([currU, currAng], match_pairs, globalSpline, medialSpline, uHigh, uLow, u1)
        oldCost = medialOverlapCostFunc_GPU([currU, currAng], input_GPU, len(match_pairs), globalSpline, medialSpline, uHigh, uLow, u1)
        
        newParam = scipy.optimize.fmin(medialOverlapCostFunc_GPU, [currU,currAng], [input_GPU, len(match_pairs), globalSpline, medialSpline, uHigh, uLow, u1], disp = 0)
        newU = newParam[0]
        newAng = newParam[1]
        
        #newCost = medialOverlapCostFunc([newU, newAng], match_pairs, globalSpline, medialSpline, uHigh, uLow, u1)
        newCost = medialOverlapCostFunc_GPU([newU, newAng], input_GPU, len(match_pairs), globalSpline, medialSpline, uHigh, uLow, u1)


        #flatMatchPairs = []
        #for pair in match_pairs:
        #    p1 = pair[0]
        #    p2 = pair[1]
        #    C1 = pair[2]
        #    C2 = pair[3]

        #    flatMatchPairs.append(p1[0])
        #    flatMatchPairs.append(p1[1])
        #    flatMatchPairs.append(C1[0][0])
        #    flatMatchPairs.append(C1[0][1])
        #    flatMatchPairs.append(C1[1][0])
        #    flatMatchPairs.append(C1[1][1])
        #    flatMatchPairs.append(p2[0])
        #    flatMatchPairs.append(p2[1])
        #    flatMatchPairs.append(C2[0][0])
        #    flatMatchPairs.append(C2[0][1])
        #    flatMatchPairs.append(C2[1][0])
        #    flatMatchPairs.append(C2[1][1])

        #saveStr += "matchPairs = " + repr(flatMatchPairs) + '\n'
        #saveStr += "numPairs = " + repr(len(match_pairs)) + '\n'

        #f = open("cudaData3.txt", 'w')
        #f.write(saveStr)
        #f.close()

        #exit()

        " set the current parameters "
        currAng = functions.normalizeAngle(newAng)
        currU = newU


        #if currU+0.02 > 1.0:
        #    pose2 = medialSpline.getUVecSet([0.98, 1.0])[0]
        #elif currU < 0.0:
        #    pose2 = medialSpline.getUVecSet([0.0, 0.02])[0]
        #else:
        #    pose2 = medialSpline.getUVecSet([currU, currU + 0.02])[0]
            
        if currU >= 1.0:
            pose2 = poses_2[-1]
        elif currU < 0.0:
            pose2 = poses_2[0]
        else:
            pose2 = poses_2[int(currU*100)]


        point1 = [pose1[0],pose1[1]]
        point2 = [pose2[0],pose2[1]]
        ang1 = pose1[2]
        ang2 = pose2[2]
    
        currPose = computeOffset(point1, point2, ang1, ang2 + currAng)

        # check for convergence condition, different between last and current cost is below threshold
        
        isTerminate = False
        if abs(lastCost - newCost) < costThresh or (numIterations - startIteration) > 10:
            isTerminate = True
            print "terminating globalOverlap:", lastCost, newCost, costThresh, numIterations-startIteration
        
        #if abs(lastCost - newCost) < costThresh or (numIterations - startIteration) > 10:
        #    break



        " update after check for termination condition "
        lastCost = newCost

        numIterations += 1
        
        if isTerminate:
            break

    " draw final position "
    if True:
        
        " set the origin of pose 1 "
        poseOrigin = Pose(currPose)
        
        pylab.clf()
        pylab.axes()
        match_global = []
        
        for pair in match_pairs:
            p1 = pair[0]
            p2 = pair[1]
            
            p1_o = dispOffset(p1, currPose)
            
            p1_g = p1_o
            p2_g = p2
            match_global.append([p1_g,p2_g])

        draw_matches(match_global, [0.0,0.0,0.0])

        
        xP = []
        yP = []
        for b in globalPoly:
            p1 = b
            xP.append(p1[0])    
            yP.append(p1[1])

        pylab.plot(xP,yP,linewidth=1, color=(0.0,0.0,1.0))

        pylab.scatter([pose1[0]],[pose1[1]],color=(0.0,0.0,1.0))

        xP = []
        yP = []
        for b in points:
            p = [b[0],b[1]]
            p1 = dispOffset(p,currPose)
            
            #p1 = poseOrigin.convertLocalToGlobal(p)
            xP.append(p1[0])    
            yP.append(p1[1])

        pylab.plot(xP,yP,linewidth=1, color=(1.0,0.0,0.0))

        #pylab.scatter([pose2[0]],[pose2[1]],color=(1.0,0.0,0.0))


        plotEnv()        
        pylab.title("%u, u1 = %1.3f, u2 = %1.3f, ang = %1.3f, cost = %f" % (n1, u1, currU, currAng, newCost))

        pylab.xlim(currPose[0]-4, currPose[0]+4)                    
        pylab.ylim(currPose[1]-3, currPose[1]+3)
        pylab.savefig("ICP_plot_%04u.png" % globalPlotCount)
        pylab.clf()
        
        " save inputs "
        saveFile = ""
        saveFile += "initGuess = " + repr(initGuess) + "\n"
        saveFile += "globalPath = " + repr(globalPath) + "\n"
        #saveFile += "hull = " + repr(hull) + "\n"
        saveFile += "medialPoints = " + repr(medialPoints) + "\n"

        f = open("icpInputSave_%04u.txt" % globalPlotCount, 'w')
        globalPlotCount += 1
        f.write(saveFile)
        f.close()        
            
        #numIterations += 1


    return currPose, newCost


def globalOverlapICP_GPU2(initGuess, globalPath, medialPoints,plotIter = False, n1 = 0, n2 = 0):

    global numIterations
    global globalPlotCount
    
    def computeOffset(point1, point2, ang1, ang2):
    
        " corner points and orientations "
        corner1Pose = [point1[0], point1[1], ang1]
        corner2Pose = [point2[0], point2[1], ang2]
        
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

    
    u1 = initGuess[0]
    u2 = initGuess[1]

    uHigh = u2 + 0.2
    uLow = u2 - 0.2

    currU = u2
    currAng = initGuess[2]
    
    globalSpline = SplineFit(globalPath, smooth=0.1)
    medialSpline = SplineFit(medialPoints, smooth=0.1)
    

    uSet = [i*0.01 for i in range(100)]
    poses_1 = globalSpline.getUVecSet(uSet)
    poses_2 = medialSpline.getUVecSet(uSet)
    
    #globalSpline = SplineFit(globalPath, smooth=0.1)
    #medialSpline = SplineFit(medialPoints, smooth=0.1)

    if u1 >= 1.0:
        pose1 = poses_1[-1]
    elif u1 < 0.0:
        pose1 = poses_1[0]
    else:
        pose1 = poses_1[int(u1*100)]

    if currU >= 1.0:
        pose2 = poses_2[-1]
    elif currU < 0.0:
        pose2 = poses_2[0]
    else:
        pose2 = poses_2[int(currU*100)]
    
    point1 = [pose1[0],pose1[1]]
    point2 = [pose2[0],pose2[1]]
    ang1 = pose1[2]
    ang2 = pose2[2]

    currPose = computeOffset(point1, point2, ang1, ang2 + currAng)

    " transform the new pose "
    poseOrigin = Pose(currPose)
    
    costThresh = 0.004
    minMatchDist = 2.0
    lastCost = 1e100
    
    startIteration = numIterations

    " set the initial guess "
    poseOrigin = Pose(currPose)
    
    " sample a series of points from the medial curves "
    " augment points with point-to-line covariances "
    " treat the points with the point-to-line constraint "
    globalPoints = addGPACVectorCovariance(globalSpline.getUniformSamples(),high_var=0.05, low_var = 0.001)
    #globalPoints = addGPACVectorCovariance(globalSpline.getUniformSamples(),high_var=0.05, low_var = 0.05)
    #localPoints = addGPACVectorCovariance(localSpline.getUniformSamples(),high_var=0.05, low_var = 0.001)
    points = addGPACVectorCovariance(medialSpline.getUniformSamples(),high_var=0.05, low_var = 0.001)
    #points = addGPACVectorCovariance(medialSpline.getUniformSamples(),high_var=0.05, low_var = 0.05)

    " transform pose 2 by initial offset guess "    
    " transform the new pose "
    points_offset = []

    globalPoly = []        
    for p in globalPoints:
        globalPoly.append([p[0],p[1]])    


    while True:
        
        " find the matching pairs "
        match_pairs = []

        " transform the target Hull with the latest offset "
        points_offset = []
        for p in points:
            result = dispPoint(p, currPose)        
            points_offset.append(result)
        
        " transformed points without associated covariance "
        #posePoly = []
        #for p in points_offset:
        #    posePoly.append([p[0],p[1]])
        
        #match_pairs = matchPairs2(points, points_offset, globalPoints, minMatchDist)

        #print points_offset

        #cProfile.run('match_pairs = matchPairs(points, points_offset, globalPoints, minMatchDist)', 'match_prof')
        #cProfile.runctx('match_pairs = matchPairs(points, points_offset, globalPoints, minMatchDist)', globals(), locals(), 'match_prof')

        match_pairs = matchPairs(points, points_offset, globalPoints, minMatchDist)       
        
        """
        " get the circles and radii "
        radius2, center2 = computeEnclosingCircle(points_offset)

        for i in range(len(globalPoints)):
            p_1 = globalPoly[i]
            #p_1 = localPoly[i]
    
            if isInCircle(p_1, radius2, center2):
        
                " for every point of B, find it's closest neighbor in transformed A "
                p_2, i_2, minDist = findClosestPointInA(points_offset, p_1)
    
                if minDist <= minMatchDist:
        
                    " add to the list of match pairs less than 1.0 distance apart "
                    " keep A points and covariances untransformed "
                    C2 = points[i_2][2]
                    
                    C1 = globalPoints[i][2]
    
                    " we store the untransformed point, but the transformed covariance of the A point "
                    match_pairs.append([points[i_2],globalPoints[i],C1,C2])
        """
        #input_GPU = convertToGPU(match_pairs)

        #oldCost = medialOverlapCostFunc([currU, currAng], match_pairs, globalSpline, medialSpline, uHigh, uLow, u1)
        #oldCost = medialOverlapCostFunc_GPU([currU, currAng], input_GPU, len(match_pairs), globalSpline, medialSpline, uHigh, uLow, u1)
        

        #newParam = scipy.optimize.fmin(medialOverlapCostFunc_GPU, [currU,currAng], [input_GPU, len(match_pairs), globalSpline, medialSpline, uHigh, uLow, u1], disp = 0)
        
        flatMatchPairs = []
        for pair in match_pairs:
            p1 = pair[0]
            p2 = pair[1]
            C1 = pair[2]
            C2 = pair[3]

            flatMatchPairs.append(p1[0])
            flatMatchPairs.append(p1[1])
            flatMatchPairs.append(C1[0][0])
            flatMatchPairs.append(C1[0][1])
            flatMatchPairs.append(C1[1][0])
            flatMatchPairs.append(C1[1][1])
            flatMatchPairs.append(p2[0])
            flatMatchPairs.append(p2[1])
            flatMatchPairs.append(C2[0][0])
            flatMatchPairs.append(C2[0][1])
            flatMatchPairs.append(C2[1][0])
            flatMatchPairs.append(C2[1][1])
        c_poses_1 = [item for sublist in poses_1 for item in sublist]
        c_poses_2 = [item for sublist in poses_2 for item in sublist]
        newParam, newCost = nelminICP.ICPmin(flatMatchPairs, len(match_pairs), initGuess, c_poses_1, c_poses_2, len(poses_1))

        newU = newParam[0]
        newAng = newParam[1]
        
        #newCost = medialOverlapCostFunc([newU, newAng], match_pairs, globalSpline, medialSpline, uHigh, uLow, u1)
        #newCost = medialOverlapCostFunc_GPU([newU, newAng], input_GPU, len(match_pairs), globalSpline, medialSpline, uHigh, uLow, u1)


        #saveStr = ""
        #saveStr += "initGuess = " + repr(initGuess) + '\n'
        #saveStr += "match_pairs = " + repr(match_pairs) + '\n'
        #saveStr += "numPairs = " + repr(len(match_pairs)) + '\n'
        #saveStr += "poses_1 = " + repr(poses_1) + '\n'
        #saveStr += "poses_2 = " + repr(poses_2) + '\n'

        #f = open("costTest.txt", 'w')
        #f.write(saveStr)
        #f.close()

        #exit()
        " set the current parameters "
        currAng = functions.normalizeAngle(newAng)
        currU = newU


        #if currU+0.02 > 1.0:
        #    pose2 = medialSpline.getUVecSet([0.98, 1.0])[0]
        #elif currU < 0.0:
        #    pose2 = medialSpline.getUVecSet([0.0, 0.02])[0]
        #else:
        #    pose2 = medialSpline.getUVecSet([currU, currU + 0.02])[0]
            
        if currU >= 1.0:
            pose2 = poses_2[-1]
        elif currU < 0.0:
            pose2 = poses_2[0]
        else:
            pose2 = poses_2[int(currU*100)]


        point1 = [pose1[0],pose1[1]]
        point2 = [pose2[0],pose2[1]]
        ang1 = pose1[2]
        ang2 = pose2[2]
    
        currPose = computeOffset(point1, point2, ang1, ang2 + currAng)

        # check for convergence condition, different between last and current cost is below threshold
        
        isTerminate = False
        if abs(lastCost - newCost) < costThresh or (numIterations - startIteration) > 10:
            isTerminate = True
            print "terminating globalOverlap_GPU2:", lastCost, newCost, costThresh, numIterations-startIteration
        
        #if abs(lastCost - newCost) < costThresh or (numIterations - startIteration) > 10:
        #    break



        " update after check for termination condition "
        lastCost = newCost

        numIterations += 1
        
        if isTerminate:
            break

    " draw final position "
    if False:
        
        " set the origin of pose 1 "
        poseOrigin = Pose(currPose)
        
        pylab.clf()
        pylab.axes()
        match_global = []
        
        for pair in match_pairs:
            p1 = pair[0]
            p2 = pair[1]
            
            p1_o = dispOffset(p1, currPose)
            
            p1_g = p1_o
            p2_g = p2
            match_global.append([p1_g,p2_g])

        draw_matches(match_global, [0.0,0.0,0.0])

        
        xP = []
        yP = []
        for b in globalPoly:
            p1 = b
            xP.append(p1[0])    
            yP.append(p1[1])

        pylab.plot(xP,yP,linewidth=1, color=(0.0,0.0,1.0))

        pylab.scatter([pose1[0]],[pose1[1]],color=(0.0,0.0,1.0))

        xP = []
        yP = []
        for b in points:
            p = [b[0],b[1]]
            p1 = dispOffset(p,currPose)
            
            #p1 = poseOrigin.convertLocalToGlobal(p)
            xP.append(p1[0])    
            yP.append(p1[1])

        pylab.plot(xP,yP,linewidth=1, color=(1.0,0.0,0.0))

        #pylab.scatter([pose2[0]],[pose2[1]],color=(1.0,0.0,0.0))


        plotEnv()        
        pylab.title("Cy: %u, u1 = %1.3f, u2 = %1.3f, ang = %1.3f, cost = %f" % (n1, u1, currU, currAng, newCost))

        pylab.xlim(currPose[0]-4, currPose[0]+4)                    
        pylab.ylim(currPose[1]-3, currPose[1]+3)
        pylab.savefig("ICP_plot_%04u.png" % globalPlotCount)
        pylab.clf()
        
        " save inputs "
        saveFile = ""
        saveFile += "initGuess = " + repr(initGuess) + "\n"
        saveFile += "globalPath = " + repr(globalPath) + "\n"
        #saveFile += "hull = " + repr(hull) + "\n"
        saveFile += "medialPoints = " + repr(medialPoints) + "\n"

        f = open("icpInputSave_%04u.txt" % globalPlotCount, 'w')
        globalPlotCount += 1
        f.write(saveFile)
        f.close()        
            
        #numIterations += 1


    return currPose, newCost


def globalOverlapTransformICP(initGuess, globalPath, hull, medialPoints,  plotIter = False, n1 = 0, n2 = 0):

    global numIterations
    
    def computeOffset(point1, point2, ang1, ang2):
    
        " corner points and orientations "
        corner1Pose = [point1[0], point1[1], ang1]
        corner2Pose = [point2[0], point2[1], ang2]
        
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

    u1 = initGuess[0]
    u2 = initGuess[1]

    uHigh = u2 + 0.2
    uLow = u2 - 0.2

    currU = u2
    currAng = initGuess[2]
    
    globalSpline = SplineFit(globalPath, smooth=0.1)
    medialSpline = SplineFit(medialPoints, smooth=0.1)

    if u1+0.02 > 1.0:
        pose1 = globalSpline.getUVecSet([0.98, 1.0])[0]
    elif u1 < 0.0:
        pose1 = globalSpline.getUVecSet([0.0, 0.02])[0]
    else:
        pose1 = globalSpline.getUVecSet([u1, u1+0.02])[0]

    if currU+0.02 > 1.0:
        pose2 = medialSpline.getUVecSet([0.98, 1.0])[0]
    elif currU < 0.0:
        pose2 = medialSpline.getUVecSet([0.0, 0.02])[0]
    else:
        pose2 = medialSpline.getUVecSet([currU, currU + 0.02])[0]
    
    point1 = [pose1[0],pose1[1]]
    point2 = [pose2[0],pose2[1]]
    ang1 = pose1[2]
    ang2 = pose2[2]

    gpacPoint1 = [pose1[0],pose1[1]]
    gpacPoint2 = [pose2[0],pose2[1]]
    
    gpacAng1 = pose1[2]
    gpacAng2 = pose2[2]

    currPose = computeOffset(gpacPoint1, gpacPoint2, gpacAng1, gpacAng2 + currAng)
    #currPose = computeOffset(point1, point2, ang1, ang2)
    #currPose[2] += currAng

    " transform the new pose "
    a_data_raw = hull
    a_data = []
    for p in a_data_raw:
        result = dispOffset(p, currPose)        
        a_data.append(result)    
    
    costThresh = 0.004
    #minMatchDist = 2.0
    minMatchDist = 1.0
    lastCost = 1e100
    
    startIteration = numIterations

    " sample a series of points from the medial curves "
    " augment points with point-to-line covariances "
    " treat the points with the point-to-line constraint "
    globalPoints = addGPACVectorCovariance(globalSpline.getUniformSamples(),high_var=0.05, low_var = 0.05)
    points = addGPACVectorCovariance(medialSpline.getUniformSamples(),high_var=0.05, low_var = 0.05)
    

    globalPoly = []        
    for p in globalPoints:
        globalPoly.append([p[0],p[1]])    


    " arc-angle transform points "
    gPoints = globalSpline.getTransformCurve()
    arcDist1 = globalSpline.dist(0.0, u1, iter = 0.005)                
    angOffset = -gPoints[0][1]
    gPointsOffset = globalSpline.getTransformCurve(angOffset)


    addPointToLineCovariance(gPoints, high_var=1.0, low_var=0.001)
    g1Covar = gPoints

    gPointHash = {}
    maxArcDist = gPoints[-1][0]
    minArcDist = 0.0

    " divide into boxes of 0.1 "
    for p in gPoints:
        thisDist = p[0]
        index = math.floor(thisDist*10)
        try:
            gPointHash[index].append(copy(p))
        except:
            gPointHash[index] = [copy(p)]

    print "gPointHash:", gPointHash

    medialPoints_trans = []
    for b in medialPoints:
        p = [b[0],b[1]]
        p = dispOffset(p,currPose)
        
        medialPoints_trans.append(p)

    medialSpline_trans = SplineFit(medialPoints_trans, smooth=0.1)
    tPoints2_trans = medialSpline_trans.getTransformCurve(compareAngle = gPoints[0][1])

    addPointToLineCovariance(tPoints2_trans, high_var=1.0, low_var=0.001)
    trans2Covar = tPoints2_trans

    arcDist2 = medialSpline_trans.dist(0.0, currU, iter = 0.005)    

    currAngOffset = 0.0
    currArcDist = arcDist2

    while True:

        """

        " transform the target Hull with the latest offset "
        points_offset = []
        for p in points:
            result = dispPoint(p, currPose)        
            points_offset.append(result)
        
        " transformed points without associated covariance "
        posePoly = []
        for p in points_offset:
            posePoly.append([p[0],p[1]])    
    
        " get the circles and radii "
        radius2, center2 = computeEnclosingCircle(points_offset)

        for i in range(len(globalPoints)):
            p_1 = globalPoly[i]
    
            if isInCircle(p_1, radius2, center2):
        
                " for every point of B, find it's closest neighbor in transformed A "
                p_2, i_2, minDist = findClosestPointInA(points_offset, p_1)
    
                if minDist <= minMatchDist:
        
                    " add to the list of match pairs less than 1.0 distance apart "
                    " keep A points and covariances untransformed "
                    C2 = points[i_2][2]
                    
                    C1 = globalPoints[i][2]
    
                    " we store the untransformed point, but the transformed covariance of the A point "
                    match_pairs.append([points[i_2],globalPoints[i],C2,C1])
        """
        #radius2, center2 = computeEnclosingCircle(points_offset)

        
        " find the matching pairs "
        match_pairs = []
                    
        tPointsTemp = []
        #for p in tPoints2_trans:
        for p in trans2Covar:
            tPointsTemp.append((p[0]+arcDist1-currArcDist, p[1] + currAngOffset, p[2]))


        #for i in range(len(gPoints)):
        #    p_1 = gPoints[i]
        #
        #    " for every point of B, find it's closest neighbor in transformed A "
        #    p_2, i_2, minDist = findClosestPointInA(tPointsTemp, p_1)
        #
        #    if minDist <= minMatchDist:
        #        " we store the untransformed point, but the transformed covariance of the A point "
        #        match_pairs.append([tPoints2_trans[i_2],gPoints[i]])

        for i in range(len(tPointsTemp)):
            p_1 = tPointsTemp[i]
        
            " for every point of B, find it's closest neighbor in transformed A "
            p_2, i_2, minDist = findClosestPointInA(gPoints, p_1)

            if minDist <= minMatchDist:
                " we store the untransformed point, but the transformed covariance of the A point "
                #match_pairs.append([tPoints2_trans[i],gPoints[i_2]])
                match_pairs.append([trans2Covar[i],gPoints[i_2]])




        if plotIter:

            fig1 = pylab.figure(2)
            fig1.clf()

            axes2 = fig1.add_subplot(122)

            #match_global = []
            
            #for pair in match_pairs:
            #    p1 = pair[0]
            #    p2 = pair[1]
                
            #    p1_o = dispOffset(p1, currPose)
                
            #    p1_g = p1_o
            #    p2_g = p2
            #    match_global.append([p1_g,p2_g])

            #draw_matches(match_global, [0.0,0.0,0.0], axes2)

            
            xP = []
            yP = []
            for b in globalPoly:
                p1 = b
                xP.append(p1[0])    
                yP.append(p1[1])
    
            axes2.plot(xP,yP,linewidth=1, color=(0.0,0.0,1.0))

            axes2.scatter([globalPoly[0][0]],[globalPoly[0][1]],color=(0.0,0.0,1.0))

            xP = []
            yP = []
            for b in points:
                p = [b[0],b[1]]
                p1 = dispOffset(p,currPose)
                
                xP.append(p1[0])    
                yP.append(p1[1])

            axes2.plot(xP,yP,linewidth=1, color=(1.0,0.0,0.0))

            u0 = points[0]
            u0L = dispOffset(u0,currPose)

            axes2.scatter([u0L[0]],[u0L[1]],color=(1.0,0.0,0.0))

            xP = []
            yP = []
            for b in a_data_raw:
                p = [b[0],b[1]]
                p1 = dispOffset(p,currPose)
                
                #p1 = poseOrigin.convertLocalToGlobal(p)
                xP.append(p1[0])    
                yP.append(p1[1])

            axes2.plot(xP,yP,linewidth=1, color=(1.0,0.0,0.0))
            
            #cost = medialOverlapCostFunc([currU, currAng], match_pairs, globalSpline, medialSpline, uHigh, uLow, u1)            
            cost = 0.0

            plotEnv(axes2)                    
            axes2.set_title("%u, u1 = %1.3f, u2 = %1.3f, ang = %1.3f, cost = %f" % (n1, u1, currU, currAng, cost))

            axes2.set_xlim(currPose[0]-4, currPose[0]+4)                    
            axes2.set_ylim(currPose[1]-3, currPose[1]+3)
        

            axes1 = fig1.add_subplot(121)
            axes1.grid(True)
            
            #gPoints = globalSpline.getTransformCurve()
            #angOffset = -gPoints[0][1]
            #gPointsOffset = globalSpline.getTransformCurve(angOffset)

            #dist1 = globalSpline.dist(0.0, u1, iter = 0.005)                

            xP = []
            yP = []
            for p in gPointsOffset:
                xP.append(p[0])
                yP.append(p[1])
            axes1.plot(xP,yP, color='b')

            #medialPoints_trans = []
            #for b in medialPoints:
            #    p = [b[0],b[1]]
            #    p = dispOffset(p,currPose)
            #    
            #    medialPoints_trans.append(p)

            medialSpline_trans = SplineFit(medialPoints_trans, smooth=0.1)
            tPointsTemp = medialSpline_trans.getTransformCurve(offset = angOffset, compareAngle = gPointsOffset[0][1])

            #dist2 = medialSpline_trans.dist(0.0, currU, iter = 0.005)    

            xP = []
            yP = []
            for p in tPointsTemp:
                xP.append(p[0]+arcDist1-currArcDist)
                yP.append(p[1]+currAngOffset)
            axes1.plot(xP,yP, color='r')


            #draw_matches(match_pairs, [0.0,0.0,0.0], axes1)
            draw_matches(match_pairs, [arcDist1-currArcDist,currAngOffset,0.0], axes1)

            thisCost = medialOverlapTransformCostFunc([currArcDist,currAngOffset], match_pairs, arcDist1)
            
            axes1.set_title("dist,ang,cost = %1.3f,%1.3f,%1.3f" % (currArcDist, currAngOffset, thisCost))
            axes1.set_xlabel("distance")
            axes1.set_ylabel("angle (radians)")
            axes1.set_xlim(-2,10)
            axes1.set_ylim(-6,6)        

            fig1.savefig("ICP_plot_%04u_0.png" % numIterations)
            pylab.clf()            
            
            pylab.figure(1)            
            
            numIterations += 1

        """
        gPoints
        arcDist1
        gPointHash
        tPoints2_trans
        arcDist2
        currAngOffset
        currArcDist
        """

        #newParam = scipy.optimize.fmin(medialOverlapCostFunc, [currU,currAng], [match_pairs, globalSpline, medialSpline, uHigh, uLow, u1], disp = 0)

        newParam = scipy.optimize.fmin(medialOverlapTransformCostFunc, [currArcDist,currAngOffset], [match_pairs, arcDist1], disp = 0)

        newArcDist = newParam[0]
        newAngOffset = newParam[1]
        
        newCost = medialOverlapTransformCostFunc([currArcDist,currAngOffset], match_pairs, arcDist1)

        currArcDist = newArcDist
        currAngOffset = newAngOffset


        dist1 = globalSpline.dist(0.0, u1, iter = 0.005)                
        dist2 = medialSpline_trans.dist(0.0, currU, iter = 0.005)    
        xP.append(p[0]+dist1-dist2)

        currU = medialSpline.getUOfDist(0.0, currArcDist, distIter = 0.001)

        currAng = initGuess[2] + currAngOffset 
        
        #if u1+0.02 > 1.0:
        #    pose1 = globalSpline.getUVecSet([0.98, 1.0])[0]
        #elif u1 < 0.0:
        #    pose1 = globalSpline.getUVecSet([0.0, 0.02])[0]
        #else:
        #    pose1 = globalSpline.getUVecSet([u1, u1+0.02])[0]
    
        if currU+0.02 > 1.0:
            pose2 = medialSpline.getUVecSet([0.98, 1.0])[0]
        elif currU < 0.0:
            pose2 = medialSpline.getUVecSet([0.0, 0.02])[0]
        else:
            pose2 = medialSpline.getUVecSet([currU, currU + 0.02])[0]

        #gpacPoint1 = [pose1[0],pose1[1]]
        gpacPoint2 = [pose2[0],pose2[1]]
        
        #gpacAng1 = pose1[2]
        #gpacAng2 = pose2[2]
    
        currPose = computeOffset(gpacPoint1, gpacPoint2, gpacAng1, gpacAng2 + currAng)
            


        #point1 = [pose1[0],pose1[1]]
        #point2 = [pose2[0],pose2[1]]
        #ang1 = pose1[2]
        #ang2 = pose2[2]    
        #currPose = computeOffset(point1, point2, ang1, ang2)
        #currPose[2] += currAng
        
        " check for convergence condition, different between last and current cost is below threshold "
        
        isTerminate = False
        if abs(lastCost - newCost) < costThresh or (numIterations - startIteration) > 10:
            isTerminate = True
            print "terminating globalOverlapTransform:", lastCost, newCost, costThresh, numIterations-startIteration
        
        " update after check for termination condition "
        lastCost = newCost

        " optionally draw the position of the points in current transform "
        if plotIter:


            fig1 = pylab.figure(2)
            fig1.clf()

            axes2 = fig1.add_subplot(122)

            match_global = []


                        
            #for pair in match_pairs:
            #    p1 = pair[0]
            #    p2 = pair[1]
                
            #    p1_o = dispOffset(p1, currPose)
                
            #    p1_g = p1_o
            #    p2_g = p2
            #    match_global.append([p1_g,p2_g])

            #draw_matches(match_global, [0.0,0.0,0.0], axes2)

            
            xP = []
            yP = []
            for b in globalPoly:
                p1 = b
                xP.append(p1[0])    
                yP.append(p1[1])
    
            axes2.plot(xP,yP,linewidth=1, color=(0.0,0.0,1.0))

            axes2.scatter([globalPoly[0][0]],[globalPoly[0][1]],color=(0.0,0.0,1.0))

            xP = []
            yP = []
            for b in points:
                p = [b[0],b[1]]
                p1 = dispOffset(p,currPose)
                
                xP.append(p1[0])    
                yP.append(p1[1])

            axes2.plot(xP,yP,linewidth=1, color=(1.0,0.0,0.0))

            u0 = points[0]
            u0L = dispOffset(u0,currPose)

            axes2.scatter([u0L[0]],[u0L[1]],color=(1.0,0.0,0.0))

            xP = []
            yP = []
            for b in a_data_raw:
                p = [b[0],b[1]]
                p1 = dispOffset(p,currPose)
                
                #p1 = poseOrigin.convertLocalToGlobal(p)
                xP.append(p1[0])    
                yP.append(p1[1])

            axes2.plot(xP,yP,linewidth=1, color=(1.0,0.0,0.0))
            
            cost = 0.0
            #cost = medialOverlapCostFunc([currU, currAng], match_pairs, globalSpline, medialSpline, uHigh, uLow, u1)            

            plotEnv(axes2)                    
            axes2.set_title("%u, u1 = %1.3f, u2 = %1.3f, ang = %1.3f, cost = %f" % (n1, u1, currU, currAng, cost))

            axes2.set_xlim(currPose[0]-4, currPose[0]+4)                    
            axes2.set_ylim(currPose[1]-3, currPose[1]+3)
        

            axes1 = fig1.add_subplot(121)
            axes1.grid(True)
            
            #gPoints = globalSpline.getTransformCurve()
            #angOffset = -gPoints[0][1]
            #gPointsOffset = globalSpline.getTransformCurve(angOffset)

            #dist1 = globalSpline.dist(0.0, u1, iter = 0.005)                

            xP = []
            yP = []
            for p in gPointsOffset:
                xP.append(p[0])
                yP.append(p[1])
            axes1.plot(xP,yP, color='b')

            #medialPoints_trans = []
            #for b in medialPoints:
            #    p = [b[0],b[1]]
            #    p = dispOffset(p,currPose)
            #    
            #    medialPoints_trans.append(p)

            medialSpline_trans = SplineFit(medialPoints_trans, smooth=0.1)
            tPointsTemp = medialSpline_trans.getTransformCurve(offset = angOffset, compareAngle = gPointsOffset[0][1])






            #dist2 = medialSpline_trans.dist(0.0, currU, iter = 0.005)    

            xP = []
            yP = []
            for p in tPointsTemp:
                xP.append(p[0]+arcDist1-currArcDist)
                yP.append(p[1]+currAngOffset)
            axes1.plot(xP,yP, color='r')

            draw_matches(match_pairs, [arcDist1-currArcDist,currAngOffset,0.0], axes1)

            
            #axes1.set_title("Transform ICP: arcDist = %1.3f, angOffset = %1.3f, cost = %1.3f" % (currArcDist, currAngOffset, newCost))
            axes1.set_title("dist,ang,cost = %1.3f,%1.3f,%1.3f" % (currArcDist, currAngOffset, thisCost))
            axes1.set_xlabel("distance")
            axes1.set_ylabel("angle (radians)")
            axes1.set_xlim(-2,10)
            axes1.set_ylim(-6,6)        

            fig1.savefig("ICP_plot_%04u_0.png" % numIterations)
            pylab.clf()            
            
            pylab.figure(1)            
            
        
        numIterations += 1
        
        if isTerminate:
            break

    " draw final position "
    if False:
        
        " set the origin of pose 1 "
        poseOrigin = Pose(currPose)
        
        pylab.clf()
        pylab.axes()
        match_global = []
        
        for pair in match_pairs:
            p1 = pair[0]
            p2 = pair[1]
            
            p1_o = dispOffset(p1, currPose)
            
            p1_g = p1_o
            p2_g = p2
            match_global.append([p1_g,p2_g])

        draw_matches(match_global, [0.0,0.0,0.0])

        
        xP = []
        yP = []
        for b in globalPoly:
            p1 = b
            xP.append(p1[0])    
            yP.append(p1[1])

        pylab.plot(xP,yP,linewidth=1, color=(0.0,0.0,1.0))

        pylab.scatter([globalPoly[0][0]],[globalPoly[0][1]],color=(0.0,0.0,1.0))

        xP = []
        yP = []
        for b in points:
            p = [b[0],b[1]]
            p1 = dispOffset(p,currPose)
            
            #p1 = poseOrigin.convertLocalToGlobal(p)
            xP.append(p1[0])    
            yP.append(p1[1])

        pylab.plot(xP,yP,linewidth=1, color=(1.0,0.0,0.0))

        u0 = points[0]
        u0L = dispOffset(u0,currPose)

        pylab.scatter([u0L[0]],[u0L[1]],color=(1.0,0.0,0.0))

        xP = []
        yP = []
        for b in a_data_raw:
            p = [b[0],b[1]]
            p1 = dispOffset(p,currPose)
            
            #p1 = poseOrigin.convertLocalToGlobal(p)
            xP.append(p1[0])    
            yP.append(p1[1])

        pylab.plot(xP,yP,linewidth=1, color=(1.0,0.0,0.0))

        plotEnv()        
        pylab.title("%u, u1 = %1.3f, u2 = %1.3f, ang = %1.3f, cost = %f" % (n1, u1, currU, currAng, newCost))

        pylab.xlim(currPose[0]-4, currPose[0]+4)                    
        pylab.ylim(currPose[1]-3, currPose[1]+3)
        pylab.savefig("ICP_plot_%04u_1.png" % numIterations)
        pylab.clf()            
        numIterations += 1
    
    print "returning Transform", currPose
    return currPose, newCost

def overlapICP2(estPose1, initGuess, hull1, hull2, medialPoints1, medialPoints2, supportLine, rootPose1, rootPose2, inPlace = False, plotIter = False, n1 = 0, n2 = 0):

    global numIterations

    def computeOffset(point1, point2, ang1, ang2):
    
        " corner points and orientations "
        corner1Pose = [point1[0], point1[1], ang1]
        corner2Pose = [point2[0], point2[1], ang2]
        
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

    u1 = initGuess[0]
    u2 = initGuess[1]
    
    #uHigh = u2 + 1.0
    #uLow = u2 - 1.0
    #uHigh = u2 + 0.08
    #uLow = u2 - 0.08
    uHigh = u2 + 0.2
    uLow = u2 - 0.2

    if inPlace:
        uHigh = u2 + 0.08
        uLow = u2 - 0.08
    
    currU = u2
    currAng = initGuess[2]
    medialSpline1 = SplineFit(medialPoints1, smooth=0.1)
    medialSpline2 = SplineFit(medialPoints2, smooth=0.1)


    ##pose1 = medialSpline1.getUVecSet([u1, u1+0.02])[0]
    #pose2 = medialSpline2.getUVecSet([currU, currU + 0.02])[0]
    if u1+0.02 > 1.0:
        pose1 = medialSpline1.getUVecSet([0.98, 1.0])[0]
    elif u1 < 0.0:
        pose1 = medialSpline1.getUVecSet([0.0, 0.02])[0]
    else:
        pose1 = medialSpline1.getUVecSet([u1, u1+0.02])[0]

    if currU+0.02 > 1.0:
        pose2 = medialSpline2.getUVecSet([0.98, 1.0])[0]
    elif currU < 0.0:
        pose2 = medialSpline2.getUVecSet([0.0, 0.02])[0]
    else:
        pose2 = medialSpline2.getUVecSet([currU, currU + 0.02])[0]
        
    point1 = [pose1[0],pose1[1]]
    point2 = [pose2[0],pose2[1]]
    ang1 = pose1[2]
    ang2 = pose2[2]

    offset = computeOffset(point1, point2, ang1, ang2 + currAng)

    if inPlace:
        offset[2] =  functions.normalizeAngle(offset[2])    
        return offset, (0,0,0)


    " transform the past poses "
    a_data_raw = hull2
    a_data = []
    for p in a_data_raw:
        result = dispOffset(p, offset)        
        a_data.append(result)    

    polyB = []        
    for p in hull1:
        polyB.append([p[0],p[1]])    
    
    costThresh = 0.004
    minMatchDist = 2.0
    #minMatchDist2 = 0.5
    minMatchDist2 = 2.0
    lastCost = 1e100
    
    startIteration = numIterations

    " set the initial guess "
    poseOrigin = Pose(estPose1)
    
    localSupport = []
    for pnt in supportLine:
        localSupport.append(poseOrigin.convertGlobalToLocal(pnt))

    supportSpline = SplineFit(localSupport, smooth=0.1)
    
    
    " sample a series of points from the medial curves "
    " augment points with point-to-line covariances "
    " treat the points with the point-to-line constraint "
    #points1 = addGPACVectorCovariance(medialSpline1.getUniformSamples(),high_var=0.1, low_var = 0.05)
    #points2 = addGPACVectorCovariance(medialSpline2.getUniformSamples(),high_var=0.1, low_var = 0.05)
    points1 = addGPACVectorCovariance(medialSpline1.getUniformSamples(),high_var=0.05, low_var = 0.001)
    points2 = addGPACVectorCovariance(medialSpline2.getUniformSamples(),high_var=0.05, low_var = 0.001)
    #points1 = addGPACVectorCovariance(medialSpline1.getUniformSamples(),high_var=0.1, low_var = 0.0005)
    #points2 = addGPACVectorCovariance(medialSpline2.getUniformSamples(),high_var=0.1, low_var = 0.0005)
    #points1 = addGPACVectorCovarianceWithAngle(medialSpline1.getUniformSamples(),high_var=0.1, low_var = 0.001)
    #points2 = addGPACVectorCovarianceWithAngle(medialSpline2.getUniformSamples(),high_var=0.1, low_var = 0.001)
    supportPoints = addGPACVectorCovariance(supportSpline.getUniformSamples(),high_var=0.05, low_var = 0.001)



    " transform pose 2 by initial offset guess "    
    " transform the new pose "
    points2_offset = []
    for p in points2:
        result = dispPoint(p, offset)        
        #result = dispPointWithAngle(p, offset)
        points2_offset.append(result)


    
    poly1 = []        
    for p in points1:
        poly1.append([p[0],p[1]])    

    while True:
        
        " find the matching pairs "
        match_pairs = []

        " transform the target Hull with the latest offset "
        points2_offset = []
        for p in points2:
            result = dispPoint(p, offset)        
            #result = dispPointWithAngle(p, offset)
            points2_offset.append(result)
        
        " transformed points without associated covariance "
        poly2 = []
        for p in points2_offset:
            poly2.append([p[0],p[1]])    
        
        " get the circles and radii "
        radius2, center2 = computeEnclosingCircle(points2_offset)
        radius1, center1 = computeEnclosingCircle(points1)
    
        for i in range(len(points2_offset)):
            p_2 = poly2[i]

            if isInCircle(p_2, radius1, center1):

                " for every transformed point of A, find it's closest neighbor in B "
                p_1, minDist = findClosestPointInB(points1, p_2, [0.0,0.0,0.0])
    
                if minDist <= minMatchDist:
        
                    " add to the list of match pairs less than 1.0 distance apart "
                    " keep A points and covariances untransformed "
                    C2 = points2[i][2]
                    C1 = p_1[2]
    
                    " we store the untransformed point, but the transformed covariance of the A point "
                    match_pairs.append([points2[i],p_1,C2,C1])

        for i in range(len(points1)):
            p_1 = poly1[i]
    
            if isInCircle(p_1, radius2, center2):
        
                #print "selected", b_p, "in circle", radiusA, centerA, "with distance"
                " for every point of B, find it's closest neighbor in transformed A "
                p_2, i_2, minDist = findClosestPointInA(points2_offset, p_1)
    
                if minDist <= minMatchDist:
        
                    " add to the list of match pairs less than 1.0 distance apart "
                    " keep A points and covariances untransformed "
                    C2 = points2[i_2][2]
                    
                    C1 = points1[i][2]
    
                    " we store the untransformed point, but the transformed covariance of the A point "
                    #match_pairs.append([points2[i_2],p_1,C2,C1])
                    match_pairs.append([points2[i_2],points1[i],C2,C1])

        
        support_pairs = []
        for i in range(len(poly2)):

            
            p_2 = poly2[i]


            " for every transformed point of A, find it's closest neighbor in B "
            #p_1, minDist = findClosestPointInB(localSupport, p_2, [0.0,0.0,0.0])
            p_1, minDist = findClosestPointInB(supportPoints, p_2, [0.0,0.0,0.0])

            if isInCircle(p_1, radius2, center2):

                if minDist <= minMatchDist2:
                    C2 = points2[i][2]
                    C1 = p_1[2]

                    " we store the untransformed point, but the transformed covariance of the A point "
                    support_pairs.append([points2[i],p_1,C2,C1])

        #if plotIter and startIteration == numIterations:
        if plotIter:
            
            " set the origin of pose 1 "
            poseOrigin = Pose(estPose1)
    
            pylab.clf()
            pylab.axes()
            match_global = []
            
            #for pair in match_pairs:
            #    p1 = pair[0]
            #    p2 = pair[1]
            
            #    p1_o = dispOffset(p1, offset)
            
            #    p1_g = poseOrigin.convertLocalToGlobal(p1_o)
            #    p2_g = poseOrigin.convertLocalToGlobal(p2)
            #    match_global.append([p1_g,p2_g])
            
            #draw_matches(match_global, [0.0,0.0,0.0])

            #localSupport = []
            #for pnt in supportLine:
            #    localSupport.append(poseOrigin.convertGlobalToLocal(pnt))

            for pair in support_pairs:
                p1 = pair[0]
                p2 = pair[1]
                
                p1_o = dispOffset(p1, offset)
                
                p1_g = poseOrigin.convertLocalToGlobal(p1_o)
                p2_g = poseOrigin.convertLocalToGlobal(p2)
                pylab.plot([p1_g[0],p2_g[0]], [p1_g[1],p2_g[1]])
                #match_global.append([p1_g,p2_g])

            #draw_matches(match_global, [0.0,0.0,0.0])

            #for pair in match_pairs:
            #    a_off = dispOffset(pair[0], offset)
            #    b = pair[1]
        
            #    pylab.plot([a_off[0],b[0]], [a_off[1],b[1]])



            xP = []
            yP = []
            for b in supportPoints:
                p1 = poseOrigin.convertLocalToGlobal(b)
    
                xP.append(p1[0])    
                yP.append(p1[1])
    
            pylab.plot(xP,yP,linewidth=1, color=(0.0,0.0,0.5))
            
            xP = []
            yP = []
            for b in poly1:
                p1 = poseOrigin.convertLocalToGlobal(b)
    
                xP.append(p1[0])    
                yP.append(p1[1])
    
            pylab.plot(xP,yP,linewidth=1, color=(0.0,0.0,1.0))

            #center1Global = poseOrigin.convertLocalToGlobal(center1)
            #cir = pylab.Circle((center1Global[0],center1Global[1]), radius=radius1,  ec='b', fc='none')
            #pylab.gca().add_patch(cir)

    
            xP = []
            yP = []
            for b in points2:
                p = [b[0],b[1]]
                p = dispOffset(p,offset)
                
                p1 = poseOrigin.convertLocalToGlobal(p)
                xP.append(p1[0])    
                yP.append(p1[1])

            pylab.plot(xP,yP,linewidth=1, color=(1.0,0.0,0.0))

            #center2Global = poseOrigin.convertLocalToGlobal(center2)
            #cir = pylab.Circle((center2Global[0],center2Global[1]), radius=radius2,  ec='r', fc='none')
            #pylab.gca().add_patch(cir)


            point2_trans = dispOffset(point2, offset)
            p1 = poseOrigin.convertLocalToGlobal(point1)
            p2 = poseOrigin.convertLocalToGlobal(point2_trans)

            xP = [p1[0],p2[0]]
            yP = [p1[1],p2[1]]
            pylab.scatter(xP,yP,color='b')


            gpac2_trans = dispOffset(rootPose1, offset)
            gpac1 = poseOrigin.convertLocalToGlobal(rootPose2)
            gpac2 = poseOrigin.convertLocalToGlobal(gpac2_trans)
            #gpac2_trans = dispOffset([0.0,0.0], offset)
            #gpac1 = poseOrigin.convertLocalToGlobal([0.0,0.0])
            #gpac2 = poseOrigin.convertLocalToGlobal(gpac2_trans)

            xP = [gpac1[0],gpac2[0]]
            yP = [gpac1[1],gpac2[1]]
            pylab.scatter(xP,yP,color='r')


            xP = []
            yP = []
            for b in polyB:
                p1 = poseOrigin.convertLocalToGlobal(b)

                xP.append(p1[0])    
                yP.append(p1[1])
            
            pylab.plot(xP,yP,linewidth=1, color=(0.0,0.0,1.0))


            xP = []
            yP = []
            for b in a_data_raw:
                p = [b[0],b[1]]
                p = dispOffset(p,offset)
                
                p1 = poseOrigin.convertLocalToGlobal(p)
                xP.append(p1[0])    
                yP.append(p1[1])

            pylab.plot(xP,yP,linewidth=1, color=(1.0,0.0,0.0))

            """
            xP = []
            yP = []
            for b in medialPoints1:
                p1 = poseOrigin.convertLocalToGlobal(b)

                xP.append(p1[0])    
                yP.append(p1[1])

            pylab.scatter(xP,yP,color=(0.0,0.0,1.0))

            xP = []
            yP = []
            for b in medialPoints2:
                p = [b[0],b[1]]
                p = dispOffset(p,offset)
                
                p1 = poseOrigin.convertLocalToGlobal(p)
                xP.append(p1[0])    
                yP.append(p1[1])
                
            pylab.scatter(xP,yP,color=(1.0,0.0,0.0))
            """
            plotEnv()        
            
            #pylab.title("%u %u, u = %f" % (n1, n2, currU))
            lowCount, midCount, highCount = overlapHistogram([currU,currAng], match_pairs, medialSpline1, medialSpline2, u1)            
            pylab.title("%u %u, u1 = %1.3f, u2 = %1.3f, ang = %1.3f, hist = %d %d %d" % (n1, n2, u1, currU, currAng, lowCount, midCount, highCount))
            pylab.xlim(estPose1[0]-4, estPose1[0]+4)                    
            pylab.ylim(estPose1[1]-3, estPose1[1]+3)
            #pylab.xlim(-3,3)
            #pylab.ylim(-3,3)
            #pylab.axis('equal')
            pylab.savefig("ICP_plot_%04u_0.png" % numIterations)
            pylab.clf()            
        

        #if inPlace:
        #    newParam = scipy.optimize.fmin(medialOverlapInPlaceCostFunc, [currAng], [match_pairs, medialSpline1, medialSpline2, currU, u1], disp = 0)
        #    newU = currU
        #    newAng = newParam[0]
        #else:
        #    pass
        #newParam = scipy.optimize.fmin(medialOverlapCostFunc, [currU,currAng], [match_pairs, medialSpline1, medialSpline2, uHigh, uLow, u1], disp = 0)
        newParam = scipy.optimize.fmin(medialOverlapCostFunc, [currU,currAng], [support_pairs, medialSpline1, medialSpline2, uHigh, uLow, u1], disp = 0)

        newU = newParam[0]
        newAng = newParam[1]
        
        # get the current cost
        #if inPlace:
        #    newCost = medialOverlapInPlaceCostFunc([newAng], match_pairs, medialSpline1, medialSpline2, newU, u1)
        #else:
        newCost = medialOverlapCostFunc([newU, newAng], match_pairs, medialSpline1, medialSpline2, uHigh, uLow, u1)
        
        " set the current parameters "
        currAng = functions.normalizeAngle(newAng)
        currU = newU
        
        #print "iteration", numIterations
        if inPlace:
            print "currU =", currU, "currAng =", currAng, "newCost =", newCost

        " compute offset from newU and newAng"
        #pose1 = medialSpline1.getUVecSet([u1, u1+0.02])[0]
        #pose2 = medialSpline2.getUVecSet([currU, currU+0.02])[0]
        if u1+0.02 > 1.0:
            pose1 = medialSpline1.getUVecSet([0.98, 1.0])[0]
        elif u1 < 0.0:
            pose1 = medialSpline1.getUVecSet([0.0, 0.02])[0]
        else:
            pose1 = medialSpline1.getUVecSet([u1, u1+0.02])[0]
        
        if currU+0.02 > 1.0:
            pose2 = medialSpline2.getUVecSet([0.98, 1.0])[0]
        elif currU < 0.0:
            pose2 = medialSpline2.getUVecSet([0.0, 0.02])[0]
        else:
            pose2 = medialSpline2.getUVecSet([currU, currU + 0.02])[0]
                
        point1 = [pose1[0],pose1[1]]
        point2 = [pose2[0],pose2[1]]
        ang1 = pose1[2]
        ang2 = pose2[2]
    
        newOffset = computeOffset(point1, point2, ang1, ang2 + currAng)

        # save the current offset and cost
        offset = newOffset

    
        # check for convergence condition, different between last and current cost is below threshold
        if abs(lastCost - newCost) < costThresh or (numIterations - startIteration) > 10:
            #print "breaking:", abs(lastCost - newCost) < costThresh, (numIterations - startIteration) > 10
            #print "lastCost =", lastCost, "newCost =", newCost, "costThresh =", costThresh
            break



        " update after check for termination condition "
        lastCost = newCost

        # optionally draw the position of the points in current transform
        if plotIter:
            
            " set the origin of pose 1 "
            poseOrigin = Pose(estPose1)
    
            pylab.clf()
            pylab.axes()
            match_global = []
            
            #for pair in match_pairs:
            #    p1 = pair[0]
            #    p2 = pair[1]
            #    
            #    p1_o = dispOffset(p1, offset)
                
            #    p1_g = poseOrigin.convertLocalToGlobal(p1_o)
            #    p2_g = poseOrigin.convertLocalToGlobal(p2)
            #    match_global.append([p1_g,p2_g])
            
            #draw_matches(match_global, [0.0,0.0,0.0])

            for pair in support_pairs:
                p1 = pair[0]
                p2 = pair[1]
                
                p1_o = dispOffset(p1, offset)
                
                p1_g = poseOrigin.convertLocalToGlobal(p1_o)
                p2_g = poseOrigin.convertLocalToGlobal(p2)
                pylab.plot([p1_g[0],p2_g[0]], [p1_g[1],p2_g[1]])

            xP = []
            yP = []
            for b in supportPoints:
                p1 = poseOrigin.convertLocalToGlobal(b)
    
                xP.append(p1[0])    
                yP.append(p1[1])
    
            pylab.plot(xP,yP,linewidth=1, color=(0.0,0.0,0.5))
                
            xP = []
            yP = []
            for b in poly1:
                p1 = poseOrigin.convertLocalToGlobal(b)

                xP.append(p1[0])    
                yP.append(p1[1])

            pylab.plot(xP,yP,linewidth=1, color=(0.0,0.0,1.0))

            #center1Global = poseOrigin.convertLocalToGlobal(center1)
            #cir = pylab.Circle((center1Global[0],center1Global[1]), radius=radius1,  ec='b', fc='none')
            #pylab.gca().add_patch(cir)

            xP = []
            yP = []
            for b in points2:
                p = [b[0],b[1]]
                p = dispOffset(p,offset)
                
                p1 = poseOrigin.convertLocalToGlobal(p)
                xP.append(p1[0])    
                yP.append(p1[1])
            
            pylab.plot(xP,yP,linewidth=1, color=(1.0,0.0,0.0))

            #center2Global = poseOrigin.convertLocalToGlobal(center2)
            #cir = pylab.Circle((center2Global[0],center2Global[1]), radius=radius2,  ec='r', fc='none')
            #pylab.gca().add_patch(cir)

            point2_trans = dispOffset(point2, offset)
            p1 = poseOrigin.convertLocalToGlobal(point1)
            p2 = poseOrigin.convertLocalToGlobal(point2_trans)

            #xP = [point1[0],point2_trans[0]]
            #yP = [point1[1],point2_trans[1]]
            xP = [p1[0],p2[0]]
            yP = [p1[1],p2[1]]

            pylab.scatter(xP,yP,color='b')

            gpac2_trans = dispOffset(rootPose1, offset)
            gpac1 = poseOrigin.convertLocalToGlobal(rootPose2)
            gpac2 = poseOrigin.convertLocalToGlobal(gpac2_trans)
            #gpac2_trans = dispOffset([0.0,0.0], offset)
            #gpac1 = poseOrigin.convertLocalToGlobal([0.0,0.0])
            #gpac2 = poseOrigin.convertLocalToGlobal(gpac2_trans)

            xP = [gpac1[0],gpac2[0]]
            yP = [gpac1[1],gpac2[1]]
            pylab.scatter(xP,yP,color='r')

            xP = []
            yP = []
            for b in polyB:
                p1 = poseOrigin.convertLocalToGlobal(b)

                xP.append(p1[0])    
                yP.append(p1[1])
            
            pylab.plot(xP,yP,linewidth=1, color=(0.0,0.0,1.0))

            xP = []
            yP = []
            for b in a_data_raw:
                p = [b[0],b[1]]
                p = dispOffset(p,offset)
                
                p1 = poseOrigin.convertLocalToGlobal(p)
                xP.append(p1[0])    
                yP.append(p1[1])
            
            pylab.plot(xP,yP,linewidth=1, color=(1.0,0.0,0.0))
            """
            xP = []
            yP = []
            for b in medialPoints1:
                p1 = poseOrigin.convertLocalToGlobal(b)

                xP.append(p1[0])    
                yP.append(p1[1])

            pylab.scatter(xP,yP,color=(0.0,0.0,1.0))

            xP = []
            yP = []
            for b in medialPoints2:
                p = [b[0],b[1]]
                p = dispOffset(p,offset)
                
                p1 = poseOrigin.convertLocalToGlobal(p)
                xP.append(p1[0])    
                yP.append(p1[1])
                
            pylab.scatter(xP,yP,color=(1.0,0.0,0.0))
            """
            plotEnv()        
            
            lowCount, midCount, highCount = overlapHistogram([currU,currAng], match_pairs, medialSpline1, medialSpline2, u1)            
            #pylab.title("%u %u, u1 = %1.3f, u2 = %1.3f, hist = %d %d %d" % (n1, n2, u1, currU, lowCount, midCount, highCount))
            pylab.title("%u %u, u1 = %1.3f, u2 = %1.3f, ang = %1.3f, hist = %d %d %d" % (n1, n2, u1, currU, currAng, lowCount, midCount, highCount))
            #pylab.title("%u %u, %u, u1 = %1.3f, u2 = %1.3f" % (n1, n2, len(match_pairs), u1, currU))
            pylab.xlim(estPose1[0]-4, estPose1[0]+4)                    
            pylab.ylim(estPose1[1]-3, estPose1[1]+3)
            #pylab.axis('equal')
            pylab.savefig("ICP_plot_%04u_1.png" % numIterations)
            pylab.clf()            

        
        numIterations += 1

        " reduce the minMatch distance for each step down to a floor value "
        #minMatchDist /= 2
        #if minMatchDist < 0.25:
        #    minMatchDist = 0.25
    
        
    histogram = overlapHistogram([currU,currAng], match_pairs, medialSpline1, medialSpline2, u1)            


    offset[2] =  functions.normalizeAngle(offset[2])    
    return offset, histogram




def cornerICP(estPose1, angGuess, point1, point2, ang1, ang2, hull1, hull2, sweepHull1, sweepHull2, pastCircles, costThresh = 0.004, minMatchDist = 2.0, plotIter = False, n1 = 0, n2 = 0):

    global numIterations

    def computeOffset(point1, point2, ang1, ang2):
    
        " corner points and orientations "
        corner1Pose = [point1[0], point1[1], ang1]
        corner2Pose = [point2[0], point2[1], ang2]
        
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

    currAng = angGuess 
    offset = computeOffset(point1, point2, ang1, ang2 + currAng)
    
    lastCost = 1e100
    
    startIteration = numIterations

    " set the initial guess "
    poseOrigin = Pose(estPose1)
    
    " transform the past poses "
    a_data_raw = hull2
    a_data = []
    for p in a_data_raw:
        result = dispPoint(p, offset)        
        a_data.append(result)    

    polyB = []        
    for p in hull1:
        polyB.append([p[0],p[1]])    

    while True:
        " find the matching pairs "
        match_pairs = []

        a_data_raw = hull2
        b_data = hull1            

        " transform the target Hull with the latest offset "
        a_data = []
        for p in a_data_raw:
            result = dispPoint(p, offset)
            a_data.append(result)

        " transformed points without associated covariance "
        polyA = []
        for p in a_data:
            polyA.append([p[0],p[1]])    
        
        " get the circles and radii "
        radiusA, centerA = computeEnclosingCircle(a_data)
        #radiusB, centerB = computeEnclosingCircle(pastHull)
        
        if True:
            for i in range(len(a_data)):
                a_p = polyA[i]
    
                #if isValid(a_p, radiusB, centerB, polyB):
                if isValidPast(a_p, pastCircles, polyB):
    
                    " for every transformed point of A, find it's closest neighbor in B "
                    b_p, minDist = findClosestPointInB(b_data, a_p, [0.0,0.0,0.0])
        
                    if minDist <= minMatchDist:
            
                        " add to the list of match pairs less than 1.0 distance apart "
                        " keep A points and covariances untransformed "
                        Ca = a_data_raw[i][2]
                        Cb = b_p[2]
        
                        " we store the untransformed point, but the transformed covariance of the A point "
                        match_pairs.append([a_data_raw[i],b_p,Ca,Cb])

        if True:
            for i in range(len(b_data)):
                b_p = polyB[i]
        
                if isValid(b_p, radiusA, centerA, polyA):
            
                    #print "selected", b_p, "in circle", radiusA, centerA, "with distance"
                    " for every point of B, find it's closest neighbor in transformed A "
                    a_p, a_i, minDist = findClosestPointInA(a_data, b_p)
        
                    if minDist <= minMatchDist:
            
                        " add to the list of match pairs less than 1.0 distance apart "
                        " keep A points and covariances untransformed "
                        Ca = a_data_raw[a_i][2]
                        
                        Cb = b_data[i][2]
        
                        " we store the untransformed point, but the transformed covariance of the A point "
                        match_pairs.append([a_data_raw[a_i],b_p,Ca,Cb])


        " plot the initial configuration "
        #if plotIter and startIteration == numIterations:
        if False and startIteration == numIterations:

            numIterations += 1

            pylab.clf()
            pylab.axes()
            match_global = []
            
            for pair in match_pairs:
                p1 = pair[0]
                p2 = pair[1]
                
                p1_o = dispOffset(p1, offset)
                #p2_o = dispOffset(p2, offset)

                
                p1_g = poseOrigin.convertLocalToGlobal(p1_o)
                p2_g = poseOrigin.convertLocalToGlobal(p2)
                match_global.append([p1_g,p2_g])
            
            draw_matches(match_global, [0.0,0.0,0.0])
            
            xP = []
            yP = []
            for b in polyB:
                p1 = poseOrigin.convertLocalToGlobal(b)

                xP.append(p1[0])    
                yP.append(p1[1])
            
            p1 = poseOrigin.convertLocalToGlobal(polyB[0])
            xP.append(p1[0])    
            yP.append(p1[1])
            
            pylab.plot(xP,yP,linewidth=1, color=(0.0,0.0,1.0))

            xP = []
            yP = []
            for b in a_data_raw:
                p = [b[0],b[1]]
                p = dispOffset(p,offset)
                
                p1 = poseOrigin.convertLocalToGlobal(p)
                xP.append(p1[0])    
                yP.append(p1[1])
            
            p = [a_data_raw[0][0],a_data_raw[0][1]]
            p = dispOffset(p,offset)

            p1 = poseOrigin.convertLocalToGlobal(p)
            xP.append(p1[0])    
            yP.append(p1[1])

            pylab.plot(xP,yP,linewidth=1, color=(0.0,0.0,1.0))

            plotEnv()        

            newCost = cornerCostFunc(currAng, match_pairs, point1, point2, ang1, ang2)
            
            pylab.xlim(-4,9)
            pylab.ylim(-3,9)
            pylab.title("%u %u, cost = %f" % (n1, n2, newCost))
            #pylab.ylim(-7,5)
            pylab.savefig("ICP_plot_%04u.png" % numIterations)
            pylab.clf()                    

            
        allCircles = [[radiusA,centerA]]

        # optimize the match error for the current list of match pairs
        #newOffset = scipy.optimize.fmin(cost_func, offset, [match_pairs, a_data_raw, polyB, allCircles], disp = 0)
        #newAng = scipy.optimize.fmin(cornerCostFunc, currAng, [match_pairs, a_data_raw, polyB, allCircles], disp = 0)
        newAng = scipy.optimize.fmin(cornerCostFunc, currAng, [match_pairs, point1, point2, ang1, ang2, a_data_raw, polyB, allCircles], disp = 0)
        newAng = newAng[0]
        
            
        # get the current cost
        #newCost = cost_func(newOffset, match_pairs)    
        newCost = cornerCostFunc(newAng, match_pairs, point1, point2, ang1, ang2)


    
        # check for convergence condition, different between last and current cost is below threshold
        if abs(lastCost - newCost) < costThresh or (numIterations - startIteration) > 10:
            currAng = newAng
            lastCost = newCost
            break
    
        numIterations += 1

        currAng = newAng
        offset = computeOffset(point1, point2, ang1, ang2 + currAng)

        " reduce the minMatch distance for each step down to a floor value "
        minMatchDist /= 2
        if minMatchDist < 0.25:
            minMatchDist = 0.25
    
        # optionally draw the position of the points in current transform
        #if plotIter:
        if False:
            pylab.clf()
            pylab.axes()
            match_global = []
            
            for pair in match_pairs:
                p1 = pair[0]
                p2 = pair[1]
                
                p1_o = dispOffset(p1, offset)
                
                p1_g = poseOrigin.convertLocalToGlobal(p1_o)
                p2_g = poseOrigin.convertLocalToGlobal(p2)
                match_global.append([p1_g,p2_g])
            
            draw_matches(match_global, [0.0,0.0,0.0])
            
            xP = []
            yP = []
            for b in polyB:
                p1 = poseOrigin.convertLocalToGlobal(b)

                xP.append(p1[0])    
                yP.append(p1[1])
            
            p1 = poseOrigin.convertLocalToGlobal(polyB[0])
            xP.append(p1[0])    
            yP.append(p1[1])
            
            pylab.plot(xP,yP,linewidth=1, color=(0.0,0.0,1.0))

            xP = []
            yP = []
            for b in a_data_raw:
                p = [b[0],b[1]]
                p = dispOffset(p,offset)
                
                p1 = poseOrigin.convertLocalToGlobal(p)
                xP.append(p1[0])    
                yP.append(p1[1])
                
            p = [a_data_raw[0][0],a_data_raw[0][1]]
            p = dispOffset(p,offset)

            p1 = poseOrigin.convertLocalToGlobal(p)
            xP.append(p1[0])    
            yP.append(p1[1])
            
            pylab.plot(xP,yP,linewidth=1, color=(0.0,0.0,1.0))

            plotEnv()        
            
            #pylab.title("%u %u" % (n1, n2))
            pylab.title("%u %u, cost = %f" % (n1, n2, newCost))
            pylab.xlim(-4,9)
            pylab.ylim(-3,9)
            pylab.savefig("ICP_plot_%04u.png" % numIterations)
            pylab.clf()            
                        
        # save the current offset and cost
        currAng = newAng
        offset = computeOffset(point1, point2, ang1, ang2 + currAng)

        lastCost = newCost
        




    currAng =  functions.normalizeAngle(currAng)
    offset = computeOffset(point1, point2, ang1, ang2 + currAng)

    " transform the new pose "
    hull2_trans = []
    for p in hull2:
        result = dispPoint(p, offset)
        hull2_trans.append(result)

    " transformed pose 2 points without associated covariance "
    poly2 = []
    for p in hull2_trans:
        poly2.append([p[0],p[1]])

    poly1 = []        
    for p in hull1:
        poly1.append([p[0],p[1]])    


    " transform the new pose "
    sweepHull2_trans = []
    for p in sweepHull2:
        result = dispPoint(p, offset)        
        sweepHull2_trans.append(result)    

    " transformed pose 2 points without associated covariance "
    polySweep2 = []
    for p in sweepHull2_trans:
        polySweep2.append([p[0],p[1]])    

    " pose 1 without associated covariance"
    polySweep1 = []        
    for p in sweepHull1:
        polySweep1.append([p[0],p[1]])    

    " find the matching pairs "
    match_pairs = []

    " get the circles and radii "
    #radius2, center2 = computeEnclosingCircle(sweepHull2_trans)
    #radius1, center1 = computeEnclosingCircle(sweepHull1)
    radius2, center2 = computeEnclosingCircle(hull2_trans)
    radius1, center1 = computeEnclosingCircle(hull1)
    
    #for i_2 in range(len(sweepHull2_trans)):
    for i_2 in range(len(hull2_trans)):
        #p_2 = polySweep2[i_2]
        p_2 = poly2[i_2]

        #if isValidPast(p_2, [(radius1, center1)], polySweep1):
        if isValidPast(p_2, [(radius1, center1)], poly1):

            " for every transformed point of hull 2, find it's closest neighbor in hull 1 "
            #p_1, minDist = findClosestPointInHull1(sweepHull1, p_2, [0.0,0.0,0.0])
            p_1, minDist = findClosestPointInHull1(hull1, p_2, [0.0,0.0,0.0])

            if minDist <= minMatchDist:
    
                " add to the list of match pairs less than 1.0 distance apart "
                " keep A points and covariances untransformed "
                #Ca = sweepHull2[i_2][2]
                Ca = hull2[i_2][2]
                Cb = p_1[2]

                " we store the untransformed point, but the transformed covariance of the A point "
                #match_pairs.append([sweepHull2[i_2],p_1,Ca,Cb])
                match_pairs.append([hull2[i_2],p_1,Ca,Cb])

    #for i_1 in range(len(sweepHull1)):
    for i_1 in range(len(hull1)):
        #p_1 = polySweep1[i_1]
        p_1 = poly1[i_1]

        #if isValid(p_1, radius2, center2, polySweep2):
        if isValid(p_1, radius2, center2, poly2):
    
            " for every point of B, find it's closest neighbor in transformed A "
            #p_2, i_2, minDist = findClosestPointInHull2(sweepHull2_trans, p_1)
            p_2, i_2, minDist = findClosestPointInHull2(hull2_trans, p_1)

            if minDist <= minMatchDist:
    
                " add to the list of match pairs less than 1.0 distance apart "
                " keep A points and covariances untransformed "
                #Ca = sweepHull2[i_2][2]
                #Cb = sweepHull1[i_1][2]
                Ca = hull2[i_2][2]
                Cb = hull1[i_1][2]

                " we store the untransformed point, but the transformed covariance of the A point "
                #match_pairs.append([sweepHull2[i_2],p_1,Ca,Cb])
                match_pairs.append([hull2[i_2],p_1,Ca,Cb])

    #sweepCost = cornerCostFunc(currAng, match_pairs, point1, point2, ang1, ang2, isPrint = True)
    threshVal = 0.01
    upCount, downCount, matchVals = cornerMatchQuality(currAng, match_pairs, point1, point2, ang1, ang2, threshVal)
    

    # optionally draw the position of the points in current transform
    if plotIter:
        pylab.clf()
        pylab.axes()
        match_global = []
        
        for pair in match_pairs:
            p1 = pair[0]
            p2 = pair[1]
            
            p1_o = dispOffset(p1, offset)
            
            p1_g = poseOrigin.convertLocalToGlobal(p1_o)
            p2_g = poseOrigin.convertLocalToGlobal(p2)
            match_global.append([p1_g,p2_g])
        
        draw_matches(match_global, [0.0,0.0,0.0])
        
        xP = []
        yP = []
        for b in polyB:
            p1 = poseOrigin.convertLocalToGlobal(b)

            xP.append(p1[0])    
            yP.append(p1[1])
        
        p1 = poseOrigin.convertLocalToGlobal(polyB[0])
        xP.append(p1[0])    
        yP.append(p1[1])
        
        pylab.plot(xP,yP,linewidth=1, color=(0.0,0.0,1.0))

        xP = []
        yP = []
        for b in a_data_raw:
            p = [b[0],b[1]]
            p = dispOffset(p,offset)
            
            p1 = poseOrigin.convertLocalToGlobal(p)
            xP.append(p1[0])    
            yP.append(p1[1])
            
        p = [a_data_raw[0][0],a_data_raw[0][1]]
        p = dispOffset(p,offset)

        p1 = poseOrigin.convertLocalToGlobal(p)
        xP.append(p1[0])    
        yP.append(p1[1])
        
        pylab.plot(xP,yP,linewidth=1, color=(0.0,0.0,1.0))        
        
        #pylab.title("%u %u" % (n1, n2))
        pylab.title("%u %u, ratio = %f, thresh = %f, over = %d, under = %d" % (n1, n2, float(upCount)/float(downCount), threshVal, upCount, downCount))
        pylab.xlim(-4,9)
        pylab.ylim(-3,9)
        pylab.savefig("ICP_plot_%04u.png" % numIterations)
        pylab.clf()            
        
    
    #return offset, lastCost

    if (downCount) > 0:
        return offset, float(upCount)/float(downCount)

    else:
        return offset, float(upCount)


def quickICP(curve1, curve2, offset, costThresh = 0.004, minMatchDist = 2.0, plotIter = False, n1 = 0, n2 = 0):

    global numIterations
    startIteration = numIterations    

    vec = curve1.getUVector(0.5)
    vecPoint = curve1.getU(0.5)
    
    angle = math.acos(vec[0])
    if math.asin(vec[1]) < 0:
        angle = -angle
    
    localAIRPose1 = [vecPoint[0], vecPoint[1], angle]    

    vec = curve2.getUVector(0.5)
    vecPoint = curve2.getU(0.5)
    
    angle = math.acos(vec[0])
    if math.asin(vec[1]) < 0:
        angle = -angle

    localAIRPose2 = [vecPoint[0], vecPoint[1], angle]    


    localAIRProfile1 = Pose(localAIRPose1)
    
    relPose2 = localAIRProfile1.convertGlobalPoseToLocal(localAIRPose2)


    #midPoint1 = points1[len(points1)/2]
    #midPoint2 = points2[len(points2)/2]

    #handlePoint1 = points1[len(points1)/2+5]
    handlePoint2 = points2[len(points2)/2+5]
    
    #points1 = [midPoint1, handlePoint1]
    #points2 = [midPoint2, handlePoint2]
    #points1 = addGPACVectorCovariance(points1,high_var=0.1)
    #points2 = addGPACVectorCovariance(points2,high_var=0.1)

    #points2_offset = []
    #for p in points2:
    #    result = dispPoint(p, offset)        
    #    points2_offset.append(result)

    " correct only the orientation "

    pass

def motionICP(points1, points2, offset, costThresh = 0.004, minMatchDist = 2.0, plotIter = False, n1 = 0, n2 = 0):

    global numIterations
    
    lastCost = 1e100
    startIteration = numIterations
    
    lastNumPairs = 1e100


    " sample a series of points from the AIRL curves "
    
    " augment points with point-to-line covariances "
    
    " treat the points with the point-to-line constraint "
    points1 = addGPACVectorCovariance(points1,high_var=0.1)
    points2 = addGPACVectorCovariance(points2,high_var=0.1)
    
    " transform pose 2 by initial offset guess "
    
    " transform the new pose "
    points2_offset = []
    for p in points2:
        result = dispPoint(p, offset)        
        points2_offset.append(result)
    
    poly1 = []        
    for p in points1:
        poly1.append([p[0],p[1]])    



    while True:
        
        " find the matching pairs "
        match_pairs = []

        " transform the target Hull with the latest offset "
        points2_offset = []
        for p in points2:
            result = dispPoint(p, offset)        
            points2_offset.append(result)
        
        " transformed points without associated covariance "
        poly2 = []
        for p in points2_offset:
            poly2.append([p[0],p[1]])    
        
        " get the circles and radii "
        radius2, center2 = computeEnclosingCircle(points2_offset)
        radius1, center1 = computeEnclosingCircle(points1)
        
        if True:
            for i in range(len(points2_offset)):
                p_2 = poly2[i]
    
                if isValid(p_2, radius1, center1, poly1):
    
                    " for every transformed point of A, find it's closest neighbor in B "
                    p_1, minDist = findClosestPointInB(points1, p_2, [0.0,0.0,0.0])
        
                    if minDist <= minMatchDist:
            
                        " add to the list of match pairs less than 1.0 distance apart "
                        " keep A points and covariances untransformed "
                        C2 = points2[i][2]
                        C1 = p_1[2]
        
                        " we store the untransformed point, but the transformed covariance of the A point "
                        match_pairs.append([points2[i],p_1,C2,C1])

        if True:
            for i in range(len(points1)):
                p_1 = poly1[i]
        
                if isValid(p_1, radius2, center2, poly2):
            
                    #print "selected", b_p, "in circle", radiusA, centerA, "with distance"
                    " for every point of B, find it's closest neighbor in transformed A "
                    p_2, i_2, minDist = findClosestPointInA(points2_offset, p_1)
        
                    if minDist <= minMatchDist:
            
                        " add to the list of match pairs less than 1.0 distance apart "
                        " keep A points and covariances untransformed "
                        C2 = points2[i_2][2]
                        
                        C1 = points1[i][2]
        
                        " we store the untransformed point, but the transformed covariance of the A point "
                        match_pairs.append([points2[i_2],p_1,C2,C1])

        if plotIter and startIteration == numIterations:
            
            " set the origin of pose 1 "
            poseOrigin = Pose([0.0,0.0,0.0])
    
            pylab.clf()
            pylab.axes()
            match_global = []
            
            for pair in match_pairs:
                p1 = pair[0]
                p2 = pair[1]
                
                p1_o = dispOffset(p1, offset)
                
                p1_g = poseOrigin.convertLocalToGlobal(p1_o)
                p2_g = poseOrigin.convertLocalToGlobal(p2)
                match_global.append([p1_g,p2_g])
            
            draw_matches(match_global, [0.0,0.0,0.0])
            
            xP = []
            yP = []
            for b in poly1:
                p1 = poseOrigin.convertLocalToGlobal(b)
    
                xP.append(p1[0])    
                yP.append(p1[1])
    
            pylab.plot(xP,yP,linewidth=1, color=(0.0,0.0,1.0))
    
            xP = []
            yP = []
            for b in points2:
                p = [b[0],b[1]]
                p = dispOffset(p,offset)
                
                p1 = poseOrigin.convertLocalToGlobal(p)
                xP.append(p1[0])    
                yP.append(p1[1])
            
            pylab.plot(xP,yP,linewidth=1, color=(0.0,0.0,1.0))
    
            #plotEnv()        
            
            pylab.title("%u %u" % (n1, n2))
            pylab.xlim(-3,3)
            pylab.ylim(-3,3)
            #pylab.axis('equal')
            pylab.savefig("ICP_plot_%04u.png" % numIterations)
            pylab.clf()            
        
        #print "pairs =", len(match_pairs)
                        
        #if len(match_pairs) < lastNumPairs:
        #    lastNumPairs = len(match_pairs)
        #else:
        #    print "breaking"
        #    break                
                        
                        
        # optimize the match error for the current list of match pairs
        newOffset = scipy.optimize.fmin(cost_func, offset, [match_pairs, [], [], []], disp = 0)
        #newOffset = scipy.optimize.leastsq(cost_func2, offset, match_pairs)
        #print "newOffset:", newOffset
        # get the current cost
        #newCost = cost_func2(newOffset[0], match_pairs)
        newCost = cost_func(newOffset, match_pairs)
    
        #print "newCost:", newCost
        
        #totalCost = 0.0
        #for i in range(len(newCost)):
        #    totalCost += newCost[i]
        #newCost = totalCost
        
        # check for convergence condition, different between last and current cost is below threshold
        if abs(lastCost - newCost) < costThresh or (numIterations - startIteration) > 10:
            #offset = newOffset[0]
            offset = newOffset
            lastCost = newCost
            break

        numIterations += 1

        " reduce the minMatch distance for each step down to a floor value "
        #minMatchDist /= 2
        if minMatchDist < 0.25:
            minMatchDist = 0.25
    
        
        # save the current offset and cost
        #offset = newOffset[0]
        offset = newOffset
        lastCost = newCost
        
        
        # optionally draw the position of the points in current transform
        if plotIter:
            
            " set the origin of pose 1 "
            poseOrigin = Pose([0.0,0.0,0.0])
    
            pylab.clf()
            pylab.axes()
            match_global = []
            
            for pair in match_pairs:
                p1 = pair[0]
                p2 = pair[1]
                
                p1_o = dispOffset(p1, offset)
                
                p1_g = poseOrigin.convertLocalToGlobal(p1_o)
                p2_g = poseOrigin.convertLocalToGlobal(p2)
                match_global.append([p1_g,p2_g])
            
            draw_matches(match_global, [0.0,0.0,0.0])
            
            xP = []
            yP = []
            for b in poly1:
                p1 = poseOrigin.convertLocalToGlobal(b)

                xP.append(p1[0])    
                yP.append(p1[1])

            pylab.plot(xP,yP,linewidth=1, color=(0.0,0.0,1.0))

            xP = []
            yP = []
            for b in points2:
                p = [b[0],b[1]]
                p = dispOffset(p,offset)
                
                p1 = poseOrigin.convertLocalToGlobal(p)
                xP.append(p1[0])    
                yP.append(p1[1])
            
            pylab.plot(xP,yP,linewidth=1, color=(0.0,0.0,1.0))

            #plotEnv()        
            
            pylab.title("%u %u" % (n1, n2))
            pylab.xlim(-3,3)
            pylab.ylim(-3,3)
            #pylab.axis('equal')
            pylab.savefig("ICP_plot_%04u.png" % numIterations)
            pylab.clf()            

    #print "offset:", offset
    offset[2] =  functions.normalizeAngle(offset[2])
    
    return offset

def plotEnv(axes = 0):

    """
    WLEN = 3.0
    WLEN2 = 5.0
    wall1 = [[-14.0, -0.2], [-4.0, -0.2], [-4.0 + WLEN*math.cos(math.pi/3), -0.2 - WLEN*math.sin(math.pi/3)]]
    wall2 = [[-4.0 + WLEN2*math.cos(math.pi/3), 0.2 + WLEN2*math.sin(math.pi/3)], [-4.0, 0.2] ,[-14.0, 0.2]]
    w1 = wall1[2]
    w2 = wall2[0]
            
    wall3 = [[w1[0] + 0.4*math.cos(math.pi/6), w1[1] + 0.4*math.sin(math.pi/6)], [0.4*math.cos(math.pi/6) - 4, 0.0], [w2[0] + 0.4*math.cos(math.pi/6), w2[1] - 0.4*math.sin(math.pi/6)], w2]
    lp = wall3[0]
    rp = wall3[2]
    
    wall6 = [lp, [lp[0] + WLEN*math.cos(math.pi/6), lp[1] + WLEN*math.sin(math.pi/6)]]
    wall6.append([wall6[1][0] + 0.4*math.cos(math.pi/3), wall6[1][1] - 0.4*math.sin(math.pi/3)])
    wall6.append([wall6[2][0] - WLEN*math.cos(math.pi/6), wall6[2][1] - WLEN*math.sin(math.pi/6)])
    wall6.append([wall6[3][0] + WLEN*math.cos(math.pi/3), wall6[3][1] - WLEN*math.sin(math.pi/3)])
    wall6.append([wall6[4][0] - 0.4*math.cos(math.pi/6), wall6[4][1] - 0.4*math.sin(math.pi/6)])
    wall6.append(w1)
    wall6.reverse()

    walls = [wall1, wall2, wall3, wall6]

    for wall in walls:
        xP = []
        yP = []
        for i in range(len(wall)):
            p = copy(wall[i])
            p[0] += 6.0
            wall[i] = p
            xP.append(p[0])
            yP.append(p[1])

        pylab.plot(xP,yP, linewidth=2, color = 'g')

    """
    #walls = [[[-8.0, -0.20000000000000001], [2.0, -0.20000000000000001], [2.0, -3.2000000000000002], [-0.12132034355964194, -5.321320343559643], [0.1615223689149774, -5.6041630560342623], [2.3999999999999999, -3.3656854249492403], [2.3999999999999999, -0.20000000000000001], [5.4000000000000004, -0.20000000000000001], [7.5213203435596423, -2.3213203435596426], [7.8041630560342616, -2.0384776310850237], [5.5656854249492387, 0.20000000000000018], [2.3999999999999999, 0.20000000000000018], [2.3999999999999999, 3.2000000000000002], [4.5213203435596423, 5.321320343559643], [4.2384776310850238, 5.6041630560342623], [2.0, 3.3656854249492389], [2.0, 0.20000000000000018], [-8.0, 0.20000000000000001]]]
    #walls = [[[-8.0, -0.20000000000000012, 0.29085552099855971, -0.41088507028870347], [2.0, -0.20000000000000012, 0.29085552099855971, -0.41088507028870347], [2.2806451181559613, -0.27638881953905159, 0.26520770061705834, 0.11321713192974797], [2.5370513151779122, -0.20863199999055054, 0.24795215590938968, -0.10981024281727976], [2.7848488018333271, -0.19987544659475237, 0.23353847966708152, 0.20374476210993273], [3.0043182007953995, -0.1200418421136879, 0.21975917568865269, 0.29657781267365446], [3.2029852426446048, -0.026097700293977899, 0.21508539406709762, 0.056088098579494067], [3.413730909394566, 0.016890561907613677, 0.20388010253525796, 0.36028613207975529], [3.5921202043788938, 0.11560406475693091, 0.19801329528243747, 0.6273076010075348], [3.7339392048689763, 0.2537946440617575, 0.1902961062461512, 0.49402567316832524], [3.8866704525983273, 0.36731016731511051, 0.18787266335561981, 0.87013002267287798], [3.9857538206767789, 0.52693041605117807, 0.18805516598324074, 0.50575791094240441], [4.1353600294183792, 0.64087218959779535, 0.18246656688688789, 1.0276187039774782], [4.2060870814468139, 0.80907365256067898, 0.18094859179440287, 1.1980967382592496], [4.2469104296909244, 0.98535708090979379, 0.18744931072062979, 0.6415632716498898], [4.3792849257343605, 1.1180757189637074, 0.17738517995107692, 0.98175418024357552], [4.4554670591833201, 1.2782686795126676, 0.17647308121096553, 1.2376680060904997], [4.4884480834825116, 1.451632460125194, 0.17466290163270304, 1.1510105858787698], [4.5358188682188327, 1.6197489017964961, 0.17559377378404828, 1.1105723460317549], [4.5902358869511843, 1.7866978804217233, 0.1926482355884373, 1.7440654185577869], [4.5299274630986961, 1.9696630046606525, 0.16930651818779779, 1.2714992386473354], [4.5559251683878728, 2.1369615892947498, 0.1784168472914523, 0.995641177482397], [4.6303055211266066, 2.2991347514513873, 0.1670242404225995, 1.2789115335675816], [4.6547287298560995, 2.4643636949331573, 0.16962507193267959, 1.6515201351823057], [4.6167426243189809, 2.6296807261619133, 0.1592835246410266, 1.4587885214924494], [4.6114672780772015, 2.7888768691556991, 0.15639033518206449, 1.2214512460096039], [4.6431825489510974, 2.9420175839632048, 0.15174813396991674, 1.5873386985585656], [4.6187553426035883, 3.0917867655507671, 0.15819212388489601, 0.92698317128697394], [4.6944133795562237, 3.2307133990797205, 0.13502839570386299, 1.483100038020851], [4.6866620448249776, 3.3655191276532617, 0.14913376918436616, 1.9058069741094452], [4.6177761636677772, 3.497790124052063, 0.12632120871674016, 0.94225945997518656], [4.6764897098202862, 3.6096371949562913, 0.11167304709277615, 1.5740992757407912], [4.6599742103711375, 3.720082238869229, 0.097529896947194508, 1.366927356104384], [4.6656994221293253, 3.8174439497665533, 0.11650119104995019, 0.33988210872246405], [4.768764282213013, 3.871759344747967, 0.089726115209923796, 1.8002207501797118], [4.7359370183749725, 3.9552647128363247, 0.085561688551505255, -0.043869843881169504], [4.8210603964578098, 3.9639141666659636, 0.087979840904632706, 1.9069193803621254], [4.780335144097684, 4.0419007438815875, 0.054159765459383231, -0.45449491884675819], [4.8319238400796536, 4.025411758421737, 0.046631346084173464, 1.3293618215083638], [4.8364075671642386, 4.0718270428553467, 0.071229361127770202, -0.99999659039555366], [4.8831569281347136, 4.0180858190553392, 0.091205748966026678, 0.95215096310744884], [4.924748110710647, 4.0992563936275896, 0.076049917897750355, -0.51346806525829836], [4.9956972198778899, 4.0718736632920409, 0.09623796606459066, -0.62208121936721161], [5.0811950363771281, 4.0276936962695133, 0.12215181295522384, -0.8204997659120179], [5.1765318097354438, 3.9513263927039599, 0.14403319063568581, -0.31964930185390095], [5.3183772304104941, 3.9263176614643767, 0.17073170501151289, -0.46074509597638574], [5.4806759142188355, 3.8733228284122925, 0.19756589258419874, -0.78896091522290712], [5.638689890309978, 3.7547316598015605, 0.22740669680793218, -0.54077936863169296], [5.8485288912848681, 3.6670880666360857, 0.2599206498104315, -0.96900556564566886], [6.0251150170105348, 3.4763629125554689, 0.28324555029563081, -0.70435785672767193], [6.2652128227196453, 3.3260928707775239, 0.3141826253410735, -0.79589300641163296], [6.5517604668270337, 3.7155076980781763, 0.32631742129069846, -0.6943075247967454], [6.3115102267736898, 3.8870252404764871, 0.29519221737038642, -0.76513510254712824], [6.0856243360704738, 4.0302523509938286, 0.26746670971530573, -0.7102269509424779], [5.9017420489637846, 4.1924242986736679, 0.24517837613839372, -0.86787655031941868], [5.7210666901139646, 4.3182692459138448, 0.22018296037930499, -0.75352318399716123], [5.5396895801414594, 4.3984769554676006, 0.1983202780701902, -0.56149431522686644], [5.3760980430666052, 4.4712086487346481, 0.17903097555452518, -0.56348166753639928], [5.2236551357450551, 4.5320148411496799, 0.16412261583542043, -0.5246722441158127], [5.0752565908243827, 4.5680581880536302, 0.15271296929406816, -0.38340138792334849], [4.9297066220320716, 4.5167265711186708, 0.15433641278907731, 0.19392148453960353], [4.7976355971529134, 4.5199451927717345, 0.13211023858118182, -0.16949864539854909], [4.6598602861442497, 4.5168327544855602, 0.13781046257675433, -0.12254624331079732], [4.5618586313324307, 4.4218483491032909, 0.13647842910760119, 0.62463180727372991], [4.4335724845084563, 4.3740996654215047, 0.13688415635231688, 0.21118450487343896], [4.3707786976247487, 4.2614968690461321, 0.12892807848855131, 0.91696001674048078], [4.2963235398023301, 4.1615774295663739, 0.12460924890356352, 0.78527216054333704], [4.2909014989818761, 4.0332682297503704, 0.12842370997631933, 1.383430753822529], [4.21976013024263, 3.9264838744375572, 0.12831209173599512, 0.83797320577251222], [4.2436972742117334, 3.7935955543596895, 0.13502700646357899, 1.6038818861998509], [4.1941760521235203, 3.6708507276947158, 0.13235801415148468, 1.0421875193893313], [4.1643019835383583, 3.5390093149246793, 0.13518364581202172, 1.2028350455872336], [4.1277479584633028, 3.4021558452858809, 0.14165122273162764, 1.1646531703620502], [4.1143378707313261, 3.2590644866120084, 0.14371836131863644, 1.3322192098904273], [4.1666223048412441, 3.1096516391997486, 0.15829674987830153, 1.7622780638642241], [4.1038458359366636, 2.9617269489371081, 0.16069411637463904, 1.0243167764197789], [4.1224922262983874, 2.8081261213640936, 0.15472847864771602, 1.5464672531734851], [4.0804234026941142, 2.6552200095837266, 0.1585877200139951, 1.1571776043902695], [4.1090228875288659, 2.497726098004212, 0.16006955587315558, 1.6052968707198705], [4.1360722515818935, 2.3386334401241258, 0.16137577850168058, 1.59407546733993], [4.139152166743834, 2.1803598261272268, 0.15830357786558111, 1.4451202226225672], [4.0620247233350728, 2.0324030907250674, 0.16685274369226427, 0.94513396529910776], [4.084465722580692, 1.8728260558812453, 0.1611472261526338, 1.5653750589779227], [4.0399984744428732, 1.7255753115732171, 0.15381845746274597, 1.1323880111091904], [3.9519494567660143, 1.5935200937984013, 0.1587173905260984, 0.8375966523047077], [3.8992836115295368, 1.4588317249962448, 0.1446189750518872, 1.0529221398022031], [3.8541877726969425, 1.3258235603835247, 0.14044503029178265, 1.0987804491268034], [3.8092897337934772, 1.1952707542613714, 0.13805748470744181, 1.094426661467458], [3.7691017100082327, 1.0624205633524242, 0.13879571492055212, 1.1319086515877685], [3.7371114230294946, 0.91966805664437978, 0.14629305052671587, 1.2052092655520044], [3.614618669201203, 0.84105662971869766, 0.14554803737505545, 0.42543045027795079], [3.5209377444026391, 0.74300639064990359, 0.13560960531084018, 0.66304988952158783], [3.4039244618215125, 0.66807361916304164, 0.13894973387204754, 0.42444846248232337], [3.3208542230365197, 0.53625754823121907, 0.15580802651885103, 0.86333469043768263], [3.1693994718595171, 0.49254877076103026, 0.15763565231896726, 0.13582597787584802], [3.0286394902609644, 0.4104173496773596, 0.16296914661630862, 0.38305495877500034], [2.8719552797225933, 0.32231329533603686, 0.17975612986324455, 0.36710696955666272], [2.671910147312957, 0.33573958465321102, 0.20049518758717727, -0.21214887903744165], [2.4727670898907519, 0.27661024647542037, 0.20773597654909609, 0.14349462469880658], [2.2476092306134596, 0.25135160988260447, 0.22657021056807342, -0.033418256283370334], [2.0, 0.20000000000000001, 0.2528780712567319, 0.059357777011700676], [-8.0, 0.20000000000000001, 0.2528780712567319, 0.059357777011700676]]]
    #walls = [[[-8.0, -0.20000000000000048, 0.31298415685345782, 0.10317334362666691], [2.0, -0.20000000000000048, 0.31298415685345782, 0.10317334362666691], [2.3071904162504442, -0.25994272770861715, 0.29327898851278633, 0.51688445216082046], [2.5933364543302773, -0.19565423733671619, 0.28087753844881347, 0.14769362818778828], [2.8711355458641687, -0.2371253881286737, 0.26272334780665285, 0.61144772257142976], [3.1208860235395708, -0.15558854747584935, 0.2435153177948646, 0.56254784555905513], [3.3557943826485523, -0.091418672756796576, 0.23181694697854999, 0.3880227689389173], [3.5866280130140673, -0.070089537782934169, 0.21872917843412332, 0.61022587331919154], [3.7946394726829471, -0.002460504484445375, 0.219272587417218, 0.34125531323021063], [4.0136864056722867, 0.0074847828721419096, 0.20308371831881264, 0.80331934627848245], [4.1911802774120819, 0.10617076030061276, 0.19378430893660895, 0.88199945092780929], [4.3526211282646869, 0.21335849311564001, 0.185269221641716, 0.95746338916155749], [4.4988027783177653, 0.32718120471348644, 0.18712474002273635, 0.62562161128164417], [4.6758466059606221, 0.38777123212842746, 0.18867760475626413, 1.3950708178977318], [4.7615667330729874, 0.55585245756541946, 0.16761874441075153, 0.99970458361550951], [4.8893551201156225, 0.66432444295740278, 0.1669907360611419, 1.0090979628973684], [5.0156440264966928, 0.77358110433470784, 0.16508856674322461, 1.3868746271803885], [5.0918500252839376, 0.92002863971229853, 0.16045234583595991, 1.3929975607110563], [5.1650430236110996, 1.0628142842293247, 0.17185166355940495, 0.99096765870188308], [5.2970251151789736, 1.1728766382542311, 0.16349301875888431, 1.2950619777546792], [5.3854738202913488, 1.31037861550124, 0.16277623868820101, 1.4341023401056761], [5.4537116792664868, 1.4581612209863495, 0.17985703828660979, 1.9961995577222997], [5.4304817967881451, 1.6365117936831404, 0.15963347717156232, 1.6267294060413413], [5.4684194727243831, 1.7915717168691165, 0.16303391527974936, 1.4504653379580683], [5.5343343171039496, 1.9406867423422789, 0.16709041486228279, 1.4260889090677713], [5.6055941039420834, 2.0918199601292149, 0.1650383871875655, 1.60133574821126], [5.6488740571654557, 2.2510823690475008, 0.166114307680051, 1.7976874560389009], [5.6603256586924093, 2.4168014809160199, 0.16691132921007992, 1.6352236805866074], [5.6986143667835432, 2.5792618350070637, 0.16770432208690567, 1.6636922970789401], [5.7324230281999373, 2.7435229504904926, 0.17261673276976497, 1.9859007891948426], [5.7118922861258472, 2.9149143886138589, 0.16926704684511745, 1.8659597796870624], [5.7120142393435156, 3.0841813915266534, 0.17363562063897819, 1.5187008096253691], [5.7712238223607475, 3.247409924674534, 0.16902067373540397, 1.775330992040727], [5.7866422722422062, 3.4157258762454497, 0.17020183346898765, 1.8613986037205381], [5.787541215296665, 3.5859253357557246, 0.17047042202585933, 1.610534379087551], [5.8307305901405098, 3.7508339246353604, 0.16872685716728622, 1.6881897734728386], [5.8606870722183775, 3.9168801878003827, 0.17012007155328246, 1.8544529346726784], [5.8627671334978828, 4.0869875424043727, 0.17686650563505055, 1.9941216468931295], [5.8402879843199837, 4.262419719235738, 0.182928928856054, 1.2588898506909498], [5.9447504591800087, 4.4125882391010629, 0.17803285200968011, 1.9845635533873958], [5.9238119338973316, 4.5893855087453845, 0.16834454812775868, 1.5187263699266922], [5.9812132211204165, 4.7476415642676679, 0.17044571470977041, 1.7482473316261276], [6.0013524487341208, 4.9168933092337035, 0.17468016773259434, 1.2394048367119672], [6.1038794079409335, 5.0583195572744826, 0.16880259792514429, 1.2963685965763856], [6.1950149870449724, 5.2004062325387368, 0.17838149723128438, 1.8039927673417662], [6.206189953036759, 5.3784373497333355, 0.17004175216199441, 1.2826008126025501], [6.2999563605589888, 5.5202895908107896, 0.1706187168286801, 1.395916958966495], [6.377343326614616, 5.6723488030176492, 0.17468176526702689, 1.1187669668215228], [6.4961462215459154, 5.8004096919445463, 0.1762897630969516, 1.3529827617329957], [6.5827751648561925, 5.9539463545164679, 0.18419938614338302, 0.92378090134499069], [6.731841227581624, 6.0621532064212191, 0.18677271514601804, 1.3194900367937827], [6.3139952781655833, 6.4274145398526201, 0.21427851231471004, 1.1565557160603011], [6.2064911530488072, 6.2466981389257263, 0.21027494972319763, 1.3300357975846544], [6.0694793274012948, 6.0954307106224945, 0.20409330031320841, 1.1306925166049124], [5.9650194275297057, 5.9252961919962841, 0.19964374548018105, 1.3160416898135956], [5.8718647861524769, 5.7528212041458957, 0.19602400017373178, 1.3714655326798053], [5.7580748537250859, 5.5934746490490843, 0.19580468161678735, 1.2465514389377454], [5.6580812004624299, 5.4274861480100656, 0.19378057996092582, 1.3244882826767079], [5.6761320177902714, 5.2173603640395729, 0.21089968491066344, 1.9523746916941911], [5.5459152090265142, 5.0654129170875262, 0.20011107895332172, 1.1581445150561049], [5.4989808724131048, 4.8852853611325342, 0.18614179640172032, 1.6117861670508931], [5.4869613582497916, 4.6983721805817398, 0.1872992412807265, 1.8024633444537721], [5.3947316556807667, 4.5324948701896739, 0.1897935724382587, 1.3592333627216537], [5.3736003985195033, 4.3514581341337815, 0.18226582189479903, 1.7504824757242203], [5.383582893545233, 4.1673315186220004, 0.18439701935433159, 1.9208426161147356], [5.3573682731168528, 3.9916332589078851, 0.17764313887896888, 1.7185703479037666], [5.3092348892421253, 3.8203270096869475, 0.1779400282830649, 1.5927643636975071], [5.3038773664988614, 3.6448493366223071, 0.17555943949023434, 1.8361586598806081], [5.291217140680617, 3.4723129346261428, 0.17300026396381957, 1.7934343976475378], [5.2215103628411548, 3.3081265154608217, 0.17837100413209064, 1.4651834283209586], [5.2775377560310393, 3.1322276709116634, 0.18460626289858834, 2.1750405820649688], [5.1669605405663148, 2.9774275499724276, 0.19023774079486372, 1.2464063882992986], [5.1454702376419181, 2.8157141998249924, 0.16313503834458548, 1.7345630290031717], [5.1387070682193432, 2.6551865916255086, 0.16067001417466009, 1.8245742788756132], [5.152787482033597, 2.4943478155425827, 0.16145392514595028, 1.9540012908977629], [5.1567605922524091, 2.3371334256362579, 0.15726458596400181, 1.8919468045415595], [5.1135111928291668, 2.1913163857736366, 0.15209575820767773, 1.5783452288228019], [5.1089131664415266, 2.0407699399652874, 0.15061664646441814, 1.8361475048274465], [5.0289468765044774, 1.9122576218012466, 0.15136057427940786, 1.3100637202119618], [4.986576297185632, 1.7789895811781979, 0.13984146968378539, 1.5588515876183222], [4.9424963509165085, 1.649945300813102, 0.13636520068563743, 1.537519209153172], [4.9301626747443184, 1.5099082071599594, 0.14057918468514274, 1.7788325862213554], [4.8012354330214517, 1.4264617966572286, 0.15357583496128818, 0.87031395761918096], [4.7829605142204752, 1.290732275143869, 0.13695428312918057, 1.7328429611406364], [4.6683639681248206, 1.2145073073862747, 0.13763216951970103, 0.88284230373731776], [4.6140639931760479, 1.0991091935526576, 0.12753514008222855, 1.4268732888950559], [4.5268976446705143, 1.0070638989026082, 0.12676872082261773, 1.1084998646351347], [4.4383883836867408, 0.91286210346483521, 0.1292589166889222, 1.1124279580686243], [4.3428021442331044, 0.81863525329976783, 0.13422156482815231, 1.0741204912591675], [4.2400881733482976, 0.72060090205048311, 0.14198906239493606, 1.0579754375844241], [4.0982563262823124, 0.67038643356515248, 0.15045851816169523, 0.6361553294904897], [3.98932168512831, 0.53776460565218931, 0.17162547970020917, 1.1790302419263774], [3.8095188637703679, 0.531994655822699, 0.1798953776207467, 0.32796335230106749], [3.6564404139892872, 0.42977788208061607, 0.18406868452188949, 0.88463001021710741], [3.4623200190877288, 0.41192882055461244, 0.19493926416732063, 0.38757452080519034], [3.269671889606057, 0.33346019537642546, 0.20801592951055811, 0.68268107153529223], [3.0428392486672049, 0.35359997380331509, 0.22772496057827946, 0.20732918961957375], [2.8145235753130073, 0.28281446850503877, 0.23903688932780309, 0.59652015575684447], [2.5612935592305224, 0.26952348155819872, 0.25357857042573051, 0.34832164285740974], [2.2901167447046831, 0.26150571448839033, 0.2712953175511269, 0.32544188573295552], [2.0, 0.20000000000000001, 0.29656479641516803, 0.5047940600366122], [-8.0, 0.20000000000000001, 0.29656479641516803, 0.5047940600366122]]]
    #walls = [[[-8.0, -0.20000000000000059, 0.27791203840248918, 1.2208669106164414], [2.0, -0.20000000000000059, 0.27791203840248918, 1.2208669106164414], [2.277911859887789, -0.20031499640799058, 0.25587222218964312, 1.287420976921964], [2.5332367288937849, -0.18358761301018775, 0.23535362996837378, 1.302551422106105], [2.7678272285294891, -0.16465012056899606, 0.2164794330037782, 1.3716173090798207], [2.9818882066335437, -0.13238183014676561, 0.19932012957350637, 1.3987731775403049], [3.178102196669828, -0.097330665155569401, 0.18400803241344268, 1.4697901307717358], [3.3564900481364646, -0.05220051512497817, 0.17036231572252442, 1.566658581363596], [3.5168335234849462, 0.0053606569575316565, 0.1584742994428473, 1.6151785368642197], [3.6632156128053257, 0.066076285177137456, 0.14861998250447794, 1.7068460643966701], [3.7947066768276212, 0.13534387703592596, 0.14038322051201521, 1.8146784442746613], [3.9111472181418669, 0.2137597487697189, 0.13371546205086593, 1.9015551525794026], [4.0151581351202221, 0.29779280143828668, 0.1283204200977201, 2.0244760539793214], [4.1043316669098804, 0.39006528618977843, 0.12381015059535978, 2.1385473996704674], [4.1796780318115179, 0.48830926150915788, 0.11998408411196582, 2.2605518515764333], [4.2405662250303777, 0.59169595488277393, 0.11683091926769172, 2.3421130569712014], [4.291455553208495, 0.69686125322074233, 0.11425995413472256, 2.4562006687127647], [4.3291928780572757, 0.80470944018800894, 0.11211751319851924, 2.5245654491405691], [4.3589069710099722, 0.91281776292876093, 0.11058256792450583, 2.5641204793440249], [4.383974754430854, 1.0205215723467276, 0.10918472839613408, 2.6514300248039255], [4.3993584748201329, 1.1286171137067766, 0.10786088601064941, 2.6940529593151186], [4.4099917603498522, 1.2359525876650435, 0.10654855995536816, 2.7202953323659722], [4.4177099083845697, 1.3422212369435436, 0.10510080201653624, 2.7503830924304333], [4.422166273745388, 1.4472275195760129, 0.10342502274374744, 2.7499551280011394], [4.4265958066303526, 1.550557643779779, 0.10152565375332427, 2.7946491255538555], [4.4264077356760954, 1.6520831233371591, 0.099173736564734119, 2.7444869426030039], [4.4311969291693529, 1.7511411550655125, 0.096530116299294699, 2.711713397570235], [4.4390153339614908, 1.8473541266403557, 0.093538110979353703, 2.6673726449808397], [4.4507165256381729, 1.9401574688157941, 0.090480928306037878, 2.6956203344362502], [4.4594952993952912, 2.0302115165784773, 0.087261251738381229, 2.5600169351483344], [4.4796250024704225, 2.1151192377177015, 0.084067618823878476, 2.506950548526456], [4.5033294940956896, 2.1957756795652354, 0.081530413608073782, 2.4514025639209458], [4.5306259641471378, 2.2726008776830007, 0.080270671386639375, 2.3963870749893976], [4.5616191818714604, 2.3466468170716981, 0.080332491892671137, 2.2463330934936012], [4.6033654944222935, 2.415280297546567, 0.082202036664384914, 2.1472631292805189], [4.6528202709486575, 2.480941548756852, 0.086395761650031516, 2.0408574349822093], [4.7118334488027678, 2.5440420437679716, 0.093091587117425767, 1.9401704500178907], [4.781932445152659, 2.6052970378727727, 0.10235609930039341, 1.8524844865486099], [4.8646097846720071, 2.6656396314033319, 0.11374269501476442, 1.7059714679085747], [4.9652896036789578, 2.7185638670043257, 0.12707592958728525, 1.6612844565127749], [5.0803004817869031, 2.7726081987821543, 0.14255225928290238, 1.587438487698559], [5.213439627447821, 2.8235504653134038, 0.15962629596725134, 1.5106271568435652], [5.3664630837783722, 2.8689858719570793, 0.17812209172977703, 1.4411595312959309], [5.5403246012187246, 2.9077112165718058, 0.19797044907284569, 1.3843755876663151], [5.7356909629006712, 2.9397156448476389, 0.21915054597446532, 1.3397442706094935], [5.9533241519657594, 2.9654597085957359, 0.24160155590154561, 1.289910000337575], [6.1943688224975819, 2.9818541778537488, 0.26528099260096843, 1.260198069051347], [6.4594563074255973, 2.9919848426836113, 0.29024598663394696, 1.2323817455052308], [6.749686653805858, 2.9949979469148218, 0.31642188798023552, 1.2005089463280847], [7.0660354700589316, 2.9881981198394589, 0.3437840226428448, 1.1775290186350424], [7.4094795980105825, 2.9729146254737744, 0.37238893504540599, 1.1646093421487571], [7.4319395986204615, 3.3793200765884617, 0.38089248055539254, 1.1602937215988161], [7.078860114925611, 3.3928205388475501, 0.35333749346401244, 1.1837826417655395], [6.7517984533166651, 3.4022837701153286, 0.32719853795582915, 1.1930743312721428], [6.4494160309940245, 3.397230218465074, 0.30242464799349522, 1.238711246164738], [6.1704868585408921, 3.3879250045551617, 0.27908434254056286, 1.255348472622059], [5.9141412434600564, 3.3665242674845972, 0.25723737270920677, 1.3052911292187819], [5.6788399462317116, 3.3389649751350925, 0.23690971924373941, 1.3385925718404832], [5.4640367039592093, 3.3007051299301553, 0.21818397889369803, 1.3982675951127672], [5.2667029809844816, 3.2610936294929127, 0.20127013983687655, 1.4201011551104139], [5.0887451611876688, 3.2068822718487855, 0.18603187072234267, 1.5176997640490331], [4.9265691520595762, 3.1481676274771102, 0.17247744026511713, 1.5693631141094921], [4.7794576380185614, 3.0829264798793199, 0.16092919220364221, 1.6394197881328278], [4.6481348224969512, 3.0081733174332062, 0.15110829617270433, 1.7394890983979268], [4.5325326328482145, 2.9243854384147445, 0.14277350917447737, 1.8491741012487657], [4.4312078031624136, 2.8340077352083375, 0.1357749989789068, 1.9503557077100673], [4.342670300854822, 2.7386951124448524, 0.13008991264860306, 2.0442331690741762], [4.2651021783472274, 2.6397819240257983, 0.1257005667154156, 2.1277619011429412], [4.2021600445133966, 2.5352186418886999, 0.12204586098200199, 2.250942718150303], [4.1527982717248886, 2.427110695025716, 0.11884406921573397, 2.3644701793516334], [4.1121526711657364, 2.3182492361784588, 0.11620190216673353, 2.4354558317491719], [4.0832876394146744, 2.207980101753515, 0.11398452554981421, 2.536772424679806], [4.0602627538314362, 2.0982744893257932, 0.11209579275897319, 2.5859205375769676], [4.0460013690800736, 1.9887557953950536, 0.11044334029411337, 2.6633065883032296], [4.0341302818719544, 1.8805215284610581, 0.10888332861481656, 2.683553766200359], [4.0280953903635268, 1.7733769288176295, 0.10731442190250613, 2.7365313886790013], [4.0225981008826945, 1.6679353974413262, 0.10558473720484438, 2.7407079345014842], [4.0186471463360451, 1.5643868713328173, 0.10362387418483381, 2.7546595923972297], [4.0143296352159226, 1.4631413786513907, 0.10133750880388354, 2.7501785147066915], [4.0143748139715942, 1.3643163806919203, 0.098825008286411539, 2.7932538360428856], [4.0044863344140929, 1.2690360481818013, 0.095792086266026713, 2.6893838913273642], [3.9965816677827348, 1.1769918361589744, 0.092383011000161955, 2.7071278536956136], [3.9819760271431406, 1.0893872250905228, 0.088813808717719267, 2.6275938880948959], [3.9650169963639481, 1.0059901977039571, 0.085103894751563339, 2.5921792240493025], [3.9431921313107843, 0.92730877968612879, 0.081652252117654783, 2.5222166305555351], [3.9190704686485276, 0.85211320586415995, 0.078969797593807484, 2.4823812028721362], [3.8854242021987613, 0.78262673813609884, 0.077203888783736296, 2.3418580636574844], [3.846146730588325, 0.71643856463077105, 0.076964888670613546, 2.2572285023744891], [3.7973325238704785, 0.65462363746439778, 0.07876491603548387, 2.1243782285481716], [3.7389927816395367, 0.59573833428782963, 0.082891522206851351, 2.012052442222056], [3.6717601052520177, 0.5360909677543596, 0.089877923366125281, 1.9476862666894645], [3.5882461773160141, 0.48236763251251835, 0.099301424505427421, 1.7936411762252422], [3.4891720356465807, 0.43265944247301258, 0.11084489029512719, 1.6870287509185984], [3.3758689099856287, 0.38049979440925197, 0.12473262272026837, 1.6534319623703835], [3.241931245707784, 0.33841831817314005, 0.14039283654950035, 1.5264214129817693], [3.0894526240973494, 0.2986931038003785, 0.15756846989540818, 1.4768645338711808], [2.9165157282037253, 0.2640836763813933, 0.17636604669718087, 1.4195184963411938], [2.7215978590380878, 0.2394365410176148, 0.19646999007918028, 1.3477816239008082], [2.5045422943062245, 0.22136659328357164, 0.21780643055745355, 1.3050591520829689], [2.2644285898815752, 0.20774261469660429, 0.24049990404377808, 1.2786792742793276], [2.0, 0.20000000000000001, 0.26454191960651186, 1.2512725370596907], [-8.0, 0.20000000000000001, 0.26454191960651186, 1.2512725370596907]]]

    #WLEN = 3.0
    #WLEN2 = 5.0
    #wall1 = [[-14.0, -0.2], [-4.0, -0.2], [-4.0 + WLEN*cos(pi/3), -0.2 - WLEN*sin(pi/3)]]
    #wall2 = [[-4.0 + WLEN2*cos(pi/3), 0.2 + WLEN2*sin(pi/3)], [-4.0, 0.2] ,[-14.0, 0.2]]
    #w1 = wall1[2]
    #w2 = wall2[0]
    
    #wall3 = [[w1[0] + 0.4*cos(pi/6), w1[1] + 0.4*sin(pi/6)], [0.4*cos(pi/6) - 4, 0.0], [w2[0] + 0.4*cos(pi/6), w2[1] - 0.4*sin(pi/6)], w2]
    #lp = wall3[0]
    #rp = wall3[2]
    
    #wall6 = [lp, [lp[0] + WLEN*cos(pi/6), lp[1] + WLEN*sin(pi/6)]]
    #wall6.append([wall6[1][0] + 0.4*cos(pi/3), wall6[1][1] - 0.4*sin(pi/3)])
    #wall6.append([wall6[2][0] - WLEN*cos(pi/6), wall6[2][1] - WLEN*sin(pi/6)])
    #wall6.append([wall6[3][0] + WLEN*cos(pi/3), wall6[3][1] - WLEN*sin(pi/3)])
    #wall6.append([wall6[4][0] - 0.4*cos(pi/6), wall6[4][1] - 0.4*sin(pi/6)])
    #wall6.append(w1)
    #wall6.reverse()
    #walls = [wall1, wall2, wall3, wall6]
        
        
    # triple-Y - multi-junction, test 1                 
    WLEN = 2.0
    WWID = 0.4
    wall1 = [[-11.0, -WWID/2.0], [-4.0, -WWID/2.0], [-4.0 + WLEN*cos(-pi/3), -WWID/2.0 + WLEN*sin(-pi/3)]]
    wall2 = [[-11.0, WWID/2.0], [-4.0, WWID/2.0], [-4.0 + WLEN*cos(pi/3), WWID/2.0 + WLEN*sin(pi/3)]]
    #wall2 = [[-4.0 + WLEN2*cos(pi/3), 0.2 + WLEN2*sin(pi/3)], [-4.0, 0.2] ,[-14.0, 0.2]]
    w1 = wall1[2]
    w2 = wall2[2]
    
    wall3 = [[w1[0] + WWID*cos(pi/6), w1[1] + WWID*sin(pi/6)], [WWID*cos(pi/6) - 4, 0.0], [w2[0] + WWID*cos(pi/6), w2[1] - 0.4*sin(pi/6)]]
    
    wall4 = [w1, [w1[0] + WLEN*cos(-2*pi/3), w1[1] + WLEN*sin(-2*pi/3)]]            
    wall5 = [w2, [w2[0] + WLEN*cos(2*pi/3), w2[1] + WLEN*sin(2*pi/3)]]
    
    wall4.append([wall4[-1][0] + WWID*cos(-2*pi/3 + pi/2), wall4[-1][1] + WWID*sin(-2*pi/3 + pi/2)])
    wall4.append([wall4[-1][0] + WLEN*cos(-2*pi/3 + pi), wall4[-1][1] + WLEN*sin(-2*pi/3 + pi)])
    wall4.append([wall4[-1][0] + WLEN*cos(-2*pi/3 + pi - pi/3), wall4[-1][1] + WLEN*sin(-2*pi/3 + pi - pi/3)])
    wall4.append([wall4[-1][0] + WWID*cos(-2*pi/3 + pi - pi/3 + pi/2.), wall4[-1][1] + WWID*sin(-2*pi/3 + pi - pi/3 + pi/2.)])
    wall4.append([wall4[-1][0] + WLEN*cos(-2*pi/3 + pi - pi/3 + pi), wall4[-1][1] + WLEN*sin(-2*pi/3 + pi - pi/3 + pi)])
    
    wall5.append([wall5[-1][0] + WWID*cos(2*pi/3 - pi/2), wall5[-1][1] + WWID*sin(2*pi/3 - pi/2)])
    wall5.append([wall5[-1][0] + WLEN*cos(2*pi/3 - pi), wall5[-1][1] + WLEN*sin(2*pi/3 - pi)])
    wall5.append([wall5[-1][0] + WLEN*cos(2*pi/3 - pi + pi/3), wall5[-1][1] + WLEN*sin(2*pi/3 - pi + pi/3)])
    wall5.append([wall5[-1][0] + WWID*cos(2*pi/3 - pi + pi/3 - pi/2.), wall5[-1][1] + WWID*sin(2*pi/3 - pi + pi/3 - pi/2.)])
    wall5.append([wall5[-1][0] + WLEN*cos(2*pi/3 - pi + pi/3 - pi), wall5[-1][1] + WLEN*sin(2*pi/3 - pi + pi/3 - pi)])
    
    wall5.reverse()
    
    wall6 = [[-11.0, WWID/2.0],[-11.0, -WWID/2.0]]
    
    walls = [wall1, wall2, wall3, wall4, wall5, wall6]
        
    for wall in walls:
        for i in range(len(wall)):
            p = copy(wall[i])
            p[0] += 6.0
            wall[i] = p


    if axes == 0:
    
        for wall in walls:
            xP = []
            yP = []
            for i in range(len(wall)):
                p = copy(wall[i])
                xP.append(p[0])
                yP.append(p[1])
    
            pylab.plot(xP,yP, linewidth=2, color = 'g')

    else:

        for wall in walls:
            xP = []
            yP = []
            for i in range(len(wall)):
                p = copy(wall[i])
                xP.append(p[0])
                yP.append(p[1])
    
            axes.plot(xP,yP, linewidth=2, color = 'g')
        
    
def draw_matches(match_pairs, offset, axes = 0):


    if axes == 0:
        for pair in match_pairs:
            a_off = dispOffset(pair[0], offset)
            b = pair[1]
    
            pylab.plot([a_off[0],b[0]], [a_off[1],b[1]])

    else:
        for pair in match_pairs:
            a_off = dispOffset(pair[0], offset)
            b = pair[1]
    
            axes.plot([a_off[0],b[0]], [a_off[1],b[1]])

def draw(a_pnts, b_pnts, filename, fileWrite = True):

    xP = []
    yP = []
    for a in a_pnts:
        xP.append(a[0])    
        yP.append(a[1])    
        #pylab.scatter([a[0]],[a[1]],linewidth=1, color=(0.6,0.6,1.0), edgecolors='none')
    if len(xP) > 0:
        pylab.plot(xP,yP, linewidth=1, color='b')
        pass
        
    xP = []
    yP = []
    for b in b_pnts:
        xP.append(b[0])    
        yP.append(b[1])    
        #pylab.scatter([b[0]],[b[1]],linewidth=1, color=(1.0,0.6,0.6), edgecolors='none')
    if len(xP) > 0:
        pylab.plot(xP,yP, linewidth=1, color='r')
        pass

    if fileWrite:
        pylab.xlim(-4.5,4.5)
        pylab.ylim(-4,4)
        pylab.savefig(filename)
        pylab.clf()

def save(filename):
    pylab.savefig(filename)
    pylab.clf()


def computeUnions(point_sets):

    currPoly = point_sets[0]

    for i in range(1,len(point_sets)):
        currPoly = computeUnion(currPoly,point_sets[i])

    return currPoly

def computeUnion(points1, points2):
    
    try:
        
        " remove any degenerate segments "
        nPoints1 = []
        nPoints2 = []
        
        for i in range(len(points1)-1):
            if not (points1[i][0] == points1[i+1][0] and points1[i][1] == points1[i+1][1]):
                nPoints1.append(points1[i])
        nPoints1.append(points1[-1])

        for i in range(len(points2)-1):
            if not (points2[i][0] == points2[i+1][0] and points2[i][1] == points2[i+1][1]):
                nPoints2.append(points2[i])
        nPoints2.append(points2[-1])

        numPoints1 = len(nPoints1)
        numPoints2 = len(nPoints2)
    
        #print numPoints1, numPoints2
    
        inputStr = str(numPoints1) + " " + str(numPoints2) + " "
        
        " Removed Gaussians because were causing self-intersections ."
        " not necessary anymore because polygons come from CGAL and are not degenerate. "
        for p in nPoints1:
            p2 = copy(p)
            #p2[0] += gauss(0.0,0.0001)
            #p2[1] += gauss(0.0,0.0001)
    
            inputStr += str(p2[0]) + " " + str(p2[1]) + " "
    
        for p in nPoints2:
            p2 = copy(p)
            #p2[0] += gauss(0.0,0.0001)
            #p2[1] += gauss(0.0,0.0001)
    
            inputStr += str(p2[0]) + " " + str(p2[1]) + " "
        
        inputStr += "\n"
        
        #print inputStr 
        
        #f = open("input2.txt", 'w')
        #f.write(inputStr)
        #f.close()
        
        if True:
            
            " start the subprocess "
            subProc = Popen(["./poly_union.exe"], stdin=PIPE, stdout=PIPE)
            
            " send input and receive output "
            sout, serr = subProc.communicate(inputStr)
        
            " convert string output to typed data "
            sArr = sout.split(" ")
        
            numVert = int(sArr[0])
            
            vertices = []
            for i in range(numVert):
                #print sArr[2*i+1], sArr[2*i+2]
                vertices.append([float(sArr[2*i + 1]), float(sArr[2*i + 2])])
        
        else:
            maxX = 0
            minX = 1e10
            maxY = 0
            minY = 1e10
            
            for p in nPoints1:
                if p[0] > maxX:
                    maxX = p[0]
                
                if p[0] < minX:
                    minX = p[0]
                
                if p[1] > maxY:
                    maxY = p[1]
                    
                if p[1] < minY:
                    minY = p[1]
            
            for p in nPoints2:
                if p[0] > maxX:
                    maxX = p[0]
                
                if p[0] < minX:
                    minX = p[0]
                
                if p[1] > maxY:
                    maxY = p[1]
                    
                if p[1] < minY:
                    minY = p[1]

            vertices = []
            vertices.append([maxX,maxY])
            vertices.append([maxX,minY])
            vertices.append([minX,minY])
            vertices.append([minX,maxY])

        return vertices
    except:
        print "polygon union failed with the following data:"
        print "sending input to poly_union:  ", inputStr
        print "rawOutput ="
        print sout
    
        print serr        
        raise

def doTest():



    f = open("icpInputSave_input.txt", 'r')        
    saveStr = f.read()
    f.close()
        
    saveStr = saveStr.replace('\r\n','\n')        
    exec(saveStr)

    uSet = [i*0.01 for i in range(100)]

    #globalSpline = SplineFit(globalPath, smooth=0.1)
    #medialSpline = SplineFit(medialPoints, smooth=0.1)

    #poses_1 = globalSpline.getUVecSet(uSet)
    #poses_2 = medialSpline.getUVecSet(uSet)

    
    """
    saveStr = ''
    saveStr += "initGuess = " + repr(initGuess) + "\n"
    saveStr += "globalPath = " + repr(globalPath) + "\n"
    saveStr += "medialPoints = " + repr(medialPoints) + "\n"
    saveStr += "poses_1 = " + repr(poses_1) + "\n"
    saveStr += "poses_2 = " + repr(poses_2) + "\n"

    f = open("cudaData.txt", 'w')
    f.write(saveStr)
    f.close()
    """


    #import nelminICP
    #param, cost = nelminICP.ICPmin(matchPairs, numPairs, initGuess, poses_1, poses_2, numPoses)

    #uSet = uSet[:20]
    
    args = []
    for i in range(len(uSet)):

        u1 = uSet[i]
        #print "doing u1 =", u1
        initGuess[0] = u1
        initGuess[1] = 0.5
        args.append(copy(initGuess))
        
    results = batchGlobalICP(globalPath, medialPoints, args)
    
    return
    
    print "initGuess:", initGuess
        
    #for u1 in uSet:
    for i in range(len(uSet)/2):

        u1 = uSet[2*i]
        #print "doing u1 =", u1
        initGuess[0] = u1
        initGuess[1] = 0.5
        globalOverlapICP_GPU2(initGuess, globalPath, medialPoints)
        #globalOverlapICP(initGuess, globalPath, medialPoints)
        #globalOverlapICP_GPU(initGuess, globalPath, medialPoints, poses_1, poses_2)

    #print "doing u1 =", u1
    #globalOverlapICP_GPU2(initGuess, globalPath, medialPoints, poses_1, poses_2)


def doTest2():


    f = open("costTest.txt", 'r')        
    saveStr = f.read()
    f.close()
        
    saveStr = saveStr.replace('\r\n','\n')        
    exec(saveStr)

    flatMatchPairs = []
    for pair in match_pairs:
        p1 = pair[0]
        p2 = pair[1]
        C1 = pair[2]
        C2 = pair[3]

        flatMatchPairs.append(p1[0])
        flatMatchPairs.append(p1[1])
        flatMatchPairs.append(C1[0][0])
        flatMatchPairs.append(C1[0][1])
        flatMatchPairs.append(C1[1][0])
        flatMatchPairs.append(C1[1][1])
        flatMatchPairs.append(p2[0])
        flatMatchPairs.append(p2[1])
        flatMatchPairs.append(C2[0][0])
        flatMatchPairs.append(C2[0][1])
        flatMatchPairs.append(C2[1][0])
        flatMatchPairs.append(C2[1][1])
    c_poses_1 = [item for sublist in poses_1 for item in sublist]
    c_poses_2 = [item for sublist in poses_2 for item in sublist]

    u1 = initGuess[0]
    currU = initGuess[1]
    currAng = initGuess[2]

    uHigh = currU + 0.2
    uLow = currU - 0.2

    #newParam, newCost = nelminICP.ICPmin(flatMatchPairs, len(match_pairs), initGuess, c_poses_1, c_poses_2, len(poses_1))

    for i in range(1):
        cost1, param, offset = nelminICP.ICPcost(flatMatchPairs, len(match_pairs), initGuess, c_poses_1, c_poses_2, len(poses_1))
        print "offset =", offset
        print "param =", param
        cost2 = medialOverlapCostFunc([currU, currAng], match_pairs, poses_1, poses_2, uHigh, uLow, u1)            
        #print cost1, cost2
        print cost1, cost2
        

    #globalOverlapICP_GPU2(initGuess, globalPath, medialPoints)
    #globalOverlapICP(initGuess, globalPath, medialPoints)
    #globalOverlapICP_GPU(initGuess, globalPath, medialPoints, poses_1, poses_2)


class ControlError(Exception):
    def __init__(self, value):
        self.value = value
        
    def __str__(self):
        return repr(self.value)
        

if __name__ == '__main__':


    numpy.set_printoptions(threshold=numpy.nan)

    #print knn_search()

    #print __num_processors()
    
    #exit()

    try:
        time1 = time.time()

        " while loop "        
        #cProfile.run('doTest()', 'test_prof')
        #doTest2()
        doTest()
        time2 = time.time()
        print time2 - time1
    
    except ControlError as inst:
        print inst.value

    except:
        traceback.print_exc()
        print "Exception:", sys.exc_info()[0]
        
    

    exit()

    inputStr = "323 175 0.782056 0.654747 0.746317 0.706635 0.702085 0.773295 0.674216 0.811786 0.636855 0.862563 0.577827 0.84605 0.543567 0.802175 0.49582 0.737985 0.457926 0.684164 0.429763 0.642696 0.384887 0.576467 0.356848 0.534852 0.328793 0.493153 0.283956 0.426898 0.256036 0.385406 0.22816 0.343921 0.183368 0.277636 0.153638 0.237144 0.105648 0.173136 0.0677733 0.119273 0.0397852 0.0777528 0.0263176 0.0578137 0.0135987 0.0611008 -0.0295399 0.0357618 -0.0598928 -0.0111072 -0.0905972 -0.0579506 -0.121016 -0.104847 -0.151364 -0.151644 -0.181771 -0.198579 -0.212376 -0.245375 -0.24257 -0.292198 -0.273265 -0.339256 -0.30371 -0.386283 -0.334194 -0.433067 -0.364533 -0.479823 -0.39506 -0.526632 -0.425435 -0.5737 -0.455891 -0.620513 -0.486375 -0.667436 -0.516476 -0.706214 -0.532903 -0.686855 -0.57036 -0.670711 -0.613448 -0.738116 -0.64888 -0.793543 -0.692269 -0.860755 -0.728033 -0.916154 -0.768024 -0.985442 -0.793122 -1.02881 -0.818296 -1.07228 -0.858302 -1.14156 -0.884855 -1.18384 -0.928024 -1.25119 -0.963537 -1.30657 -1.00673 -1.3739 -1.04062 -1.43014 -1.06582 -1.47371 -1.1058 -1.543 -1.13087 -1.58638 -1.15599 -1.62983 -1.1575 -1.63247 -1.17686 -1.65925 -1.20714 -1.70648 -1.24877 -1.7748 -1.28131 -1.83188 -1.29765 -1.86184 -1.29955 -1.86499 -1.31375 -1.88969 -1.34695 -1.94408 -1.38121 -2.0002 -1.39092 -2.01617 -1.41818 -2.05639 -1.44636 -2.09791 -1.49108 -2.16424 -1.51926 -2.20581 -1.54738 -2.24723 -1.59207 -2.31359 -1.62019 -2.35527 -1.64828 -2.39689 -1.69321 -2.46308 -1.72261 -2.50339 -1.77028 -2.56764 -1.80805 -2.62132 -1.83637 -2.66299 -1.88125 -2.72922 -1.90929 -2.77077 -1.9373 -2.8123 -1.99986 -2.86216 -2.06242 -2.91201 -2.12499 -2.96187 -2.18755 -3.01173 -2.25011 -3.06159 -2.31267 -3.11145 -2.37524 -3.16131 -2.42599 -3.20175 -2.42762 -3.21243 -2.43229 -3.21474 -2.504 -3.25021 -2.57571 -3.28568 -2.64742 -3.32115 -2.71913 -3.35661 -2.79083 -3.39208 -2.86254 -3.42755 -2.87831 -3.43535 -2.9253 -3.43372 -2.98099 -3.43349 -3.03695 -3.43319 -3.09278 -3.42989 -3.17252 -3.42347 -3.23807 -3.42122 -3.29407 -3.42095 -3.34996 -3.42083 -3.40583 -3.42057 -3.46168 -3.42026 -3.5175 -3.42003 -3.57351 -3.41985 -3.62938 -3.4196 -3.68545 -3.41938 -3.74128 -3.41917 -3.79723 -3.41892 -3.85313 -3.41867 -3.89148 -3.4297 -3.89873 -3.43928 -3.95947 -3.43567 -4.0153 -3.4369 -4.09516 -3.43229 -4.1608 -3.42851 -4.21679 -3.42965 -4.27264 -3.43075 -4.35251 -3.4262 -4.41826 -3.42246 -4.47385 -3.42332 -4.52982 -3.42437 -4.60969 -3.41979 -4.67559 -3.41602 -4.73121 -3.41707 -4.78723 -3.4178 -4.84316 -3.41922 -4.88737 -3.44257 -4.88737 -3.44257 -4.89015 -3.44238 -4.97015 -3.44187 -5.01499 -3.44168 -5.09499 -3.44136 -5.13999 -3.44114 -5.21999 -3.44075 -5.26516 -3.44279 -5.3152 -3.44537 -5.3952 -3.44492 -5.44008 -3.44467 -5.52008 -3.44426 -5.56516 -3.44403 -5.64516 -3.44366 -5.68997 -3.44344 -5.76997 -3.44302 -5.81526 -3.44278 -5.89516 -3.4468 -5.94495 -3.44927 -5.99548 -3.45177 -6.03537 -3.48146 -6.07342 -3.53214 -6.09594 -3.58336 -6.12648 -3.6573 -6.14148 -3.71626 -6.15172 -3.78604 -6.14205 -3.84119 -6.06205 -3.84151 -6.01715 -3.84169 -5.93715 -3.8421 -5.89206 -3.84232 -5.81207 -3.84273 -5.76699 -3.84295 -5.68785 -3.85468 -5.60871 -3.8664 -5.52958 -3.87812 -5.45234 -3.88956 -5.39736 -3.87987 -5.34216 -3.87004 -5.28712 -3.86027 -5.23209 -3.85062 -5.15219 -3.84654 -5.0723 -3.84247 -5.05231 -3.84145 -4.97241 -3.8376 -4.92197 -3.83724 -4.84197 -3.83753 -4.7969 -3.83777 -4.71691 -3.83823 -4.67183 -3.83843 -4.59183 -3.83879 -4.54711 -3.839 -4.46711 -3.83937 -4.43143 -3.83379 -4.39914 -3.83567 -4.34327 -3.83463 -4.28733 -3.83347 -4.20747 -3.8381 -4.14184 -3.84191 -4.08594 -3.84088 -4.03013 -3.83987 -3.97418 -3.83873 -3.89431 -3.84334 -3.82854 -3.84715 -3.77267 -3.84609 -3.71672 -3.84493 -3.66095 -3.84413 -3.58108 -3.84874 -3.51534 -3.85253 -3.45951 -3.85129 -3.4035 -3.85035 -3.34754 -3.84918 -3.29182 -3.84822 -3.23607 -3.8471 -3.15609 -3.84869 -3.07611 -3.85029 -2.99612 -3.85188 -2.91614 -3.85347 -2.83615 -3.85506 -2.77004 -3.85638 -2.7588 -3.85766 -2.67893 -3.86218 -2.6287 -3.86503 -2.57858 -3.86787 -2.49871 -3.87233 -2.44876 -3.8702 -2.39361 -3.86089 -2.31765 -3.83972 -2.31481 -3.83873 -2.31037 -3.83875 -2.25433 -3.83906 -2.19845 -3.83598 -2.13808 -3.83116 -2.10701 -3.83127 -2.06513 -3.81597 -2.06033 -3.81398 -1.99629 -3.83316 -1.94736 -3.84302 -1.86817 -3.83793 -1.82208 -3.81877 -1.76936 -3.7586 -1.74054 -3.69552 -1.72471 -3.6171 -1.70888 -3.53868 -1.69306 -3.46026 -1.67723 -3.38184 -1.66141 -3.30343 -1.64779 -3.23598 -1.63821 -3.21955 -1.59622 -3.15146 -1.57041 -3.10316 -1.54224 -3.05478 -1.50037 -2.98662 -1.49527 -2.97711 -1.48681 -2.9639 -1.44506 -2.89566 -1.41077 -2.83947 -1.36915 -2.77115 -1.33485 -2.7148 -1.29325 -2.64646 -1.2565 -2.5919 -1.22808 -2.54401 -1.20403 -2.50012 -1.16229 -2.43188 -1.12803 -2.37586 -1.08638 -2.30755 -1.05204 -2.25122 -1.01035 -2.18295 -0.976081 -2.12683 -0.943624 -2.08139 -0.905314 -2.01115 -0.881258 -1.96736 -0.864868 -1.93757 -0.852933 -1.91983 -0.808275 -1.85346 -0.763617 -1.78708 -0.718959 -1.72071 -0.674301 -1.65433 -0.645067 -1.60963 -0.606052 -1.54938 -0.575654 -1.50266 -0.545161 -1.45578 -0.514778 -1.40875 -0.484331 -1.36195 -0.435266 -1.29876 -0.401757 -1.2556 -0.371163 -1.20875 -0.345223 -1.16836 -0.321475 -1.15082 -0.276606 -1.08458 -0.248607 -1.04315 -0.220534 -1.00159 -0.175639 -0.935372 -0.1459 -0.894693 -0.0980938 -0.830549 -0.0603663 -0.776854 -0.032353 -0.735423 0.0123763 -0.669096 0.0404848 -0.627549 0.0686508 -0.585952 0.113469 -0.519685 0.141469 -0.478052 0.169376 -0.436498 0.214248 -0.370267 0.24237 -0.328761 0.270532 -0.287198 0.315167 -0.220807 0.343199 -0.179228 0.371258 -0.13764 0.419252 -0.0736363 0.458505 -0.0212898 0.503151 0.0450937 0.531337 0.086775 0.559558 0.128448 0.604252 0.194799 0.632253 0.236319 0.660345 0.277961 0.705247 0.344171 0.731318 0.382803 0.734514 0.385505 0.734732 0.387866 0.761258 0.427221 0.806076 0.493488 0.825424 0.534521 0.829394 0.590255 -3.84460529197 -3.82900135089 -3.76460540964 -3.82913856068 -3.71495033922 -3.83117007124 -3.63505189559 -3.83519979889 -3.56902055147 -3.83658476609 -3.51892676239 -3.83667500634 -3.43892697064 -3.83685754034 -3.38882622381 -3.83703303216 -3.33869715289 -3.83722408104 -3.2586973238 -3.83738944761 -3.20839117506 -3.8375155713 -3.1584067702 -3.8376465555 -3.0784070987 -3.8378758163 -3.02813688751 -3.8379083647 -2.97814373788 -3.83791225956 -2.89814407845 -3.83814569522 -2.84786524893 -3.83831920567 -2.79809505061 -3.83849777794 -2.71809539593 -3.83873283143 -2.66900492785 -3.84935510983 -2.59186031575 -3.87053785809 -2.51471570365 -3.89172060636 -2.43757109155 -3.91290335462 -2.36042647945 -3.93408610288 -2.28328186736 -3.95526885115 -2.21592011365 -3.95770077147 -2.16307593811 -3.91521348968 -2.12708025114 -3.84376903209 -2.09108456417 -3.77232457451 -2.07815332484 -3.73292888855 -2.05362342379 -3.65678242471 -2.03460052889 -3.59773097665 -2.00748227285 -3.5224674384 -1.99950994572 -3.46065816842 -1.99936675208 -3.41048352814 -2.02696047509 -3.36191822725 -2.05460686723 -3.31334170307 -2.0962044665 -3.28529373996 -2.13762567848 -3.25770106524 -2.19308664517 -3.25055352901 -2.26810166159 -3.2783518619 -2.34311667802 -3.30615019478 -2.41813169444 -3.33394852767 -2.49314671086 -3.36174686055 -2.56816172728 -3.38954519344 -2.6431767437 -3.41734352633 -2.71354302609 -3.44341917671 -2.79354244643 -3.44311463518 -2.84375602929 -3.4429899247 -2.89370887589 -3.44288273927 -2.97370853605 -3.44264955437 -3.02399973154 -3.44253246532 -3.07419109443 -3.44242318371 -3.15419098987 -3.4422938382 -3.20429594883 -3.44225834405 -3.25444544011 -3.44223432718 -3.33444506631 -3.44198976963 -3.38460735199 -3.44179940065 -3.43478855318 -3.44159954985 -3.51478852488 -3.44153225349 -3.56489678094 -3.44142140448 -3.61505675509 -3.44129306251 -3.69505639907 -3.44105439283 -3.74513095058 -3.43904044074 -3.82503250852 -3.43507294278 -3.88545835061 -3.44071462704 -3.95486234214 -3.45437801521 -4.03479214956 -3.45102750987 -4.11472195697 -3.44767700452 -4.19465176439 -3.44432649918 -4.27458157181 -3.44097599383 -4.35451137923 -3.43762548849 -4.43444118665 -3.43427498315 -4.4818553633 -3.43228747114 -4.56185471357 -3.43196504842 -4.60654200111 -3.43883254971 -4.68652435287 -3.43715224976 -4.76650670462 -3.4354719498 -4.84648905638 -3.43379164985 -4.92647140813 -3.43211134989 -5.00645375988 -3.43043104994 -5.08643611164 -3.42875074998 -5.16641846339 -3.42707045003 -5.24640081514 -3.42539015007 -5.3263831669 -3.42370985012 -5.38853899809 -3.42309471871 -5.43857377542 -3.42290750312 -5.51857374387 -3.4228364482 -5.56855078121 -3.42254646646 -5.61889810939 -3.42219261746 -5.69889787103 -3.4223879059 -5.74904069333 -3.42224832631 -5.79906185715 -3.42204281239 -5.87906120845 -3.42172064468 -5.92933357241 -3.42157766291 -5.97934683101 -3.42145061243 -6.05925094485 -3.41753492705 -6.12481820616 -3.41432181631 -6.20481719041 -3.41391867807 -6.2553648537 -3.41377550739 -6.30526332591 -3.41366321382 -6.38526331167 -3.41361548099 -6.43545380999 -3.41356517604 -6.48533692631 -3.41351002441 -6.56533689772 -3.41344238596 -6.6154585591 -3.41135902962 -6.6953508716 -3.40720952335 -6.7613936331 -3.40577722896 -6.81134372749 -3.40567775068 -6.89134351581 -3.40549371962 -6.94158754726 -3.40537322491 -6.99155977241 -3.40525213044 -7.07155968949 -3.40513694929 -7.12179525101 -3.40496005024 -7.17193373625 -3.40475678492 -7.24140803634 -3.41858973932 -7.26930020062 -3.45999080177 -7.29682111845 -3.50155641234 -7.32636971025 -3.55928081314 -7.35699748917 -3.63318574638 -7.38762526808 -3.70709067963 -7.40861204281 -3.75773183745 -7.36169827837 -3.79448404428 -7.28180082377 -3.79853333409 -7.21620573024 -3.79992253968 -7.1659977426 -3.80003803477 -7.08599802119 -3.80024915972 -7.03605066192 -3.80031331068 -6.98584523617 -3.80036087924 -6.90584606208 -3.80072439656 -6.85565374297 -3.80275372855 -6.77575209027 -3.80671931756 -6.70998286143 -3.80815431498 -6.65988319544 -3.8083497464 -6.57997688387 -3.81222032408 -6.51434672589 -3.8153993799 -6.4344547492 -3.81124341307 -6.38945352434 -3.80890245713 -6.30945370867 -3.80907419415 -6.2592619437 -3.81108788911 -6.17936117792 -3.81507130831 -6.11832188645 -3.83480137123 -6.04734670986 -3.87171375815 -5.97637153327 -3.90862614506 -5.90539635668 -3.94553853198 -5.83442118009 -3.9824509189 -5.7634460035 -4.01936330582 -5.6924708269 -4.05627569273 -5.62149565031 -4.09318807965 -5.55116256828 -4.1297665295 -5.47488661224 -4.15089592967 -5.42647739228 -4.12327464324 -5.35842361949 -4.0812189182 -5.2903698467 -4.03916319316 -5.2223160739 -3.99710746812 -5.15426230111 -3.95505174308 -5.08620852832 -3.91299601804 -5.01815475553 -3.87094029299 -4.95010098274 -3.82888456795 -4.91247492605 -3.80563249787 -4.84311373223 -3.79183448947 -4.76843992134 -3.78557932998 -4.68845043522 -3.78687629199 -4.6084609491 -3.788173254 -4.52847146298 -3.78947021601 -4.44848197686 -3.79076717802 -4.36849249075 -3.79206414003 -4.30242177961 -3.79313542085 -4.22713693455 -3.82019446933 -4.16544060725 -3.8282156323 -4.11519277984 -3.82849551347 -4.03519296498 -3.82866762414 -3.9850251929 -3.82869939508 -3.93478107715 -3.82871182562 -3.85478148949 -3.82896867842"
    points = [0.782056, 0.654747, 0.746317, 0.706635, 0.702085, 0.773295, 0.674216, 0.811786, 0.636855, 0.862563, 0.577827, 0.84605, 0.543567, 0.802175, 0.49582, 0.737985, 0.457926, 0.684164, 0.429763, 0.642696, 0.384887, 0.576467, 0.356848, 0.534852, 0.328793, 0.493153, 0.283956, 0.426898, 0.256036, 0.385406, 0.22816, 0.343921, 0.183368, 0.277636, 0.153638, 0.237144, 0.105648, 0.173136, 0.0677733, 0.119273, 0.0397852, 0.0777528, 0.0263176, 0.0578137, 0.0135987, 0.0611008, -0.0295399, 0.0357618, -0.0598928, -0.0111072, -0.0905972, -0.0579506, -0.121016, -0.104847, -0.151364, -0.151644, -0.181771, -0.198579, -0.212376, -0.245375, -0.24257, -0.292198, -0.273265, -0.339256, -0.30371, -0.386283, -0.334194, -0.433067, -0.364533, -0.479823, -0.39506, -0.526632, -0.425435, -0.5737, -0.455891, -0.620513, -0.486375, -0.667436, -0.516476, -0.706214, -0.532903, -0.686855, -0.57036, -0.670711, -0.613448, -0.738116, -0.64888, -0.793543, -0.692269, -0.860755, -0.728033, -0.916154, -0.768024, -0.985442, -0.793122, -1.02881, -0.818296, -1.07228, -0.858302, -1.14156, -0.884855, -1.18384, -0.928024, -1.25119, -0.963537, -1.30657, -1.00673, -1.3739, -1.04062, -1.43014, -1.06582, -1.47371, -1.1058, -1.543, -1.13087, -1.58638, -1.15599, -1.62983, -1.1575, -1.63247, -1.17686, -1.65925, -1.20714, -1.70648, -1.24877, -1.7748, -1.28131, -1.83188, -1.29765, -1.86184, -1.29955, -1.86499, -1.31375, -1.88969, -1.34695, -1.94408, -1.38121, -2.0002, -1.39092, -2.01617, -1.41818, -2.05639, -1.44636, -2.09791, -1.49108, -2.16424, -1.51926, -2.20581, -1.54738, -2.24723, -1.59207, -2.31359, -1.62019, -2.35527, -1.64828, -2.39689, -1.69321, -2.46308, -1.72261, -2.50339, -1.77028, -2.56764, -1.80805, -2.62132, -1.83637, -2.66299, -1.88125, -2.72922, -1.90929, -2.77077, -1.9373, -2.8123, -1.99986, -2.86216, -2.06242, -2.91201, -2.12499, -2.96187, -2.18755, -3.01173, -2.25011, -3.06159, -2.31267, -3.11145, -2.37524, -3.16131, -2.42599, -3.20175, -2.42762, -3.21243, -2.43229, -3.21474, -2.504, -3.25021, -2.57571, -3.28568, -2.64742, -3.32115, -2.71913, -3.35661, -2.79083, -3.39208, -2.86254, -3.42755, -2.87831, -3.43535, -2.9253, -3.43372, -2.98099, -3.43349, -3.03695, -3.43319, -3.09278, -3.42989, -3.17252, -3.42347, -3.23807, -3.42122, -3.29407, -3.42095, -3.34996, -3.42083, -3.40583, -3.42057, -3.46168, -3.42026, -3.5175, -3.42003, -3.57351, -3.41985, -3.62938, -3.4196, -3.68545, -3.41938, -3.74128, -3.41917, -3.79723, -3.41892, -3.85313, -3.41867, -3.89148, -3.4297, -3.89873, -3.43928, -3.95947, -3.43567, -4.0153, -3.4369, -4.09516, -3.43229, -4.1608, -3.42851, -4.21679, -3.42965, -4.27264, -3.43075, -4.35251, -3.4262, -4.41826, -3.42246, -4.47385, -3.42332, -4.52982, -3.42437, -4.60969, -3.41979, -4.67559, -3.41602, -4.73121, -3.41707, -4.78723, -3.4178, -4.84316, -3.41922, -4.88737, -3.44257, -4.88737, -3.44257, -4.89015, -3.44238, -4.97015, -3.44187, -5.01499, -3.44168, -5.09499, -3.44136, -5.13999, -3.44114, -5.21999, -3.44075, -5.26516, -3.44279, -5.3152, -3.44537, -5.3952, -3.44492, -5.44008, -3.44467, -5.52008, -3.44426, -5.56516, -3.44403, -5.64516, -3.44366, -5.68997, -3.44344, -5.76997, -3.44302, -5.81526, -3.44278, -5.89516, -3.4468, -5.94495, -3.44927, -5.99548, -3.45177, -6.03537, -3.48146, -6.07342, -3.53214, -6.09594, -3.58336, -6.12648, -3.6573, -6.14148, -3.71626, -6.15172, -3.78604, -6.14205, -3.84119, -6.06205, -3.84151, -6.01715, -3.84169, -5.93715, -3.8421, -5.89206, -3.84232, -5.81207, -3.84273, -5.76699, -3.84295, -5.68785, -3.85468, -5.60871, -3.8664, -5.52958, -3.87812, -5.45234, -3.88956, -5.39736, -3.87987, -5.34216, -3.87004, -5.28712, -3.86027, -5.23209, -3.85062, -5.15219, -3.84654, -5.0723, -3.84247, -5.05231, -3.84145, -4.97241, -3.8376, -4.92197, -3.83724, -4.84197, -3.83753, -4.7969, -3.83777, -4.71691, -3.83823, -4.67183, -3.83843, -4.59183, -3.83879, -4.54711, -3.839, -4.46711, -3.83937, -4.43143, -3.83379, -4.39914, -3.83567, -4.34327, -3.83463, -4.28733, -3.83347, -4.20747, -3.8381, -4.14184, -3.84191, -4.08594, -3.84088, -4.03013, -3.83987, -3.97418, -3.83873, -3.89431, -3.84334, -3.82854, -3.84715, -3.77267, -3.84609, -3.71672, -3.84493, -3.66095, -3.84413, -3.58108, -3.84874, -3.51534, -3.85253, -3.45951, -3.85129, -3.4035, -3.85035, -3.34754, -3.84918, -3.29182, -3.84822, -3.23607, -3.8471, -3.15609, -3.84869, -3.07611, -3.85029, -2.99612, -3.85188, -2.91614, -3.85347, -2.83615, -3.85506, -2.77004, -3.85638, -2.7588, -3.85766, -2.67893, -3.86218, -2.6287, -3.86503, -2.57858, -3.86787, -2.49871, -3.87233, -2.44876, -3.8702, -2.39361, -3.86089, -2.31765, -3.83972, -2.31481, -3.83873, -2.31037, -3.83875, -2.25433, -3.83906, -2.19845, -3.83598, -2.13808, -3.83116, -2.10701, -3.83127, -2.06513, -3.81597, -2.06033, -3.81398, -1.99629, -3.83316, -1.94736, -3.84302, -1.86817, -3.83793, -1.82208, -3.81877, -1.76936, -3.7586, -1.74054, -3.69552, -1.72471, -3.6171, -1.70888, -3.53868, -1.69306, -3.46026, -1.67723, -3.38184, -1.66141, -3.30343, -1.64779, -3.23598, -1.63821, -3.21955, -1.59622, -3.15146, -1.57041, -3.10316, -1.54224, -3.05478, -1.50037, -2.98662, -1.49527, -2.97711, -1.48681, -2.9639, -1.44506, -2.89566, -1.41077, -2.83947, -1.36915, -2.77115, -1.33485, -2.7148, -1.29325, -2.64646, -1.2565, -2.5919, -1.22808, -2.54401, -1.20403, -2.50012, -1.16229, -2.43188, -1.12803, -2.37586, -1.08638, -2.30755, -1.05204, -2.25122, -1.01035, -2.18295, -0.976081, -2.12683, -0.943624, -2.08139, -0.905314, -2.01115, -0.881258, -1.96736, -0.864868, -1.93757, -0.852933, -1.91983, -0.808275, -1.85346, -0.763617, -1.78708, -0.718959, -1.72071, -0.674301, -1.65433, -0.645067, -1.60963, -0.606052, -1.54938, -0.575654, -1.50266, -0.545161, -1.45578, -0.514778, -1.40875, -0.484331, -1.36195, -0.435266, -1.29876, -0.401757, -1.2556, -0.371163, -1.20875, -0.345223, -1.16836, -0.321475, -1.15082, -0.276606, -1.08458, -0.248607, -1.04315, -0.220534, -1.00159, -0.175639, -0.935372, -0.1459, -0.894693, -0.0980938, -0.830549, -0.0603663, -0.776854, -0.032353, -0.735423, 0.0123763, -0.669096, 0.0404848, -0.627549, 0.0686508, -0.585952, 0.113469, -0.519685, 0.141469, -0.478052, 0.169376, -0.436498, 0.214248, -0.370267, 0.24237, -0.328761, 0.270532, -0.287198, 0.315167, -0.220807, 0.343199, -0.179228, 0.371258, -0.13764, 0.419252, -0.0736363, 0.458505, -0.0212898, 0.503151, 0.0450937, 0.531337, 0.086775, 0.559558, 0.128448, 0.604252, 0.194799, 0.632253, 0.236319, 0.660345, 0.277961, 0.705247, 0.344171, 0.731318, 0.382803, 0.734514, 0.385505, 0.734732, 0.387866, 0.761258, 0.427221, 0.806076, 0.493488, 0.825424, 0.534521, 0.829394, 0.590255, -3.84460529197, -3.82900135089, -3.76460540964, -3.82913856068, -3.71495033922, -3.83117007124, -3.63505189559, -3.83519979889, -3.56902055147, -3.83658476609, -3.51892676239, -3.83667500634, -3.43892697064, -3.83685754034, -3.38882622381, -3.83703303216, -3.33869715289, -3.83722408104, -3.2586973238, -3.83738944761, -3.20839117506, -3.8375155713, -3.1584067702, -3.8376465555, -3.0784070987, -3.8378758163, -3.02813688751, -3.8379083647, -2.97814373788, -3.83791225956, -2.89814407845, -3.83814569522, -2.84786524893, -3.83831920567, -2.79809505061, -3.83849777794, -2.71809539593, -3.83873283143, -2.66900492785, -3.84935510983, -2.59186031575, -3.87053785809, -2.51471570365, -3.89172060636, -2.43757109155, -3.91290335462, -2.36042647945, -3.93408610288, -2.28328186736, -3.95526885115, -2.21592011365, -3.95770077147, -2.16307593811, -3.91521348968, -2.12708025114, -3.84376903209, -2.09108456417, -3.77232457451, -2.07815332484, -3.73292888855, -2.05362342379, -3.65678242471, -2.03460052889, -3.59773097665, -2.00748227285, -3.5224674384, -1.99950994572, -3.46065816842, -1.99936675208, -3.41048352814, -2.02696047509, -3.36191822725, -2.05460686723, -3.31334170307, -2.0962044665, -3.28529373996, -2.13762567848, -3.25770106524, -2.19308664517, -3.25055352901, -2.26810166159, -3.2783518619, -2.34311667802, -3.30615019478, -2.41813169444, -3.33394852767, -2.49314671086, -3.36174686055, -2.56816172728, -3.38954519344, -2.6431767437, -3.41734352633, -2.71354302609, -3.44341917671, -2.79354244643, -3.44311463518, -2.84375602929, -3.4429899247, -2.89370887589, -3.44288273927, -2.97370853605, -3.44264955437, -3.02399973154, -3.44253246532, -3.07419109443, -3.44242318371, -3.15419098987, -3.4422938382, -3.20429594883, -3.44225834405, -3.25444544011, -3.44223432718, -3.33444506631, -3.44198976963, -3.38460735199, -3.44179940065, -3.43478855318, -3.44159954985, -3.51478852488, -3.44153225349, -3.56489678094, -3.44142140448, -3.61505675509, -3.44129306251, -3.69505639907, -3.44105439283, -3.74513095058, -3.43904044074, -3.82503250852, -3.43507294278, -3.88545835061, -3.44071462704, -3.95486234214, -3.45437801521, -4.03479214956, -3.45102750987, -4.11472195697, -3.44767700452, -4.19465176439, -3.44432649918, -4.27458157181, -3.44097599383, -4.35451137923, -3.43762548849, -4.43444118665, -3.43427498315, -4.4818553633, -3.43228747114, -4.56185471357, -3.43196504842, -4.60654200111, -3.43883254971, -4.68652435287, -3.43715224976, -4.76650670462, -3.4354719498, -4.84648905638, -3.43379164985, -4.92647140813, -3.43211134989, -5.00645375988, -3.43043104994, -5.08643611164, -3.42875074998, -5.16641846339, -3.42707045003, -5.24640081514, -3.42539015007, -5.3263831669, -3.42370985012, -5.38853899809, -3.42309471871, -5.43857377542, -3.42290750312, -5.51857374387, -3.4228364482, -5.56855078121, -3.42254646646, -5.61889810939, -3.42219261746, -5.69889787103, -3.4223879059, -5.74904069333, -3.42224832631, -5.79906185715, -3.42204281239, -5.87906120845, -3.42172064468, -5.92933357241, -3.42157766291, -5.97934683101, -3.42145061243, -6.05925094485, -3.41753492705, -6.12481820616, -3.41432181631, -6.20481719041, -3.41391867807, -6.2553648537, -3.41377550739, -6.30526332591, -3.41366321382, -6.38526331167, -3.41361548099, -6.43545380999, -3.41356517604, -6.48533692631, -3.41351002441, -6.56533689772, -3.41344238596, -6.6154585591, -3.41135902962, -6.6953508716, -3.40720952335, -6.7613936331, -3.40577722896, -6.81134372749, -3.40567775068, -6.89134351581, -3.40549371962, -6.94158754726, -3.40537322491, -6.99155977241, -3.40525213044, -7.07155968949, -3.40513694929, -7.12179525101, -3.40496005024, -7.17193373625, -3.40475678492, -7.24140803634, -3.41858973932, -7.26930020062, -3.45999080177, -7.29682111845, -3.50155641234, -7.32636971025, -3.55928081314, -7.35699748917, -3.63318574638, -7.38762526808, -3.70709067963, -7.40861204281, -3.75773183745, -7.36169827837, -3.79448404428, -7.28180082377, -3.79853333409, -7.21620573024, -3.79992253968, -7.1659977426, -3.80003803477, -7.08599802119, -3.80024915972, -7.03605066192, -3.80031331068, -6.98584523617, -3.80036087924, -6.90584606208, -3.80072439656, -6.85565374297, -3.80275372855, -6.77575209027, -3.80671931756, -6.70998286143, -3.80815431498, -6.65988319544, -3.8083497464, -6.57997688387, -3.81222032408, -6.51434672589, -3.8153993799, -6.4344547492, -3.81124341307, -6.38945352434, -3.80890245713, -6.30945370867, -3.80907419415, -6.2592619437, -3.81108788911, -6.17936117792, -3.81507130831, -6.11832188645, -3.83480137123, -6.04734670986, -3.87171375815, -5.97637153327, -3.90862614506, -5.90539635668, -3.94553853198, -5.83442118009, -3.9824509189, -5.7634460035, -4.01936330582, -5.6924708269, -4.05627569273, -5.62149565031, -4.09318807965, -5.55116256828, -4.1297665295, -5.47488661224, -4.15089592967, -5.42647739228, -4.12327464324, -5.35842361949, -4.0812189182, -5.2903698467, -4.03916319316, -5.2223160739, -3.99710746812, -5.15426230111, -3.95505174308, -5.08620852832, -3.91299601804, -5.01815475553, -3.87094029299, -4.95010098274, -3.82888456795, -4.91247492605, -3.80563249787, -4.84311373223, -3.79183448947, -4.76843992134, -3.78557932998, -4.68845043522, -3.78687629199, -4.6084609491, -3.788173254, -4.52847146298, -3.78947021601, -4.44848197686, -3.79076717802, -4.36849249075, -3.79206414003, -4.30242177961, -3.79313542085, -4.22713693455, -3.82019446933, -4.16544060725, -3.8282156323, -4.11519277984, -3.82849551347, -4.03519296498, -3.82866762414, -3.9850251929, -3.82869939508, -3.93478107715, -3.82871182562, -3.85478148949, -3.82896867842]
    points1 = points[:323*2]
    points2 = points[323*2:]
    
    nPoints1 = []
    nPoints2 = []
    
    print len(points1)
    print len(points2)


    for i in range(len(points1)/2-1):
        if False and points1[2*i] == points1[2*i+2] and points1[2*i+1] == points1[2*i+3]:
            print "points1 ", i
            print points1[2*i], points1[2*i+2], points1[2*i+1], points1[2*i+3]
        else:
            nPoints1.append([points1[2*i], points1[2*i+1]])

    nPoints1.append([points1[-2], points1[-1]])


    for i in range(len(points2)/2-1):
        if False and points2[2*i] == points2[2*i+2] and points2[2*i+1] == points2[2*i+3]:
            print "points2 ", i
        else:
            nPoints2.append([points2[2*i],points2[2*i+1]])

    nPoints2.append([points2[-2], points2[-1]])
            
    print len(nPoints1)
    print len(nPoints2)
    
    xP = []
    yP = []
    for i in range(len(points1)/2):
        
        xP.append(points1[i*2])
        yP.append(points1[i*2+1])

    pylab.plot(xP,yP)

    xP = []
    yP = []
    for i in range(len(points2)/2):
        xP.append(points2[i*2])
        yP.append(points2[i*2+1])

    pylab.plot(xP,yP)
    pylab.show()
    
    #computeUnion(nPoints1, nPoints2)
    
    exit()
    
    # TUNE ME:  threshold cost difference between iterations to determine if converged
    costThresh = 0.004

    # TUNE ME:   minimum match distance before point is discarded from consideration
    minMatchDist = 2.0

    # plot the best fit at each iteration of the algorithm?
    plotIteration = True

    # initial guess for x, y, theta parameters
    offset = [0.0,0.0,-math.pi/4]

    # sample data
    a_data = [[1.0,1.0],[1.1,1.1],[1.2,1.2],[1.3,1.31],[1.4,1.4],[1.51,1.5],[1.6,1.6]]
    b_data = [[0.3,1.0],[0.3,1.1],[0.3,1.2],[0.31,1.3],[0.3,1.4],[0.3,1.5],[0.3,1.6]]

    # treat the points with the point-to-line constraint
    addPointToLineCovariance(a_data, high_var=1.0, low_var=0.001)
    addPointToLineCovariance(b_data, high_var=1.0, low_var=0.001)

    #addDistanceFromOriginCovariance(a_data, tan_var=0.1, perp_var=0.01)
    #addDistanceFromOriginCovariance(b_data, tan_var=0.1, perp_var=0.01)

    # plot the data without A transformed, plot 997
    draw(a_data, b_data, "rawData.png")

    # transform the points in A by 'offset'
    a_trans = []
    for p in a_data:
        a_trans.append(dispOffset(p, offset))

    # plot the data with A transformed, plot 998
    draw(a_trans, b_data, "initialGuess.png") 

    # run generalized ICP (a plot is made for each iteration of the algorithm)
    offset = gen_ICP(offset, a_data, b_data, costThresh, minMatchDist, plotIteration)

    # transform the points of A with parameters determined by algorithm
    a_trans = []
    for p in a_data:
        a_trans.append(dispPoint(p, offset))

    # plot the final result, plot 999
    draw(a_trans, b_data, "finalOutput.png") 
