
import nelminICP
import time

f = open('cudaData3.txt' , 'r')
saveStr = f.read()
f.close()
exec(saveStr)


param, cost = nelminICP.ICPmin(matchPairs, numPairs, initGuess, initGuess[1]+0.2, initGuess[1]-0.2, poses_1, poses_2, numPoses)


print param
print cost

param, cost = nelminICP.ICPmin(matchPairs, numPairs, initGuess, initGuess[1]+0.2, initGuess[1]-0.2, poses_1, poses_2, numPoses)

print param
print cost



f = open('icpInput_0977.txt' , 'r')
saveStr = f.read()
f.close()
exec(saveStr)

time1 = time.time()
newParam, newCost = nelminICP.ICPmin(flatMatchPairs, numPairs, [u1,currU,currAng], uHigh, uLow, c_poses_1, c_poses_2, numPoses)

time2 = time.time()
print "time:", time2 - time1


print newParam
print newCost
