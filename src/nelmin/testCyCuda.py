
import nelminICP


f = open('cudaData3.txt' , 'r')
saveStr = f.read()
f.close()

exec(saveStr)

param, cost = nelminICP.ICPmin(matchPairs, numPairs, initGuess, poses_1, poses_2, numPoses)


print param
print cost

param, cost = nelminICP.ICPmin(matchPairs, numPairs, initGuess, poses_1, poses_2, numPoses)

print param
print cost
