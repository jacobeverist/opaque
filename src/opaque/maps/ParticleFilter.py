import random

class ParticleFilter:

    def __init__(self, initPose = [0.2,0.1], numParticles = 0):
    
        self.numParticles = 20

        " pathID, distMean, distVar "
        self.initPose = [0, 0.1, 0.02]


        self.particles = []
        
        " each particle:  [pathID, distMean, distVar] "
        pathID = self.initPose[0]
        distMean = self.initPose[1]
        distVar = self.initPose[2]
        
        for i in range(numParticles):
            self.particles.append([pathID, random.gauss(distMean,distVar), 0.02])

        
        

    def move(self, moveDist = 0.0):
        
        newParticles = []
        
        
        for part in self.particles:    
            pathID = part[0]
            distMean = part[1]
            distVar = part[2]
            newPart = [pathID, random.gauss(distMean + moveDist, distVar), distVar + 0.01]
            newParticles.append(newPart)
        



    
        " move forward or backward "
        
        


