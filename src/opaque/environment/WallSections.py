#!/usr/bin/python

import ogre.renderer.OGRE as ogre
import ogre.physics.OgreOde as OgreOde
from math import *

class WallSections:

	def __init__(self, world, lpnts, mapFile = ""):

		self._world = world
		self._mgr = world.getSceneManager()
		self._space = world.getDefaultSpace()

		self.trimeshes =  []

		self.npnts = []

		for pnts in lpnts:
			self.npnts.append(pnts)

		self.setupMyWorld()


	def setupMyWorld(self):

		# create the static geometry
		s = self._mgr.createStaticGeometry("corridor")

		# region dimensions (FIXME: what does this mean?)
		s.setRegionDimensions((160.0, 100.0, 160.0))

		## Set the region origin so the center is at 0 world
		s.setOrigin(ogre.Vector3().ZERO)


		for j, points in enumerate(self.npnts):

			# start with a ManualObject for the mesh
			# create the wall1 mesh
			mo1_int = ogre.ManualObject("WallObject_Interior" + str(j))
			mo1_int.begin("WallMaterial_Int", ogre.RenderOperation.OT_TRIANGLE_LIST)

			mo1_ext = ogre.ManualObject("WallObject_Exterior" + str(j))
			mo1_ext.begin("WallMaterial_Ext", ogre.RenderOperation.OT_TRIANGLE_LIST)

			mo1_top = ogre.ManualObject("WallObject_Top" + str(j))
			mo1_top.begin("WallMaterial_Top", ogre.RenderOperation.OT_TRIANGLE_LIST)

			topPoints = []
			for i in range(len(points)-1):
				vec = [points[i][0]-points[i+1][0], points[i][1]-points[i+1][1]]
				vecMag = sqrt(vec[0]**2 + vec[1]**2)
				vec[0] /= vecMag
				vec[1] /= vecMag
				
				orthoVec = [vec[0]*cos(pi/2) - vec[1]*sin(pi/2), vec[0]*sin(pi/2) + vec[1]*cos(pi/2)]
				
				topPoints.append(ogre.Vector3(points[i][0], 0.2, points[i][1]))
				topPoints.append(ogre.Vector3(points[i][0] + 0.05*orthoVec[0], 0.2, points[i][1] + 0.05*orthoVec[1]))
				topPoints.append(ogre.Vector3(points[i+1][0], 0.2, points[i+1][1]))
				topPoints.append(ogre.Vector3(points[i+1][0] + 0.05*orthoVec[0], 0.2, points[i+1][1] + 0.05*orthoVec[1]))
			
			for p in topPoints:
				mo1_top.position(p[0], p[1], p[2])

			for i in range(len(points) - 1):
				mo1_top.triangle(4*i+2, 4*i+1, 4*i)
				mo1_top.triangle(4*i+1, 4*i+2, 4*i+3)
				if i < len(points)-2:
					mo1_top.triangle(4*i+2, 4*i+5, 4*i+3)
					
			wall1Points = []
			for i in range(len(points)):
				wall1Points.append(ogre.Vector3(points[i][0], -100.0, points[i][1]))
				wall1Points.append(ogre.Vector3(points[i][0], 1.5, points[i][1]))

			# build the triangle meshes of the wall
			for p in wall1Points:
				#mo1_int.position(p[0],p[1],p[2])
				#mo1_ext.position(p[0],p[1],p[2])
				if p[1] < 0.0:
					mo1_int.position(p[0],0.0,p[2])
					mo1_ext.position(p[0],0.0,p[2])
				else:
					mo1_int.position(p[0],0.2,p[2])
					mo1_ext.position(p[0],0.2,p[2])

			# indices for ODE geometry
			indices1 = []
			indices2 = []

			for i in range(len(points) - 1):
				indices1 += [2*i+2, 2*i+1, 2*i]
				indices1 += [2*i+1, 2*i+2, 2*i+3]
				mo1_int.triangle(2*i+2, 2*i+1, 2*i)
				mo1_int.triangle(2*i+1, 2*i+2, 2*i+3)

				indices2 += [2*i, 2*i+1, 2*i+2]
				indices2 += [2*i+3, 2*i+2, 2*i+1]
				mo1_ext.triangle(2*i, 2*i+1, 2*i+2)
				mo1_ext.triangle(2*i+3, 2*i+2, 2*i+1)


			mo1_int.end()
			mo1_ext.end()
			mo1_top.end()

			# convert to mesh
			mo1_int.convertToMesh("WallMesh_Interior" + str(j))
			mo1_ext.convertToMesh("WallMesh_Exterior" + str(j))
			mo1_top.convertToMesh("WallMesh_Top" + str(j))
			# create Ogre entity
			print "creating entity " + str(j)
			entity1_int = self._mgr.createEntity("WallEntity_Interior" + str(j), "WallMesh_Interior" + str(j))
			entity1_ext = self._mgr.createEntity("WallEntity_Exterior" + str(j), "WallMesh_Exterior" + str(j))
			entity1_top = self._mgr.createEntity("WallEntity_Top" + str(j), "WallMesh_Top" + str(j))

			# add ODE geometry to the entity
			entity1_int.setCastShadows(False)
			entity1_ext.setCastShadows(False)
			entity1_top.setCastShadows(False)

			# create the ODE geometry
			self.trimeshes.append(OgreOde.makeTriangleMeshGeometry(wall1Points, len(wall1Points), indices1, len(indices1), self._world, self._space))
			#self.trimeshes.append(OgreOde.makeTriangleMeshGeometry(wall1Points, len(wall1Points), indices2, len(indices2), self._world, self._space))

			# add the entity to the Ogre static geometry
			s.addEntity(entity1_int, ogre.Vector3(0.0,0.0,0.0))
			s.addEntity(entity1_ext, ogre.Vector3(0.0,0.0,0.0))
			s.addEntity(entity1_top, ogre.Vector3(0.0,0.0,0.0))


		# now build the StaticGeometry so it will be rendered
		s.build()


