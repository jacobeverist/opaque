
from numpy import *
from functions import normalizeAngle
import graph

def mahab_dist(c1, c2, r1, r2, E):

	" distance between centroids "
	d = c2 - c1
	
	
	s = max(0, sqrt(d[0,0]**2 + d[1,0]**2) - r1 - r2) * d / sqrt(d[0,0]**2 + d[1,0]**2)
	E_inv = linalg.inv(E)
	dist = s.T * E_inv * s
	return dist[0,0]

def dijkstra_proj(initNode, numNodes, poseGraph):
	
	identTransform = matrix([[0.], [0.], [0.]],dtype=float)
	zeroCov = matrix([[0.,0.,0.],[0.,0.,0.],[0.,0.,0.]],dtype=float)
	
	paths = {initNode : [identTransform, zeroCov] }
	
	
	visited = [False for i in range(numNodes)]
	distances = [Inf for i in range(numNodes)]
	optimals = {}
	
	" initial node 0 "
	distances[initNode] = 0.0
	visited[initNode] = True
	
	" for each edge out of 0, add to path "
	neighbors = poseGraph.neighbors(initNode)
	incidents = poseGraph.incidents(initNode)

	for neigh in neighbors:
		if not visited[neigh]:
			transform, covE = poseGraph.get_edge_attributes(initNode,neigh)
			
			dist = linalg.det(covE)
			if dist < distances[neigh]:
				paths[neigh] = [transform, covE]
				distances[neigh] = dist


	for incid in incidents:
		if not visited[incid]:
			transform, covE = poseGraph.get_edge_attributes(incid, initNode)
			
			xA = transform[0,0]
			yA = transform[1,0]
			pA = transform[2,0]
			
			x1 = -cos(pA)*xA - sin(pA)*yA
			y1 = sin(pA)*xA - cos(pA)*yA
			p1 = normalizeAngle(-pA)

			newTransform = matrix([[x1], [y1], [p1]],dtype=float)

			paths[incid] = [newTransform, covE]
			dist = linalg.det(covE)
			if dist < distances[incid]:
				paths[incid] = [newTransform, covE]
				distances[incid] = dist


	" while some nodes unvisited "
	while visited.count(False) > 0:

		" find the minimum uncertainty path p "

		minDist = Inf
		minDest = -1

		for i in range(numNodes):
			if not visited[i]:
				if distances[i] < minDist:
					minDist = distances[i]
					minDest = i

		
		" an unvisited node is unreachable "
		if minDest == -1:
			break
		else:
			dest = minDest
			minUnc = Inf
			minTrans = 0
			minCov = 0
			
			p = paths[dest]
			uncertainty = linalg.det(p[1]) 
			if uncertainty < minUnc:
				minUnc = uncertainty
				minTrans = p[0]
				minCov = p[1]
					
			" mark this as visited and record the optimal path "
			visited[dest] = True
			optimals[dest] = [minTrans, minCov]
		
			" for all edges leaving dest, add the composed path "

			" for each edge out of dest, add to path "
			neighbors = poseGraph.neighbors(dest)
			incidents = poseGraph.incidents(dest)

			for neigh in neighbors:
				if not visited[neigh]:
					#print dest, neigh
					transform, covE = poseGraph.get_edge_attributes(dest,neigh)

					T_b_a = optimals[dest][0]
					T_c_b = transform

					
					E_a_b = optimals[dest][1]
					E_b_c = covE
					
					x1 = T_b_a[0,0]
					y1 = T_b_a[1,0]
					p1 = T_b_a[2,0]
					
					x2 = T_c_b[0,0]
					y2 = T_c_b[1,0]
					p2 = T_c_b[2,0]
					
					T_c_a = matrix([[x1 + x2*cos(p1) - y2*sin(p1)],
									[y1 + x2*sin(p1) + y2*cos(p1)],
									[p1+p2]],dtype=float)
					
					J1 = matrix([[1,0,-x2*sin(p1) - y2*cos(p1)],[0,1,x2*cos(p1) - y2*sin(p1)],[0,0,1]],dtype=float)
					J2 = matrix([[cos(p1), -sin(p1), 0], [sin(p1), cos(p1), 0], [0, 0, 1]],dtype=float)
					
					E_a_c = J1 * E_a_b * J1.T + J2 * E_b_c * J2.T

					dist = linalg.det(E_a_c)
					if dist < distances[neigh]:
						paths[neigh] = [T_c_a, E_a_c]
						distances[neigh] = dist
	
			for incid in incidents:
				if not visited[incid]:
					transform, covE = poseGraph.get_edge_attributes(incid, dest)

					T_b_a = optimals[dest][0]
					T_c_b = transform

					E_a_b = optimals[dest][1]
					E_b_c = covE
					
					x1 = T_b_a[0,0]
					y1 = T_b_a[1,0]
					p1 = T_b_a[2,0]

					xA = T_c_b[0,0]
					yA = T_c_b[1,0]
					pA = T_c_b[2,0]
					
					x2 = -cos(pA)*xA - sin(pA)*yA
					y2 = sin(pA)*xA - cos(pA)*yA
					p2 = -pA

					T_c_a = matrix([[x1 + x2*cos(p1) - y2*sin(p1)],
									[y1 + x2*sin(p1) + y2*cos(p1)],
									[p1+p2]],dtype=float)
			
					J1 = matrix([[1,0,-x2*sin(p1) - y2*cos(p1)],[0,1,x2*cos(p1) - y2*sin(p1)],[0,0,1]],dtype=float)
					J2 = matrix([[cos(p1), -sin(p1), 0], [sin(p1), cos(p1), 0], [0, 0, 1]],dtype=float)
	
					E_a_c = J1 * E_a_b * J1.T + J2 * E_b_c * J2.T

					dist = linalg.det(E_a_c)
					if dist < distances[incid]:
						paths[incid] = [T_c_a, E_a_c]
						distances[incid] = dist
		
			" compute the set of pairwise poses to consider using djikstra projection "
	
	for key, value in paths.iteritems():
		
		offset = value[0]
		covE = value[1]
		
	return paths


if __name__ == '__main__':

	fname = "../../testData/clustering/motion_graph_16.txt"
	#fname = "../../graph_out.txt"
	fp = open(fname)
	
	nodes, edges, edge_attrs = eval(fp.read())	

	newGraph = graph.digraph()
	
	for i in range(len(nodes)):
		newGraph.add_node(nodes[i])
		
	for i in range(len(edges)):
		newGraph.add_edge(edges[i][0], edges[i][1], attrs = edge_attrs[i])

	numNodes = len(nodes)
	paths = []
	for i in range(numNodes):
		paths.append(dijkstra_proj(i, numNodes, newGraph))	
	
	print paths[2][10]
	
	#nodes = eval(fp.read())
	#edges = eval(fp.read())
	#edge_attrs = eval(fp.read())
	
	#nodes = self.poseGraph.nodes()
	#edges = self.poseGraph.edges()
		
	#edge_attrs = []
	#for i in range(len(edges)):
	#	edge_attrs.append( self.poseGraph.get_edge_attributes(edges[i][0],edges[i][1]) )
	
	