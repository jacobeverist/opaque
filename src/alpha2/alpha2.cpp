/*! \file simple_join_intersect.cpp
 * Computing the union and the intersection of two simple polygons.
 */

#include <iostream>
#include <fstream>
#include "cgal_types.h"
#include <CGAL/Unique_hash_map.h>

int main ()
{

	//Alpha_shape A(1000000.0, Alpha_shape::REGULARIZED);
	Alpha_shape A;
	CGALPointlist L;

	// read points from stdin
	//double points[20] = {1.0,0.0, 0.0,1.0, -1.0, -0.5, 0.0, 0.0 };

	// read number of points
	int numPoints = 0;
	std::cin >> numPoints;

	double radius = 0.0;
	std::cin >> radius;

	// read points from stdin
	for ( int i = 0 ; i < numPoints ; i++ ) {
		double x, y;

		std::cin >> x;
		std::cin >> y;

		//std::cout << x << " " << y << std::endl;
		L.push_back(Point_2(x,y));
	}


	//for ( int i = 0 ; i < 4120 ; i++ ) {
	//	L.push_back(Point_2(points[i*2],points[i*2+1]));
	//}

	// input the points for the shape
	A.make_alpha_shape(L.begin(), L.end());

	// find the optimal alpha for a single connected component
	Alpha_iterator at = A.find_optimal_alpha(1);

	// set the optimal alpha
	double optAlpha = *at;
	//A.set_alpha(0.05);
	//A.set_alpha(0.5);
	//A.set_alpha(2.0);
	A.set_alpha(radius);

	// the number of solid components
	int n = A.number_of_solid_components();

	int numVertices = 0;

	CGAL::Unique_hash_map< Vertex_handle , int > V;

	// cycle through the vertices of the interior of the alpha shape
	Alpha_shape::Alpha_shape_vertices_iterator vt = A.alpha_shape_vertices_begin();
	while ( vt != A.alpha_shape_vertices_end()) {
		V[*vt] = numVertices;
		vt++;
		numVertices ++;
	}

	// edges on the boundary?
	Alpha_shape::Alpha_shape_edges_iterator et = A.alpha_shape_edges_begin();

	// Edge is of type pair(pair(v1,v2), int)
	int numEdges = 0;

	std::vector<int> Indices;
	std::vector<Face_handle> Faces;
	Indices.resize(numVertices);
	Faces.resize(numVertices);

	while ( et != A.alpha_shape_edges_end() ) {

		Face_handle f=(*et).first; int i=(*et).second;
		Faces[V[f->vertex(f->ccw(i))]] = (*et).first;
		Indices[V[f->vertex(f->ccw(i))]] = (*et).second;
		et++;
		numEdges++;
	}

	// let's cycle through the edges until we return
	int init = 0;
	int curr = -1;

	Face_handle fh = Faces[init];
	int ind = Indices[init]; 
	int next = V[fh->vertex(fh->cw(ind))];
	curr = next;

	std::cout << numVertices << " ";
	std::cout << fh->vertex(fh->ccw(ind))->point().x() << " " << fh->vertex(fh->ccw(ind))->point().y() << " ";

	fh = Faces[curr];
	ind = Indices[curr]; 

	while ( curr != init ) {

		Face_handle fh = Faces[curr];
		int ind = Indices[curr]; 
		std::cout << fh->vertex(fh->ccw(ind))->point().x() << " " << fh->vertex(fh->ccw(ind))->point().y() << " ";

		int next = V[fh->vertex(fh->cw(ind))];
		curr = next;
	}

	fh = Faces[curr];
	ind = Indices[curr]; 
	std::cout << fh->vertex(fh->ccw(ind))->point().x() << " " << fh->vertex(fh->ccw(ind))->point().y();
	std::cout << std::endl;

  return 0;
}

