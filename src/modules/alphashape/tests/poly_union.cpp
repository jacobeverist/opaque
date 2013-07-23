/*! \file simple_join_intersect.cpp
 * Computing the union and the intersection of two simple polygons.
 */

#include <iostream>
#include <fstream>
#include "cgal_types.h"
#include <CGAL/Unique_hash_map.h>

template<class Kernel, class Container>
void print_polygon (const CGAL::Polygon_2<Kernel, Container>& P)
{
  typename CGAL::Polygon_2<Kernel, Container>::Vertex_const_iterator  vit;

  for (vit = P.vertices_begin(); vit != P.vertices_end(); ++vit) {
    std::cout << (*vit)[0] << " " << (*vit)[1] << " ";
	}
  std::cout << std::endl;

}

int main ()
{

  Polygon_2 P;
  Polygon_2 Q;

	// read number of points
	int numPoints1 = 0;
	std::cin >> numPoints1;

	int numPoints2 = 0;
	std::cin >> numPoints2;

	// read points from stdin
	for ( int i = 0 ; i < numPoints1 ; i++ ) {
		double x, y;

		std::cin >> x;
		std::cin >> y;

		P.push_back(Point_2(x,y));
	}

	// read points from stdin
	for ( int i = 0 ; i < numPoints2 ; i++ ) {
		double x, y;

		std::cin >> x;
		std::cin >> y;

		Q.push_back(Point_2(x,y));
	}

	// Compute the union of P and Q.
  Polygon_with_holes_2 unionR;

	// perform the union operation
  CGAL::join (P, Q, unionR);

	//CGAL::Polygon_2 result = unionR.outer_boundary();

	int numVertices = unionR.outer_boundary().size();

	std::cout << numVertices << " ";

	print_polygon(unionR.outer_boundary());

	/*
  typename CGAL::Polygon_2<Rep, Container>::Vertex_const_iterator  vit;
  for (vit = unionR.outer_boundary().vertices_begin(); vit != unionR.outer_boundary().vertices_end(); ++vit) {
    std::cout << *vit.x() << " " << *vit.y() << " ";
	}
  std::cout << std::endl;
	*/


	return 1;

}

