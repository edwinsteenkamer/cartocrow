//
// Created by s165558 on 22-7-2024.
//

#ifndef CARTOCROW_GEOPHYLOGENY_H
#define CARTOCROW_GEOPHYLOGENY_H

#include "../../core/core.h"
#include "site.h"
#include "tree.h"
#include "node.h"
#include "../boundary/boundary_shape.h"
#include "../boundary/linear_boundary.h"
#include "../boundary/circular_boundary.h"

namespace cartocrow {
namespace geophylogeny {

struct Geophylogeny {

	/// The type of boundary of the geophylogeny.
	enum class BoundaryType {
		circular,
		linear
	};

	/// The type of positioning of the leaves of the geophylogeny.
	enum class PositionType {
		fix,
		sliding
	};

	/// Constructor of a geophylogeny.
	//Geophylogeny(std::shared_ptr<Tree> tree, BoundaryType boundary_type, PositionType position_type);

	/// Constructor of a geophylogeny with separate sites.
	Geophylogeny(std::shared_ptr<Tree> tree, std::vector<std::shared_ptr<Site>> sites, BoundaryType boundary_type, PositionType position_type);

	/// Creates the boundary box/circle of the geophylogeny.
	//void createBoundary(std::vector<std::shared_ptr<Site>> sites, BoundaryType boundary_type);

	void createBoundary1();

	///returns a list containing the position of each site.
	static std::vector<Point<Inexact>> getSitePositions(std::vector<std::shared_ptr<Site>> sites);

	///returns a list containing either the x or y coordinate of each site (if i = 0, x-coordinates are returned, if i=1, y-coordinates).
	static std::vector<Number<Inexact>> getSiteCoordinates(std::vector<Point<Inexact>> positions, int i);

	///normalizes the positions of the sites.
	static void normalizeSitePositions(std::vector<std::shared_ptr<Site>> sites);

	///The tree of this geophylogeny.
	std::shared_ptr<Tree> m_tree;

	///The sites of this geophylogeny.
	const std::vector<std::shared_ptr<Site>> m_sites;

	///The type of boundary of the geophylogeny.
	const BoundaryType m_boundary_type;

	/// The type of positioning of the leaves of the geophylogeny.
	const PositionType m_position_type;

	/// The linear bounding box of the geophylogeny.
	Polygon<Inexact> linear_boundary;

	/// The circular bounding box of the geophylogeny.
	Circle<Inexact> circular_boundary;

	/// The boundary of the geophylogeny
	std::shared_ptr<Boundary> boundary;




};
}
}


#endif //CARTOCROW_GEOPHYLOGENY_H
