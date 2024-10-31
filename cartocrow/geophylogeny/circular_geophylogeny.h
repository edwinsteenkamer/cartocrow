//
// Created by s165558 on 1-8-2024.
//

#ifndef CARTOCROW_CIRCULAR_GEOPHYLOGENY_H
#define CARTOCROW_CIRCULAR_GEOPHYLOGENY_H

#include "../core/core.h"
#include "phylo_tree/tree.h"
#include "phylo_tree/site.h"
#include "phylo_tree/geophylogeny.h"
#include "../flow_map/polar_point.h"
#include "color_sites.h"

namespace cartocrow {
namespace geophylogeny {

struct CircularGeophylogeny {

	/// The type of positioning of the leaves of the geophylogeny.
	enum class PositionType {
		fixed,
		sliding
	};

	/// Constructor of a geophylogeny with a circular tree and boundary.
	CircularGeophylogeny(std::shared_ptr<Tree> tree, std::vector<std::shared_ptr<Site>> sites, const PositionType position_type);

	/// Constructor of a geophylogeny with a circular tree and boundary with parameters.
	CircularGeophylogeny(std::shared_ptr<Tree> tree, std::vector<std::shared_ptr<Site>> sites, const PositionType position_type, Number<Inexact> color_difference, std::string color_distance);

	/// The tree of this geophylogeny.
	std::shared_ptr<Tree> m_tree;

	/// The sites of this geophylogeny.
	const std::vector<std::shared_ptr<Site>> m_sites;

	/// The type of positioning of the leaves of this geophylogeny.
	const PositionType m_position_type;

	/// The boundary of the geophylogeny.
	Circle<Inexact> boundary;

	/// The space between each leaf when positions are fixed.
	Number<Inexact> leaf_step;

	/// List of site positions with type polar point.
	std::vector<flow_map::PolarPoint> polar_sites;

	/// The radius of the circular boundary
	Number<Inexact> radius;

	/// The parameter to determine color difference in the tree.
	Number<Inexact> m_color_difference = 0.0;

	/// creates a circular boundary on which the leaves are placed.
	void createBoundary();

	/// returns the centroid of the original site positions.
	Point<Inexact> computeCentroid(std::vector<Point<Inexact>> positions);

	/// Sets the centroid of the sites to the origin by translating the sites.
	void setOriginAsCentroid();

	/// returns the smallest enclosing circle of the original site positions.
	Circle<Inexact> computeSmallestEnclosingCircle(std::vector<Point<Inexact>> positions);

	/// Sets the center of the smallest enclosing circle of the sites to the origin by translating the sites.
	void setOriginAsCenter();

	/// Computes the radius of the circular boundary.
	void computeRadius();

	/// Sets the positions of the leaves of the tree
	void setPositionsOfLeaves();

	/// Sets the positions of the inner vertices of the tree.
	void setInnerPositions(std::shared_ptr<Node> node);

	/// returns the polar position of the parent of two children in the tree.
	flow_map::PolarPoint determineParentPositions(flow_map::PolarPoint first_child, flow_map::PolarPoint second_child);

	///Returns phi based on the position.
	flow_map::PolarPoint getPolarPointByPosition(int position);


};

}
}

#endif //CARTOCROW_CIRCULAR_GEOPHYLOGENY_H
