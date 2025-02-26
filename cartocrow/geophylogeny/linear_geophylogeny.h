//
// Created by s165558 on 1-8-2024.
//

#ifndef CARTOCROW_LINEAR_GEOPHYLOGENY_H
#define CARTOCROW_LINEAR_GEOPHYLOGENY_H

#include "../core/core.h"
#include "phylo_tree/tree.h"
#include "phylo_tree/site.h"
#include "phylo_tree/geophylogeny.h"
#include "color_sites.h"


namespace cartocrow {
namespace geophylogeny {

struct RectangularGeophylogeny {

	/// The type of positioning of the leaves of the geophylogeny.
	enum class PositionType {
		fixed,
		sliding
	};

	/// Constructor of a geophylogeny with a rectangular tree and boundary.
	RectangularGeophylogeny(std::shared_ptr<Tree> tree, std::vector<std::shared_ptr<Site>> sites, const PositionType position_type);

	/// Constructor of a geophylogeny with a rectangular tree and boundary with parameters.
	RectangularGeophylogeny(std::shared_ptr<Tree> tree, std::vector<std::shared_ptr<Site>> sites, const PositionType position_type, Number<Inexact> color_difference, std::string color_distance);

	/// The tree of this geophylogeny.
	std::shared_ptr<Tree> m_tree;

	/// The sites of this geophylogeny.
	std::vector<std::shared_ptr<Site>> m_sites;

	/// The type of positioning of the leaves of this geophylogeny.
	const PositionType m_position_type;

	/// The boundary of the geophylogeny.
	Polygon<Inexact> boundary;

	/// The boundary coordinate of the leaves of the geophylogeny.
	Number<Inexact> boundary_coordinate;

	/// The bounding box of the geophylogeny.
	Box boundary_box;

	/// The space between each leaf when positions are fixed.
	Number<Inexact> leaf_step;

	/// The parameter to determine color difference in the tree.
	Number<Inexact> m_color_difference = 0.0;

	/// creates a rectangular boundary on which the leaves are placed.
	void createBoundary();

	/// sets the y_coordinate of the leaves of the tree to the boundary coordinate.
	void setYCoordinates();

	/// set the x_coordinate of the leaves of the tree.
	void setXByID();

	/// set the position of the leaves of the tree.
	void setPositionOfLeaves();

	/// Determines the position of parent node of two children.
	static Point<Inexact> determineParentPosition(const Point<Inexact>& first_child, const Point<Inexact>& second_child);

	/// Sets the positions of the inner nodes of the tree.
	void setInnerPositions(std::shared_ptr<Node> node);

	/// Returns the point corresponding to the given fixed position.
	Point<Inexact> getPointByPosition(int position);

	/// returns the smallest enclosing circle of the original site positions.
	Circle<Inexact> computeSmallestEnclosingCircle(std::vector<Point<Inexact>> positions);

	/// Sets the center of the smallest enclosing circle of the sites to the origin by translating the sites.
	void setOriginAsCenter();

	void rotateSites();







};

}
}
#endif //CARTOCROW_LINEAR_GEOPHYLOGENY_H
