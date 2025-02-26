//
// Created by s165558 on 1-8-2024.
//

#include "linear_geophylogeny.h"
#include <CGAL/Min_circle_2.h>
#include <CGAL/Min_circle_2_traits_2.h>

namespace cartocrow {
namespace geophylogeny {

RectangularGeophylogeny::RectangularGeophylogeny(std::shared_ptr<Tree> tree, std::vector<std::shared_ptr<Site>> sites, const PositionType position_type)
	: m_sites(sites), m_tree(tree), m_position_type(position_type) {
	Geophylogeny::normalizeSitePositions(m_sites);
	createBoundary();
	setYCoordinates();
	leaf_step = (boundary_box.xmax() - boundary_box.xmin()) / (m_tree->m_leaves.size() - 1);
	//m_tree->setColors(m_color_difference);
}

RectangularGeophylogeny::RectangularGeophylogeny(std::shared_ptr<Tree> tree, std::vector<std::shared_ptr<Site>> sites, const PositionType position_type, Number<Inexact> color_difference, std::string color_distance)
    : m_sites(sites), m_tree(tree), m_position_type(position_type), m_color_difference(color_difference) {
	Geophylogeny::normalizeSitePositions(m_sites);
	//setOriginAsCenter();
	//rotateSites();
	createBoundary();
	setYCoordinates();
	leaf_step = (boundary_box.xmax() - boundary_box.xmin()) / (m_tree->m_leaves.size() - 1);
	if (color_distance == "leaves") {
		m_tree->setColors(m_color_difference);
	} else {
		ColorSites::colorSitesRectangular(m_sites, boundary_box);
	}

}


void RectangularGeophylogeny::createBoundary() {
	std::vector<Point<Inexact>> positions = Geophylogeny::getSitePositions(m_sites);
	auto x_coordinates = Geophylogeny::getSiteCoordinates(positions, 0);
	auto y_coordinates = Geophylogeny::getSiteCoordinates(positions, 1);
	auto min_x = *std::min_element(x_coordinates.begin(), x_coordinates.end()) - 2;
	auto max_x = *std::max_element(x_coordinates.begin(), x_coordinates.end()) + 2;
	auto min_y = *std::min_element(y_coordinates.begin(), y_coordinates.end()) - 1;
	auto max_y = *std::max_element(y_coordinates.begin(), y_coordinates.end()) + 1;

	boundary_box = Box(min_x, min_y, max_x, max_y);

	boundary_coordinate = max_y;


	boundary.push_back(Point<Inexact>(min_x, min_y));
	boundary.push_back(Point<Inexact>(max_x, min_y));
	boundary.push_back(Point<Inexact>(max_x, max_y));
	boundary.push_back(Point<Inexact>(min_x, max_y));
}
void RectangularGeophylogeny::setYCoordinates() {
	for (auto& leaf: m_tree->m_leaves) {
		leaf->m_y = boundary_box.ymax();
	}
}

void RectangularGeophylogeny::setXByID() {
	for (auto& leaf: m_tree->m_leaves) {
		leaf->m_x = boundary_box.xmin() + ((leaf->m_ID + 1) * leaf_step);
	}
}

void RectangularGeophylogeny::setPositionOfLeaves() {
	for (auto& leaf: m_tree->m_leaves) {
		leaf->m_position = Point<Inexact>(leaf->m_x, leaf->m_y);
	}
}

Point<Inexact> RectangularGeophylogeny::determineParentPosition(const Point<Inexact>& first_child, const Point<Inexact>& second_child) {
	Number<Inexact> x;

	if (first_child.x() <= second_child.x()) {
		x = first_child.x() + (second_child.x() - first_child.x()) / 2;
	} else {
		x = second_child.x() + (first_child.x() - second_child.x()) / 2;
	}

	return Point<Inexact>(x, std::max(first_child.y(), second_child.y()) + 0.3);
}

void RectangularGeophylogeny::setInnerPositions(std::shared_ptr<Node> node) {
	if (node->getType() != Node::ConnectionType::kLeaf) {
		if (node->m_first_child) {
			setInnerPositions(node->m_first_child);
		}
		if (node->m_second_child) {
			setInnerPositions(node->m_second_child);
		}
		node->setPosition(determineParentPosition(node->m_first_child->m_position,
		                                          node->m_second_child->m_position));
	}
}

Point<Inexact> RectangularGeophylogeny::getPointByPosition(int position) {
	return Point<Inexact>(boundary_box.xmin() + position * leaf_step, boundary_box.ymax());
}

Circle<Inexact> RectangularGeophylogeny::computeSmallestEnclosingCircle(std::vector<Point<Inexact>> positions) {
	CGAL::Min_circle_2<CGAL::Min_circle_2_traits_2<Inexact>> enclosing_circle(positions.begin(), positions.end());
	auto optimization_circle = enclosing_circle.circle();
	auto center = optimization_circle.center();
	auto squared_radius = optimization_circle.squared_radius();
	return Circle<Inexact>(center, squared_radius);
}

void RectangularGeophylogeny::setOriginAsCenter() {
	std::vector<Point<Inexact>> positions = Geophylogeny::getSitePositions(m_sites);
	auto center = computeSmallestEnclosingCircle(positions).center();
	for (auto& site: m_sites) {
		site->m_x -= center.x();
		site->m_y -= center.y();
		site->m_position = Point<Inexact>(site->m_x, site->m_y);
	}
}

void RectangularGeophylogeny::rotateSites() {
	for (auto& site: m_sites) {
		auto copy_x = site->m_x;
		auto copy_y = site->m_y;
		site->m_x = copy_y;
		site->m_y = -copy_x;
		if (site->m_y > 0) {
			site->m_y = site->m_y - 2;
		} else {
			site->m_y = site->m_y + 2;
		}
		site->m_position = Point<Inexact>(site->m_x, site->m_y);
	}
}




}
}