//
// Created by s165558 on 1-8-2024.
//

#include "circular_geophylogeny.h"

#include <CGAL/centroid.h>
#include <CGAL/Origin.h>
#include <CGAL/number_utils.h>
#include <CGAL/Homogeneous.h>
#include <CGAL/Min_circle_2.h>
#include <CGAL/Min_circle_2_traits_2.h>
#include <CGAL/Gmpz.h>

namespace cartocrow {
namespace geophylogeny {


CircularGeophylogeny::CircularGeophylogeny(std::shared_ptr<Tree> tree, std::vector<std::shared_ptr<Site>> sites, const PositionType position_type)
	: m_sites(sites), m_tree(tree), m_position_type(position_type) {
	Geophylogeny::normalizeSitePositions(m_sites);
	setOriginAsCenter();
	//rotateSites();
	for (auto& site: m_sites) {
		polar_sites.push_back(flow_map::PolarPoint(site->m_position));
	}
	computeRadius();
	createBoundary();
	leaf_step = M_2xPI / m_tree->m_leaves.size();
	//m_tree->setColors(m_color_difference);
	//ColorSites::colorSitesCircular(m_sites, radius);
}

CircularGeophylogeny::CircularGeophylogeny(std::shared_ptr<Tree> tree, std::vector<std::shared_ptr<Site>> sites, const PositionType position_type, Number<Inexact> color_difference, std::string color_distance)
    : m_sites(sites), m_tree(tree), m_position_type(position_type), m_color_difference(color_difference) {
	Geophylogeny::normalizeSitePositions(m_sites);
	setOriginAsCenter();
	//rotateSites();
	for (auto& site: m_sites) {
		polar_sites.push_back(flow_map::PolarPoint(site->m_position));
	}
	computeRadius();
	createBoundary();
	leaf_step = M_2xPI / m_tree->m_leaves.size();
	if (color_distance == "leaves") {
		m_tree->setColors(m_color_difference);
	} else {
		ColorSites::colorSitesCircular(m_sites, radius);
	}
}


void CircularGeophylogeny::createBoundary() {
	boundary = Circle<Inexact>(Point<Inexact>(CGAL::ORIGIN), CGAL::square(radius));
}

Point<Inexact> CircularGeophylogeny::computeCentroid(std::vector<Point<Inexact>> positions) {
	return CGAL::centroid(positions.begin(), positions.end());
}

Circle<Inexact> CircularGeophylogeny::computeSmallestEnclosingCircle(std::vector<Point<Inexact>> positions) {
	CGAL::Min_circle_2<CGAL::Min_circle_2_traits_2<Inexact>> enclosing_circle(positions.begin(), positions.end());
	auto optimization_circle = enclosing_circle.circle();
	auto center = optimization_circle.center();
	auto squared_radius = optimization_circle.squared_radius();
	return Circle<Inexact>(center, squared_radius);
}

void CircularGeophylogeny::rotateSites() {
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

void CircularGeophylogeny::setOriginAsCenter() {
	std::vector<Point<Inexact>> positions = Geophylogeny::getSitePositions(m_sites);
	auto center = computeSmallestEnclosingCircle(positions).center();
	for (auto& site: m_sites) {
		site->m_x -= center.x();
		site->m_y -= center.y();
		site->m_position = Point<Inexact>(site->m_x, site->m_y);
	}
}


void CircularGeophylogeny::setOriginAsCentroid() {
	std::vector<Point<Inexact>> positions = Geophylogeny::getSitePositions(m_sites);
	auto centroid = computeCentroid(positions);
	for (auto& site: m_sites) {
		site->m_x -= centroid.x();
		site->m_y -= centroid.y();
		site->m_position = Point<Inexact>(site->m_x, site->m_y);
	}
}

void CircularGeophylogeny::computeRadius() {
	Number<Inexact> max_distance = 0;
	for (auto& polar: polar_sites) {
		max_distance = std::max(max_distance, polar.r());
	}
	radius = max_distance + 1;
}

void CircularGeophylogeny::setPositionsOfLeaves() {
	for (auto& leaf: m_tree->m_leaves) {
		leaf->polar_position = flow_map::PolarPoint(radius, leaf->m_ID * leaf_step);
		leaf->m_position = leaf->polar_position.toCartesian();
	}
}


flow_map::PolarPoint CircularGeophylogeny::determineParentPositions(flow_map::PolarPoint first_child, flow_map::PolarPoint second_child) {
	Number<Inexact> phi;

	if (second_child.phi() < first_child.phi()) {
		phi = first_child.phi() + (second_child.phi() + M_2xPI - first_child.phi()) / 2;
	} else {
		phi = first_child.phi() + (second_child.phi() - first_child.phi()) / 2;
	}
	return flow_map::PolarPoint(std::max(first_child.r(), second_child.r()) + 0.3, phi);
}

void CircularGeophylogeny::setInnerPositions(std::shared_ptr<Node> node) {
	if (node->getType() != Node::ConnectionType::kLeaf) {
		if (node->m_first_child) {
			setInnerPositions(node->m_first_child);
		}
		if (node->m_second_child) {
			setInnerPositions(node->m_second_child);
		}
		if (node->firstAsRightChild) {
			node->setPolarPosition(determineParentPositions(node->m_first_child->polar_position,
			                                                node->m_second_child->polar_position));
		} else {
			node->setPolarPosition(determineParentPositions(node->m_second_child->polar_position,
			                                                node->m_first_child->polar_position));
		}
	}
}

flow_map::PolarPoint CircularGeophylogeny::getPolarPointByPosition(int position) {
	return flow_map::PolarPoint(radius, position * leaf_step);
}






}
}
