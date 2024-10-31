//
// Created by s165558 on 22-7-2024.
//

#include <CGAL/bounding_box.h>

#include "geophylogeny.h"

namespace cartocrow {
namespace geophylogeny {

Geophylogeny::Geophylogeny(std::shared_ptr<Tree> tree, std::vector<std::shared_ptr<Site>> sites, BoundaryType boundary_type, PositionType position_type)
	: m_tree(tree), m_sites(sites), m_boundary_type(boundary_type), m_position_type(position_type) {
	//createBoundary(sites, boundary_type);
	createBoundary1();
}

void Geophylogeny::createBoundary1() {
	auto site_locations = getSitePositions(m_sites);
	if (m_boundary_type == Geophylogeny::BoundaryType::linear) {
		boundary = std::make_unique<LinearBoundary>(site_locations);
		//boundary->computeBoundary();
	}
	if (m_boundary_type == Geophylogeny::BoundaryType::circular) {
		boundary = std::make_unique<CircularBoundary>(site_locations);
		//boundary->computeBoundary();
	}
	boundary->computeBoundary();
}

std::vector<Point<Inexact>> Geophylogeny::getSitePositions(std::vector<std::shared_ptr<Site>> sites) {
	std::vector<Point<Inexact>> positions;
	for (const auto& site: sites) {
		const Point<Inexact>& pos = site->m_position;
		positions.push_back(pos);
	}
	return positions;
}

std::vector<Number<Inexact>> Geophylogeny::getSiteCoordinates(std::vector<Point<Inexact>> positions, int i) {
	std::vector<Number<Inexact>> coordinates;
	for (const auto& position: positions) {
		Number<Inexact> cor;
		if (i == 0) {
			cor = position.x();
		}
		if (i == 1) {
			cor = position.y();
		}
		coordinates.push_back(cor);
	}
	return coordinates;
}

void Geophylogeny::normalizeSitePositions(std::vector<std::shared_ptr<Site>> sites) {
	std::vector<Point<Inexact>> positions = getSitePositions(sites);
	auto x_coordinates = getSiteCoordinates(positions, 0);
	auto y_coordinates = getSiteCoordinates(positions, 1);
	auto min_x = *std::min_element(x_coordinates.begin(), x_coordinates.end());
	auto max_x = *std::max_element(x_coordinates.begin(), x_coordinates.end());
	auto min_y = *std::min_element(y_coordinates.begin(), y_coordinates.end());
	auto max_y = *std::max_element(y_coordinates.begin(), y_coordinates.end());

	auto length_x = max_x - min_x;
	auto length_y = max_y - min_y;

	auto x_ratio = 10 / length_x;
	auto y_ratio = (length_y / length_x * 10) / length_y;

	for (auto& site: sites) {
		site->m_x = site->m_x * x_ratio;
		site->m_y = site->m_y * y_ratio;
		site->m_position = Point<Inexact>(site->m_x, site->m_y);
	}
}


}
}

