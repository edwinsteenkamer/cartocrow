//
// Created by s165558 on 1-8-2024.
//

#include "linear_boundary.h"

namespace cartocrow {
namespace geophylogeny {


LinearBoundary::LinearBoundary(const std::vector<Point<Inexact>>& site_locations)
	: Boundary(site_locations) {
}

void LinearBoundary::computeBoundary() {
	auto x_coordinates = Geophylogeny::getSiteCoordinates(m_site_locations, 0);
	auto y_coordinates = Geophylogeny::getSiteCoordinates(m_site_locations, 1);
	auto min_x = *std::min_element(x_coordinates.begin(), x_coordinates.end());
	auto max_x = *std::max_element(x_coordinates.begin(), x_coordinates.end());
	auto min_y = *std::min_element(y_coordinates.begin(), y_coordinates.end());
	auto max_y = *std::max_element(y_coordinates.begin(), y_coordinates.end());

	computed_boundary.push_back(Point<Inexact>(min_x-1, min_y-1));
	computed_boundary.push_back(Point<Inexact>(max_x+1, min_y-1));
	computed_boundary.push_back(Point<Inexact>(max_x+1, max_y+1));
	computed_boundary.push_back(Point<Inexact>(min_x-1, max_y+1));

	boundary_coordinate = max_y + 1;
}

std::variant<Polygon<Inexact>, Circle<Inexact>> LinearBoundary::getBoundary() {
	return computed_boundary;
}


}
}
