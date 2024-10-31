//
// Created by s165558 on 2-8-2024.
//

#include "circular_boundary.h"
#include <CGAL/centroid.h>
#include <CGAL/squared_distance_2.h>

namespace cartocrow {
namespace geophylogeny {

CircularBoundary::CircularBoundary(const std::vector<Point<Inexact>>& site_locations)
	: Boundary(site_locations) {}

void CircularBoundary::computeBoundary() {
	double max_distance = 0;
	Point<Inexact> centroid = CGAL::centroid(m_site_locations.begin(), m_site_locations.end());
	for (const auto& site: m_site_locations) {
		double site_distance = CGAL::squared_distance(site, centroid);
		max_distance = std::max(max_distance, site_distance);
	}
	Number<Inexact> radius = max_distance + 10;

	boundary_coordinate = radius;

	computed_boundary = Circle<Inexact>(centroid, radius);
}

std::variant<Polygon<Inexact>, Circle<Inexact>> CircularBoundary::getBoundary() {
	return computed_boundary;
}

}
}
