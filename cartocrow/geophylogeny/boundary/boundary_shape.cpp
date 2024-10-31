//
// Created by s165558 on 1-8-2024.
//

#include "boundary_shape.h"

namespace cartocrow {
namespace geophylogeny {

Boundary::Boundary(const std::vector<Point<Inexact>> site_locations)
	: m_site_locations(site_locations) {}

}
}
