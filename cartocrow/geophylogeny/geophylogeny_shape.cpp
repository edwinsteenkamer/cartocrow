//
// Created by s165558 on 1-8-2024.
//

#include "geophylogeny_shape.h"
#include <vector>

namespace cartocrow {
namespace geophylogeny {

GeophyShape::GeophyShape(const std::vector<Point<Inexact>> site_locations, std::shared_ptr<Tree> tree)
	: m_site_locations(site_locations), m_tree(tree) {}

}
}
