//
// Created by s165558 on 1-8-2024.
//

#ifndef CARTOCROW_GEOPHYLOGENY_SHAPE_H
#define CARTOCROW_GEOPHYLOGENY_SHAPE_H

#include "../core/core.h"
#include "phylo_tree/tree.h"

namespace cartocrow {
namespace geophylogeny {

class GeophyShape {
  public:

	virtual ~GeophyShape() = default;

  protected:
	GeophyShape(const std::vector<Point<Inexact>> site_locations, std::shared_ptr<Tree> tree);

	///Locations of the sites.
	std::vector<Point<Inexact>> m_site_locations;

	///Tree of the geophylogeny
	std::shared_ptr<Tree> m_tree;
};
}
}

#endif //CARTOCROW_GEOPHYLOGENY_SHAPE_H
