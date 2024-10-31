//
// Created by s165558 on 1-8-2024.
//

#ifndef CARTOCROW_LINEAR_BOUNDARY_H
#define CARTOCROW_LINEAR_BOUNDARY_H

#include "boundary_shape.h"
#include "../../core/core.h"
#include <vector>
#include "../phylo_tree/geophylogeny.h"
#include <variant>

namespace cartocrow {
namespace geophylogeny {

class LinearBoundary : public Boundary {
  public:

	LinearBoundary(const std::vector<Point<Inexact>>& site_locations);

	void computeBoundary() override;

	std::variant<Polygon<Inexact>, Circle<Inexact>> getBoundary() override;

	//Number<Inexact> getBoundaryCoordinate();

	Polygon<Inexact> computed_boundary;

};

}
}



#endif //CARTOCROW_LINEAR_BOUNDARY_H
