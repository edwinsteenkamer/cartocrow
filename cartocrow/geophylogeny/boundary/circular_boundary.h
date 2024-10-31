//
// Created by s165558 on 2-8-2024.
//

#ifndef CARTOCROW_CIRCULAR_BOUNDARY_H
#define CARTOCROW_CIRCULAR_BOUNDARY_H

#include "boundary_shape.h"
#include "../../core/core.h"
#include <vector>

namespace cartocrow {
namespace geophylogeny {

class CircularBoundary : public Boundary {
  public:

	CircularBoundary(const std::vector<Point<Inexact>>& site_locations);

	void computeBoundary() override;

	std::variant<Polygon<Inexact>, Circle<Inexact>> getBoundary() override;

	//Number<Inexact> getBoundaryCoordinate() override;

	Circle<Inexact> computed_boundary;

};
}
}

#endif //CARTOCROW_CIRCULAR_BOUNDARY_H
