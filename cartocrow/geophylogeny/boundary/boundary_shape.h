//
// Created by s165558 on 1-8-2024.
//

#ifndef CARTOCROW_BOUNDARY_SHAPE_H
#define CARTOCROW_BOUNDARY_SHAPE_H

#include "../../core/core.h"
#include <vector>
#include <variant>

namespace cartocrow {
namespace geophylogeny {

class Boundary {
  public:
	virtual ~Boundary() = default;

	virtual void computeBoundary() = 0;

	virtual std::variant<Polygon<Inexact>, Circle<Inexact>> getBoundary() = 0;

	//virtual Number<Inexact> getBoundaryCoordinate() = 0;

	Number<Inexact> boundary_coordinate;

  protected:

	Boundary(std::vector<Point<Inexact>> site_locations);

	std::vector<Point<Inexact>> m_site_locations;



};

}
}

#endif //CARTOCROW_BOUNDARY_SHAPE_H
