//
// Created by s165558 on 19-7-2024.
//

#ifndef CARTOCROW_SITE_H
#define CARTOCROW_SITE_H

#include "../../core/core.h"
#include "../../necklace_map/range.h"
#include "../../necklace_map/circular_range.h"

namespace cartocrow {
namespace geophylogeny{

struct Site {

	/// Constructor of a site
	//Site(const Point<Inexact>& position);

	/// Constructor of a site with an associated leaf
	Site(Number<Inexact> x, Number<Inexact> y, Color color);

	/// The x-coordinate of the site.
	Number<Inexact> m_x;

	/// The y-coordinate of the site.
	Number<Inexact> m_y;

	/// The position of this site.
	Point<Inexact> m_position;

	/// The color of this site.
	Color m_color;

	/// The corresponding interval of this site.
	necklace_map::Range m_interval;

	/// The corresponding circular interval of this site.
	necklace_map::CircularRange m_circular_interval;

	/// Translated phi.
	Number<Inexact> m_trans_phi;

};
}
}



#endif //CARTOCROW_SITE_H
