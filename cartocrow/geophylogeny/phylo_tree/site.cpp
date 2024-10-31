//
// Created by s165558 on 19-7-2024.
//

#include "site.h"
#include <iostream>

namespace cartocrow {
namespace geophylogeny {

Site::Site(Number<Inexact> x, Number<Inexact> y, Color color)
	: m_x(x), m_y(y), m_color(color), m_interval(m_x, m_x), m_circular_interval(m_x, m_x) {
	m_position = Point<Inexact>(m_x, m_y);
}

}
}


