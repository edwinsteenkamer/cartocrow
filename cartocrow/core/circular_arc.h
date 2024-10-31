//
// Created by s165558 on 12-8-2024.
//

#ifndef CARTOCROW_CIRCULAR_ARC_H
#define CARTOCROW_CIRCULAR_ARC_H

#include "core.h"

namespace cartocrow {

class CircularArc {

  public:
	CircularArc(const Circle<Inexact>& circle, const Number<Inexact>& start_angle, const Number<Inexact>& span_angle);

	Circle<Inexact> circle() const;

	Number<Inexact> startAngle() const;

	Number<Inexact> spanAngle() const;

	Circle<Inexact> m_circle;
	Number<Inexact> m_start_angle;
	Number<Inexact> m_span_angle;

};

}

#endif //CARTOCROW_CIRCULAR_ARC_H
