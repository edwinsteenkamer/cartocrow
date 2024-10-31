//
// Created by s165558 on 12-8-2024.
//

#include "circular_arc.h"

namespace cartocrow {

CircularArc::CircularArc(const Circle<Inexact>& circle, const Number<Inexact>& start_angle, const Number<Inexact>& span_angle)
	: m_circle(circle), m_start_angle(start_angle), m_span_angle(span_angle) {}

Circle<Inexact> CircularArc::circle() const {
	return m_circle;
}

Number<Inexact> CircularArc::startAngle() const {
	return m_start_angle;
}

Number<Inexact> CircularArc::spanAngle() const {
	return m_span_angle;
}

}
