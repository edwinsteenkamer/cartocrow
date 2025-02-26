//
// Created by s165558 on 15-8-2024.
//

#ifndef CARTOCROW_RECTANGULAR_SLIDING_OPTIMIZATION_H
#define CARTOCROW_RECTANGULAR_SLIDING_OPTIMIZATION_H

#include "../../core/core.h"
#include "../phylo_tree/site.h"
#include "../linear_geophylogeny.h"

namespace cartocrow {
namespace geophylogeny {

struct RectangularSlideOrdener {

	/// Constructor.
	RectangularSlideOrdener(std::shared_ptr<RectangularGeophylogeny> geophylogeny);

	/// Constructor with parameters.
	RectangularSlideOrdener(std::shared_ptr<RectangularGeophylogeny> geophylogeny, int num_cycles, Number<Inexact> aversion_centroid_ratio, Number<Inexact> interval_margin, bool allowed_outside_interval);

	/// Computes feasible intervals
	void computeIntervals(Number<Inexact> min_angle);

	/// Widens the intervals such that the leaves can adhere to the given margin.
	void computeIntervalsWithCorrectMargins();

	/// Sets the positions of the leaves of the tree.
	void setPositionsOfLeaves();

	/// Rrecursively places leaves that cannot be placed at the left end of their interval.
	void recursiveLeafPlacement(std::vector<std::shared_ptr<Node>> leaves, Number<Inexact> left_boundary, Number<Inexact> right_boundary);

	/// returns the center of an interval.
	Number<Inexact> middle(necklace_map::Range interval);

	/// Returns the intersecting angle of the intervals of two sites.
	Number<Inexact> intersectingAngle(std::shared_ptr<Site> first_site, std::shared_ptr<Site> second_site);

	/// Returns the intersecting angle of the intervals of two sites where there is a margin between the intersection on the boundary.
	Number<Inexact> intersectingAngleWithMargin(std::shared_ptr<Site> first_site, std::shared_ptr<Site> second_site, Number<Inexact> margin);

	/// Returns the minimum angle needed to make the tree suitable for the set of intervals.
	Number<Inexact> minAngleSuitable();

	/// Applies a force-based model to optimize the positions of the leaves.
	void forceBasedPositioning(Number<Inexact> aversion_ratio, std::vector<std::shared_ptr<Node>> leaves);

	/// The geophylogeny being optimized.
	std::shared_ptr<RectangularGeophylogeny> m_geophylogeny;

	/// The number of cycles used for the force-based position computing.
	int m_num_cycles;

	///The ratio between the aversion and centroid force.
	Number<Inexact> m_aversion_centroid_ratio;

	/// The margin used to enlarge the intervals.
	Number<Inexact> m_interval_margin;

	/// Tells if leaves are allowed to be pushed outside their intervals.
	bool m_allowed_outside_interval = false;

	/// The number of nodes that are a conflict
	int m_num_conflicts = 0;

	Number<Inexact> m_angle;

};

}
}
#endif //CARTOCROW_RECTANGULAR_SLIDING_OPTIMIZATION_H
