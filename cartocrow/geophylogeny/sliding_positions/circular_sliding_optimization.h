//
// Created by s165558 on 23-8-2024.
//

#ifndef CARTOCROW_CIRCULAR_SLIDING_OPTIMIZATION_H
#define CARTOCROW_CIRCULAR_SLIDING_OPTIMIZATION_H

#include "../../core/core.h"
#include "../../core/timer.h"
#include "../phylo_tree/site.h"
#include "../circular_geophylogeny.h"

namespace cartocrow {
namespace geophylogeny {

struct CircularSlideOrdener {

	/// Constructor.
	CircularSlideOrdener(std::shared_ptr<CircularGeophylogeny> geophylogeny);

	/// Constructor with parameters.
	CircularSlideOrdener(std::shared_ptr<CircularGeophylogeny> geophylogeny, int num_cycles, Number<Inexact> aversion_centroid_ratio,
	                     Number<Inexact> interval_margin, bool allowed_outside_interval);

	/// Computes feasible intervals.
	void computeIntervals(Number<Inexact> min_value);

	/// Sets the positions of the leaves of the tree
	void setPositionsOfLeaves();

	/// Recursively places leaves that cannot be placed at the left end of their interval.
	void recursiveLeafPlacement(std::vector<std::shared_ptr<Node>> leaves, Number<Inexact> left_boundary, Number<Inexact> right_boundary);

	/// Translate the intervals of the sites.
	void translate(Number<Inexact> phi);

	/// Computes translated intervals;
	void translatedIntervals(Number<Inexact> min_value);

	/// Widens the intervals such that the leaves can adhere to the given margin.
	void computeIntervalsWithCorrectMargins();

	/// Returns the intersecting value of the intervals of two sites.
	Number<Inexact> intersectingValue(std::shared_ptr<Site> first_site, std::shared_ptr<Site> second_site);

	/// Returns the intersecting value of the intervals of two sites.
	Number<Inexact> intersectingValueWithMargin(std::shared_ptr<Site> first_site, std::shared_ptr<Site> second_site, Number<Inexact> margin);

	/// Returns the minimum value needed to make the tree suitable for the set of intervals.
	Number<Inexact> minValueSuitable();

	/// Applies a force-based model to optimize the positions of the leaves.
	void forceBasedPositioning(Number<Inexact> aversion_ratio, std::vector<std::shared_ptr<Node>> leaves);

	/// The geophylogeny being optimized.
	std::shared_ptr<CircularGeophylogeny> m_geophylogeny;

	/// The first site of the geophylogeny (to determine an order of the sites).
	std::shared_ptr<Site> first_site;

	/// The number of cycles used for the force-based position computing.
	int m_num_cycles;

	///The ratio between the aversion and centroid force.
	Number<Inexact> m_aversion_centroid_ratio;

	/// The margin used to enlarge the intervals.
	Number<Inexact> m_interval_margin;

	/// Tells if leaves are allowed to be pushed outside their intervals.
	bool m_allowed_outside_interval = false;

	Number<Inexact> m_min_value;

	/// The number of conflicts;
	int m_num_conflicts = 10000;

};

}
}

#endif //CARTOCROW_CIRCULAR_SLIDING_OPTIMIZATION_H
