//
// Created by s165558 on 19-9-2024.
//

#ifndef CARTOCROW_DP_ORDENER_RECTANGULAR_H
#define CARTOCROW_DP_ORDENER_RECTANGULAR_H

#include "../../core/core.h"
#include "../linear_geophylogeny.h"

namespace cartocrow {
namespace geophylogeny {

struct DPOrdenerRectangular {

	/// The strategy the DP uses for the quality measure.
	enum class DPStrategy {
		/// The Euclidean distance.
		kEuclidean,
		/// The horizontal distance.
		kHorizontal
	};

	/// Constructor.
	DPOrdenerRectangular(std::shared_ptr<RectangularGeophylogeny> geophylogeny, DPStrategy strategy);

	/// orders the leaves according to the quality measure.
	void orderLeaves();

	///Recovers the order of leaves from the computed minimum value of the quality measure.
	void recoverOrder(std::shared_ptr<Node> node, int position);

	void forceBasedPositioning(Number<Inexact> aversion_ratio, std::vector<std::shared_ptr<Node>> leaves);

	/// The geophylogeny being optimized.
	std::shared_ptr<RectangularGeophylogeny> m_geophylogeny;

	/// The strategy to optimize the geophylogeny.
	DPStrategy m_strategy;

	/// List storing the value of the quality measure of each node at each position.
	std::vector<std::vector<Number<Inexact>>> value_of_node_at_position;

	/// List of booleans storing if the first child should be the left child or not.
	std::vector<std::vector<bool>> first_as_left_child_at_position;
};
}
}

#endif //CARTOCROW_DP_ORDENER_RECTANGULAR_H
