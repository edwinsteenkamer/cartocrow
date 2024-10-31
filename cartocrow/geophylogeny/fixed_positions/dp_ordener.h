//
// Created by s165558 on 5-8-2024.
//

#ifndef CARTOCROW_DP_ORDENER_H
#define CARTOCROW_DP_ORDENER_H

#include "../../core/core.h"
#include "../circular_geophylogeny.h"

namespace cartocrow {
namespace geophylogeny {

struct DPOrdener {

	/// The strategy the DP uses for the quality measure.
	enum class DPStrategy {
		/// The Euclidean distance.
		kEuclidean,
		/// The radial distance.
		kRadial
	};

	/// Constructor.
	DPOrdener(std::shared_ptr<CircularGeophylogeny> geophylogeny, DPStrategy strategy);

	/// orders the leaves according to the quality measure.
	void orderLeaves();

	///Recovers the order of leaves from the computed minimum value of the quality measure.
	void recoverOrder(std::shared_ptr<Node> node, int position);

	/// Applies a force-based model to optimize the positions of the leaves.
	void forceBasedPositioning(Number<Inexact> aversion_ratio, std::vector<std::shared_ptr<Node>> leaves);

	/// The geophylogeny being optimized.
	std::shared_ptr<CircularGeophylogeny> m_geophylogeny;

	/// The strategy to optimize the geophylogeny.
	DPStrategy m_strategy;

	/// List storing the value of the quality measure of each node at each position.
	std::vector<std::vector<Number<Inexact>>> value_of_node_at_position;

	/// List of booleans storing if the first child should be the left child or not.
	std::vector<std::vector<bool>> first_as_right_child_at_position;
};
}
}

#endif //CARTOCROW_DP_ORDENER_H
