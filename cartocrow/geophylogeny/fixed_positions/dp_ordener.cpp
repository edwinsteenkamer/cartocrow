//
// Created by s165558 on 5-8-2024.
//

#include "dp_ordener.h"

#include <CGAL/squared_distance_2.h>
#include <cmath>

namespace cartocrow {
namespace geophylogeny {

DPOrdener::DPOrdener(std::shared_ptr<CircularGeophylogeny> geophylogeny, DPStrategy strategy)
	: m_geophylogeny(geophylogeny), m_strategy(strategy) {}

void DPOrdener::orderLeaves() {
	size_t numNodes = m_geophylogeny->m_tree->nodes().size();
	size_t numPositions = m_geophylogeny->m_sites.size();

	// Resize the 2D vector to the correct dimensions
	value_of_node_at_position.resize(numNodes);
	first_as_right_child_at_position.resize(numNodes);
	for (size_t i = 0; i < numNodes; ++i) {
		value_of_node_at_position[i].resize(numPositions);
		first_as_right_child_at_position[i].resize(numPositions);
	}

	for (auto& leaf: m_geophylogeny->m_tree->leaves()) {
		for (int position = 0; position < numPositions; position++) {
			if (m_strategy == DPOrdener::DPStrategy::kEuclidean) {
				value_of_node_at_position[leaf->m_ID][position] = sqrt(squared_distance(leaf->m_site->m_position,
				                                                                        m_geophylogeny->getPolarPointByPosition(position).toCartesian()));
				//std::cout << leaf->m_ID << " " << value_of_node_at_position[leaf->m_ID][position] << std::endl;
			}
			if (m_strategy == DPOrdener::DPStrategy::kRadial) {
				Number<Inexact> phi_difference = std::abs(flow_map::PolarPoint(leaf->m_site->m_position).phi() -
				                                         m_geophylogeny->getPolarPointByPosition(position).phi());
				value_of_node_at_position[leaf->m_ID][position] = std::min(phi_difference, M_2xPI - phi_difference) * flow_map::PolarPoint(leaf->m_site->m_position).r();

			}
		}
	}

	for (auto& node: m_geophylogeny->m_tree->innerNodes()) {
		for (int position = 0; position < numPositions; position++) {
			Number<Inexact> first_child_left = value_of_node_at_position[node->m_first_child->m_ID][position] +
			                                   value_of_node_at_position[node->m_second_child->m_ID][(position + node->m_first_child->m_clade_size) % numPositions];
			Number<Inexact> second_child_left = value_of_node_at_position[node->m_second_child->m_ID][position] +
			                                   value_of_node_at_position[node->m_first_child->m_ID][(position + node->m_second_child->m_clade_size) % numPositions];
			first_as_right_child_at_position[node->m_ID][position] = first_child_left <= second_child_left;
			value_of_node_at_position[node->m_ID][position] = std::min(first_child_left, second_child_left);
		}
	}

	Number<Inexact> min_value = 100000;
	int optimal_position;
	for (int position = 0; position < numPositions; position++) {
		if (min_value > value_of_node_at_position[m_geophylogeny->m_tree->m_root->m_ID][position]) {
			min_value = value_of_node_at_position[m_geophylogeny->m_tree->m_root->m_ID][position];
			optimal_position = position;
		}
	}

	recoverOrder(m_geophylogeny->m_tree->m_root, optimal_position);
	std::vector<std::shared_ptr<Node>> leaves = m_geophylogeny->m_tree->leavesByTreeOrderRight(m_geophylogeny->m_tree->m_root);
	//forceBasedPositioning(0.1, leaves);
}

void DPOrdener::recoverOrder(std::shared_ptr<Node> node, int position) {
	size_t numPositions = m_geophylogeny->m_sites.size();
	if (node->getType() == Node::ConnectionType::kLeaf) {
		node->setPolarPosition(m_geophylogeny->getPolarPointByPosition(position));
		return;
	}

	if (first_as_right_child_at_position[node->m_ID][position]) {
		node->firstAsRightChild = true;
		recoverOrder(node->m_first_child, position);
		recoverOrder(node->m_second_child, (position + node->m_first_child->m_clade_size) % numPositions);
	} else {
		node->firstAsRightChild = false;
		recoverOrder(node->m_second_child, position);
		recoverOrder(node->m_first_child, (position + node->m_second_child->m_clade_size) % numPositions);
	}
}

void DPOrdener::forceBasedPositioning(Number<Inexact> aversion_ratio, std::vector<std::shared_ptr<Node>> leaves) {
	Number<Inexact> centroid_ratio = 1 - aversion_ratio;
	Number<Inexact> m_interval_margin = std::min(0.2 / m_geophylogeny->radius, M_2xPI / (m_geophylogeny->m_tree->leaves().size() - 1));
	for (int cycles = 0; cycles < 100; cycles++) {
		for (int i = leaves.size() - 1; i >= 0; i--) {
			std::shared_ptr<Node>& leaf = leaves[i];

			const size_t index_prev = (i + leaves.size() - 1) % leaves.size();
			const size_t index_next = (i + 1) % leaves.size();

			const std::shared_ptr<Node>& left_neighbor = leaves[index_prev];
			const std::shared_ptr<Node>& right_neighbor = leaves[index_next];

			double left_boundary = left_neighbor->polar_position.phi();
			double right_boundary = right_neighbor->polar_position.phi();



			Number<Inexact> leaf_interval_centroid = flow_map::PolarPoint(leaf->m_site->m_position).phi();

			if (left_boundary > right_boundary) {
				if (left_boundary > leaf_interval_centroid && leaf_interval_centroid < leaves[0]->polar_position.phi()) {
					leaf_interval_centroid += M_2xPI;
				}
				right_boundary += M_2xPI;
			}

			necklace_map::CircularRange possible_range = necklace_map::CircularRange(left_boundary, right_boundary);

			//get right value of the centroid of the interval.
			if (!possible_range.contains(leaf_interval_centroid)) {
				Number<Inexact> to_centroid = necklace_map::CircularRange(right_boundary, leaf_interval_centroid).length();
				Number<Inexact> centroid_from = necklace_map::CircularRange(leaf_interval_centroid, left_boundary).length();

				if (to_centroid < centroid_from) {
					leaf_interval_centroid = right_boundary + to_centroid;
				} else {
					leaf_interval_centroid = left_boundary - centroid_from;
				}
			}

			Number<Inexact> margin = (right_boundary - left_boundary) / 2.0 - 0.001;
			left_boundary = left_boundary + std::min(margin, m_interval_margin);
			right_boundary = right_boundary - std::min(margin, m_interval_margin);
			if (right_boundary <= left_boundary) {
				std::cout << "Wrong computation of boundaries with margins!" << std::endl;
			}

			const double w_3 = centroid_ratio;
			const double w_2 =
			    -centroid_ratio * (leaf_interval_centroid + left_boundary + right_boundary);
			const double w_1 =
			    centroid_ratio * (leaf_interval_centroid * (left_boundary + right_boundary) +
			                      left_boundary * right_boundary) -
			    2 * aversion_ratio;
			const double w_0 =
			    aversion_ratio * (left_boundary + right_boundary) -
			    centroid_ratio * leaf_interval_centroid * left_boundary * right_boundary;

			const double q = (3 * w_3 * w_1 - w_2 * w_2) / (9 * w_3 * w_3);
			const double r = (9 * w_3 * w_2 * w_1 - 27 * w_3 * w_3 * w_0 - 2 * w_2 * w_2 * w_2) /
			                 (54 * w_3 * w_3 * w_3);

			double d = q * q * q + r * r;

			double cos_theta = r / std::sqrt(-q * q * q);

			// Ensure that cos_theta is within the valid range [-1, 1] to avoid domain errors
			if (cos_theta > 1.0) {
				cos_theta = 1.0;
				std::cout << "cos_theta is changed here" << std::endl;
			}
			if (cos_theta < -1.0) {
				cos_theta = -1.0;
				std::cout << "cos_theta is changed here" << std::endl;
			}

			// Step 2: Calculate theta
			double theta = std::acos(cos_theta);

			// Step 3: Calculate the square root of 2 * sqrt(Q)
			double sqrt2Q = 2 * std::sqrt(-q);

			// Step 4: Calculate the roots x1, x2, x3 using the provided formulas
			double x1 = sqrt2Q * std::cos(theta / 3.0) - w_2 / (3.0 * w_3);
			double x2 = sqrt2Q * std::cos((theta + 2.0 * M_PI) / 3.0) - w_2 / (3.0 * w_3);
			double x3 = sqrt2Q * std::cos((theta + 4.0 * M_PI) / 3.0) - w_2 / (3.0 * w_3);

			Number<Inexact> computed_position;



			if (possible_range.contains(x1)) {
				computed_position = x1;
			}
			if (possible_range.contains(x2)) {
				computed_position = x2;
			}
			if (possible_range.contains(x3)) {
				computed_position = x3;
			}

			//std::cout << aversion_ratio * (1 / (computed_position - left_boundary) - 1 / (right_boundary - computed_position)) +
			//                 centroid_ratio * (leaf_interval_centroid - computed_position) << std::endl;

			if (left_boundary < computed_position && computed_position < right_boundary) {
				leaf->setPolarPosition(flow_map::PolarPoint(m_geophylogeny->radius, wrapAngle(computed_position)));
			}
		}
	}
}

}
}
