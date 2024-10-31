//
// Created by s165558 on 19-9-2024.
//

#include "dp_ordener_rectangular.h"

#include <CGAL/squared_distance_2.h>
#include <cmath>

namespace cartocrow {
namespace geophylogeny {

DPOrdenerRectangular::DPOrdenerRectangular(std::shared_ptr<RectangularGeophylogeny> geophylogeny, DPStrategy strategy)
    : m_geophylogeny(geophylogeny), m_strategy(strategy) {}

void DPOrdenerRectangular::orderLeaves() {
	size_t numNodes = m_geophylogeny->m_tree->nodes().size();
	size_t numPositions = m_geophylogeny->m_sites.size();

	// Resize the 2D vector to the correct dimensions
	value_of_node_at_position.resize(numNodes);
	first_as_left_child_at_position.resize(numNodes);
	for (size_t i = 0; i < numNodes; ++i) {
		value_of_node_at_position[i].resize(numPositions);
		first_as_left_child_at_position[i].resize(numPositions);
	}

	for (auto& leaf: m_geophylogeny->m_tree->leaves()) {
		for (int position = 0; position < numPositions; position++) {
			if (m_strategy == DPOrdenerRectangular::DPStrategy::kEuclidean) {
				value_of_node_at_position[leaf->m_ID][position] = sqrt(squared_distance(leaf->m_site->m_position,
				                                                                        m_geophylogeny->getPointByPosition(position)));
			}
			if (m_strategy == DPOrdenerRectangular::DPStrategy::kHorizontal) {
				value_of_node_at_position[leaf->m_ID][position] = std::abs(leaf->m_site->m_position.x() - m_geophylogeny->getPointByPosition(position).x());
			}
		}
	}

	for (auto& node: m_geophylogeny->m_tree->innerNodes()) {
		for (int position = 0; position < numPositions; position++) {
			Number<Inexact> first_child_left = value_of_node_at_position[node->m_first_child->m_ID][position] +
			                                   value_of_node_at_position[node->m_second_child->m_ID][position + node->m_first_child->m_clade_size];
			Number<Inexact> second_child_left = value_of_node_at_position[node->m_second_child->m_ID][position] +
			                                    value_of_node_at_position[node->m_first_child->m_ID][position + node->m_second_child->m_clade_size];
			first_as_left_child_at_position[node->m_ID][position] = first_child_left <= second_child_left;
			value_of_node_at_position[node->m_ID][position] = std::min(first_child_left, second_child_left);
		}
	}

	recoverOrder(m_geophylogeny->m_tree->m_root, 0);
	std::vector<std::shared_ptr<Node>> leaves = m_geophylogeny->m_tree->leavesByTreeOrder(m_geophylogeny->m_tree->m_root);
	//forceBasedPositioning(0.1, leaves);
}

void DPOrdenerRectangular::recoverOrder(std::shared_ptr<Node> node, int position) {
	if (node->getType() == Node::ConnectionType::kLeaf) {
		node->setPosition(m_geophylogeny->getPointByPosition(position));
		return;
	}

	if (first_as_left_child_at_position[node->m_ID][position]) {
		node->firstAsLeftChild = true;
		recoverOrder(node->m_first_child, position);
		recoverOrder(node->m_second_child, position + node->m_first_child->m_clade_size);
	} else {
		node->firstAsLeftChild = false;
		recoverOrder(node->m_second_child, position);
		recoverOrder(node->m_first_child, position + node->m_second_child->m_clade_size);
	}
}

void DPOrdenerRectangular::forceBasedPositioning(Number<Inexact> aversion_ratio, std::vector<std::shared_ptr<Node>> leaves) {
	Number<Inexact> centroid_ratio = 1 - aversion_ratio;
	Number<Inexact> m_interval_margin = std::min(0.2, (m_geophylogeny->boundary_box.xmax() - m_geophylogeny->boundary_box.xmin()) / (m_geophylogeny->m_tree->leaves().size() - 1));
	for (int cycles = 0; cycles < 100; cycles++) {
		for (int i = leaves.size() - 1; i >= 0; i--) {
			double left_boundary;
			double right_boundary;
			std::shared_ptr<Node>& leaf = leaves[i];

			if (i == leaves.size() - 1) {
				left_boundary = leaves[i-1]->m_position.x();
				right_boundary = m_geophylogeny->boundary_box.xmax();
			} else if (i == 0) {
				left_boundary = m_geophylogeny->boundary_box.xmin();
				right_boundary = leaves[i+1]->m_position.x();
			} else {
				std::shared_ptr<Node>& left_neighbor = leaves[i - 1];
				std::shared_ptr<Node>& right_neighbor = leaves[i + 1];

				left_boundary = left_neighbor->m_position.x();
				right_boundary = right_neighbor->m_position.x();
			}

			Number<Inexact> margin = (right_boundary - left_boundary) / 2.0 - 0.001;
			left_boundary = left_boundary + std::min(margin, m_interval_margin);
			right_boundary = right_boundary - std::min(margin, m_interval_margin);
			if (right_boundary <= left_boundary) {
				std::cout << "Wrong computation of boundaries with margins!" << std::endl;
			}

			Number<Inexact> leaf_interval_centroid = leaf->m_site->m_position.x();

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

			std::complex<double> sqrtD = sqrt(std::complex<double>(d));
			std::complex<double> s = std::pow(r + sqrtD, 1.0 / 3.0);
			std::complex<double> t = std::pow(r - sqrtD, 1.0 / 3.0);

			double common_real_part = -w_2 / (3.0 * w_3);

			std::complex<double> x1 = s + t + common_real_part;
			std::complex<double> x2 = -0.5 * (s + t) + common_real_part + std::complex<double>(0, sqrt(3) / 2.0) * (s - t);
			std::complex<double> x3 = -0.5 * (s + t) + common_real_part - std::complex<double>(0, sqrt(3) / 2.0) * (s - t);

			//std::cout << aversion_ratio * (1 / (std::real(x3) - left_boundary) - 1 / (right_boundary - std::real(x3))) +
			//                 centroid_ratio * (leaf_interval_centroid - std::real(x3)) << std::endl;
			if (left_boundary < std::real(x3) && std::real(x3) < right_boundary) {
				if (std::real(x3) > m_geophylogeny->boundary_box.xmax()) {
					leaf->setPosition(Point<Inexact>(m_geophylogeny->boundary_box.xmax(), m_geophylogeny->boundary_box.ymax()));
				} else if (std::real(x3) < m_geophylogeny->boundary_box.xmin()) {
					leaf->setPosition(Point<Inexact>(m_geophylogeny->boundary_box.xmin(), m_geophylogeny->boundary_box.ymax()));
				} else {
					leaf->setPosition(Point<Inexact>(std::real(x3), m_geophylogeny->boundary_box.ymax()));
				}
			}
		}
	}
}

}
}