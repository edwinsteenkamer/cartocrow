//
// Created by s165558 on 23-8-2024.
//

#include "circular_sliding_optimization.h"

namespace cartocrow {
namespace geophylogeny {

CircularSlideOrdener::CircularSlideOrdener(std::shared_ptr<CircularGeophylogeny> geophylogeny)
	: m_geophylogeny(geophylogeny) {
	m_num_cycles = 100;
	m_aversion_centroid_ratio = 0.5;
	m_interval_margin = 0.5;
	m_allowed_outside_interval = false;
}

CircularSlideOrdener::CircularSlideOrdener(std::shared_ptr<CircularGeophylogeny> geophylogeny, int num_cycles, Number<Inexact> aversion_centroid_ratio,
                                           Number<Inexact> interval_margin, bool allowed_outside_interval)
	: m_geophylogeny(geophylogeny), m_num_cycles(num_cycles), m_aversion_centroid_ratio(aversion_centroid_ratio),
      m_interval_margin(interval_margin), m_allowed_outside_interval(allowed_outside_interval) {
	m_interval_margin = std::min(m_interval_margin / m_geophylogeny->radius, M_2xPI / (m_geophylogeny->m_tree->leaves().size() - 1));
}

Number<Inexact> CircularSlideOrdener::intersectingValue(std::shared_ptr<Site> first_site, std::shared_ptr<Site> second_site) {
	flow_map::PolarPoint first_site_polar = flow_map::PolarPoint(first_site->m_position);
	flow_map::PolarPoint second_site_polar = flow_map::PolarPoint(second_site->m_position);

	Number<Inexact> a = m_geophylogeny->radius - first_site_polar.r();
	Number<Inexact> b = m_geophylogeny->radius - second_site_polar.r();
	Number<Inexact> beta = std::abs(first_site->m_trans_phi - second_site->m_trans_phi);

	return beta * m_geophylogeny->radius / (a + b);
}

Number<Inexact> CircularSlideOrdener::intersectingValueWithMargin(std::shared_ptr<Site> first_site, std::shared_ptr<Site> second_site, Number<Inexact> margin) {
	flow_map::PolarPoint first_site_polar = flow_map::PolarPoint(first_site->m_position);
	flow_map::PolarPoint second_site_polar = flow_map::PolarPoint(second_site->m_position);

	Number<Inexact> a = m_geophylogeny->radius - first_site_polar.r();
	Number<Inexact> b = m_geophylogeny->radius - second_site_polar.r();
	Number<Inexact> beta = std::abs(first_site->m_trans_phi - second_site->m_trans_phi) + margin;

	return beta * m_geophylogeny->radius / (a + b);
}

Number<Inexact> CircularSlideOrdener::minValueSuitable() {
	Number<Inexact> min_value = 10000;
	for (auto& site: m_geophylogeny->m_sites) {
		Number<Inexact> min_value_site = 0;
		std::vector<bool> current_firstAsRightChild;

		for (auto& site1: m_geophylogeny->m_sites) {
			site1->m_trans_phi = wrapAngle(flow_map::PolarPoint(site1->m_position).phi() - flow_map::PolarPoint(site->m_position).phi());
		}

		for (auto& node: m_geophylogeny->m_tree->innerNodes()) {
			Tree subtree_first = Tree(node->m_first_child);
			Tree subtree_second = Tree(node->m_second_child);

			Number<Inexact> min_value_first_second = 0;
			Number<Inexact> min_value_second_first = 0;

			for (auto& leaf_first: subtree_first.leaves()) {
				for (auto& leaf_second: subtree_second.leaves()) {
					if (leaf_first->m_site->m_trans_phi < leaf_second->m_site->m_trans_phi) {
						min_value_first_second = std::max(min_value_first_second, intersectingValue(leaf_first->m_site, leaf_second->m_site));
					}
					if (leaf_second->m_site->m_trans_phi < leaf_first->m_site->m_trans_phi) {
						min_value_second_first = std::max(min_value_second_first, intersectingValue(leaf_first->m_site, leaf_second->m_site));
					}
				}
			}

			if (min_value_second_first <= min_value_first_second) {
				current_firstAsRightChild.push_back(true);
			} else {
				current_firstAsRightChild.push_back(false);
			}

			Number<Inexact> min_value_node = std::min(min_value_first_second, min_value_second_first);
			min_value_site = std::max(min_value_site, min_value_node);
		}

		if (min_value_site < min_value) {
			min_value = min_value_site;
			first_site = site;

			int node_index = 0;
			for (auto& node: m_geophylogeny->m_tree->innerNodes()) {
				node->firstAsRightChild = current_firstAsRightChild[node_index];
				node_index++;
			}
		}
	}
	return min_value;
}

void CircularSlideOrdener::computeIntervals(Number<Inexact> min_value) {
	translate(flow_map::PolarPoint(first_site->m_position).phi());

	for (auto& site: m_geophylogeny->m_sites) {
		flow_map::PolarPoint site_polar = flow_map::PolarPoint(site->m_position);
		Number<Inexact> a = m_geophylogeny->radius - site_polar.r();
		Number<Inexact> half_interval_size = a * (min_value / m_geophylogeny->radius);
		site->m_circular_interval = necklace_map::CircularRange(site_polar.phi() - std::min(M_PI - 0.0001, half_interval_size),
		                                                        site_polar.phi() + std::min(M_PI - 0.0001, half_interval_size));
	}
}

void CircularSlideOrdener::translatedIntervals(Number<Inexact> min_value) {
	translate(flow_map::PolarPoint(first_site->m_position).phi());

	for (auto& site: m_geophylogeny->m_sites) {
		flow_map::PolarPoint site_polar = flow_map::PolarPoint(site->m_position);
		Number<Inexact> a = m_geophylogeny->radius - site_polar.r();
		Number<Inexact> half_interval_size = a * (min_value / m_geophylogeny->radius);
		site->m_circular_interval = necklace_map::CircularRange(std::max(site->m_trans_phi - half_interval_size, Number<Inexact>(0)),
		                                                        std::min(site->m_trans_phi + half_interval_size, Number<Inexact>(M_2xPI)));
	}
}

void CircularSlideOrdener::translate(Number<Inexact> phi) {
	for (auto& site: m_geophylogeny->m_sites) {
		site->m_trans_phi = wrapAngle(flow_map::PolarPoint(site->m_position).phi() - phi);
	}
}

void CircularSlideOrdener::computeIntervalsWithCorrectMargins() {
	Number<Inexact> min_value = minValueSuitable();
	translatedIntervals(min_value);
	std::vector<std::shared_ptr<Node>> leaves_by_order = m_geophylogeny->m_tree->leavesByTreeOrderRight(m_geophylogeny->m_tree->m_root);
	for (int i = 0; i < leaves_by_order.size(); i++) {
		auto& leaf = leaves_by_order[i];
		Number<Inexact> left_boundary = leaf->m_site->m_circular_interval.from();
		Number<Inexact> num_leaves = 1;
		for (int j = i + 1; j < leaves_by_order.size(); j++) {
			auto& comparing_leaf = leaves_by_order[j];
			Number<Inexact> right_boundary = comparing_leaf->m_site->m_circular_interval.to();
			Number<Inexact> margin = (right_boundary - left_boundary) / num_leaves;
			if (margin < m_interval_margin) {
				Number<Inexact> min_distance = m_interval_margin * num_leaves;
				min_value = std::max(min_value, intersectingValueWithMargin(leaf->m_site, comparing_leaf->m_site, min_distance));
			}
			num_leaves++;
		}
	}
	m_min_value = min_value;
	translatedIntervals(min_value);
}

void CircularSlideOrdener::setPositionsOfLeaves() {
	computeIntervalsWithCorrectMargins();
	std::vector<std::shared_ptr<Node>> leaves = m_geophylogeny->m_tree->leavesByTreeOrderRight(m_geophylogeny->m_tree->m_root);

	std::vector<std::shared_ptr<Node>> leaves_recursive;
	std::vector<Number<Inexact>> boundaries;
	Number<Inexact> left_boundary;
	Number<Inexact> right_boundary;
	for (int i = 0; i < leaves.size(); i++) {
		if (m_allowed_outside_interval) {
			if (i == 0) {
				leaves[i]->setPolarPosition(flow_map::PolarPoint(m_geophylogeny->radius, leaves[i]->m_site->m_circular_interval.from()));
			} else {
				leaves[i]->setPolarPosition(flow_map::PolarPoint(m_geophylogeny->radius, leaves[0]->polar_position.phi() + i * m_geophylogeny->leaf_step));
			}
		} else {
			if (i == 0) {
				leaves[i]->setPolarPosition(flow_map::PolarPoint(m_geophylogeny->radius, leaves[i]->m_site->m_circular_interval.from()));
				boundaries.push_back(leaves[i]->polar_position.phi());
			}
			if (i == leaves.size() - 1) {
				leaves[i]->setPolarPosition(flow_map::PolarPoint(m_geophylogeny->radius, leaves[i]->m_site->m_circular_interval.to()));
				boundaries.push_back(leaves[i]->polar_position.phi());
				recursiveLeafPlacement(leaves_recursive, boundaries[0], boundaries[1]);
			}
			if (i > 0 && i < leaves.size() - 1 && leaves[i]->m_site->m_circular_interval.from() > boundaries[0]) {
				Number<Inexact> min_distance = boundaries[0] + m_interval_margin * (leaves_recursive.size() + 1);
				leaves[i]->setPolarPosition(flow_map::PolarPoint(m_geophylogeny->radius, std::max(leaves[i]->m_site->m_circular_interval.from(), min_distance)));
				boundaries.push_back(leaves[i]->polar_position.phi());
				recursiveLeafPlacement(leaves_recursive, boundaries[0], boundaries[1]);
				leaves_recursive.clear();
				boundaries.erase(boundaries.begin());
			}	else if (i > 0 && i < leaves.size() - 1) {
				leaves_recursive.push_back(leaves[i]);
			}
		}
	}
	for (auto& leaf: leaves) {
		leaf->setPolarPosition(flow_map::PolarPoint(m_geophylogeny->radius, leaf->polar_position.phi() + flow_map::PolarPoint(first_site->m_position).phi()));
	}
	computeIntervals(m_min_value);
	forceBasedPositioning(m_aversion_centroid_ratio, leaves);

}

void CircularSlideOrdener::recursiveLeafPlacement(std::vector<std::shared_ptr<Node>> leaves, Number<Inexact> left_boundary, Number<Inexact> right_boundary) {
	if (leaves.empty()) {
		return;
	}
	int k = leaves.size() - 1;
	right_boundary = wrapAngle(right_boundary, left_boundary);
	Number<Inexact> step_size = std::abs(right_boundary - left_boundary) / (leaves.size() + 1);
	if (leaves[k]->m_site->m_circular_interval.to() < right_boundary - step_size) {
		leaves[k]->setPolarPosition(flow_map::PolarPoint(m_geophylogeny->radius, leaves[k]->m_site->m_circular_interval.to()));
	} else {
		leaves[k]->setPolarPosition(flow_map::PolarPoint(m_geophylogeny->radius, right_boundary - step_size));
	}
	Number<Inexact> new_boundary = leaves[k]->polar_position.phi();
	leaves.pop_back();
	recursiveLeafPlacement(leaves, left_boundary, new_boundary);
}

void CircularSlideOrdener::forceBasedPositioning(Number<Inexact> aversion_ratio, std::vector<std::shared_ptr<Node>> leaves) {
	Number<Inexact> centroid_ratio = 1 - aversion_ratio;
	for (int cycles = 0; cycles < m_num_cycles; cycles++) {
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
				if (m_allowed_outside_interval) {
					leaf->setPolarPosition(flow_map::PolarPoint(m_geophylogeny->radius, wrapAngle(computed_position)));
				} else {
					if (leaf->m_site->m_circular_interval.contains(wrapAngle(computed_position))) {
						leaf->setPolarPosition(flow_map::PolarPoint(m_geophylogeny->radius,
						                                            wrapAngle(computed_position)));
					} else {
						if (!possible_range.contains(leaf->m_site->m_circular_interval.from())) {
							leaf->setPolarPosition(flow_map::PolarPoint(
							    m_geophylogeny->radius, leaf->m_site->m_circular_interval.to()));
						} else if (!possible_range.contains(leaf->m_site->m_circular_interval.to())) {
							leaf->setPolarPosition(flow_map::PolarPoint(
							    m_geophylogeny->radius, leaf->m_site->m_circular_interval.from()));
						} else if (necklace_map::CircularRange(
						               left_boundary, leaf->m_site->m_circular_interval.from())
						               .contains(wrapAngle(computed_position))) {
							leaf->setPolarPosition(flow_map::PolarPoint(
							    m_geophylogeny->radius, leaf->m_site->m_circular_interval.from()));
						} else if (necklace_map::CircularRange(
						               leaf->m_site->m_circular_interval.to(), right_boundary)
						               .contains(wrapAngle(computed_position))) {
							leaf->setPolarPosition(flow_map::PolarPoint(
							    m_geophylogeny->radius, leaf->m_site->m_circular_interval.to()));
						}
					}
				}
			}
		}
	}
}

}
}
