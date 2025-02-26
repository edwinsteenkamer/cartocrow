//
// Created by s165558 on 15-8-2024.
//

#include <cmath>
#include <complex>

#include "rectangular_sliding_optimization.h"

namespace cartocrow {
namespace geophylogeny {

RectangularSlideOrdener::RectangularSlideOrdener(std::shared_ptr<RectangularGeophylogeny> geophylogeny)
	: m_geophylogeny(geophylogeny) {
	m_num_cycles = 100;
	m_aversion_centroid_ratio = 0.5;
	m_allowed_outside_interval = false;
	m_interval_margin = 0.5;
}

RectangularSlideOrdener::RectangularSlideOrdener(std::shared_ptr<RectangularGeophylogeny> geophylogeny, int num_cycles,
                                                 Number<Inexact> aversion_centroid_ratio, Number<Inexact> interval_margin, bool allowed_outside_interval)
	: m_geophylogeny(geophylogeny), m_num_cycles(num_cycles), m_aversion_centroid_ratio(aversion_centroid_ratio),
      m_interval_margin(interval_margin), m_allowed_outside_interval(allowed_outside_interval) {
	m_interval_margin = std::min(m_interval_margin, (m_geophylogeny->boundary_box.xmax() - m_geophylogeny->boundary_box.xmin()) / (m_geophylogeny->m_tree->leaves().size() - 1));
}

Number<Inexact> RectangularSlideOrdener::intersectingAngle(std::shared_ptr<Site> first_site, std::shared_ptr<Site> second_site) {
	Number<Inexact> a = std::abs(m_geophylogeny->boundary_box.ymax() - first_site->m_position.y());
	Number<Inexact> b = std::abs(m_geophylogeny->boundary_box.ymax() - second_site->m_position.y());
	Number<Inexact> l = std::abs(first_site->m_position.x() - second_site->m_position.x());

	return std::atan(l / (a + b));
}

Number<Inexact> RectangularSlideOrdener::intersectingAngleWithMargin(std::shared_ptr<Site> first_site, std::shared_ptr<Site> second_site, Number<Inexact> margin) {
	Number<Inexact> a = std::abs(m_geophylogeny->boundary_box.ymax() - first_site->m_position.y());
	Number<Inexact> b = std::abs(m_geophylogeny->boundary_box.ymax() - second_site->m_position.y());
	Number<Inexact> l = std::abs(first_site->m_position.x() - second_site->m_position.x()) + margin;

	return std::atan(l / (a + b));
}

Number<Inexact> RectangularSlideOrdener::minAngleSuitable() {
	Number<Inexact> min_angle = 0;
	std::vector<Number<Inexact>> conflict_values;
	for (auto& node: m_geophylogeny->m_tree->innerNodes()) {
		Tree subtree_first = Tree(node->m_first_child);
		Tree subtree_second = Tree(node->m_second_child);

		Number<Inexact> min_angle_first_second = 0;
		Number<Inexact> min_angle_second_first = 0;

		for (auto& leaf_first: subtree_first.leaves()) {
			for (auto& leaf_second: subtree_second.leaves()) {
				if (leaf_first->m_site->m_position.x() < leaf_second->m_site->m_position.x()) {
					min_angle_first_second = std::max(min_angle_first_second, intersectingAngle(leaf_first->m_site, leaf_second->m_site));
				}
				if (leaf_second->m_site->m_position.x() < leaf_first->m_site->m_position.x()) {
					min_angle_second_first = std::max(min_angle_second_first, intersectingAngle(leaf_first->m_site, leaf_second->m_site));
				}
			}
		}

		if (min_angle_second_first <= min_angle_first_second) {
			node->firstAsLeftChild = true;
		} else {
			node->firstAsLeftChild = false;
		}

		Number<Inexact> min_angle_node = std::min(min_angle_first_second, min_angle_second_first);
		if (min_angle_node > 0) {
			m_num_conflicts++;
		}
		conflict_values.push_back(min_angle_node);
		min_angle = std::max(min_angle, min_angle_node);
	}
	std::sort(conflict_values.begin(), conflict_values.end());
	for (auto& value: conflict_values) {
		if (value > 0) {
			//std::cout << value << std::endl;
		}
	}
	return min_angle;
}

void RectangularSlideOrdener::computeIntervals(Number<Inexact> min_angle) {
	for (auto& site: m_geophylogeny->m_sites) {
		Number<Inexact> a = std::abs(m_geophylogeny->boundary_box.ymax() - site->m_position.y());
		Number<Inexact> half_interval_size = std::abs(std::tan(min_angle) * a);
		site->m_interval = necklace_map::Range(site->m_position.x() - half_interval_size, site->m_position.x() + half_interval_size);
	}


}

void RectangularSlideOrdener::computeIntervalsWithCorrectMargins() {
	Number<Inexact> min_angle = minAngleSuitable();
	computeIntervals(min_angle);
	std::vector<std::shared_ptr<Node>> leaves_by_order = m_geophylogeny->m_tree->leavesByTreeOrder(m_geophylogeny->m_tree->m_root);
	for (int i = 0; i < leaves_by_order.size(); i++) {
		auto& leaf = leaves_by_order[i];
		Number<Inexact> left_boundary = std::max(m_geophylogeny->boundary_box.xmin(), leaf->m_site->m_interval.from());
		Number<Inexact> num_leaves = 1;
		for (int j = i + 1; j < leaves_by_order.size(); j++) {
			auto& comparing_leaf = leaves_by_order[j];
			Number<Inexact> right_boundary = std::min(m_geophylogeny->boundary_box.xmax(), comparing_leaf->m_site->m_interval.to());
			Number<Inexact> margin = (right_boundary - left_boundary) / num_leaves;
			if (margin < m_interval_margin) {
				Number<Inexact> min_distance = m_interval_margin * num_leaves;
				min_angle = std::max(min_angle, intersectingAngleWithMargin(leaf->m_site, comparing_leaf->m_site, min_distance));
			}
			num_leaves++;
		}
	}
	m_angle = min_angle * 180 / M_PI;
	computeIntervals(min_angle);
}

void RectangularSlideOrdener::setPositionsOfLeaves() {
	computeIntervalsWithCorrectMargins();
	std::vector<std::shared_ptr<Node>> leaves = m_geophylogeny->m_tree->leavesByTreeOrder(m_geophylogeny->m_tree->m_root);

	std::vector<std::shared_ptr<Node>> leaves_recursive;
	std::vector<Number<Inexact>> boundaries;
	Number<Inexact> left_boundary;
	Number<Inexact> right_boundary;
	for (int i = 0; i < leaves.size(); i++) {
		if (m_allowed_outside_interval) {
			leaves[i]->setPosition(Point<Inexact>(m_geophylogeny->boundary_box.xmin() + i * (m_geophylogeny->leaf_step), m_geophylogeny->boundary_box.ymax()));
		} else {
			if (i == 0) {
				leaves[i]->setPosition(Point<Inexact>(std::max(m_geophylogeny->boundary_box.xmin(), leaves[i]->m_site->m_interval.from()), m_geophylogeny->boundary_box.ymax()));
				boundaries.push_back(leaves[i]->m_position.x());
			}
			if (i == leaves.size() - 1) {
				leaves[i]->setPosition(Point<Inexact>(std::min(m_geophylogeny->boundary_box.xmax(), leaves[i]->m_site->m_interval.to()), m_geophylogeny->boundary_box.ymax()));
				boundaries.push_back(leaves[i]->m_position.x());
				recursiveLeafPlacement(leaves_recursive, boundaries[0], boundaries[1]);
			}
			if (i > 0 && i < leaves.size() - 1 && leaves[i]->m_site->m_interval.from() > boundaries[0]) {
				Number<Inexact> min_distance = boundaries[0] + m_interval_margin * (leaves_recursive.size() + 1);
				leaves[i]->setPosition(Point<Inexact>(std::max(std::max(m_geophylogeny->boundary_box.xmin(), leaves[i]->m_site->m_interval.from()), min_distance), m_geophylogeny->boundary_box.ymax()));
				boundaries.push_back(leaves[i]->m_position.x());
				if (boundaries.size() == 2) {
					recursiveLeafPlacement(leaves_recursive, boundaries[0], boundaries[1]);
					leaves_recursive.clear();
					boundaries.erase(boundaries.begin());
				}
			}	else if (i > 0 && i < leaves.size() - 1) {
				leaves_recursive.push_back(leaves[i]);
			}
		}
	}
	forceBasedPositioning(m_aversion_centroid_ratio, leaves);
}

void RectangularSlideOrdener::recursiveLeafPlacement(std::vector<std::shared_ptr<Node>> leaves, Number<Inexact> left_boundary, Number<Inexact> right_boundary) {
	if (leaves.empty()) {
		return;
	}
	int k = leaves.size() - 1;
	Number<Inexact> step_size = std::abs(right_boundary - left_boundary) / (leaves.size() + 1);
	if (leaves[k]->m_site->m_interval.to() < right_boundary - step_size) {
		leaves[k]->setPosition(Point<Inexact>(leaves[k]->m_site->m_interval.to(), m_geophylogeny->boundary_box.ymax()));
	} else {
		leaves[k]->setPosition(Point<Inexact>(right_boundary - step_size, m_geophylogeny->boundary_box.ymax()));
	}
	if (leaves[k]->m_position.y() < m_geophylogeny->boundary_box.ymax()) {
		std::cout << leaves[k]->m_ID << ": ";
		std::cout << "Something goes wrong here!" << std::endl;
	}
	Number<Inexact> new_boundary = leaves[k]->m_position.x();
	leaves.pop_back();
	recursiveLeafPlacement(leaves, left_boundary, new_boundary);
}

Number<Inexact> RectangularSlideOrdener::middle(necklace_map::Range interval) {
	return (interval.from() + interval.to()) * 0.5;
}

void RectangularSlideOrdener::forceBasedPositioning(Number<Inexact> aversion_ratio, std::vector<std::shared_ptr<Node>> leaves) {
	Number<Inexact> centroid_ratio = 1 - aversion_ratio;
	for (int cycles = 0; cycles < m_num_cycles; cycles++) {
		for (int i = leaves.size() - 1; i >= 0; i--) {
			double left_boundary;
			double right_boundary;
			std::shared_ptr<Node>& leaf = leaves[i];

			if (i == leaves.size() - 1) {
				left_boundary = leaves[i-1]->m_position.x();
				right_boundary = std::max(leaf->m_site->m_interval.to(), m_geophylogeny->boundary_box.xmax());
			} else if (i == 0) {
				left_boundary = std::min(leaf->m_site->m_interval.from(), m_geophylogeny->boundary_box.xmin());
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

			Number<Inexact> leaf_interval_centroid = middle(leaf->m_site->m_interval);

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
				if (m_allowed_outside_interval) {
					if (std::real(x3) > m_geophylogeny->boundary_box.xmax()) {
						leaf->setPosition(Point<Inexact>(m_geophylogeny->boundary_box.xmax(), m_geophylogeny->boundary_box.ymax()));
					} else if (std::real(x3) < m_geophylogeny->boundary_box.xmin()) {
						leaf->setPosition(Point<Inexact>(m_geophylogeny->boundary_box.xmin(), m_geophylogeny->boundary_box.ymax()));
					} else {
						leaf->setPosition(Point<Inexact>(std::real(x3), m_geophylogeny->boundary_box.ymax()));
					}
				} else {
					if (std::real(x3) < std::max(m_geophylogeny->boundary_box.xmin(), leaf->m_site->m_interval.from())) {
						leaf->setPosition(Point<Inexact>(std::max(m_geophylogeny->boundary_box.xmin(), leaf->m_site->m_interval.from()), m_geophylogeny->boundary_box.ymax()));
					} else if (std::real(x3) > std::min(m_geophylogeny->boundary_box.xmax(), leaf->m_site->m_interval.to())) {
						leaf->setPosition(Point<Inexact>(std::min(m_geophylogeny->boundary_box.xmax(), leaf->m_site->m_interval.to()), m_geophylogeny->boundary_box.ymax()));
					} else {
						leaf->setPosition(Point<Inexact>(std::real(x3), m_geophylogeny->boundary_box.ymax()));
					}
				}



			}
		}
	}
}

}
}