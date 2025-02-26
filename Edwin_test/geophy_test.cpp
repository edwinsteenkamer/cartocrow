//
// Created by s165558 on 30-7-2024.
//

#include <iostream>
#include <memory>

#include "cartocrow/core/core.h"
#include "cartocrow/flow_map/polar_point.h"
#include "cartocrow/geophylogeny/circular_geophylogeny.h"
#include "cartocrow/geophylogeny/linear_geophylogeny.h"
#include "cartocrow/geophylogeny/phylo_tree/geophylogeny.h"
#include "cartocrow/geophylogeny/phylo_tree/node.h"
#include "cartocrow/geophylogeny/phylo_tree/site.h"
#include "cartocrow/geophylogeny/phylo_tree/tree.h"
#include "cartocrow/geophylogeny/fixed_positions/dp_ordener.h"

using namespace cartocrow;
using namespace cartocrow::geophylogeny;

int main(int argc, char* argv[]) {
	std::shared_ptr<Site> s_1 = std::make_shared<Site>(0, 3, Color{2, 100, 8});
	std::shared_ptr<Site> s_2 = std::make_shared<Site>(3, 0, Color{1,0,1});
	std::shared_ptr<Node> v_2 = std::make_shared<Node>(0, s_1);
	std::shared_ptr<Node> v_3 = std::make_shared<Node>(1, s_2);
	std::shared_ptr<Node> v_1 = std::make_shared<Node>(2, v_2, v_3);

	std::vector<std::shared_ptr<Site>> sites;
	sites.push_back(s_1);
	sites.push_back(s_2);

	std::cout << "Hello" << std::endl;

	std::shared_ptr<Tree> tree = std::make_shared<Tree>(v_1);
	/*Geophylogeny geophy = Geophylogeny(tree, sites, Geophylogeny::BoundaryType::linear, Geophylogeny::PositionType::fix);
	for (const auto leaf: geophy.m_tree->leaves()) {
		std::shared_ptr<Node> parent = leaf->m_parent.lock();
		std::cout << parent->m_ID << std::endl;
	}

	RectangularGeophylogeny geophy1 = RectangularGeophylogeny(tree, sites, RectangularGeophylogeny::PositionType::fixed);
	for (const auto leaf: geophy1.m_tree->leaves()) {
		std::shared_ptr<Node> parent = leaf->m_parent.lock();
		std::cout << parent->m_ID << std::endl;
	}

	for (const auto node: geophy1.m_tree->nodes()) {
		std::cout << node->m_position.x() << std::endl;
	}*/

	std::shared_ptr<CircularGeophylogeny> geophylo = std::make_shared<CircularGeophylogeny>(tree, sites, CircularGeophylogeny::PositionType::fixed);
	for (auto& site: geophylo->polar_sites) {
		std::cout << site.r() << std::endl;
	}

	DPOrdener order = DPOrdener(geophylo, DPOrdener::DPStrategy::kEuclidean);
	order.orderLeaves();
	for (size_t i = 0; i < order.value_of_node_at_position.size(); ++i) {
		std::cout << "Leaf " << i << ": ";
		for (size_t j = 0; j < order.value_of_node_at_position[i].size(); ++j) {
			std::cout << "Position " << j << ": ";
			std::cout << order.value_of_node_at_position[i][j] << " ";
		}
		std::cout << std::endl;
	}










	return 0;
}

