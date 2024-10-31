//
// Created by s165558 on 20-8-2024.
//

#include <memory>
#include <iostream>


#include <QApplication>


#include "cartocrow/geophylogeny/phylo_tree/node.h"
#include "cartocrow/geophylogeny/phylo_tree/tree.h"
#include "cartocrow/geophylogeny/painting.h"
#include "cartocrow/geophylogeny/circular_painting.h"
#include "cartocrow/core/core.h"
#include "cartocrow/renderer/geometry_painting.h"
#include "cartocrow/renderer/geometry_widget.h"
#include "cartocrow/geophylogeny/linear_geophylogeny.h"
#include "cartocrow/geophylogeny/circular_geophylogeny.h"
#include "cartocrow/flow_map/polar_point.h"
#include "cartocrow/geophylogeny/fixed_positions/dp_ordener.h"
#include "cartocrow/geophylogeny/sliding_positions/rectangular_sliding_optimization.h"
#include "cartocrow/geophylogeny/sliding_positions/circular_sliding_optimization.h"

using namespace cartocrow;
using namespace cartocrow::geophylogeny;



int main(int argc, char* argv[]) {
	std::shared_ptr<renderer::GeometryPainting> painting;
	std::shared_ptr<renderer::GeometryPainting> painting1;

	std::shared_ptr<Site> s_1 = std::make_shared<Site>(2, -5, Color{0, 172, 130});
	std::shared_ptr<Site> s_2 = std::make_shared<Site>(5.5, -8, Color{0, 102, 204});
	std::shared_ptr<Site> s_3 = std::make_shared<Site>(5, -1, Color{214, 0, 74});
	std::shared_ptr<Site> s_4 = std::make_shared<Site>(9, -2, Color{132, 210, 0});
	std::shared_ptr<Site> s_5 = std::make_shared<Site>(12, -7, Color{255, 154, 0});
	std::shared_ptr<Site> s_6 = std::make_shared<Site>(9.5, -4.5, Color{255, 221, 0});
	std::vector<std::shared_ptr<Site>> sites;
	sites.push_back(s_1);
	sites.push_back(s_2);
	sites.push_back(s_3);
	sites.push_back(s_4);
	sites.push_back(s_5);
	sites.push_back(s_6);

	std::vector<std::shared_ptr<Site>> sites_one;
	sites_one.push_back(s_1);
	sites_one.push_back(s_2);
	sites_one.push_back(s_3);
	sites_one.push_back(s_4);
	sites_one.push_back(s_5);
	sites_one.push_back(s_6);



	std::shared_ptr<Node> v_1 = std::make_shared<Node>(0, s_1);
	std::shared_ptr<Node> v_2 = std::make_shared<Node>(1, s_2);
	std::shared_ptr<Node> v_3 = std::make_shared<Node>(2, s_3);
	std::shared_ptr<Node> v_4 = std::make_shared<Node>(3, s_4);
	std::shared_ptr<Node> v_5 = std::make_shared<Node>(4, s_5);
	std::shared_ptr<Node> v_6 = std::make_shared<Node>(5, s_6);
	std::shared_ptr<Node> v_7 = std::make_shared<Node>(6, v_1, v_2);
	std::shared_ptr<Node> v_8 = std::make_shared<Node>(7, v_3, v_4);
	std::shared_ptr<Node> v_9 = std::make_shared<Node>(8, v_7, v_8);
	std::shared_ptr<Node> v_10 = std::make_shared<Node>(9, v_5, v_6);
	std::shared_ptr<Node> v_11 = std::make_shared<Node>(10, v_9, v_10);

	/*v_6->m_position = Point<Inexact>(12, 0);
	v_5->m_position = Point<Inexact>(10, 0);
	v_4->m_position = Point<Inexact>(8, 0);
	v_3->m_position = Point<Inexact>(6, 0);
	v_2->m_position = Point<Inexact>(4, 0);
	v_1->m_position = Point<Inexact>(2, 0);
	*/


	std::shared_ptr<Tree> tree = std::make_shared<Tree>(v_11);
	std::shared_ptr<Tree> tree_one = std::make_shared<Tree>(v_11);

	//std::shared_ptr<Geophylogeny> geophy = std::make_shared<Geophylogeny>(tree, sites, Geophylogeny::BoundaryType::linear, Geophylogeny::PositionType::fix);
	/*std::shared_ptr<RectangularGeophylogeny> geophy1 = std::make_shared<RectangularGeophylogeny>(tree, sites, RectangularGeophylogeny::PositionType::fixed);



	RectangularSlideOrdener slide = RectangularSlideOrdener(geophy1);
	slide.setPositionsOfLeaves();

	geophy1->setInnerPositions(geophy1->m_tree->m_root);



	for (auto& node: geophy1->m_tree->innerNodes()) {
		if (node->firstAsLeftChild) {
			std::cout << "Yes" << std::endl;
		} else {
			std::cout << "No" << std::endl;
		}
	}



	for (auto& leaf: geophy1->m_tree->leavesByTreeOrder(geophy1->m_tree->m_root)) {
		std::cout << leaf->m_ID << std::endl;
	}*/


	std::shared_ptr<CircularGeophylogeny> geophy2 = std::make_shared<CircularGeophylogeny>(tree_one, sites_one, CircularGeophylogeny::PositionType::fixed);

	DPOrdener order = DPOrdener(geophy2, DPOrdener::DPStrategy::kEuclidean);
	order.orderLeaves();

	//CircularSlideOrdener slide = CircularSlideOrdener(geophy2);

	//slide.setPositionsOfLeaves();
	geophy2->setInnerPositions(geophy2->m_tree->m_root);



	for (auto& site: geophy2->m_sites) {
		std::cout << site->m_circular_interval.from() <<" ";
		//std::cout << site->m_position.x() << " ";
		std::cout << site->m_circular_interval.to() << std::endl;
	}



	/*for (size_t i = 0; i < order.value_of_node_at_position.size(); ++i) {
		std::cout << "Leaf " << i << ": ";
		for (size_t j = 0; j < order.value_of_node_at_position[i].size(); ++j) {
			std::cout << "Position " << j << ": ";
			std::cout << order.value_of_node_at_position[i][j] << " ";
		}
		std::cout << std::endl;
	}


	 */

	/*for (auto& leaf: geophy2->m_tree->m_nodes) {
		std::cout << leaf->m_ID << std::endl;
		std::cout << leaf->polar_position.phi() << std::endl;
	}*/


	painting1 = std::make_shared<CircularPainting>(geophy2);

	QApplication a(argc, argv);
	a.setApplicationName("CartoCrow GUI");

	cartocrow::renderer::GeometryWidget widget(painting1);


	widget.show();
	return a.exec();
}



