//
// Created by s165558 on 26-7-2024.
//
#include <memory>
#include <iostream>
#include <vector>
#include <algorithm>
#include <string>
#include <cstdlib>  // For rand and srand
#include <ctime>    // For time


#include <QApplication>

#include <nlohmann/json.hpp>

//#include "geophylogeny_demo.h"
#include "cartocrow/geophylogeny/phylo_tree/node.h"
#include "cartocrow/geophylogeny/phylo_tree/tree.h"
#include "cartocrow/geophylogeny/painting.h"
#include "cartocrow/geophylogeny/circular_painting.h"
#include "cartocrow/core/core.h"
#include "cartocrow/renderer/geometry_painting.h"
#include "cartocrow/renderer/geometry_widget.h"
#include "cartocrow/renderer/ipe_renderer.h"
#include "cartocrow/geophylogeny/linear_geophylogeny.h"
#include "cartocrow/geophylogeny/circular_geophylogeny.h"
#include "cartocrow/flow_map/polar_point.h"
#include "cartocrow/geophylogeny/fixed_positions/dp_ordener.h"
#include "cartocrow/geophylogeny/fixed_positions/dp_ordener_rectangular.h"
#include "cartocrow/geophylogeny/sliding_positions/rectangular_sliding_optimization.h"
#include "cartocrow/geophylogeny\sliding_positions/circular_sliding_optimization.h"
#include "JSON_to_geophylogeny.h"
#include "cartocrow/geophylogeny/colormap.hpp"
#include "cartocrow/core/timer.h"

using namespace cartocrow;
using namespace cartocrow::geophylogeny;
using json = nlohmann::json;

void seedRandom() {
	std::srand(std::time(nullptr)); // Seed with current time
}


int main(int argc, char* argv[]) {
	Timer timer;
	std::shared_ptr<renderer::GeometryPainting> painting;
	std::shared_ptr<renderer::GeometryPainting> painting1;
	seedRandom();

	const std::filesystem::path projectFilename = argv[1];

	std::string outputFilename = "";
	if (argc == 3) {
		outputFilename = argv[2];
	}

	//std::cout << projectFilename.string() << std::endl;

	std::ifstream file(projectFilename);
	json jsonGeophylogeny = json::parse(file);

	std::vector<std::shared_ptr<Site>> sites = JSONToGeophylogeny::readSitesFromJSON(projectFilename);

	std::shared_ptr<Tree> tree = JSONToGeophylogeny::readTreeFromJSON(projectFilename, sites);

	json jsonParameters = jsonGeophylogeny["parameters"];

	int num_cycles = jsonParameters["num_cycles"].get<int>();
	Number<Inexact> aversion_centroid_ratio = jsonParameters["aversion_centroid_ratio"].get<double>();
	Number<Inexact> interval_margin_rectangular = jsonParameters["interval_margin_rectangular"].get<double>();
	Number<Inexact> interval_margin_circular = jsonParameters["interval_margin_circular"].get<double>();
	Number<Inexact> color_difference = jsonParameters["color_difference"].get<double>();
	bool allowed_outside_interval = jsonParameters["allowed_outside_interval"].get<bool>();


	if (jsonParameters["tree_type"] == "rectangular") {
		std::shared_ptr<RectangularGeophylogeny> geophy = std::make_shared<RectangularGeophylogeny>(tree, sites, RectangularGeophylogeny::PositionType::fixed, color_difference, jsonParameters["color_distance"]);
		if (jsonParameters["position_type"] == "fixed") {
			DPOrdenerRectangular order = DPOrdenerRectangular(geophy, DPOrdenerRectangular::DPStrategy::kEuclidean);
			order.orderLeaves();

		} else if (jsonParameters["position_type"] == "sliding") {
			RectangularSlideOrdener slide = RectangularSlideOrdener(geophy, num_cycles, aversion_centroid_ratio, interval_margin_rectangular, allowed_outside_interval);
			slide.setPositionsOfLeaves();


		}
		geophy->setInnerPositions(geophy->m_tree->m_root);
		/*std::vector<std::shared_ptr<Node>> leavesInOrder = geophy->m_tree->leavesByTreeOrder(geophy->m_tree->m_root);
		int num_colors = leavesInOrder.size();
		int num_colors_first_half = num_colors / 2;
		int num_colors_second_half = num_colors - num_colors_first_half;
		std::vector<unsigned char> colormap(num_colors * 3);
		std::vector<unsigned char> colormap_first(num_colors_first_half * 3);
		std::vector<unsigned char> colormap_second(num_colors_second_half * 3);
		int colors1 = ColorMap::BrewerQualitative(num_colors_first_half, colormap_first.data(), 0.0f, 6.283185307f, 0.5f, 1.0f, 1.0f);
		int colors2 = ColorMap::BrewerQualitative(num_colors_second_half, colormap_second.data(), 0.0f, 6.283185307f, 0.5f, 1.0f, 0.7f);
		colormap_first.insert(colormap_first.end(), colormap_second.begin(), colormap_second.end());
		//std::merge(colormap_first.begin(), colormap_first.end(), colormap_second.begin(), colormap_second.end(), colormap.begin());
		std::vector<unsigned char> colormap_alt(num_colors * 3);
		int colors = ColorMap::BrewerQualitative(num_colors, colormap_alt.data(), 0.0f, 6.283185307f, 0.4f, 2.0f, 0.8f);
		int i = 0;
		int r;
		int g;
		int b;
		for (auto& leaf: leavesInOrder) {
			r = int(colormap_alt[i]);
			i = i + 1;
			g = int(colormap_alt[i]);
			i = i + 1;
			b = int(colormap_alt[i]);
			i = i + 1;
			leaf->m_site->m_color = Color(r, g ,b);
		}*/
		/*for (auto& site: geophy->m_sites) {
			int r = std::rand() % 256;
			int g = std::rand() % 256;
			int b = std::rand() % 256;
			site->m_color = Color(r, g, b);
		}*/
		timer.stamp("Order tree");
		painting = std::make_shared<Painting>(geophy);

	} else if (jsonParameters["tree_type"] == "circular") {
		std::shared_ptr<CircularGeophylogeny> geophy1 = std::make_shared<CircularGeophylogeny>(tree, sites, CircularGeophylogeny::PositionType::fixed, color_difference, jsonParameters["color_distance"]);
		if (jsonParameters["position_type"] == "fixed") {
			DPOrdener order = DPOrdener(geophy1, DPOrdener::DPStrategy::kEuclidean);
			order.orderLeaves();
		} else if (jsonParameters["position_type"] == "sliding") {
			CircularSlideOrdener slide = CircularSlideOrdener(geophy1, num_cycles, aversion_centroid_ratio, interval_margin_circular, allowed_outside_interval);
			slide.setPositionsOfLeaves();
			for (auto& leaf: tree->leaves()) {
				//std::cout << std::sqrt(squared_distance(leaf->m_site->m_position, leaf->m_position)) << std::endl;
			}
		}
		geophy1->setInnerPositions(geophy1->m_tree->m_root);
		timer.stamp("Order tree");
		painting = std::make_shared<CircularPainting>(geophy1);
	}



	if (outputFilename == "") {
		QApplication a(argc, argv);
		a.setApplicationName("CartoCrow GUI");

		cartocrow::renderer::GeometryWidget widget(painting);
		widget.show();
		timer.output();
		return a.exec();

	} else {
		cartocrow::renderer::IpeRenderer renderer(painting);
		renderer.save(outputFilename);
		timer.stamp("drawing");
		timer.output();
		return 0;
	}

	/*

	//std::shared_ptr<Geophylogeny> geophy = std::make_shared<Geophylogeny>(tree, sites, Geophylogeny::BoundaryType::linear, Geophylogeny::PositionType::fix);
	std::shared_ptr<RectangularGeophylogeny> geophy1 = std::make_shared<RectangularGeophylogeny>(tree, sites, RectangularGeophylogeny::PositionType::fixed);


	RectangularSlideOrdener slide = RectangularSlideOrdener(geophy1);
	slide.setPositionsOfLeaves();
	geophy1->setInnerPositions(geophy1->m_tree->m_root);


	std::shared_ptr<CircularGeophylogeny> geophy2 = std::make_shared<CircularGeophylogeny>(tree, sites, CircularGeophylogeny::PositionType::fixed);

	//DPOrdener order = DPOrdener(geophy2, DPOrdener::DPStrategy::kEuclidean);
	//order.orderLeaves();

	CircularSlideOrdener slide = CircularSlideOrdener(geophy2);
	slide.setPositionsOfLeaves();
	geophy2->setInnerPositions(geophy2->m_tree->m_root);


	painting = std::make_shared<Painting>(geophy);
	//painting1 = std::make_shared<CircularPainting>(geophy2);



	 */
}


