//
// Created by s165558 on 20-9-2024.
//
#include <memory>
#include <iostream>
#include <vector>
#include <algorithm>
#include <string>
#include <filesystem>
#include <numeric>
#include <cmath>
#include <fstream>
#include <regex>
#include <CGAL/squared_distance_2.h>


#include <QApplication>

#include <nlohmann/json.hpp>

#include "demos/geophylogeny_demo/geophylogeny_JSON/JSON_to_geophylogeny.h"
#include "cartocrow/core/core.h"
#include "cartocrow/flow_map/polar_point.h"
#include "cartocrow/geophylogeny/circular_geophylogeny.h"
#include "cartocrow/geophylogeny/circular_painting.h"
#include "cartocrow/geophylogeny/colormap.hpp"
#include "cartocrow/geophylogeny/fixed_positions/dp_ordener.h"
#include "cartocrow/geophylogeny/fixed_positions/dp_ordener_rectangular.h"
#include "cartocrow/geophylogeny/linear_geophylogeny.h"
#include "cartocrow/geophylogeny/painting.h"
#include "cartocrow/geophylogeny/phylo_tree/node.h"
#include "cartocrow/geophylogeny/phylo_tree/tree.h"
#include "cartocrow/geophylogeny/sliding_positions/rectangular_sliding_optimization.h"
#include "cartocrow/geophylogeny\sliding_positions/circular_sliding_optimization.h"
#include "cartocrow/renderer/geometry_painting.h"
#include "cartocrow/renderer/geometry_widget.h"
#include "data_collection.h"

using namespace cartocrow;
using namespace cartocrow::geophylogeny;
using json = nlohmann::json;

static double calculateStandardDeviation(const std::vector<double>& data) {
	if (data.empty()) {
		return 0.0; // Edge case for empty input
	}

	// Calculate the mean
	double mean = std::accumulate(data.begin(), data.end(), 0.0) / data.size();
	std::cout << mean << std::endl;

	// Calculate the variance
	double variance = std::accumulate(data.begin(), data.end(), 0.0,
	                                  [mean](double accumulator, double value) {
		                                  return accumulator + (value - mean) * (value - mean);
	                                  }) / data.size();

	// Return the standard deviation (square root of variance)
	return std::sqrt(variance);
}

Number<Inexact> computeAverageQuality(std::vector<std::shared_ptr<Node>> leaves, char quality_measure) {
	Number<Inexact> quality_measure_cumulative = 0;
	for (auto& leaf : leaves) {
		flow_map::PolarPoint polar_site = flow_map::PolarPoint(leaf->m_site->m_position);
		switch (quality_measure) {
			case 'E': {
				quality_measure_cumulative += sqrt(squared_distance(leaf->m_site->m_position, leaf->m_position));
				break;
			}
			case 'H': {
				quality_measure_cumulative += std::abs(leaf->m_site->m_position.x() - leaf->m_position.x());
				break;
			}
			case 'R': {
				Number<Inexact> phi_difference = std::abs(polar_site.phi() - leaf->polar_position.phi());
				quality_measure_cumulative += std::min(phi_difference, M_2xPI - phi_difference) * polar_site.r();
				break;
			}
			case 'A': {
				quality_measure_cumulative += std::abs(leaf->m_site->m_position.x() - leaf->m_position.x()) / std::abs(leaf->m_position.y() - leaf->m_site->m_position.y());
				break;
			}
			case 'B': {
				Number<Inexact> phi_difference_1 = std::abs(polar_site.phi() - leaf->polar_position.phi());
				quality_measure_cumulative += std::min(phi_difference_1, M_2xPI - phi_difference_1) * leaf->polar_position.r() / std::abs(leaf->polar_position.r() - polar_site.r());
				break;
			}
		}
	}
	return quality_measure_cumulative / leaves.size();
}

Number<Inexact> computeMaxQuality(std::vector<std::shared_ptr<Node>> leaves, char distance_type) {
	Number<Inexact> max_value = 0.0;
	for (auto& leaf : leaves) {
		flow_map::PolarPoint polar_site = flow_map::PolarPoint(leaf->m_site->m_position);
		switch (distance_type) {
			case 'E': {
			    Number<Inexact> euclidean_distance = std::sqrt(squared_distance(leaf->m_site->m_position, leaf->m_position));
			    max_value = std::max(max_value, euclidean_distance);
			    break;
			}
		    case 'H': {
			    Number<Inexact> horizontal_distance = std::abs(leaf->m_site->m_position.x() - leaf->m_position.x());
			    max_value = std::max(max_value, horizontal_distance);
			    break;
		    }
		    case 'R': {
			    Number<Inexact> phi_difference = std::abs(polar_site.phi() - leaf->polar_position.phi());
			    Number<Inexact> radial_distance = std::min(phi_difference, M_2xPI - phi_difference) * polar_site.r();
			    max_value = std::max(max_value, radial_distance);
			    break;
		    }
			case 'A': {
					max_value = std::max(
						max_value, std::abs(leaf->m_site->m_position.x() - leaf->m_position.x()) /
									   std::abs(leaf->m_position.y() - leaf->m_site->m_position.y()));
			        break;
			}
			case 'B': {
					Number<Inexact> phi_difference_1 =
						std::abs(polar_site.phi() - leaf->polar_position.phi());
					max_value =
						std::max(max_value, std::min(phi_difference_1, M_2xPI - phi_difference_1) *
												leaf->polar_position.r() /
												std::abs(leaf->polar_position.r() - polar_site.r()));
			        break;
			}
		}
	}
	return max_value;
}

Number<Inexact> computeAverageSpanningQuality(std::vector<std::shared_ptr<Node>> inner_nodes, char distance_type) {
	Number<Inexact> spanning_quality_cumulative = 0;
	for(auto& node: inner_nodes) {
		Tree subtree = Tree(node);
		std::vector<std::shared_ptr<Site>> sites = subtree.sites();
		if (distance_type == 'H') {
			std::vector<std::shared_ptr<Node>> leaves_in_order = subtree.leavesByTreeOrder(node);
			std::vector<Number<Inexact>> sites_x_coordinates = Geophylogeny::getSiteCoordinates(Geophylogeny::getSitePositions(sites), 0);
			auto min_site = *std::min_element(sites_x_coordinates.begin(), sites_x_coordinates.end());
			auto max_site = *std::max_element(sites_x_coordinates.begin(), sites_x_coordinates.end());
			spanning_quality_cumulative += std::abs(min_site - leaves_in_order.front()->m_position.x()) +
			                               std::abs(max_site - leaves_in_order.back()->m_position.x());
		}
	}
	Number<Inexact> spanning_quality = spanning_quality_cumulative / inner_nodes.size();
	return spanning_quality;
}

Number<Inexact> computeMaxSpanningQuality(std::vector<std::shared_ptr<Node>> inner_nodes, char distance_type) {
	Number<Inexact> max_spanning_quality = 0;
	for(auto& node: inner_nodes) {
		Tree subtree = Tree(node);
		std::vector<std::shared_ptr<Site>> sites = subtree.sites();
		if (distance_type == 'H') {
			std::vector<std::shared_ptr<Node>> leaves_in_order = subtree.leavesByTreeOrder(node);
			std::vector<Number<Inexact>> sites_x_coordinates = Geophylogeny::getSiteCoordinates(Geophylogeny::getSitePositions(sites), 0);
			auto min_site = *std::min_element(sites_x_coordinates.begin(), sites_x_coordinates.end());
			auto max_site = *std::max_element(sites_x_coordinates.begin(), sites_x_coordinates.end());
			max_spanning_quality = std::max(max_spanning_quality, std::abs(min_site - leaves_in_order.front()->m_position.x()) +
			                               std::abs(max_site - leaves_in_order.back()->m_position.x()));
		}
	}
	return max_spanning_quality;
}

std::tuple<double, double, double, double>  computeAverageDistanceForDirectory(const std::filesystem::path& directory) {
	std::vector<double> results_rect_fixed;
	std::vector<double> results_rect_sliding;
	std::vector<double> results_circ_fixed;
	std::vector<double> results_circ_sliding;
	for (const auto& entry: std::filesystem::directory_iterator(directory)) {
		std::filesystem::path projectFilename = entry.path();

		std::ifstream file(projectFilename);
		json jsonGeophylogeny = json::parse(file);

		std::vector<std::shared_ptr<Site>> sites =
		    JSONToGeophylogeny::readSitesFromJSON(projectFilename);

		std::shared_ptr<Tree> tree = JSONToGeophylogeny::readTreeFromJSON(projectFilename, sites);

		json jsonParameters = jsonGeophylogeny["parameters"];

		int num_cycles = jsonParameters["num_cycles"].get<int>();
		Number<Inexact> aversion_centroid_ratio =
		    jsonParameters["aversion_centroid_ratio"].get<double>();
		Number<Inexact> interval_margin_rectangular =
		    jsonParameters["interval_margin_rectangular"].get<double>();
		Number<Inexact> interval_margin_circular =
		    jsonParameters["interval_margin_circular"].get<double>();
		Number<Inexact> color_difference = jsonParameters["color_difference"].get<double>();
		bool allowed_outside_interval = false;

		std::shared_ptr<RectangularGeophylogeny> geophy = std::make_shared<RectangularGeophylogeny>(tree, sites, RectangularGeophylogeny::PositionType::fixed,
		                                                                                            color_difference, jsonParameters["color_distance"]);

		DPOrdenerRectangular order = DPOrdenerRectangular(geophy, DPOrdenerRectangular::DPStrategy::kHorizontal);
		order.orderLeaves();
		results_rect_fixed.push_back(computeAverageQuality(tree->leaves(), 'A'));

		RectangularSlideOrdener slide = RectangularSlideOrdener(geophy, num_cycles, aversion_centroid_ratio,
																interval_margin_rectangular, allowed_outside_interval);
		slide.setPositionsOfLeaves();
		results_rect_sliding.push_back(computeAverageQuality(tree->leaves(), 'A'));

		std::shared_ptr<CircularGeophylogeny> geophy1 = std::make_shared<CircularGeophylogeny>(tree, sites, CircularGeophylogeny::PositionType::fixed, color_difference,
															 jsonParameters["color_distance"]);
		DPOrdener order1 = DPOrdener(geophy1, DPOrdener::DPStrategy::kRadial);
		order1.orderLeaves();
		results_circ_fixed.push_back(computeAverageQuality(tree->leaves(), 'B'));

		CircularSlideOrdener slide1 = CircularSlideOrdener(geophy1, num_cycles, aversion_centroid_ratio, interval_margin_circular, allowed_outside_interval);
		slide1.setPositionsOfLeaves();
		results_circ_sliding.push_back(computeAverageQuality(tree->leaves(), 'B'));
	}
	double avg_rect_fixed = std::accumulate(results_rect_fixed.begin(), results_rect_fixed.end(), 0.0) / results_rect_fixed.size();
	double avg_rect_sliding = std::accumulate(results_rect_sliding.begin(), results_rect_sliding.end(), 0.0) / results_rect_sliding.size();
	double avg_circ_fixed = std::accumulate(results_circ_fixed.begin(), results_circ_fixed.end(), 0.0) / results_circ_fixed.size();
	double avg_circ_sliding = std::accumulate(results_circ_sliding.begin(), results_circ_sliding.end(), 0.0) / results_circ_sliding.size();

	return {avg_rect_fixed, avg_circ_fixed, avg_rect_sliding, avg_circ_sliding};
}

std::string extractNumber(const std::string& directory) {
	std::regex re("\\d+");
	std::smatch match;
	if (std::regex_search(directory, match, re)) {
		return match.str();
	}
	return "";
}

void generateCSV(const std::vector<std::tuple<std::string, double, double, double, double>>& data, const std::filesystem::path& csv_filename) {
	std::ofstream csv_file(csv_filename);
	if (!csv_file.is_open()) {
		std::cerr << "Failed to open or create CSV file: " << csv_filename << std::endl;
		return;
	}


	// Write CSV header
	csv_file << "n, Fixed rect., Fixed circ., Sliding rect., Sliding circ.\n";

	// Write data for each directory
	for (const auto& [dir_name, avg_1, avg_2, avg_3, avg_4] : data) {
		csv_file << dir_name << ", " << avg_1 << ", " << avg_2 << ", " << avg_3 << ", " << avg_4 << "\n";
	}

	csv_file.close();
	std::cout << "CSV file created: " << csv_filename << std::endl;
}


int main(int argc, char* argv[]) {

	const std::filesystem::path parent_directory = argv[1];
	std::vector<std::tuple<std::string, double, double, double, double>> directory_averages;

	for (const auto& dir_entry: std::filesystem::directory_iterator(parent_directory)) {
		std::filesystem::path directory = dir_entry.path();

		auto [avg_1, avg_2, avg_3, avg_4] = computeAverageDistanceForDirectory(directory);

		std::string directory_number = extractNumber(directory.filename().string());
		directory_averages.emplace_back(directory_number, avg_1, avg_2, avg_3, avg_4);
	}
	std::string correlation = parent_directory.parent_path().filename().string();

	std::string datatype = parent_directory.parent_path().parent_path().filename().string();

	std::cout << datatype << std::endl;

	generateCSV(directory_averages, "metric-results/" + datatype + "/force_fixed/average_angle_" + correlation + ".csv");


}




