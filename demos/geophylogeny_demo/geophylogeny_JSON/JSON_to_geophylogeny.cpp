//
// Created by s165558 on 3-9-2024.
//

#include "JSON_to_geophylogeny.h"



using namespace cartocrow;
using namespace cartocrow::geophylogeny;
using json = nlohmann::json;

std::shared_ptr<Tree> JSONToGeophylogeny::readTreeFromJSON(std::filesystem::path filepath, std::vector<std::shared_ptr<Site>> sites) {
	std::ifstream file(filepath);
	json jsonGeophylogeny = json::parse(file);

	json jsonTree = jsonGeophylogeny["tree"];
	std::shared_ptr<Node> root = parseJsonVertex(jsonTree, sites);
	std::shared_ptr<Tree> tree = std::make_shared<Tree>(root);


	return tree;
}

std::shared_ptr<Node> JSONToGeophylogeny::parseJsonVertex(json& jsonVertex, std::vector<std::shared_ptr<Site>> sites) {
	int id = jsonVertex["id"].get<int>();
	if (jsonVertex["leaf"].get<bool>()) {
		std::shared_ptr<Node> leaf = std::make_shared<Node>(id, sites[jsonVertex["site_id"].get<int>()]);
		return leaf;
	} else {
		std::shared_ptr<Node> firstChild = parseJsonVertex(jsonVertex["left"], sites);
		std::shared_ptr<Node> secondChild = parseJsonVertex(jsonVertex["right"], sites);
		std::shared_ptr<Node> inner_node = std::make_shared<Node>(id, firstChild, secondChild);
		return inner_node;
	}
}

std::vector<std::shared_ptr<Site>> JSONToGeophylogeny::readSitesFromJSON(std::filesystem::path filepath) {
	std::ifstream f(filepath);
	json jsonGeophylogeny = json::parse(f);
	int numSites = jsonGeophylogeny["num_sites"].get<int>();
	json jsonSites = jsonGeophylogeny["sites"];

	std::vector<std::shared_ptr<Site>> sites;

	for (int i = 0; i < numSites; i++) {
		double x = jsonSites[i]["x"].get<double>();
		double y = jsonSites[i]["y"].get<double>();
		std::shared_ptr<Site> site = std::make_shared<Site>(x, y, Color{0, 0, 0});
		sites.push_back(site);
	}
	return sites;
}


