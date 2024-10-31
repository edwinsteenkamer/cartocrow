//
// Created by s165558 on 3-9-2024.
//

#ifndef CARTOCROW_JSONTOGEOPHYLOGENY_H
#define CARTOCROW_JSONTOGEOPHYLOGENY_H

#include <iostream>
#include <string>
#include <vector>
#include <cstdlib>
#include <nlohmann/json.hpp>

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
#include "cartocrow/geophylogeny/colormap.hpp"

using namespace cartocrow;
using namespace cartocrow::geophylogeny;
using json = nlohmann::json;

struct JSONToGeophylogeny {

	static std::shared_ptr<Tree> readTreeFromJSON(std::filesystem::path filepath, std::vector<std::shared_ptr<Site>> sites);

	static std::vector<std::shared_ptr<Site>> readSitesFromJSON(std::filesystem::path filepath);

	static std::shared_ptr<Node> parseJsonVertex(json& jsonVertex, std::vector<std::shared_ptr<Site>> sites);
};

#endif //CARTOCROW_JSONTOGEOPHYLOGENY_H
