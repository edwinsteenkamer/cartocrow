#include "../catch.hpp"

#include <cmath>
#include <ctime>

#include <QApplication>

#include "../../cartocrow/flow_map/painting.h"
#include "../../cartocrow/flow_map/spiral_tree.h"
#include "../../cartocrow/flow_map/spiral_tree_obstructed_algorithm.h"
#include "../../cartocrow/renderer/geometry_painting.h"
#include "../../cartocrow/renderer/geometry_widget.h"
#include "../../cartocrow/renderer/ipe_renderer.h"
#include "cartocrow/renderer/geometry_widget.h"

using namespace cartocrow;
using namespace cartocrow::flow_map;

TEST_CASE("Computing a spiral tree with one obstacle") {
	auto tree = std::make_shared<SpiralTree>(Point<Inexact>(0, 0), 0.5061454830783556);
	tree->addPlace("p1", Point<Inexact>(0, 400), 1);
	Polygon<Inexact> obstacle;
	obstacle.push_back(Point<Inexact>(0, 50));
	obstacle.push_back(Point<Inexact>(8, 95));
	obstacle.push_back(Point<Inexact>(50, 140));
	obstacle.push_back(Point<Inexact>(-43, 134));
	obstacle.push_back(Point<Inexact>(-50, 100));
	tree->addObstacle(obstacle);
	SpiralTreeObstructedAlgorithm algorithm(tree);
	algorithm.run();

	// TODO add assertions
}