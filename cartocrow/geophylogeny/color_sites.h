//
// Created by s165558 on 16-9-2024.
//

#include <iostream>
#include <vector>
#include <tuple>
#include <cmath>
#include "phylo_tree/site.h"
#include "phylo_tree/geophylogeny.h"
#include "../core/core.h"

#ifndef CARTOCROW_COLOR_SITES_H
#define CARTOCROW_COLOR_SITES_H

namespace cartocrow {
namespace geophylogeny {

struct ColorSites {

	static void colorSitesRectangular(std::vector<std::shared_ptr<Site>> sites, Box boundary_box);

	static void colorSitesCircular(std::vector<std::shared_ptr<Site>> sites, Number<Inexact> geophylogeny_radius);

	static std::tuple<float, float> normalize(std::shared_ptr<Site> site, Box boundary_box);

	static Number<Inexact> normalizeRadius(std::shared_ptr<Site> site, Number<Inexact> geophylogeny_radius);

	static std::tuple<float, float, float> positionToHSV(float xNorm, float yNorm);

	static std::tuple<int, int ,int> HSVToRGB(float hue, float saturation, float brightness);

	/// The sites being colored.
	//std::vector<std::shared_ptr<Site>> m_sites;





};

}
}

#endif //CARTOCROW_COLOR_SITES_H
