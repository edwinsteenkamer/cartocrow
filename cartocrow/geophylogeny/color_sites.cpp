//
// Created by s165558 on 16-9-2024.
//

#include "color_sites.h"

namespace cartocrow {
namespace geophylogeny {

void ColorSites::colorSitesRectangular(std::vector<std::shared_ptr<Site>> sites, Box boundary_box) {
	for (auto& site: sites) {
		auto [x_norm, y_norm] = normalize(site, boundary_box);
		auto [hue, saturation, brightness] = positionToHSV(x_norm, y_norm);
		auto [r, g, b] = HSVToRGB(hue, saturation, brightness);

		site->m_color = Color{r, g, b};
	}
}

void ColorSites::colorSitesCircular(std::vector<std::shared_ptr<Site>> sites, Number<Inexact> geophylogeny_radius) {
	for (auto& site: sites) {
		auto polar_site = flow_map::PolarPoint(site->m_position);
		auto radius_norm = normalizeRadius(site, geophylogeny_radius - 1);
		float hue = polar_site.phi() * (180 / M_PI);
		float saturation = 1 - 0.5 * radius_norm;
		float brightness = 0.8 + 0.2 * radius_norm;
		auto [r, g, b] = HSVToRGB(hue, saturation, brightness);

		site->m_color = Color{r, g, b};
	}
}

std::tuple<float, float> ColorSites::normalize(std::shared_ptr<Site> site, Box boundary_box) {
	float x_norm = (site->m_position.x() - boundary_box.xmin() - 2) / (boundary_box.xmax() - boundary_box.xmin() - 4);
	float y_norm = ((site->m_position.y() - boundary_box.ymin() - 1) / (boundary_box.ymax() - boundary_box.ymin() - 1.5));

	return {x_norm, y_norm};
}

Number<Inexact> ColorSites::normalizeRadius(std::shared_ptr<Site> site, Number<Inexact> geophylogeny_radius) {
	return flow_map::PolarPoint(site->m_position).r() / geophylogeny_radius;
}

std::tuple<float, float, float> ColorSites::positionToHSV(float x_norm, float y_norm) {
	float hue = x_norm * 350;
	float saturation = 1 - 0.5 * y_norm;
	float brightness = 0.8 + 0.2 * y_norm;

	return {hue, saturation, brightness};
}

std::tuple<int, int ,int> ColorSites::HSVToRGB(float hue, float saturation, float brightness) {
	float c = brightness * saturation;
	float x = c * (1 - std::fabs(std::fmod(hue / 60.0, 2) - 1));
	float m = brightness - c;

	float rPrime = 0, gPrime = 0, bPrime = 0;

	if (0 <= hue && hue < 60) {
		rPrime = c; gPrime = x; bPrime = 0;
	} else if (60 <= hue && hue < 120) {
		rPrime = x; gPrime = c; bPrime = 0;
	} else if (120 <= hue && hue < 180) {
		rPrime = 0; gPrime = c; bPrime = x;
	} else if (180 <= hue && hue < 240) {
		rPrime = 0; gPrime = x; bPrime = c;
	} else if (240 <= hue && hue < 300) {
		rPrime = x; gPrime = 0; bPrime = c;
	} else {
		rPrime = c; gPrime = 0; bPrime = x;
	}

	int r = (rPrime + m) * 255;
	int g = (gPrime + m) * 255;
	int b = (bPrime + m) * 255;

	return {r, g, b};
}

}
}
