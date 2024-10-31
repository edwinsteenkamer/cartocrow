//
// Created by s165558 on 9-8-2024.
//

#include "circular_painting.h"
#include "../core/circular_arc.h"

namespace cartocrow {
namespace geophylogeny {


CircularPainting::CircularPainting(std::shared_ptr<CircularGeophylogeny> circular_geophylogeny)
    : m_circular_geophylogeny(circular_geophylogeny) {}

void CircularPainting::paint(renderer::GeometryRenderer& renderer) const {
	paintBoundary(renderer);
	paintEdges(renderer);
	//paintIntervals(renderer);
	//paintLeaders(renderer);
	paintNodes(renderer);
	paintSites(renderer);
}

void CircularPainting::paintNodes(renderer::GeometryRenderer& renderer) const {
	renderer.setMode(renderer::GeometryRenderer::vertices);
	for (const auto& node: m_circular_geophylogeny->m_tree->nodes()) {
		if (node->getType() == Node::ConnectionType::kLeaf) {
			renderer.setStroke(Color{0, 0, 0}, 3);
			renderer.setStrokeOpacity(255);
			renderer.setFill(node->m_site->m_color);
			auto leaf_drawing = Circle<Inexact>(node->m_position, 0.0075);
			renderer.setMode(renderer::GeometryRenderer::stroke);
			renderer.draw(leaf_drawing);
			renderer.setMode(renderer::GeometryRenderer::fill);
			renderer.draw(leaf_drawing);
		} else {
			renderer.setStroke(Color{0,0,0}, 3);
			renderer.draw(node->m_position);
		}

	}
}

void CircularPainting::paintEdges(renderer::GeometryRenderer& renderer) const {
	renderer.setMode(renderer::GeometryRenderer::stroke);
	renderer.setStroke(Color{0, 0, 0}, 2);
	for (const auto& node: m_circular_geophylogeny->m_tree->nodes()) {
		if (node->getType() != Node::ConnectionType::kLeaf) {
			flow_map::PolarPoint corner_first = flow_map::PolarPoint(node->polar_position.r(), node->m_first_child->polar_position.phi());
			flow_map::PolarPoint corner_second = flow_map::PolarPoint(node->polar_position.r(), node->m_second_child->polar_position.phi());
			Circle<Inexact> supporting_circle = Circle<Inexact>(CGAL::ORIGIN, corner_first.rSquared());
			Number<Inexact> span_angle;

			CircularArc corner_arc = node->firstAsRightChild
				? (corner_second.phi() < corner_first.phi()
					? CircularArc(supporting_circle, corner_first.phi(), corner_second.phi() + M_2xPI - corner_first.phi())
					: CircularArc(supporting_circle, corner_first.phi(), corner_second.phi() - corner_first.phi()))
				: (corner_first.phi() < corner_second.phi()
					? CircularArc(supporting_circle, corner_second.phi(), corner_first.phi() + M_2xPI - corner_second.phi())
					: CircularArc(supporting_circle, corner_second.phi(), corner_first.phi() - corner_second.phi()));

			Segment<Inexact> first_vertical = Segment<Inexact>(node->m_first_child->m_position, corner_first.toCartesian());
			Segment<Inexact> second_vertical = Segment<Inexact>(node->m_second_child->m_position, corner_second.toCartesian());
			renderer.draw(corner_arc);
			renderer.draw(first_vertical);
			renderer.draw(second_vertical);

		}
		flow_map::PolarPoint top = flow_map::PolarPoint(m_circular_geophylogeny->m_tree->m_root->polar_position.r() + 0.5,
		                                                m_circular_geophylogeny->m_tree->m_root->polar_position.phi());
		renderer.draw(Segment<Inexact>(m_circular_geophylogeny->m_tree->m_root->m_position, top.toCartesian()));
	}
}

void CircularPainting::paintSites(renderer::GeometryRenderer& renderer) const {
	for (const auto& site: m_circular_geophylogeny->m_sites) {
		renderer.setStroke(Color{0, 0, 0}, 3);
		renderer.setStrokeOpacity(255);
		renderer.setFill(site->m_color);
		auto site_drawing = Circle<Inexact>(site->m_position, 0.0125);
		renderer.setMode(renderer::GeometryRenderer::stroke);
		renderer.draw(site_drawing);
		renderer.setMode(renderer::GeometryRenderer::fill);
		renderer.draw(site_drawing);
	}
}

void CircularPainting::paintBoundary(renderer::GeometryRenderer& renderer) const {
	renderer.setMode(renderer::GeometryRenderer::stroke);
	renderer.setStroke(Color{0, 0, 0}, 2);
	Circle<Inexact> circular_boundary = m_circular_geophylogeny->boundary;
	renderer.draw(circular_boundary);
}

void CircularPainting::paintIntervals(renderer::GeometryRenderer& renderer) const {
	renderer.setMode(renderer::GeometryRenderer::stroke);
	renderer.setStroke(Color{0, 0, 0}, 0.5);
	for (const auto& site: m_circular_geophylogeny->m_sites) {
		Point<Inexact> left_end = flow_map::PolarPoint(m_circular_geophylogeny->radius, site->m_circular_interval.from()).toCartesian();
		Point<Inexact> right_end = flow_map::PolarPoint(m_circular_geophylogeny->radius, site->m_circular_interval.to()).toCartesian();
		Segment<Inexact> left_side = Segment<Inexact>(site->m_position, left_end);
		Segment<Inexact> right_side = Segment<Inexact>(site->m_position, right_end);
		renderer.draw(left_side);
		renderer.draw(right_side);
	}
}

void CircularPainting::paintLeaders(renderer::GeometryRenderer& renderer) const {
	renderer.setMode(renderer::GeometryRenderer::stroke);
	renderer.setStroke(Color{0, 0, 0}, 0.5);
	for (const auto& leaf: m_circular_geophylogeny->m_tree->m_leaves) {
		Segment<Inexact> leader = Segment<Inexact>(leaf->m_position, leaf->m_site->m_position);
		renderer.draw(leader);
	}
}

}
}
