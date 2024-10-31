//
// Created by s165558 on 24-7-2024.
//

#include "../core/core.h"
#include "painting.h"



namespace cartocrow {
namespace geophylogeny {

Painting::Painting(std::shared_ptr<Geophylogeny> geophylogeny)
	: m_geophylogeny(geophylogeny) {}

Painting::Painting(std::shared_ptr<RectangularGeophylogeny> rectangular_geophylogeny)
	: m_rectangular_geophylogeny(rectangular_geophylogeny) {}

void Painting::paint(renderer::GeometryRenderer& renderer) const {
	paintBoundary(renderer);
	paintEdges(renderer);
	//paintLeaders(renderer);
	//paintIntervals(renderer);
	paintNodes(renderer);
	paintSites(renderer);
}

void Painting::paintNodes(renderer::GeometryRenderer& renderer) const {
	renderer.setMode(renderer::GeometryRenderer::vertices);
	for (const auto& node: m_rectangular_geophylogeny->m_tree->nodes()) {
		if (node->getType() == Node::ConnectionType::kLeaf) {
			renderer.setStroke(Color{0, 0, 0}, 3);
			renderer.setStrokeOpacity(255);
			renderer.setFill(node->m_site->m_color);
			auto leaf_drawing = Circle<Inexact>(node->m_position, 0.0075);
			renderer.setMode(renderer::GeometryRenderer::stroke);
			renderer.draw(leaf_drawing);
			renderer.setMode(renderer::GeometryRenderer::fill);
			renderer.draw(leaf_drawing);
			/*renderer.setStroke(node->m_site->m_color, 2);
			renderer.draw(node->m_position);*/
		} else {
			renderer.setStroke(Color{0,0,0}, 3);
			renderer.draw(node->m_position);
		}

	}
}

void Painting::paintEdges(renderer::GeometryRenderer& renderer) const {
	renderer.setMode(renderer::GeometryRenderer::stroke);
	renderer.setStroke(Color{0, 0, 0}, 2);
	for (const auto& node: m_rectangular_geophylogeny->m_tree->nodes()) {
		if(node->getType() != Node::ConnectionType::kLeaf) {
			Point<Inexact> corner_first = Point<Inexact>(node->m_first_child->m_position.x(), node->m_position.y());
			Point<Inexact> corner_second = Point<Inexact>(node->m_second_child->m_position.x(), node->m_position.y());
			Segment<Inexact> horizontal = Segment<Inexact>(corner_first, corner_second);
			Segment<Inexact> first_vertical = Segment<Inexact>(node->m_first_child->m_position, corner_first);
			Segment<Inexact> second_vertical = Segment<Inexact>(node->m_second_child->m_position, corner_second);
			renderer.draw(horizontal);
			renderer.draw(first_vertical);
			renderer.draw(second_vertical);
		}
		renderer.draw(Segment<Inexact>(m_rectangular_geophylogeny->m_tree->m_root->m_position, Point<Inexact>(m_rectangular_geophylogeny->m_tree->m_root->m_position.x(),
		                                                                                          m_rectangular_geophylogeny->m_tree->m_root->m_position.y() + 0.5)));
	}
}

void Painting::paintSites(renderer::GeometryRenderer& renderer) const {
	for (const auto& site: m_rectangular_geophylogeny->m_sites) {
		renderer.setStroke(Color{0, 0, 0}, 3);
		renderer.setStrokeOpacity(255);
		renderer.setFill(site->m_color);
		auto site_drawing = Circle<Inexact>(site->m_position,  0.0125);
		renderer.setMode(renderer::GeometryRenderer::stroke);
		renderer.draw(site_drawing);
		renderer.setMode(renderer::GeometryRenderer::fill);
		renderer.draw(site_drawing);
		/*renderer.setMode(renderer::GeometryRenderer::vertices);
		renderer.setStroke(site->m_color, 2);
		renderer.draw(site->m_position);*/
	}
}

void Painting::paintBoundary(renderer::GeometryRenderer& renderer) const {
	renderer.setMode(renderer::GeometryRenderer::stroke);
	renderer.setStroke(Color{0, 0, 0}, 2);

	Polygon<Inexact> linear_boundary = m_rectangular_geophylogeny->boundary;
	renderer.draw(linear_boundary);
}

void Painting::paintIntervals(renderer::GeometryRenderer& renderer) const {
	renderer.setMode(renderer::GeometryRenderer::stroke);
	renderer.setStroke(Color{0, 0, 0}, 0.5);
	for (const auto& site: m_rectangular_geophylogeny->m_sites) {
		Point<Inexact> left_end = Point<Inexact>(site->m_interval.from(), m_rectangular_geophylogeny->boundary_box.ymax());
		Point<Inexact> right_end = Point<Inexact>(site->m_interval.to(), m_rectangular_geophylogeny->boundary_box.ymax());
		Segment<Inexact> left_side = Segment<Inexact>(site->m_position, left_end);
		Segment<Inexact> right_side = Segment<Inexact>(site->m_position, right_end);
		renderer.draw(left_side);
		renderer.draw(right_side);
	}
}

void Painting::paintLeaders(renderer::GeometryRenderer& renderer) const {
	renderer.setMode(renderer::GeometryRenderer::stroke);
	renderer.setStroke(Color{0, 0, 0}, 0.5);
	for (const auto& leaf: m_rectangular_geophylogeny->m_tree->m_leaves) {
		Segment<Inexact> leader = Segment<Inexact>(leaf->m_position, leaf->m_site->m_position);
		renderer.draw(leader);
	}
}




}
}

