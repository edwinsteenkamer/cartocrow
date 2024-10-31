//
// Created by s165558 on 24-7-2024.
//

#ifndef CARTOCROW_PAINTING_H
#define CARTOCROW_PAINTING_H

#include "../renderer/geometry_painting.h"
#include "../renderer/geometry_renderer.h"
#include "phylo_tree/node.h"
#include "phylo_tree/tree.h"
#include "phylo_tree/geophylogeny.h"
#include "phylo_tree/site.h"
#include "../core/core.h"
#include "linear_geophylogeny.h"
#include "circular_geophylogeny.h"

namespace cartocrow {
namespace geophylogeny {

/// The \ref renderer::GeometryPainting "GeometryPainting" for a geophylogeny.
class Painting : public renderer::GeometryPainting {


  public:

	/// Creates a new painting with the given geophylogeny.
	Painting(std::shared_ptr<Geophylogeny>);

	Painting(std::shared_ptr<RectangularGeophylogeny>);

  protected:
	void paint(renderer::GeometryRenderer& renderer) const override;

  private:
	void paintEdges(renderer::GeometryRenderer& renderer) const;
	void paintNodes(renderer::GeometryRenderer& renderer) const;
	void paintSites(renderer::GeometryRenderer& renderer) const;
	void paintBoundary(renderer::GeometryRenderer& renderer) const;
	void paintIntervals(renderer::GeometryRenderer& renderer) const;
	void paintLeaders(renderer::GeometryRenderer& renderer) const;

	std::shared_ptr<Geophylogeny> m_geophylogeny;

	std::shared_ptr<RectangularGeophylogeny> m_rectangular_geophylogeny;
};

}
}

#endif //CARTOCROW_PAINTING_H
