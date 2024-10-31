//
// Created by s165558 on 9-8-2024.
//

#ifndef CARTOCROW_CIRCULAR_PAINTING_H
#define CARTOCROW_CIRCULAR_PAINTING_H

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
class CircularPainting : public renderer::GeometryPainting {


  public:


	CircularPainting(std::shared_ptr<CircularGeophylogeny>);

  protected:
	void paint(renderer::GeometryRenderer& renderer) const override;

  private:
	void paintEdges(renderer::GeometryRenderer& renderer) const;
	void paintNodes(renderer::GeometryRenderer& renderer) const;
	void paintSites(renderer::GeometryRenderer& renderer) const;
	void paintBoundary(renderer::GeometryRenderer& renderer) const;
	void paintIntervals(renderer::GeometryRenderer& renderer) const;
	void paintLeaders(renderer::GeometryRenderer& renderer) const;

	std::shared_ptr<CircularGeophylogeny> m_circular_geophylogeny;
};

}
}

#endif //CARTOCROW_CIRCULAR_PAINTING_H
