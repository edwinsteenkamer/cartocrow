//
// Created by s165558 on 10-7-2024.
//

#ifndef CARTOCROW_NODE_H
#define CARTOCROW_NODE_H


#include <vector>
#include <iostream>
#include <string>
#include <cmath>

#include  "../../core/core.h"
#include "site.h"
#include "../../flow_map/polar_point.h"


namespace cartocrow {
namespace geophylogeny {


	struct Node {

		/// The type of node, as defined by its connections.
		enum class ConnectionType {
			/// The root node, the only node without a parent.
			kRoot,
			/// A leaf node, a node without any children.
			kLeaf,
			/// A join node, a node with two children.
			kInner
		};

		/// Constructor of a leaf node.
		Node(const int ID);

		///Constructor of a leaf node with a corresponding site.
		Node(const int ID, const std::shared_ptr<Site>& site);


		///Constructor of a non-leaf node.
		Node(const int ID, const std::shared_ptr<Node>& first_child, const std::shared_ptr<Node>& second_child);

		///The ID of the node
		const int m_ID;

		/// Determines the type of this node, based on its whether it has a parent and/or children.
		///
		/// Each node is either the root, a leaf, or a n inner node
		/// (see \ref ConnectionType).
		ConnectionType getType() const;

		/// Determines the position of this node.
		void setPosition(Point<Inexact> position);

	    /// Sets the position of this node based on a polar position.
	    void setPolarPosition(flow_map::PolarPoint polarpoint);

		void setYCoordinate(Number<Inexact> y);

	    bool firstAsRightChild = true;

	    bool firstAsLeftChild = true;

	    /// Number of colors needed to draw the subtree of this node.
	    int num_colors;

		/// The parent of this node, or \c nullptr if this is the root.
		std::weak_ptr<Node> m_parent;

		/// The first child of this node, or \c nullptr if this is a leaf.
		std::shared_ptr<Node> m_first_child;

		/// The second child of this node, or \c nullptr if this is a leaf.
		std::shared_ptr<Node> m_second_child;

		std::shared_ptr<Site> m_site;

		/// The number of leaves in the subtree of this node, also called the clade of the node
		int m_clade_size;

		/// The position of this node.
		Point<Inexact> m_position;

		/// x-coordinate of this node.
		Number<Inexact> m_x;

		/// y-coordinate of this node.
		Number<Inexact> m_y;

		/// The polar position of this node.
		flow_map::PolarPoint polar_position;
	};


}};

 // namespace cartocrow::geophylogeny_demo

#endif //CARTOCROW_NODE_H
