//
// Created by s165558 on 18-7-2024.
//

#ifndef CARTOCROW_TREE_H
#define CARTOCROW_TREE_H

#include <vector>
#include <iostream>
#include <string>

#include "node.h"
#include "../colormap.hpp"

namespace cartocrow {
namespace geophylogeny {

struct Tree {

	///Constructor of a binary tree.
	Tree(const std::shared_ptr<Node> root);

	/// Returns a list of the nodes of this binary tree.
	const std::vector<std::shared_ptr<Node>>& nodes() const;

	/// Returns a list of the leaves of this binary tree.
	const std::vector<std::shared_ptr<Node>>& leaves() const;

	/// Returns a list of the inner nodes of this binary tree.
	const std::vector<std::shared_ptr<Node>>& innerNodes() const;

	/// Returns a list of the sites of the leaves of this binary tree.
	std::vector<std::shared_ptr<Site>> sites();

	/// Returns a list of the leaves of the tree ordered from left to right.
	std::vector<std::shared_ptr<Node>> leavesByTreeOrder(std::shared_ptr<Node> node);

	/// Adds the leaves to a list in order of the tree from left to right.
	void addLeavesInOrder(std::shared_ptr<Node> node, std::vector<std::shared_ptr<Node>>& leaves);

	/// Returns a list of the leaves of the tree ordered from right to left.
	std::vector<std::shared_ptr<Node>> leavesByTreeOrderRight(std::shared_ptr<Node> node);

	/// Adds the leaves to a list in order of the tree from right to left.
	void addLeavesInOrderRight(std::shared_ptr<Node> node, std::vector<std::shared_ptr<Node>>& leaves);

	/// List of all nodes in this binary tree.
	std::vector<std::shared_ptr<Node>> m_nodes;

	/// List of all leaves in this binary tree.
	std::vector<std::shared_ptr<Node>> m_leaves;

	/// List of all inner nodes in this binary tree.
	std::vector<std::shared_ptr<Node>> m_inner_nodes;

	/// The root node of this binary tree.
	std::shared_ptr<Node> m_root;

	/// Populates the list of nodes using recursion.
	void addNodes(const std::shared_ptr<Node>& node);

	/// Sets the parent of the nodes in the tree.
	void addParents(const std::shared_ptr<Node>& node);

	/// Populates the list of leaves and inner nodes.
	void addLeavesAndInnerNodes();

	/// Sets the colors of the leaves and their site.
	void setColors(Number<Inexact> color_difference);

	void recursiveColoring(std::shared_ptr<Node> root, std::vector<unsigned char> colormap);

	/// Returns the depth of the tree.
	int depthOfTree(std::shared_ptr<Node> node);

	/// Returns the number of colors needed.
	Number<Inexact> numberOfColors(std::shared_ptr<Node> node, Number<Inexact> color_difference);








};

}
}

#endif //CARTOCROW_TREE_H
