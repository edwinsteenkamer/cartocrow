//
// Created by s165558 on 18-7-2024.
//

#include "tree.h"
#include <iostream>


namespace cartocrow {
namespace geophylogeny {

Tree::Tree(const std::shared_ptr<Node> root)
	: m_root(root) {
	addNodes(root);
	addParents(root);
	addLeavesAndInnerNodes();
}

const std::vector<std::shared_ptr<Node>>& Tree::nodes() const {
	return m_nodes;
}

const std::vector<std::shared_ptr<Node>>& Tree::leaves() const {
	return m_leaves;
}

const std::vector<std::shared_ptr<Node>>& Tree::innerNodes() const {
	return m_inner_nodes;
}

std::vector<std::shared_ptr<Site>> Tree::sites() {
	std::vector<std::shared_ptr<Site>> sites;
	for (auto& leaf: leaves()) {
		sites.push_back(leaf->m_site);
	}
	return sites;
}

std::vector<std::shared_ptr<Node>> Tree::leavesByTreeOrder(const std::shared_ptr<Node> node) {
	std::vector<std::shared_ptr<Node>> leaves;
	addLeavesInOrder(node, leaves);

	return leaves;
}

void Tree::addLeavesInOrder(std::shared_ptr<Node> node, std::vector<std::shared_ptr<Node>>& leaves) {
	if (node->getType() == Node::ConnectionType::kLeaf) {
		leaves.push_back(node);
		return;
	}

	if (node->firstAsLeftChild) {
		addLeavesInOrder(node->m_first_child, leaves);
		addLeavesInOrder(node->m_second_child, leaves);
	} else {
		addLeavesInOrder(node->m_second_child, leaves);
		addLeavesInOrder(node->m_first_child, leaves);
	}
}

std::vector<std::shared_ptr<Node>> Tree::leavesByTreeOrderRight(const std::shared_ptr<Node> node) {
	std::vector<std::shared_ptr<Node>> leaves;
	addLeavesInOrderRight(node, leaves);

	return leaves;
}

void Tree::addLeavesInOrderRight(std::shared_ptr<Node> node, std::vector<std::shared_ptr<Node>>& leaves) {
	if (node->getType() == Node::ConnectionType::kLeaf) {
		leaves.push_back(node);
		return;
	}

	if (node->firstAsRightChild) {
		addLeavesInOrderRight(node->m_first_child, leaves);
		addLeavesInOrderRight(node->m_second_child, leaves);
	} else {
		addLeavesInOrderRight(node->m_second_child, leaves);
		addLeavesInOrderRight(node->m_first_child, leaves);
	}
}

void Tree::addNodes(const std::shared_ptr<Node>& node) {
	if (node == nullptr) {
		return;
	}

	// Recursively add first and second child
	if (node->getType() != Node::ConnectionType::kLeaf) {
		addNodes(node->m_first_child);
		addNodes(node->m_second_child);
	}

	// Add the current node to the list
	m_nodes.push_back(node);
}

void Tree::addParents(const std::shared_ptr<Node>& node) {
	if (node == nullptr) {
		return;
	}

	if (node->getType() == Node::ConnectionType::kLeaf) {
		return;
	} else {
		node->m_first_child->m_parent = node;
		node->m_second_child->m_parent = node;
		addParents(node->m_first_child);
		addParents(node->m_second_child);
	}
}

void Tree::addLeavesAndInnerNodes() {
	for (const auto& node: m_nodes) {
		if (node->getType() == Node::ConnectionType::kLeaf) {
			m_leaves.push_back(node);
		} else {
			m_inner_nodes.push_back(node);
		}
	}
}

void Tree::setColors(Number<Inexact> color_difference) {
	int num_colors = numberOfColors(m_root, color_difference);
	int num_colors_first_half = num_colors / 2;
	int num_colors_second_half = num_colors - num_colors_first_half;
	std::vector<unsigned char> colormap(num_colors_first_half * 3);
	std::vector<unsigned char> colormap_second(num_colors_second_half * 3);
	int colors1 = ColorMap::BrewerQualitative(num_colors_first_half, colormap.data(), 0.0f, 6.283185307f, 0.4f, 2.0f, 1.0f);
	int colors2 = ColorMap::BrewerQualitative(num_colors_second_half, colormap_second.data(), 0.0f, 6.283185307f, 0.4f, 2.0f, 0.75f);
	colormap.insert(colormap.end(), colormap_second.begin(), colormap_second.end());

	std::vector<unsigned char> colormap_alt(num_colors * 3);
	int colors = ColorMap::BrewerQualitative(num_colors, colormap_alt.data(), 0.0f, 6.283185307f, 0.4f, 2.0f, 0.85f);

	// not necesarry at the moment
	std::vector<std::vector<unsigned char>> colorMap2D(num_colors, std::vector<unsigned char>(3));

	for (int i = 0; i < num_colors; i++) {
		colorMap2D[i][0] = colormap[3 * i];
		colorMap2D[i][1] = colormap[3 * i + 1];
		colorMap2D[i][2] = colormap[3 * i + 2];
	}
	recursiveColoring(m_root, colormap);
	/*for (auto& leaf: m_leaves) {
		leaf->m_site->m_color = Color(rand() % 256, rand() % 256, rand() % 256);
	}*/
}

void Tree::recursiveColoring(std::shared_ptr<Node> root, std::vector<unsigned char> colormap) {
	if (root->getType() == Node::ConnectionType::kLeaf) {
		if (colormap.size() == 3) {
			root->m_site->m_color = Color(colormap[0], colormap[1], colormap[2]);
		} else if (colormap.size() > 3) {
			int index = (int(colormap.size()) / 3) / 2 * 3;
			root->m_site->m_color = Color(int(colormap[index]), int(colormap[index + 1]), int(colormap[index + 2]));
		}
	} else {
		std::vector<unsigned char> colors_first_child(root->m_first_child->num_colors * 3);
		std::vector<unsigned char> colors_second_child(root->m_second_child->num_colors * 3);


		std::copy(colormap.begin(), colormap.begin() + root->m_first_child->num_colors * 3, colors_first_child.begin());

		if (root == m_root) {
			int start = (int(colormap.size()) / 3 - root->m_first_child->num_colors - root->m_second_child->num_colors) / 2;
			std::copy(colormap.end() - (root->m_second_child->num_colors + start) * 3 , colormap.end() - start * 3, colors_second_child.begin());
		} else {
			std::copy(colormap.end() - root->m_second_child->num_colors * 3 , colormap.end(), colors_second_child.begin());
		}

		recursiveColoring(root->m_first_child, colors_first_child);
		recursiveColoring(root->m_second_child, colors_second_child);
	}
}

int Tree::depthOfTree(std::shared_ptr<Node> node) {
	if (node->getType() == Node::ConnectionType::kLeaf) {
		return 0;
	}
	int depth = std::max(depthOfTree(node->m_first_child), depthOfTree(node->m_second_child)) + 1;

	return depth;
}

Number<Inexact> Tree::numberOfColors(std::shared_ptr<Node> node, Number<Inexact> color_difference) {
	if (node->getType() == Node::ConnectionType::kLeaf) {
		node->num_colors = 1;
		return node->num_colors;
	}
	int num_colors_first = numberOfColors(node->m_first_child, color_difference);
	int num_colors_second = numberOfColors(node->m_second_child, color_difference);
	int num_colors_total = num_colors_first + num_colors_second + std::max(num_colors_first, num_colors_second) * color_difference;
	node->num_colors = num_colors_total;

	return num_colors_total;
}



}
}



