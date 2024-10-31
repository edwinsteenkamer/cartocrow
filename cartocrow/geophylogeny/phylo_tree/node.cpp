//
// Created by s165558 on 10-7-2024.
//

#include "node.h"

namespace cartocrow {
namespace geophylogeny {

Node::Node(const int ID)
	: m_ID(ID), m_clade_size(1) {}

Node::Node(const int ID, const std::shared_ptr<Site>& site)
	: m_ID(ID), m_site(site), m_clade_size(1) {}

Node::Node(const int ID, const std::shared_ptr<Node>& first_child, const std::shared_ptr<Node>& second_child)
	: m_ID(ID), m_first_child(first_child), m_second_child(second_child) {
	m_clade_size = m_first_child->m_clade_size + m_second_child->m_clade_size;
}

Node::ConnectionType Node::getType() const {
	if(m_parent.expired()) {
		return ConnectionType::kRoot;
	} else if (m_first_child == nullptr && m_second_child == nullptr) {
		return ConnectionType::kLeaf;
	} else {    															//if (m_first_child != nullptr && m_second_child != nullptr) {
		return ConnectionType::kInner;
	}
}

void Node::setPolarPosition(flow_map::PolarPoint polarpoint) {
	polar_position = polarpoint;
	setPosition(polar_position.toCartesian());
}

void Node::setPosition(Point<Inexact> position) {
	m_position = position;
}




void Node::setYCoordinate(Number<Inexact> y) {
	 m_y = y;
}




}
}

/*using namespace cartocrow;
using namespace cartocrow::geophylogeny;
*/
/*int main(int argc, char* argv[]) {
	std::cout <<"Hello" << std::endl;

	Point<Inexact> p = Point<Inexact>(5.6565656, 5.2525252);
	Point<Inexact> q = Point<Inexact>(7.6565656, 4.2525252);
	Point<Inexact> r = Node::determineParentPosition(p, q);
	std::cout << r << std::endl;


	std::shared_ptr<Node> v_1 = std::make_shared<Node>(1);
	std::shared_ptr<Node> v_2 = std::make_shared<Node>(2);
	std::cout << "Hello" << std::endl;
	std::shared_ptr<Node> v_3 = std::make_shared<Node>(3, v_1, v_2);
	std::cout << "Hello" << std::endl;
	v_2->m_position = p;
	v_1->m_position = q;
	//Node::addRelation(v_1, v_2, v_3);
	std::cout << v_3->m_position << std::endl;
	v_3->setPosition();
	std::cout << v_3->m_position << std::endl;
	//v_2->cladeSize();
	//v_3->cladeSize();
	//v_1->cladeSize();
	std::cout << v_3->m_clade_size << std::endl;
	//if (auto parentShared = v_1->m_parent.lock()) {
	//	std::cout << parentShared->m_ID << std::endl;
	//};
	std::cout << v_1->m_parent->m_ID << std::endl;


	return 0;
}*/

