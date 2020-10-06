/*
The Flow Map library implements the algorithmic geo-visualization
method by the same name, developed by Kevin Verbeek, Kevin Buchin,
and Bettina Speckmann at TU Eindhoven
(DOI: 10.1007/s00453-013-9867-z & 10.1109/TVCG.2011.202).
Copyright (C) 2019  Netherlands eScience Center and TU Eindhoven

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <https://www.gnu.org/licenses/>.

Created by tvl (t.vanlankveld@esciencecenter.nl) on 04-09-2020
*/

#ifndef GEOVIZ_FLOW_MAP_NODE_H
#define GEOVIZ_FLOW_MAP_NODE_H

#include <memory>

#include "geoviz/common/core_types.h"
#include "geoviz/common/polar_point.h"


namespace geoviz
{
namespace flow_map
{

struct Node
{
  using Ptr = std::shared_ptr<Node>;

  Node(const std::string& id, const PolarPoint& position);

  std::string id;

  PolarPoint position;

  Number flow_in;
}; // struct Node

} // namespace flow_map
} // namespace geoviz

#endif //GEOVIZ_FLOW_MAP_NODE_H
