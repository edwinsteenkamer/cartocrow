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

Created by tvl (t.vanlankveld@esciencecenter.nl) on 10-09-2019
*/

#ifndef GEOVIZ_FLOW_MAP_FLOW_MAP_H
#define GEOVIZ_FLOW_MAP_FLOW_MAP_H

#include <string>

#include <geoviz/common/core_types.h>
#include <geoviz/common/region.h>

#include "geoviz/flow_map/parameters.h"
#include "geoviz/flow_map/io/data_reader.h"
#include "geoviz/flow_map/io/svg_reader.h"
#include "geoviz/flow_map/io/svg_writer.h"
//#include "geoviz/flow_map/io/type_parsers.h"


namespace geoviz
{
namespace flow_map
{

/**@brief Dummy method for running the flow map algorithm.
 * @return a dummy return string.
 */
std::string proc_flow_map();

} // namespace flow_map
} // namespace geoviz

#endif //GEOVIZ_FLOW_MAP_FLOW_MAP_H
