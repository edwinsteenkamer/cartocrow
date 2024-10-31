set(SOURCES
        phylo_tree/geophylogeny.cpp
        phylo_tree/node.cpp
        phylo_tree/site.cpp
        phylo_tree/tree.cpp
        boundary/boundary_shape.cpp
        boundary/circular_boundary.cpp
        boundary/linear_boundary.cpp
        painting.cpp
        circular_painting.cpp
        linear_geophylogeny.cpp
        circular_geophylogeny.cpp
        fixed_positions/dp_ordener.cpp
        fixed_positions/dp_ordener_rectangular.cpp
        sliding_positions/rectangular_sliding_optimization.cpp
        sliding_positions/circular_sliding_optimization.cpp
        colormap.cpp
        color_sites.cpp
        )
set(HEADERS
        phylo_tree/geophylogeny.h
        phylo_tree/node.h
        phylo_tree/site.h
        phylo_tree/tree.h
        boundary/boundary_shape.h
        boundary/circular_boundary.h
        boundary/linear_boundary.h
        painting.h
        circular_painting.h
        linear_geophylogeny.h
        circular_geophylogeny.h
        fixed_positions/dp_ordener.h
        fixed_positions/dp_ordener_rectangular.h
        sliding_positions/rectangular_sliding_optimization.h
        sliding_positions/circular_sliding_optimization.h
        colormap.hpp
        color_sites.h
        )


add_library(geophylogeny ${SOURCES})
target_link_libraries(geophylogeny
        PUBLIC core
        PUBLIC flow_map
        PUBLIC necklace_map
        PRIVATE glog::glog
        )

cartocrow_install_module(geophylogeny)
install(FILES ${HEADERS} DESTINATION ${CARTOCROW_INSTALL_DIR}/geophylogeny)
