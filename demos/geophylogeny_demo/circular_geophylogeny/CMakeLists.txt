set(SOURCES
        circular_geophylogeny.cpp
        )

add_executable(circular_geophylogeny ${SOURCES})

#add_library(geophylogeny_demo ${SOURCES})

target_link_libraries(
        circular_geophylogeny
        PRIVATE
        ${COMMON_CLA_TARGET}
        core
        #flow_map
        geophylogeny
        renderer
        CGAL::CGAL
        Qt5::Widgets
        glog::glog
)

install(TARGETS circular_geophylogeny DESTINATION ${INSTALL_BINARY_DIR})
