set(SOURCES
        rectangular_geophylogeny.cpp
        )


add_executable(rectangular_geophylogeny ${SOURCES})

#add_library(geophylogeny_demo ${SOURCES})

target_link_libraries(
        rectangular_geophylogeny
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

install(TARGETS rectangular_geophylogeny DESTINATION ${INSTALL_BINARY_DIR})
