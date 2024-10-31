set(SOURCES
    geophy_test.cpp
)

add_executable(geophy_test ${SOURCES})

target_link_libraries(
        geophy_test
        PRIVATE
        ${COMMON_CLA_TARGET}
        core
        flow_map
        geophylogeny
        #renderer
        #CGAL::CGAL
        #Qt5::Widgets
        glog::glog
)

install(TARGETS geophy_test DESTINATION ${INSTALL_BINARY_DIR})