set(SOURCES
        geophylogeny_demo.cpp
        JSON_to_geophylogeny.cpp
        )

set(HEADERS
        JSON_to_geophylogeny.h
        )

add_executable(geophylogeny_demo ${SOURCES})

find_package(nlohmann_json REQUIRED)

#add_library(geophylogeny_demo ${SOURCES})

target_link_libraries(
        geophylogeny_demo
        PRIVATE
        ${COMMON_CLA_TARGET}
        core
        #flow_map
        geophylogeny
        renderer
        CGAL::CGAL
        Qt5::Widgets
        glog::glog
        nlohmann_json::nlohmann_json
)

install(TARGETS geophylogeny_demo DESTINATION ${INSTALL_BINARY_DIR})