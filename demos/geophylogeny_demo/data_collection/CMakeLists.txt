set(SOURCES
        data_collection.cpp
        ../geophylogeny_JSON/JSON_to_geophylogeny.cpp
        )

add_executable(data_collection ${SOURCES})

#add_library(geophylogeny_demo ${SOURCES})

target_link_libraries(
        data_collection
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

install(TARGETS data_collection DESTINATION ${INSTALL_BINARY_DIR})
