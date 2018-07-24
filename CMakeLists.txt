cmake_minimum_required(VERSION 3.9)
project (dualmc)

set (CMAKE_CXX_STANDARD 11)

include_directories("${CMAKE_SOURCE_DIR}/include/")

set(EXAMPLE_APP_SOURCES
    include/dualmc.h
    include/dualmc.cpp
    include/vertex.h
    include/quad.h
    include/edges.h
    include/tables.h
    apps/example/example.cpp
    apps/example/main.cpp
)

set(GENTABLES_APP_SOURCES
    apps/gentables/gentables.cpp
    apps/gentables/main.cpp
)

# build application
add_executable(dmc ${EXAMPLE_APP_SOURCES})
add_executable(gentables ${GENTABLES_APP_SOURCES})
