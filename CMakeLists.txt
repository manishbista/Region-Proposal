cmake_minimum_required(VERSION 2.6)
project (Region-Proposal)
# version number
set (Region-Proposal VERSION_MAJOR 1)
set (Region-Proposal VERSION_MINOR 0)

# add the binary tree to the search path for include files
include_directories("header")

include(FindPkgConfig)
find_package(OpenGL REQUIRED)
find_package(GLEW REQUIRED)
find_package(SDL2 REQUIRED)
# find_package(SDL2_image REQUIRED)
pkg_search_module(SDL2IMAGE REQUIRED SDL2_image>=2.0.0)

include_directories( ${OPENGL_INCLUDE_DIRS} ${SDL2_INCLUDE_DIRS} ${GLEW_INCLUDE_DIRS} ${SDL2IMAGE_INCLUDES})
add_library(Graphics source/disjointSet.cpp source/mesh.cpp source/meshLoader.cpp source/region.cpp source/shader.cpp source/vector3d.cpp source/segment.cpp source/matrices.cpp)

# add executable
add_executable(main source/main.cpp)
target_link_libraries(main ${OPENGL_LIBRARIES} ${SDL2_LIBRARY} ${GLEW_LIBRARY} ${SDL2IMAGE_LIBRARIES} Graphics)

#install
#install(TARGETS Transform DESTINATION bin)
