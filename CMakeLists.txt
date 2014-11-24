cmake_minimum_required(VERSION 2.6)
project(rnnlib)

option(RNNLIB_USE_FLOAT "Use single-precision real numbers" ON)

# Non-standard find modules
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH}
  "${CMAKE_SOURCE_DIR}/CMake/Modules/")

# Main executable
add_executable(rnnlib src/Main.cpp src/Random.cpp src/DataExporter.cpp)

# Find NetCDF C++ library
set (NETCDF_CXX "YES")
find_package(NetCDF REQUIRED)
target_link_libraries(rnnlib ${NETCDF_LIBRARIES})

# Find Boost
find_package(Boost REQUIRED)
include_directories(${Boost_INCLUDE_DIRS})

# Should we use
if (RNNLIB_USE_FLOAT)
  add_definitions(-DREAL_FLOATS)
endif (RNNLIB_USE_FLOAT)

# make install
install(TARGETS rnnlib DESTINATION bin)

