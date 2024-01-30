cmake_minimum_required(VERSION 3.10)
project(YourProjectName)

# Set C++ standard
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Compiler flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -lm -fopenmp -O3")

# Include directories
include_directories(includes)
include_directories(external/eigen)

# Find source files
file(GLOB SOURCES "src/*.cpp")

# Add executable
add_executable(main ${SOURCES})

# Find and link Boost libraries
find_package(Boost REQUIRED COMPONENTS system filesystem)
if(Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIRS})
    target_link_libraries(main ${Boost_LIBRARIES})
endif()


