cmake_minimum_required(VERSION 3.10)
project(NLDKit)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Include directories
include_directories(include)

# Use system-installed Eigen
find_package(Eigen3 REQUIRED)
message(STATUS "Found Eigen: ${EIGEN3_INCLUDE_DIR}")

# Select the active boilerplate here:
# set(BOILERPLATE src/boilerplates/example_rhs.cpp)

add_executable(simulation
    src/main.cpp
    src/core/DynamicalSystem.cpp
    src/core/Integrator.cpp
)
target_link_libraries(simulation Eigen3::Eigen)

add_executable(test_damped_oscillator
    tests/test_damped_oscillator.cpp
    src/core/Integrator.cpp
    src/core/DynamicalSystem.cpp
)
target_link_libraries(test_damped_oscillator Eigen3::Eigen)