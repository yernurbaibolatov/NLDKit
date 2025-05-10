cmake_minimum_required(VERSION 3.10)

project(NLDKit)

# Optionally specify a Python interpreter manually (e.g., from a virtualenv)
set(PYTHON_VENV_PATH "" CACHE PATH "Path to the Python interpreter (e.g., from virtualenv)")

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

# Use system-installed Eigen
find_package(Eigen3 REQUIRED)
message(STATUS "Found Eigen: ${EIGEN3_INCLUDE_DIR}")

if (PYTHON_VENV_PATH)
    set(Python3_EXECUTABLE ${PYTHON_VENV_PATH})
endif()
find_package(Python3 REQUIRED COMPONENTS Interpreter)

if (NOT PYTHON_VENV_PATH)
    message(WARNING "PYTHON_VENV_PATH is not set. Make sure CMake is using the correct Python interpreter that has all required packages (e.g., pandas, matplotlib). You can specify it with -DPYTHON_VENV_PATH=/path/to/venv/bin/python.")
endif()

# Select the active boilerplate here:
# set(BOILERPLATE src/boilerplates/example_rhs.cpp)

add_executable(simulation
    src/main.cpp
    src/core/DynamicalSystem.cpp
    src/core/Integrator.cpp
)
target_include_directories(simulation PRIVATE include)
target_link_libraries(simulation Eigen3::Eigen)

add_executable(test_damped_oscillator
    tests/test_damped_oscillator.cpp
    src/core/Integrator.cpp
    src/core/DynamicalSystem.cpp
)
target_include_directories(test_damped_oscillator PRIVATE include)
target_link_libraries(test_damped_oscillator Eigen3::Eigen)

add_custom_target(run_all
    COMMAND ${CMAKE_COMMAND} --build . --target test_damped_oscillator
    COMMAND ./test_damped_oscillator
    COMMAND ${Python3_EXECUTABLE} ../../python/print_python_info.py
    COMMAND ${Python3_EXECUTABLE} ../../python/plot_damped_oscillator.py
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
    COMMENT "Building, running the C++ test, and plotting with Python."
)