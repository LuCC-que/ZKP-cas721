cmake_minimum_required(VERSION 3.15)

# Set your project name
project(ZKP)

# Get the OpenMP
#find_package(OpenMP REQUIRED)

# Specify the C++ standard
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Compiler options for debugging
#set(CMAKE_BUILD_TYPE Debug)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)  
endif()
set(CMAKE_CXX_FLAGS "-Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -march=native -mtune=native")

# Add the executable
add_executable(${PROJECT_NAME} main.cpp Prover.cpp Verifier.cpp)
#target_link_libraries(${PROJECT_NAME} PRIVATE OpenMP::OpenMP_CXX)

# Include directories
target_include_directories(${PROJECT_NAME} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})

# Output directory for binaries
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/../build)
