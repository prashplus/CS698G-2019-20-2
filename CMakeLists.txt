cmake_minimum_required(VERSION 2.8.9)

project(CloudMPI)

set(CMAKE_C_COMPILER /usr/bin/mpicc)
set(CMAKE_CXX_COMPILER /usr/bin/mpicxx)
set(SOURCES src/main.cpp src/cloudmpi.cpp)
############################################################
# Create a library
############################################################

#Generate the static library from the library sources
add_library(CloudMPI STATIC
        src/cloudmpi.cpp
        )

target_include_directories(CloudMPI
        PUBLIC
        ${PROJECT_SOURCE_DIR}/include
        )


############################################################
# Create an executable
############################################################

# Add an executable with the above sources
add_executable(main
            src/main.cpp
        )

# link the new hello_library target with the hello_binary target
target_link_libraries(main
        PRIVATE
            CloudMPI
        )