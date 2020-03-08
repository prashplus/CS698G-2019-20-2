
cmake_minimum_required(VERSION 2.8.9)
project(CloudMPI)

#Bring the headers, such as Student.h into the project
include_directories(include)

# Find MPI
find_package(MPI REQUIRED)

#Can manually add the sources using the set command as follows:
set(SOURCES src/main.cpp src/cloudmpi.cpp)

#However, the file(GLOB...) allows for wildcard additions:
file(GLOB SOURCES "src/*.cpp")

#............................................................................................
add_library(cloudmpi ${SOURCES})


# Add executables
#add_executable(mainExec main.cpp)

# Link against MPI
target_link_libraries(cloudmpi ${MPI_LIBRARIES})

# Include MPI includes
include_directories(${MPI_INCLUDE_PATH})
if(MPI_COMPILE_FLAGS)
    set_target_properties(cloudmpi PROPERTIES
            COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
    #    set_target_properties(cloudmpi PROPERTIES
    #            COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
endif()
if(MPI_LINK_FLAGS)
    set_target_properties(cloudmpi PROPERTIES
            LINK_FLAGS "${MPI_LINK_FLAGS}")
    #    set_target_properties(minisat-lib-static PROPERTIES
    #            LINK_FLAGS "${MPI_LINK_FLAGS}")
endif()

#............................................................................................

add_executable(testmpi ${SOURCES})