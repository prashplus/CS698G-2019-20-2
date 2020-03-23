cmake_minimum_required(VERSION 2.8.9)

project(CloudMPI)

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
##Bring the headers, such as Student.h into the project
#include_directories(include)
#
## Find MPI
#find_package(MPI REQUIRED)
#
##Can manually add the sources using the set command as follows:
#set(SOURCES src/main.cpp src/cloudmpi.cpp)
#
##However, the file(GLOB...) allows for wildcard additions:
#file(GLOB SOURCES "src/*.cpp")
#
##............................................................................................
#add_library(cloudmpi src/cloudmpi.cpp)
#
#
## Add executables
##add_executable(mainExec main.cpp)
#
## Link against MPI
##target_link_libraries(cloudmpi ${MPI_LIBRARIES})
#
### Include MPI includes
##include_directories(${MPI_INCLUDE_PATH})
##if(MPI_COMPILE_FLAGS)
##    set_target_properties(cloudmpi PROPERTIES
##            COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
##    #    set_target_properties(cloudmpi PROPERTIES
##    #            COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
##endif()
##if(MPI_LINK_FLAGS)
##    set_target_properties(cloudmpi PROPERTIES
##            LINK_FLAGS "${MPI_LINK_FLAGS}")
##    #    set_target_properties(minisat-lib-static PROPERTIES
##    #            LINK_FLAGS "${MPI_LINK_FLAGS}")
##endif()
#
##............................................................................................
#
#add_executable(testmpi src/main.cpp)