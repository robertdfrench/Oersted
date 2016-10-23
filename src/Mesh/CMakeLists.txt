PROJECT(Oersted_Mesh)

set(SOURCE_FILES
        ./include/Mesh.hpp

        ./src/Mesh.h    ./src/Mesh.cpp
        ./src/Edge.h    ./src/Edge.cpp)

add_library(mesh STATIC ${SOURCE_FILES})

target_link_libraries(mesh ${Boost_LIBRARIES})