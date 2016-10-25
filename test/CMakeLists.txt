project(Oersted_Test)

set(SOURCE_FILES
        main.cpp

        Sketch/test_Sketch.hpp

        Sketch/test_Vertex.cpp
        Sketch/test_LineSegment.cpp
        Sketch/test_CircularArc.cpp
        Sketch/test_Star.cpp
        Sketch/test_Contour.cpp
        Sketch/test_Sketch.cpp
        Sketch/test_Constraint.cpp
        Sketch/test_Pattern.cpp

        Mesh/test_Mesh.hpp
        Mesh/test_Mesh.cpp

        UseCases/test_UseCases.hpp
        UseCases/test_Stator.cpp
        UseCases/test_Rotor.cpp)

add_executable(run_tests ${SOURCE_FILES})

target_link_libraries(run_tests gtest gtest_main sketch mesh)