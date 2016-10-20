project(Oersted_Test)

set(SOURCE_FILES
        main.cpp

        Sketch/test_Sketch.h
        Sketch/test_Vertex.cpp
        Sketch/test_LineSegment.cpp
        Sketch/test_CircularArc.cpp
        Sketch/test_Star.cpp
        Sketch/test_Contour.cpp)

add_executable(run_tests ${SOURCE_FILES})

target_link_libraries(run_tests gtest gtest_main sketch)