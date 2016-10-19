PROJECT(Oersted_Sketch)

set(SOURCE_FILES
        ./include/Sketch.hpp ./src/Sketch.h ./src/Sketch.cpp

        ./src/Vertex.h ./src/Vertex.cpp

        ./src/Curve.h ./src/Curve.cpp
            ./src/LineSegment.h ./src/LineSegment.cpp
            ./src/CircularArc.h ./src/CircularArc.cpp

        ./src/Constraint.h ./src/Constraint.cpp
        ./src/Pattern.h ./src/Pattern.cpp

        ./src/Branch.h ./src/Branch.cpp
        ./src/Star.h ./src/Star.cpp
        ./src/Constellation.h ./src/Constellation.cpp

        ./src/Contour.h ./src/Contour.cpp)

add_library(sketch STATIC ${SOURCE_FILES})