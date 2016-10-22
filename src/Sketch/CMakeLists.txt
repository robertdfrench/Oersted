PROJECT(Oersted_Sketch)

set(SOURCE_FILES
        ./include/Sketch.hpp

        ./src/Sketch.h  ./src/Sketch.cpp

        ./src/Vertex.h  ./src/Vertex.cpp

        ./src/Curve.h   ./src/Curve.cpp
            ./src/LineSegment.h ./src/LineSegment.cpp
            ./src/CircularArc.h ./src/CircularArc.cpp

        ./src/Constraint.h  ./src/Constraint.cpp
            ./src/Angle.h           ./src/Angle.cpp
            ./src/Coincident.cpp    ./src/Coincident.h
            ./src/Distance.cpp      ./src/Distance.h
            ./src/Fixation.cpp      ./src/Fixation.h
            ./src/Horizontal.cpp    ./src/Horizontal.h
            ./src/Length.cpp        ./src/Length.h
            ./src/Radius.cpp        ./src/Radius.h
            ./src/Rotation.cpp      ./src/Rotation.h
            ./src/Symmetry.cpp      ./src/Symmetry.h
            ./src/Tangency.cpp      ./src/Tangency.h
            ./src/Vertical.cpp      ./src/Vertical.h

        ./src/Pattern.h ./src/Pattern.cpp
            ./src/MirrorCopy.cpp    ./src/MirrorCopy.h
            ./src/RotateCopy.cpp    ./src/RotateCopy.h

        ./src/Branch.h  ./src/Branch.cpp
        ./src/Star.h    ./src/Star.cpp
        ./src/Constellation.h   ./src/Constellation.cpp

        ./src/Contour.h ./src/Contour.cpp)

add_library(sketch STATIC ${SOURCE_FILES})

target_link_libraries(sketch ${Boost_LIBRARIES})