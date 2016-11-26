PROJECT(Oersted_Sketch)

set(SOURCE_FILES
        ./include/Sketch.hpp

        ./src/Sketch.h  ./src/Sketch.cpp

        ./src/Variable.h ./src/Variable.cpp

        ./src/SketchElement.h ./src/SketchElement.cpp

        ./src/Vertex.h  ./src/Vertex.cpp

        ./src/Curve.h   ./src/Curve.cpp
            ./src/LineSegment.h ./src/LineSegment.cpp
            ./src/CircularArc.h ./src/CircularArc.cpp

        ./src/Constraint.h  ./src/Constraint.cpp
            ./src/Angle.h           ./src/Angle.cpp
            ./src/Coincident.h    ./src/Coincident.cpp
            ./src/Distance.h      ./src/Distance.cpp
            ./src/Fixation.h      ./src/Fixation.cpp
            ./src/Horizontal.h    ./src/Horizontal.cpp
            ./src/Length.h        ./src/Length.cpp
            ./src/Radius.h        ./src/Radius.cpp
            ./src/Rotation.h      ./src/Rotation.cpp
            ./src/Symmetry.h      ./src/Symmetry.cpp
            ./src/Tangency.h      ./src/Tangency.cpp
            ./src/Vertical.h      ./src/Vertical.cpp

        ./src/Pattern.h ./src/Pattern.cpp
            ./src/MirrorCopy.h    ./src/MirrorCopy.cpp
            ./src/RotateCopy.h    ./src/RotateCopy.cpp

        ./src/Branch.h  ./src/Branch.cpp
        ./src/Star.h    ./src/Star.cpp
        ./src/Constellation.h   ./src/Constellation.cpp

        ./src/Contour.h ./src/Contour.cpp src/sPoint.cpp src/sPoint.h)

add_library(sketch STATIC ${SOURCE_FILES})

target_link_libraries(sketch ${Boost_LIBRARIES})