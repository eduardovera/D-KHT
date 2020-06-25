project(D_KHT)

add_definitions(-DDLIB_PNG_SUPPORT)

find_library(DLIB_RELEASE 
    NAMES dlib 
    HINTS ${PROJECT_SOURCE_DIR}/libs/dlib/lib)

find_library(DLIB_DEBUG 
    NAMES dlibd
    HINTS ${PROJECT_SOURCE_DIR}/libs/dlib/lib)

find_library(LIBPNG 
    NAMES libpng16 
    HINTS ${PROJECT_SOURCE_DIR}/libs/libpng/lib)

add_executable(${PROJECT_NAME} 
    # headers
    accum_ball_cell_t.h
    accum_cell_t.h
    accumulatorball_t.h
    bin_t.h
    hough.h
    kernel_t.h
    peak_detection.h
    plane_t.h 
    reader_file.h
    settings.h
    voting.h
    quadtree_t.h
    sat.h
    logger.h

    # sources
    accumulatorball_t.cpp
    hough.cpp
    voting.cpp
    quadtree_t.cpp
    logger.cpp
    main.cpp)

# cmake alternative to Qt for setting current location
set(WORKSPACE_LOCATION "${PROJECT_SOURCE_DIR}")
configure_file(configuration/workspace_location.h.in ${PROJECT_BINARY_DIR}/configuration/workspace_location.h @ONLY)

target_link_libraries(${PROJECT_NAME} 
    optimized ${DLIB_RELEASE}
    debug ${DLIB_DEBUG}
    ${LIBPNG})

target_include_directories(${PROJECT_NAME} 
    PUBLIC libs/dlib/
    PUBLIC .
    PUBLIC ${PROJECT_BINARY_DIR}/configuration)