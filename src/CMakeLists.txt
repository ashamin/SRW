cmake_minimum_required(VERSION 2.8)
project(src)

FIND_PACKAGE( OpenMP REQUIRED)
if(OPENMP_FOUND)
message("OPENMP FOUND")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()

#link_directories(${PROJECT_SOURCE_DIR}/lsolvers)

add_subdirectory(${PROJECT_SOURCE_DIR}/areas)
add_subdirectory(${PROJECT_SOURCE_DIR}/graphics)
add_subdirectory(${PROJECT_SOURCE_DIR}/precond)
add_subdirectory(${PROJECT_SOURCE_DIR}/lsolvers)

include_directories(${PROJECT_SOURCE_DIR}/areas)
include_directories(${PROJECT_SOURCE_DIR}/precond)
include_directories(${PROJECT_SOURCE_DIR}/lsolvers)


set(D_H
    ${AREAS_H}
    ${PRECOND_H}
    ${LSOLVERS_H}
)

set(D_SRC
    ${AREAS_SRC}
    ${PRECOND_SRC}
    ${GRAPHICS_SRC}
    ${LSOLVERS_SRC}
)

add_executable(${PROJECT_NAME}
    main.cpp
    ${D_H}
    ${D_SRC}
)

#target_link_libraries(${PROJECT_NAME} lsolvers)

