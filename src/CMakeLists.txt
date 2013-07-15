cmake_minimum_required(VERSION 2.8)
project(src)

include_directories(lsolvers)

add_executable(${PROJECT_NAME} main.cpp lsolvers/lsolver.h lsolvers/lsolver.cpp
defs.h lsolvers/minres5d.h lsolvers/minres5d.cpp lsolvers/minres5dseq.h lsolvers/minres5dseq.cpp)

#install(TARGETS src RUNTIME DESTINATION bin)

add_subdirectory(lsolvers)

FIND_PACKAGE( OpenMP REQUIRED)
if(OPENMP_FOUND)
message("OPENMP FOUND")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()


