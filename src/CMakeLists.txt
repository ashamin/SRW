cmake_minimum_required(VERSION 2.8)
project(src)

include_directories(lsolvers)
include_directories(areas)
include_directories(precond)

#install(TARGETS src RUNTIME DESTINATION bin)

add_subdirectory(lsolvers)
add_subdirectory(areas)
add_subdirectory(precond)

add_executable(${PROJECT_NAME} main.cpp lsolvers/lsolver.h lsolvers/lsolver.cpp lsolvers/minres5d.h lsolvers/minres5d.cpp 
lsolvers/minres5dOmpSSOR.h lsolvers/minres5dOmpSSOR.cpp defs.h lsolvers/implexplBouss.h lsolvers/implexplBouss.cpp)

target_link_libraries(${PROJECT_NAME} MathArea2d sampleArea2d BArea)
target_link_libraries(${PROJECT_NAME} SSORpar)
#target_link_libraries(${PROJECT_NAME} lsolver minres5d minres5dOmpSSOR)

FIND_PACKAGE( OpenMP REQUIRED)
if(OPENMP_FOUND)
message("OPENMP FOUND")
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
endif()