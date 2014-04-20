project(lsolvers)

set(LSOLVERS_H
    ${CMAKE_CURRENT_SOURCE_DIR}/lsolver.h
    ${CMAKE_CURRENT_SOURCE_DIR}/minres5d.h
    ${CMAKE_CURRENT_SOURCE_DIR}/minres5dOmpSSOR.h
    ${CMAKE_CURRENT_SOURCE_DIR}/implexplBouss.h
PARENT_SCOPE
)

set(LSOLVERS_SRC
    ${CMAKE_CURRENT_SOURCE_DIR}/lsolver.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/minres5d.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/minres5dOmpSSOR.cpp
    ${CMAKE_CURRENT_SOURCE_DIR}/implexplBouss.cpp
PARENT_SCOPE
)

#add_library(${PROJECT_NAME} ${LSOLVERS_H} ${LSOLVERS_CPP})
