include_directories(${GSL_INCLUDE_DIRS})

include_directories(${ROOT_INCLUDE_DIRS})

include_directories(${Boost_INCLUDE_DIR})

include_directories(${PROJECT_SOURCE_DIR}/include)

file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.h)

add_library(bsg ${headers} ${sources})

target_link_libraries(bsg ${GSL_LIBRARIES} ${ROOT_LIBRARIES} ${Boost_LIBRARIES} nme)

install(TARGETS bsg EXPORT BSGTargets RUNTIME "${INSTALL_BIN_DIR}" COMPONENT bin)
