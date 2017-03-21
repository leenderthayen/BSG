file(GLOB sources ${PROJECT_SOURCE_DIR}/BSG/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/BSG/include/*.h)

add_library(bsg ${headers} ${sources})

target_link_libraries(bsg ${GSL_LIBRARIES} ${ROOT_LIBRARIES} ${Boost_LIBRARIES} nme)

#install(TARGETS bsg EXPORT BSGTargets RUNTIME "${INSTALL_BIN_DIR}" COMPONENT bin)
