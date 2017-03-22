file(GLOB sources ${PROJECT_SOURCE_DIR}/BSG/src/*.cpp)
#file(GLOB headers ${PROJECT_SOURCE_DIR}/BSG/include/*.h)

add_library(bsg ${sources})

target_link_libraries(bsg nme ${Boost_LIBRARIES})

#install(TARGETS bsg EXPORT BSGTargets RUNTIME "${INSTALL_BIN_DIR}" COMPONENT bin)
