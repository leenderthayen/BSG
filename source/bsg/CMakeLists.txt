file(GLOB sources ${PROJECT_SOURCE_DIR}/BSG/src/*.cpp)

add_library(bsg ${sources})

target_link_libraries(bsg nme ${Boost_LIBRARIES} ${CMAKE_THREAD_LIBS_INIT})

#install(TARGETS bsg EXPORT BSGTargets RUNTIME "${INSTALL_BIN_DIR}" COMPONENT bin)
