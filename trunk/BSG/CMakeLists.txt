include_directories(${GSL_INCLUDE_DIRS})

include_directories(${ROOT_INCLUDE_DIRS})

include_directories(${Boost_INCLUDE_DIR})

include_directories(${PROJECT_SOURCE_DIR}/include)

file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cpp)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.h)

add_executable(BDSM BDSM.cc ${sources} ${headers})

target_link_libraries(BDSM ${GSL_LIBRARIES} ${ROOT_LIBRARIES} ${Boost_LIBRARIES} NME)

install(TARGETS BDSM EXPORT BDSMTargets RUNTIME "${INSTALL_BIN_DIR}" COMPONENT bin)
