file(GLOB NMEheaders ${PROJECT_SOURCE_DIR}/NME/include/*.h)
file(GLOB NMEsrc ${PROJECT_SOURCE_DIR}/NME/src/*.cpp)

add_library(nme STATIC ${NMEsrc})

target_link_libraries(nme ${GSL_LIBRARIES})
 
set_target_properties(nme PROPERTIES PUBLIC_HEADER ${NMEheaders} LINKER_LANGUAGE CXX)
 
install(TARGETS nme
  # IMPORTANT: Add the foo library to the "export-set"
  EXPORT BSGTargets
  RUNTIME DESTINATION "${INSTALL_BIN_DIR}" COMPONENT bin
  LIBRARY DESTINATION "${INSTALL_LIB_DIR}" COMPONENT shlib
  PUBLIC_HEADER DESTINATION "${INSTALL_INCLUDE_DIR}/NME"
    COMPONENT dev)
