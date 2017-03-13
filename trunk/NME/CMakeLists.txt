file(GLOB NMEheaders ./include/*.h)
file(GLOB NMEsrc ./src/*.cpp)

add_library(nme STATIC ${NMEsrc})
 
set_target_properties(nme PROPERTIES
  PUBLIC_HEADER ${NMEheaders})
 
install(TARGETS nme
  # IMPORTANT: Add the foo library to the "export-set"
  EXPORT BSGTargets
  RUNTIME DESTINATION "${INSTALL_BIN_DIR}" COMPONENT bin
  LIBRARY DESTINATION "${INSTALL_LIB_DIR}" COMPONENT shlib
  PUBLIC_HEADER DESTINATION "${INSTALL_INCLUDE_DIR}/NME"
    COMPONENT dev)
