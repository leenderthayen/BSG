file(GLOB NMEheaders ./include/*.h)
file(GLOB NMEsrc ./src/*.cpp)

add_library(NME SHARED ${NMEheaders} ${NMEsrc})
 
set_target_properties(NME PROPERTIES
  PUBLIC_HEADER ${NMEheaders})
 
install(TARGETS NME
  # IMPORTANT: Add the foo library to the "export-set"
  EXPORT BDSMTargets
  RUNTIME DESTINATION "${INSTALL_BIN_DIR}" COMPONENT bin
  LIBRARY DESTINATION "${INSTALL_LIB_DIR}" COMPONENT shlib
  PUBLIC_HEADER DESTINATION "${INSTALL_INCLUDE_DIR}/_NME"
    COMPONENT dev)
