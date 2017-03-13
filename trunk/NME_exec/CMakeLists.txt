add_executable(NME NME.cc)

target_link_libraries(NME nme)

install(TARGETS NME EXPORT BSGTargets RUNTIME DESTINATION "${INSTALL_BIN_DIR}" COMPONENT bin)
