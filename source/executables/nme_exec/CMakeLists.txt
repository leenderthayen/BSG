add_executable(nme_exec NME.cc)

target_link_libraries(nme_exec nme)

install(TARGETS nme_exec EXPORT BSGTargets
        RUNTIME DESTINATION bin)
